#!/usr/bin/env python
# coding: utf-8

# # 4. Interferon Response Capacity (IRC) Scoring and Cell-Level Differential Expression
# 
# This section calculates an Interferon Response Capacity (IRC) score using the six-gene interferon-stimulated protein signature and examines how IRC scores vary across cells and samples.
# 
# Patient cells are then divided into IRC-High and IRC-Low groups using the median **cell-level** IRC score, followed by exploratory differential-expression analysis.
# 
# ## Input
# 
# `results/preprocessing_tables/clustering_res1.0.h5ad`
# 
# ## Main outputs
# 
# ### Figures
# - `results/figures/15_irc_umap.png`
# - `results/figures/16_irc_violin_by_sample.png`
# - `results/figures/17_irc_violin_by_donor_type.png`
# 
# ### Tables
# - `results/preprocessing_tables/irc_signature_availability.csv`
# - `results/preprocessing_tables/irc_score_summary.csv`
# - `results/preprocessing_tables/irc_score_by_sample.csv`
# - `results/preprocessing_tables/patient_irc_group_counts.csv`
# - `results/preprocessing_tables/IRC_High_vs_Low_DE.csv`
# 
# ### Processed AnnData
# - `results/preprocessing_tables/ifn_irc_scored.h5ad`
# - `results/preprocessing_tables/patient_irc_DE.h5ad`
# 
# > **Important interpretation note:** The current stratification is performed at the **cell level**. Cells from patient samples are divided into IRC-High and IRC-Low based on each cell's IRC score. Therefore, the resulting differential-expression comparison is not equivalent to assigning whole patients to IRC-High or IRC-Low groups.
# 

# ## 4.1 Import libraries and locate the project
# 
# No Google Colab, Google Drive, or machine-specific absolute paths are used. The notebook automatically searches for the repository root.
# 

# In[1]:


from pathlib import Path
from importlib.metadata import version

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.set_figure_params(
    dpi=120,
    dpi_save=300,
    facecolor="white",
    frameon=False,
    figsize=(6, 4)
)

def find_project_root():
    current = Path.cwd().resolve()

    for path in [current, *current.parents]:
        if (
            (path / "data").is_dir()
            and (path / "notebooks").is_dir()
            and (path / "results").is_dir()
        ):
            return path

    raise FileNotFoundError(
        "Could not locate the project root. "
        "Run this notebook from inside the sc-project3-IFN-I repository."
    )

project_root = find_project_root()

results_dir = project_root / "results"
figures_dir = results_dir / "figures"
preprocessing_dir = results_dir / "preprocessing_tables"

figures_dir.mkdir(parents=True, exist_ok=True)
preprocessing_dir.mkdir(parents=True, exist_ok=True)

input_file = preprocessing_dir / "clustering_res1.0.h5ad"

print("Project root:", project_root)
print("Scanpy:", sc.__version__)
print("AnnData:", version("anndata"))
print("Pandas:", pd.__version__)
print("Input exists:", input_file.exists())


# ## 4.2 Load the Section 3 clustered dataset
# 
# The input dataset contains PCA, neighborhood graph, UMAP coordinates, Leiden clusters, and `adata.raw` containing the complete normalized/log-transformed expression matrix.
# 

# In[2]:


if not input_file.exists():
    raise FileNotFoundError(
        f"Required input file was not found: {input_file}\n"
        "Run notebooks/03_dimensionality_reduction.ipynb first."
    )

adata = sc.read_h5ad(input_file)

print(adata)
print(f"Cells: {adata.n_obs:,}")
print(f"Genes in HVG matrix: {adata.n_vars:,}")
print("adata.raw available:", adata.raw is not None)


# ## 4.3 Validate required metadata and expression data
# 
# The IRC analysis requires sample metadata, donor type, Leiden clusters, UMAP coordinates, and the full normalized/log-transformed expression matrix stored in `adata.raw`.
# 

# In[3]:


required_obs = ["sample", "donor_type", "leiden_res_1.0"]

missing_obs = [
    col for col in required_obs
    if col not in adata.obs.columns
]

if missing_obs:
    raise KeyError(
        "Required metadata columns are missing: "
        + ", ".join(missing_obs)
    )

if adata.raw is None:
    raise ValueError(
        "adata.raw is missing. Section 2 should preserve the full normalized/log-transformed "
        "matrix before HVG subsetting."
    )

if "X_umap" not in adata.obsm:
    raise KeyError(
        "UMAP coordinates are missing. Run Section 3 before Section 4."
    )

print("Required metadata, UMAP coordinates, and adata.raw are available.")


# ## 4.4 Define the IRC signature
# 
# The six core interferon-stimulated genes used for the IRC signature are:
# 
# - BST2
# - EIF2AK2 (PKR)
# - ISG15
# - MX1
# - IFIT3
# - IRF7
# 

# In[4]:


isg_signature = [
    "BST2",
    "EIF2AK2",
    "ISG15",
    "MX1",
    "IFIT3",
    "IRF7"
]

valid_isgs = [
    gene for gene in isg_signature
    if gene in adata.raw.var_names
]

signature_availability = pd.DataFrame({
    "gene": isg_signature,
    "available_in_raw": [
        gene in adata.raw.var_names
        for gene in isg_signature
    ]
})

signature_availability.to_csv(
    preprocessing_dir / "irc_signature_availability.csv",
    index=False
)

print(
    f"IRC genes available: "
    f"{len(valid_isgs)}/{len(isg_signature)}"
)
print(valid_isgs)

display(signature_availability)


# ## 4.5 Calculate the IRC module score
# 
# `sc.tl.score_genes()` calculates a per-cell module score for the available IRC-signature genes.
# 
# `use_raw=True` ensures that the score is calculated from the full normalized/log-transformed expression matrix rather than only the 2,500 HVGs.
# 

# In[5]:


if len(valid_isgs) == 0:
    raise ValueError(
        "None of the IRC signature genes were found in adata.raw."
    )

sc.tl.score_genes(
    adata,
    gene_list=valid_isgs,
    score_name="IRC_score",
    use_raw=True,
    random_state=537
)

print("IRC score calculated successfully.")


# ## 4.6 Summarize IRC scores
# 
# Overall and sample-level descriptive statistics are saved so that the IRC-score distribution can be reviewed independently of the notebook.
# 

# In[6]:


irc_global_summary = (
    adata.obs["IRC_score"]
    .describe()
    .rename("value")
    .reset_index()
    .rename(columns={"index": "metric"})
)

irc_global_summary.to_csv(
    preprocessing_dir / "irc_score_summary.csv",
    index=False
)

irc_by_sample = (
    adata.obs
    .groupby("sample", observed=True)["IRC_score"]
    .agg(
        cells="size",
        mean="mean",
        median="median",
        std="std",
        minimum="min",
        maximum="max"
    )
    .reset_index()
)

irc_by_sample.to_csv(
    preprocessing_dir / "irc_score_by_sample.csv",
    index=False
)

print("Saved:", preprocessing_dir / "irc_score_summary.csv")
print("Saved:", preprocessing_dir / "irc_score_by_sample.csv")

display(irc_global_summary)
display(irc_by_sample)


# ## 4.7 Visualize IRC scores on UMAP
# 
# The existing UMAP coordinates are not recalculated. Cells are simply colored by IRC score, donor type, and Leiden cluster.
# 

# In[7]:


sc.pl.umap(
    adata,
    color=[
        "IRC_score",
        "donor_type",
        "leiden_res_1.0"
    ],
    cmap="viridis",
    ncols=1,
    show=False
)

plt.savefig(
    figures_dir / "15_irc_umap.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Saved:", figures_dir / "15_irc_umap.png")


# ## 4.8 Compare IRC-score distributions across samples
# 
# This violin plot shows the distribution of cell-level IRC scores within each sample.
# 

# In[8]:


sc.pl.violin(
    adata,
    keys="IRC_score",
    groupby="sample",
    rotation=45,
    show=False
)

plt.savefig(
    figures_dir / "16_irc_violin_by_sample.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Saved:", figures_dir / "16_irc_violin_by_sample.png")


# ## 4.9 Compare IRC scores by donor type
# 
# This visualization compares cell-level IRC-score distributions between healthy-donor and patient cells. It should not be interpreted as an independent patient-level statistical test.
# 

# In[9]:


sc.pl.violin(
    adata,
    keys="IRC_score",
    groupby="donor_type",
    show=False
)

plt.savefig(
    figures_dir / "17_irc_violin_by_donor_type.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Saved:", figures_dir / "17_irc_violin_by_donor_type.png")


# ## 4.10 Save the IRC-scored AnnData object
# 
# The cell-level IRC score is stored in:
# 
# `adata.obs["IRC_score"]`
# 

# In[10]:


irc_file = preprocessing_dir / "ifn_irc_scored.h5ad"

adata.write_h5ad(irc_file)

print("Saved:", irc_file)


# ## 4.11 Restrict the analysis to patient cells
# 
# Only cells with `donor_type == "Patient"` are retained for the IRC-High versus IRC-Low comparison.
# 

# In[11]:


patient_adata = adata[
    adata.obs["donor_type"] == "Patient"
].copy()

if patient_adata.n_obs == 0:
    raise ValueError(
        "No patient cells were found in adata.obs['donor_type']."
    )

print(f"Patient cells: {patient_adata.n_obs:,}")


# ## 4.12 Cell-level IRC-High and IRC-Low stratification
# 
# The median IRC score across **patient cells** is used as the cutoff.
# 
# - Cells with IRC score greater than or equal to the median are labeled `IRC_High`.
# - Cells below the median are labeled `IRC_Low`.
# 
# This is a **cell-level stratification**. It does not classify entire patients.
# 

# In[12]:


median_score = patient_adata.obs["IRC_score"].median()

patient_adata.obs["IRC_group"] = np.where(
    patient_adata.obs["IRC_score"] >= median_score,
    "IRC_High",
    "IRC_Low"
)

irc_group_counts = (
    patient_adata.obs["IRC_group"]
    .value_counts()
    .rename_axis("IRC_group")
    .reset_index(name="cells")
)

irc_group_counts["fraction_of_patient_cells"] = (
    irc_group_counts["cells"] / patient_adata.n_obs
)

irc_group_counts["median_cutoff"] = median_score

irc_group_counts.to_csv(
    preprocessing_dir / "patient_irc_group_counts.csv",
    index=False
)

print(f"Median patient-cell IRC score: {median_score:.4f}")
display(irc_group_counts)


# ## 4.13 Differential expression: IRC-High vs IRC-Low patient cells
# 
# The Wilcoxon rank-sum test is used to compare gene expression between the two cell-level IRC groups.
# 
# Because multiple cells come from the same patient, cells are not independent biological replicates. Therefore, these results should be interpreted as an **exploratory cell-level expression comparison**, not as formal patient-level differential expression.
# 
# A patient-level analysis would require patient/sample identity to be incorporated explicitly, for example through a pseudobulk or other replicate-aware framework.
# 

# In[13]:


sc.tl.rank_genes_groups(
    patient_adata,
    groupby="IRC_group",
    groups=["IRC_High"],
    reference="IRC_Low",
    method="wilcoxon",
    use_raw=True
)

print("Differential-expression analysis complete.")


# ## 4.14 Extract and save the differential-expression table
# 
# Positive scores and log-fold changes indicate genes enriched in IRC-High cells relative to IRC-Low cells.
# 

# In[14]:


de_df = sc.get.rank_genes_groups_df(
    patient_adata,
    group="IRC_High"
)

de_file = preprocessing_dir / "IRC_High_vs_Low_DE.csv"

de_df.to_csv(
    de_file,
    index=False
)

print("Saved:", de_file)
print("\nTop genes associated with IRC-High patient cells:")

display(de_df.head(15))


# ## 4.15 Save the patient IRC/DE AnnData object
# 
# This file contains:
# 
# - patient cells only
# - `IRC_score`
# - `IRC_group`
# - differential-expression results in `patient_adata.uns["rank_genes_groups"]`
# 

# In[15]:


patient_file = preprocessing_dir / "patient_irc_DE.h5ad"

patient_adata.write_h5ad(patient_file)

print("Saved:", patient_file)

print("\n==============================")
print("SECTION 4 FINAL DATASETS")
print("==============================")
print(f"All scored cells: {adata.n_obs:,}")
print(f"Patient cells: {patient_adata.n_obs:,}")
print(f"IRC genes used: {len(valid_isgs)}/{len(isg_signature)}")
print(f"Median patient-cell IRC cutoff: {median_score:.4f}")


# # Section 4 output summary
# 
# | Output | Meaning |
# |---|---|
# | `irc_signature_availability.csv` | Records which of the six IRC-signature genes are available in `adata.raw` |
# | `irc_score_summary.csv` | Overall descriptive statistics for cell-level IRC scores |
# | `irc_score_by_sample.csv` | Sample-level descriptive statistics for cell-level IRC scores |
# | `15_irc_umap.png` | UMAP colored by IRC score, donor type, and Leiden cluster |
# | `16_irc_violin_by_sample.png` | Distribution of cell-level IRC scores across samples |
# | `17_irc_violin_by_donor_type.png` | Distribution of cell-level IRC scores across donor types |
# | `ifn_irc_scored.h5ad` | Clustered AnnData object with `IRC_score` added to cell metadata |
# | `patient_irc_group_counts.csv` | Cell counts and proportions in the IRC-High and IRC-Low patient-cell groups |
# | `IRC_High_vs_Low_DE.csv` | Exploratory cell-level Wilcoxon differential-expression results for IRC-High vs IRC-Low patient cells |
# | `patient_irc_DE.h5ad` | Patient-only AnnData object containing IRC groups and DE results |
# 
# ## Interpretation
# 
# Section 3 established the major transcriptional structure of the dataset using PCA, UMAP, and Leiden clustering.
# 
# Section 4 adds a per-cell IRC score based on the six-gene interferon-response signature and visualizes its distribution across the PBMC dataset.
# 
# Patient cells are then divided into IRC-High and IRC-Low groups according to the median **patient-cell IRC score**.
# 
# The resulting differential-expression analysis therefore compares **cells**, not independent patients. These DE results should be interpreted as exploratory evidence of gene-expression patterns associated with high versus low IRC states.
# 
# The main outputs are:
# 
# `ifn_irc_scored.h5ad`
# 
# and
# 
# `patient_irc_DE.h5ad`
# 
# with the differential-expression table saved as:
# 
# `IRC_High_vs_Low_DE.csv`
# 
