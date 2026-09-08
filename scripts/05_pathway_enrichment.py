#!/usr/bin/env python
# coding: utf-8

# # 5. Differential-Expression Visualization and Pathway Enrichment
# 
# This section uses the IRC-High versus IRC-Low differential-expression results generated in Section 4.
# 
# The workflow is:
# 
# **DE table → significance filtering → volcano plots → Enrichr pathway enrichment → top enriched pathways**
# 
# ## Input
# 
# `results/preprocessing_tables/IRC_High_vs_Low_DE.csv`
# 
# ## Main outputs
# 
# ### Figures
# - `results/figures/18_volcano_plot.png`
# - `results/figures/19_volcano_plot_labeled.png`
# - `results/figures/20_pathway_enrichment_barplot.png`
# 
# ### Tables
# - `results/preprocessing_tables/significant_DE_genes.csv`
# - `results/preprocessing_tables/enrichment_input_genes.csv`
# - `results/preprocessing_tables/pathway_enrichment_results.csv`
# - `results/preprocessing_tables/top_pathways.csv`
# 
# > **Internet requirement:** The Enrichr step uses GSEApy to query online Enrichr gene-set libraries. Therefore, an internet connection is required when that cell is executed.
# 
# > **Interpretation note:** The DE table comes from the exploratory **cell-level** IRC-High versus IRC-Low comparison in Section 4. Pathway-enrichment results inherit that interpretation and should not be treated as patient-level differential-expression evidence.
# 

# ## 5.1 Import libraries and locate the project
# 
# The notebook automatically detects the repository root and does not rely on Google Colab, Google Drive, or machine-specific absolute paths.
# 

# In[1]:


from pathlib import Path
from importlib.metadata import version

import gseapy as gp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text

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

de_file = preprocessing_dir / "IRC_High_vs_Low_DE.csv"

print("Project root:", project_root)
print("GSEApy:", version("gseapy"))
print("Pandas:", pd.__version__)
print("Input exists:", de_file.exists())


# ## 5.2 Load the differential-expression results
# 
# The input table was produced in Section 4 from the cell-level comparison of IRC-High versus IRC-Low patient cells.
# 

# In[2]:


if not de_file.exists():
    raise FileNotFoundError(
        f"Required input file was not found: {de_file}\n"
        "Run notebooks/04_IRCscore.ipynb first."
    )

de_df = pd.read_csv(de_file)

required_columns = [
    "names",
    "logfoldchanges",
    "pvals_adj"
]

missing_columns = [
    col for col in required_columns
    if col not in de_df.columns
]

if missing_columns:
    raise KeyError(
        "The DE table is missing required columns: "
        + ", ".join(missing_columns)
    )

print("DE results:", de_df.shape)
print(de_df.columns.tolist())
display(de_df.head())


# ## 5.3 Define significant differential-expression genes
# 
# Genes are considered significant when:
# 
# - adjusted P-value < 0.05
# - absolute log2 fold change > 0.5
# 
# The complete significant-gene table is saved.
# 
# For pathway enrichment, only genes upregulated in IRC-High cells are used:
# 
# - log2 fold change > 0.5
# - adjusted P-value < 0.05
# 
# The first 100 genes in the ranked DE table that meet these criteria are passed to Enrichr, matching the original analysis.
# 

# In[3]:


de_df["-log10_padj"] = -np.log10(
    de_df["pvals_adj"].clip(lower=1e-300)
)

de_df["significant"] = (
    (de_df["pvals_adj"] < 0.05)
    & (de_df["logfoldchanges"].abs() > 0.5)
)

significant_de = de_df[
    de_df["significant"]
].copy()

significant_de.to_csv(
    preprocessing_dir / "significant_DE_genes.csv",
    index=False
)

sig_genes = de_df[
    (de_df["logfoldchanges"] > 0.5)
    & (de_df["pvals_adj"] < 0.05)
]["names"].dropna().astype(str).tolist()[:100]

enrichment_input = pd.DataFrame({
    "gene": sig_genes
})

enrichment_input.to_csv(
    preprocessing_dir / "enrichment_input_genes.csv",
    index=False
)

print("Significant genes (both directions):", len(significant_de))
print("Upregulated genes used for enrichment:", len(sig_genes))
print("\nFirst 20 enrichment genes:")
print(sig_genes[:20])

print("\nSaved:", preprocessing_dir / "significant_DE_genes.csv")
print("Saved:", preprocessing_dir / "enrichment_input_genes.csv")


# ## 5.4 Volcano plot
# 
# The volcano plot shows effect size on the x-axis and statistical significance on the y-axis.
# 
# Threshold lines indicate:
# 
# - log2 fold change = ±0.5
# - adjusted P-value = 0.05
# 

# In[4]:


plt.figure(figsize=(10, 7))

plt.scatter(
    de_df["logfoldchanges"],
    de_df["-log10_padj"],
    alpha=0.4,
    s=12
)

sig = de_df[de_df["significant"]]

plt.scatter(
    sig["logfoldchanges"],
    sig["-log10_padj"],
    alpha=0.8,
    s=20
)

plt.axvline(
    0.5,
    linestyle="--",
    linewidth=1
)

plt.axvline(
    -0.5,
    linestyle="--",
    linewidth=1
)

plt.axhline(
    -np.log10(0.05),
    linestyle="--",
    linewidth=1
)

plt.xlabel("log2 fold change")
plt.ylabel("-log10 adjusted P-value")
plt.title("Differential Expression: IRC-High vs IRC-Low")

plt.tight_layout()

volcano_file = figures_dir / "18_volcano_plot.png"

plt.savefig(
    volcano_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved:", volcano_file)


# ## 5.5 Labeled volcano plot
# 
# The ten most statistically significant genes among those passing the significance thresholds are labeled.
# 
# `adjustText` is used to reduce overlap between gene labels.
# 

# In[5]:


top_genes = (
    de_df[
        de_df["significant"]
    ]
    .sort_values("pvals_adj")
    .head(10)
)

fig, ax = plt.subplots(figsize=(12, 9))

ax.scatter(
    de_df["logfoldchanges"],
    de_df["-log10_padj"],
    s=10,
    alpha=0.35
)

ax.scatter(
    sig["logfoldchanges"],
    sig["-log10_padj"],
    s=18,
    alpha=0.7
)

texts = []

for _, row in top_genes.iterrows():
    texts.append(
        ax.text(
            row["logfoldchanges"],
            row["-log10_padj"],
            row["names"],
            fontsize=9
        )
    )

if texts:
    adjust_text(
        texts,
        ax=ax,
        force_text=(1.5, 2.0),
        force_static=(1.0, 1.0),
        only_move={
            "text": "xy",
            "static": "xy",
            "explode": "xy",
            "pull": "xy"
        },
        arrowprops=dict(
            arrowstyle="-",
            linewidth=0.8
        )
    )

ax.axvline(
    0.5,
    linestyle="--",
    linewidth=1
)

ax.axvline(
    -0.5,
    linestyle="--",
    linewidth=1
)

ax.axhline(
    -np.log10(0.05),
    linestyle="--",
    linewidth=1
)

ax.set_xlabel("log2 Fold Change")
ax.set_ylabel("-log10 Adjusted P-value")
ax.set_title(
    "Differential Expression: IRC-High vs IRC-Low"
)

plt.tight_layout()

labeled_volcano_file = (
    figures_dir / "19_volcano_plot_labeled.png"
)

plt.savefig(
    labeled_volcano_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved:", labeled_volcano_file)


# ## 5.6 Pathway enrichment using Enrichr
# 
# The upregulated IRC-High genes are tested against three gene-set libraries:
# 
# - Reactome 2022
# - KEGG 2021 Human
# - GO Biological Process 2023
# 
# If no genes pass the significance thresholds, this step is skipped cleanly.
# 
# Because Enrichr is an online service, this step requires internet access.
# 

# In[6]:


gene_sets = [
    "Reactome_2022",
    "KEGG_2021_Human",
    "GO_Biological_Process_2023"
]

pathway_results_file = (
    preprocessing_dir /
    "pathway_enrichment_results.csv"
)

if len(sig_genes) > 0:
    try:
        enr = gp.enrichr(
            gene_list=sig_genes,
            gene_sets=gene_sets,
            organism="human",
            outdir=None
        )

        pathway_df = enr.results.copy()

        pathway_df.to_csv(
            pathway_results_file,
            index=False
        )

        print("Pathway enrichment completed.")
        print("Pathways returned:", len(pathway_df))
        print("Saved:", pathway_results_file)

    except Exception as exc:
        raise RuntimeError(
            "Enrichr pathway enrichment failed. "
            "Check the internet connection and availability of the requested Enrichr libraries."
        ) from exc

else:
    pathway_df = pd.DataFrame()
    print(
        "No genes passed the enrichment thresholds. "
        "Pathway enrichment was skipped."
    )


# ## 5.7 Select the top enriched pathways
# 
# The 15 pathways with the smallest adjusted P-values are retained for visualization and saved separately.
# 

# In[7]:


top_pathways_file = (
    preprocessing_dir /
    "top_pathways.csv"
)

if not pathway_df.empty:

    required_pathway_columns = [
        "Gene_set",
        "Term",
        "Adjusted P-value"
    ]

    missing_pathway_columns = [
        col for col in required_pathway_columns
        if col not in pathway_df.columns
    ]

    if missing_pathway_columns:
        raise KeyError(
            "Enrichr results are missing expected columns: "
            + ", ".join(missing_pathway_columns)
        )

    top_pathways = (
        pathway_df
        .sort_values("Adjusted P-value")
        .head(15)
        .copy()
    )

    top_pathways["-log10_padj"] = -np.log10(
        top_pathways[
            "Adjusted P-value"
        ].clip(lower=1e-300)
    )

    top_pathways.to_csv(
        top_pathways_file,
        index=False
    )

    display(
        top_pathways[
            [
                "Gene_set",
                "Term",
                "Adjusted P-value"
            ]
        ]
    )

    print("Saved:", top_pathways_file)

else:
    top_pathways = pd.DataFrame()
    print("No pathway results are available to rank.")


# ## 5.8 Pathway-enrichment bar plot
# 
# The top enriched pathways are visualized using `-log10(adjusted P-value)`.
# 
# Larger values indicate stronger statistical enrichment.
# 

# In[8]:


if not top_pathways.empty:

    plot_df = top_pathways.sort_values(
        "-log10_padj",
        ascending=True
    )

    plt.figure(figsize=(11, 8))

    plt.barh(
        plot_df["Term"],
        plot_df["-log10_padj"]
    )

    plt.xlabel("-log10 adjusted P-value")
    plt.ylabel("Enriched pathway")
    plt.title(
        "Pathway Enrichment of Genes Upregulated in IRC-High Cells"
    )

    plt.tight_layout()

    pathway_plot_file = (
        figures_dir /
        "20_pathway_enrichment_barplot.png"
    )

    plt.savefig(
        pathway_plot_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("Saved:", pathway_plot_file)

else:
    print(
        "Pathway bar plot was not generated because "
        "no enrichment results were available."
    )


# # Section 5 output summary
# 
# | Output | Meaning |
# |---|---|
# | `significant_DE_genes.csv` | All genes passing adjusted P-value < 0.05 and absolute log2 fold change > 0.5 |
# | `enrichment_input_genes.csv` | Up to 100 IRC-High-upregulated genes submitted to Enrichr |
# | `18_volcano_plot.png` | Volcano plot showing all DE genes and the significance thresholds |
# | `19_volcano_plot_labeled.png` | Volcano plot with the ten most significant threshold-passing genes labeled |
# | `pathway_enrichment_results.csv` | Complete Enrichr results from Reactome, KEGG, and GO Biological Process |
# | `top_pathways.csv` | Fifteen pathways with the smallest adjusted P-values |
# | `20_pathway_enrichment_barplot.png` | Bar plot of the top enriched pathways |
# 
# ## Interpretation
# 
# Section 4 identified genes associated with the cell-level IRC-High versus IRC-Low comparison.
# 
# Section 5 first visualizes those differential-expression results using volcano plots. Genes upregulated in IRC-High cells with:
# 
# - log2 fold change > 0.5
# - adjusted P-value < 0.05
# 
# are then used as the input for pathway enrichment.
# 
# Enrichr is queried using Reactome, KEGG, and GO Biological Process gene-set libraries. The resulting pathways provide a functional summary of biological processes represented among genes associated with the IRC-High cell state.
# 
# Because the source differential-expression analysis is a **cell-level exploratory comparison**, the enrichment analysis should also be interpreted as exploratory rather than as independent patient-level evidence.
# 
# The main pathway output is:
# 
# `pathway_enrichment_results.csv`
# 
