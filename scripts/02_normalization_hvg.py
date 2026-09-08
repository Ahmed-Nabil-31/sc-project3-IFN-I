#!/usr/bin/env python
# coding: utf-8

# # 2. Normalization and Highly Variable Genes
# 
# The quality-controlled AnnData object generated in Section 1 is prepared for
# downstream single-cell RNA-seq analysis.
# 
# ## Workflow
# 
# Final QC-filtered data → preserve counts → normalization → log transformation
# → highly variable gene (HVG) selection → HVG subsetting → regression →
# scaling → final dataset for dimensionality reduction.
# 
# ## Input
# 
# `results/preprocessing_tables/merged_qc_filtered.h5ad`
# 
# ## Main outputs
# 
# - `results/preprocessing_tables/normalization_summary.csv`
# - `results/figures/08_highly_variable_genes.png`
# - `results/preprocessing_tables/hvg_summary.csv`
# - `results/preprocessing_tables/top_hvgs.csv`
# - `results/preprocessing_tables/hvg_subset.h5ad`
# - `results/preprocessing_tables/normalized_hvg_regressed_scaled.h5ad`
# 
# > **Reproducibility:** No Google Colab or Google Drive paths are used. The
# > project root is detected automatically, so the notebook can be run after
# > cloning the repository either from the repository root or from the
# > `notebooks/` directory.
# 

# ## 2.1 Import libraries and locate the project
# 
# The notebook searches the current directory and its parent directories for
# the project structure. This avoids hard-coded user-specific paths.
# 

# In[1]:


from pathlib import Path
from importlib.metadata import version

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
sc.set_figure_params(
    dpi=100,
    dpi_save=300,
    figsize=(6, 4)
)

def find_project_root():
    current = Path.cwd().resolve()

    candidates = [current, *current.parents]

    for path in candidates:
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

data_dir = project_root / "data"
results_dir = project_root / "results"
preprocessing_dir = results_dir / "preprocessing_tables"
figures_dir = results_dir / "figures"

preprocessing_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

input_file = preprocessing_dir / "merged_qc_filtered.h5ad"

print("Project root:", project_root)
print("Pandas:", pd.__version__)
print("Scanpy:", sc.__version__)
print("AnnData:", version("anndata"))
print("Input exists:", input_file.exists())


# ## 2.2 Load the final QC-filtered dataset
# 
# This is the output of Section 1. It has already undergone:
# 
# - basic cell QC
# - Scrublet-predicted doublet removal
# - gene filtering for genes detected in at least five cells
# 
# This dataset is the starting point for expression preprocessing.
# 

# In[2]:


if not input_file.exists():
    raise FileNotFoundError(
        f"Required input file was not found: {input_file}\n"
        "Run notebooks/01_data_loading_qc.ipynb first."
    )

adata = sc.read_h5ad(input_file)

print(adata)
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")


# ## 2.3 Preserve the raw count matrix
# 
# The original count matrix is copied to `adata.layers["counts"]` before any
# normalization or transformation.
# 
# This preserves the count data for later analyses that require the original
# expression counts.
# 

# In[3]:


adata.layers["counts"] = adata.X.copy()

print("Count layer created:", "counts" in adata.layers)


# ## 2.4 Normalize and log-transform expression
# 
# Total-count normalization scales each cell to a target library size of
# 10,000 counts. This reduces the effect of differences in sequencing depth
# between cells.
# 
# The normalized expression values are then log-transformed using `log1p`.
# 
# The complete normalized/log-transformed matrix is stored in `adata.raw`
# before the AnnData object is subset to HVGs.
# 

# In[4]:


# Record sample-level information before normalization
normalization_summary = (
    adata.obs
    .groupby("sample", observed=True)
    .agg(
        cells=("sample", "size"),
        median_total_counts=("total_counts", "median"),
        median_n_genes=("n_genes_by_counts", "median")
    )
)

# Normalize each cell to a common target library size
sc.pp.normalize_total(
    adata,
    target_sum=1e4
)

# Log-transform
sc.pp.log1p(adata)

# Preserve the complete normalized/log-transformed expression matrix
adata.raw = adata

print("Normalization and log transformation complete.")
print("Target sum per cell: 10,000")
print("Full normalized/log-transformed matrix stored in adata.raw.")


# In[5]:


normalization_summary["target_sum"] = 10_000

normalization_summary.to_csv(
    preprocessing_dir / "normalization_summary.csv",
    index=True
)

display(normalization_summary)
print(
    "Saved:",
    preprocessing_dir / "normalization_summary.csv"
)


# ## 2.5 Identify highly variable genes
# 
# Highly variable genes (HVGs) show substantial expression variation across
# cells and are informative for downstream dimensionality reduction and
# clustering.
# 
# The top 2,500 HVGs are selected. `sample` is supplied as the batch key so
# that sample-specific effects are considered during HVG selection.
# 

# In[6]:


sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2500,
    batch_key="sample"
)

n_hvgs = int(adata.var["highly_variable"].sum())

print("Highly variable genes:", n_hvgs)


# ## 2.6 HVG diagnostic plot
# 
# This plot provides a visual quality check of the highly variable gene
# selection.
# 
# It is a preprocessing diagnostic rather than a biological result figure.
# 

# In[7]:


sc.pl.highly_variable_genes(
    adata,
    show=False
)

hvg_figure = figures_dir / "08_highly_variable_genes.png"

plt.savefig(
    hvg_figure,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved:", hvg_figure)


# ## 2.7 Save HVG summary tables
# 
# The summary table records the number of genes available and the number
# selected as HVGs.
# 
# The selected-gene table provides the gene names and the HVG statistics
# calculated by Scanpy.
# 

# In[8]:


hvg_summary = pd.DataFrame({
    "metric": [
        "total_genes",
        "selected_hvgs",
        "requested_hvgs"
    ],
    "value": [
        adata.n_vars,
        int(adata.var["highly_variable"].sum()),
        2500
    ]
})

hvg_summary.to_csv(
    preprocessing_dir / "hvg_summary.csv",
    index=False
)

hvg_columns = ["highly_variable"]

for col in [
    "highly_variable_rank",
    "means",
    "dispersions",
    "dispersions_norm",
    "highly_variable_nbatches"
]:
    if col in adata.var.columns:
        hvg_columns.append(col)

top_hvgs = adata.var.loc[
    adata.var["highly_variable"],
    hvg_columns
].copy()

top_hvgs.insert(0, "gene", top_hvgs.index)

top_hvgs.to_csv(
    preprocessing_dir / "top_hvgs.csv",
    index=False
)

print("Saved:", preprocessing_dir / "hvg_summary.csv")
print("Saved:", preprocessing_dir / "top_hvgs.csv")

display(hvg_summary)
display(top_hvgs.head(20))


# ## 2.8 Subset the dataset to HVGs
# 
# The AnnData object is restricted to the selected 2,500 HVGs.
# 
# The number of cells does not change. The number of genes decreases to the
# selected HVG set.
# 
# The full normalized/log-transformed expression matrix remains available
# through `adata.raw`.
# 

# In[9]:


adata = adata[
    :,
    adata.var["highly_variable"]
].copy()

print("After HVG subsetting:")
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")
print("Full normalized matrix retained in adata.raw:", adata.raw is not None)


# In[10]:


hvg_subset_file = preprocessing_dir / "hvg_subset.h5ad"

adata.write_h5ad(hvg_subset_file)

print("Saved:", hvg_subset_file)


# ## 2.9 Regress out technical effects
# 
# Variation associated with total UMI counts and mitochondrial RNA percentage
# is regressed out before scaling and PCA.
# 
# This is intended to reduce technical/quality-associated variation that could
# otherwise influence downstream dimensionality reduction.
# 

# In[11]:


sc.pp.regress_out(
    adata,
    [
        "total_counts",
        "pct_counts_mt"
    ]
)

print("Regression complete.")


# ## 2.10 Scale the expression matrix
# 
# Genes are standardized so that they have comparable scales during PCA.
# 
# Values are clipped at ±10 to limit the influence of extreme values.
# 

# In[12]:


sc.pp.scale(
    adata,
    max_value=10
)

print("Scaling complete.")


# ## 2.11 Save the final Section 2 dataset
# 
# This file is the handoff point to the dimensionality-reduction section.
# 
# The next section can start from this dataset and perform:
# 
# PCA → batch correction → neighbors → UMAP → clustering.
# 

# In[13]:


final_file = (
    preprocessing_dir /
    "normalized_hvg_regressed_scaled.h5ad"
)

adata.write_h5ad(final_file)

print("Saved:", final_file)

print("\n==============================")
print("SECTION 2 FINAL DATASET")
print("==============================")
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")
print(f"HVGs represented in adata: {adata.n_vars:,}")


# # Section 2 output summary
# 
# | Output | Meaning |
# |---|---|
# | `normalization_summary.csv` | Sample-level record of cell counts and QC metrics entering normalization |
# | `08_highly_variable_genes.png` | Diagnostic visualization of HVG selection |
# | `hvg_summary.csv` | Number of genes available and number selected as HVGs |
# | `top_hvgs.csv` | Selected HVGs and their Scanpy variability statistics |
# | `hvg_subset.h5ad` | Normalized/log-transformed dataset restricted to selected HVGs |
# | `normalized_hvg_regressed_scaled.h5ad` | Final Section 2 dataset after HVG selection, regression, and scaling |
# 
# ## Interpretation
# 
# Section 1 established which cells and genes should be retained.
# 
# Section 2 transforms that trusted expression matrix into a form suitable for
# dimensionality reduction and clustering.
# 
# The final output is:
# 
# `normalized_hvg_regressed_scaled.h5ad`
# 
# which becomes the input to Section 3.
# 
