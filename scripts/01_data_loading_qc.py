#!/usr/bin/env python
# coding: utf-8

# # Type I Interferon Response Capacity in PBMCs
# 
# Single-cell RNA-seq analysis of pre-therapy PBMC samples to investigate
# Type I interferon response capacity and its association with downstream
# transcriptional and pathway-level changes.
# 
# ## Analysis workflow
# 
# Raw 10X data → Quality control → Doublet detection → Normalization →
# Highly variable genes → Dimensionality reduction → Clustering →
# IRC score → Differential expression → Pathway enrichment

# # 1. Data Loading and Quality Control
# 
# Raw 10X Genomics count matrices from ten samples were loaded individually
# and combined into a single AnnData object.
# 
# The dataset contains two healthy donor samples (HD1 and HD2) and eight
# patient samples (P1–P8).
# 
# Quality-control metrics were calculated for:
# - number of detected genes
# - total UMI counts
# - mitochondrial transcript percentage
# - hemoglobin transcript percentage
# 
# Cells with fewer than 300 detected genes, mitochondrial content ≥15%,
# or hemoglobin content ≥5% were removed.
# 
# Because unusually high gene and UMI counts can indicate potential
# doublets, Scrublet was subsequently applied independently to each sample.
# Cells predicted as doublets were removed before downstream analysis.

# In[45]:


# ============================================================
# 1. DATA LOADING, QUALITY CONTROL & DOUBLET DETECTION
# ============================================================

import scanpy as sc
import anndata as ad
import scrublet as scr

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from importlib.metadata import version


# ------------------------------------------------------------
# Reproducibility and plotting settings
# ------------------------------------------------------------

sc.settings.verbosity = 1

sc.set_figure_params(
    dpi=100,
    dpi_save=300,
    figsize=(6, 4)
)


# ------------------------------------------------------------
# Define directories
# ------------------------------------------------------------

data_dir = Path("../data")

results_dir = Path("../results")
figures_dir = results_dir / "figures"
preprocessing_dir = results_dir / "preprocessing_tables"

figures_dir.mkdir(parents=True, exist_ok=True)
preprocessing_dir.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Software versions
# ------------------------------------------------------------

print("Scanpy version:", sc.__version__)
print("AnnData version:", version("anndata"))
print("Scrublet version:", version("scrublet"))


# ## 1.1 Load raw 10X data
# 
# Each sample is loaded separately so that sample-level metadata can be
# retained and Scrublet can later be applied independently to each sample.

# In[47]:


# ------------------------------------------------------------
# Identify sample directories
# ------------------------------------------------------------

sample_dirs = sorted([
    path for path in data_dir.iterdir()
    if path.is_dir()
])

print("Samples found:")

for path in sample_dirs:
    print(" -", path.name)


# In[48]:


# ------------------------------------------------------------
# Load each 10X dataset
# ------------------------------------------------------------

adatas = []

for sample_path in sample_dirs:

    sample_id = sample_path.name

    print(f"\nLoading {sample_id}...")

    sample_adata = sc.read_10x_mtx(
        sample_path,
        var_names="gene_symbols",
        cache=False
    )

    sample_adata.var_names_make_unique()

    # Store sample information
    sample_adata.obs["sample"] = sample_id

    # Classify donor
    if sample_id.startswith("HD"):
        sample_adata.obs["donor_type"] = "Healthy donor"
    else:
        sample_adata.obs["donor_type"] = "Patient"

    adatas.append(sample_adata)

    print(
        f"{sample_id}: "
        f"{sample_adata.n_obs} cells × "
        f"{sample_adata.n_vars} genes"
    )


# In[49]:


# ------------------------------------------------------------
# Merge all samples
# ------------------------------------------------------------

adata = ad.concat(
    adatas,
    join="outer",
    label="batch",
    keys=[x.obs["sample"].iloc[0] for x in adatas],
    index_unique="-"
)

adata.var_names_make_unique()

print("\nMerged dataset:")
print(adata)


# In[50]:


n_cells_raw = adata.n_obs
n_genes_raw = adata.n_vars


# In[51]:


# ------------------------------------------------------------
# Save merged raw dataset
# ------------------------------------------------------------

adata.write(
    preprocessing_dir / "merged_raw.h5ad"
)

print(
    "Saved:",
    preprocessing_dir / "merged_raw.h5ad"
)


# ## 1.2 Calculate quality-control metrics
# 
# Mitochondrial, ribosomal, and hemoglobin genes were annotated to allow
# calculation of their relative contribution to each cell's transcript
# counts.

# In[52]:


# ------------------------------------------------------------
# Define QC gene categories
# ------------------------------------------------------------

adata.var["mt"] = (
    adata.var_names.str.upper().str.startswith("MT-")
)

adata.var["ribo"] = (
    adata.var_names.str.upper().str.startswith("RPS")
    | adata.var_names.str.upper().str.startswith("RPL")
)

adata.var["hb"] = (
    adata.var_names.str.upper().str.startswith("HBA")
    | adata.var_names.str.upper().str.startswith("HBB")
    | adata.var_names.str.upper().str.startswith("HBD")
    | adata.var_names.str.upper().str.startswith("HBE")
    | adata.var_names.str.upper().str.startswith("HBG")
    | adata.var_names.str.upper().str.startswith("HBM")
    | adata.var_names.str.upper().str.startswith("HBQ")
    | adata.var_names.str.upper().str.startswith("HBZ")
)


# ------------------------------------------------------------
# Calculate QC metrics
# ------------------------------------------------------------

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    inplace=True
)


# In[53]:


adata


# In[54]:


adata.obs[
    [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb"
    ]
].head()


# ## 1.3 QC before filtering
# 
# QC distributions were examined across samples before filtering.
# This provides a baseline for assessing the effect of the subsequent
# quality-control and doublet-removal steps.

# In[55]:


# ------------------------------------------------------------
# QC violin plots before filtering by sample
# ------------------------------------------------------------

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"],
    groupby="sample",
    jitter=0.4,
    multi_panel=True,
    show=False
)

plt.savefig(
    figures_dir / "01_qc_before_filtering_by_sample.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()


# In[56]:


# ------------------------------------------------------------
# QC violin plots before filtering by donor type
# ------------------------------------------------------------

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"],
    groupby="donor_type",
    jitter=0.4,
    multi_panel=True,
    show=False
)

plt.savefig(
    figures_dir / "02_qc_before_filtering_by_donor_type.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()


# In[57]:


# ------------------------------------------------------------
# QC summary by sample
# ------------------------------------------------------------

qc_summary = (
    adata.obs
    .groupby("sample")
    [
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_hb"
        ]
    ]
    .median()
)

qc_summary.to_csv(
    preprocessing_dir /
    "qc_summary_before_filtering_by_sample.csv"
)

display(qc_summary)


# In[58]:


# ------------------------------------------------------------
# QC summary by donor type
# ------------------------------------------------------------

qc_summary = (
    adata.obs
    .groupby("donor_type")
    [
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_hb"
        ]
    ]
    .median()
)

qc_summary.to_csv(
    preprocessing_dir /
    "qc_summary_before_filtering_by_donor_type.csv"
)

display(qc_summary)


# ## 1.4 Basic QC filtering
# 
# Cells with fewer than 300 detected genes were removed as low-complexity
# cells. Cells with mitochondrial transcript percentages ≥15% or
# hemoglobin transcript percentages ≥5% were also excluded.
# 
# No arbitrary upper threshold was imposed on total UMI counts or number
# of detected genes because unusually high values were investigated using
# the subsequent doublet-detection step.

# In[59]:


# ------------------------------------------------------------
# Record initial cell number
# ------------------------------------------------------------

n_cells_before_basic_qc = adata.n_obs

print(
    f"Cells before basic QC: "
    f"{n_cells_before_basic_qc:,}"
)


# ------------------------------------------------------------
# Minimum detected genes
# ------------------------------------------------------------

sc.pp.filter_cells(
    adata,
    min_genes=300
)

print(
    f"After minimum gene filter: "
    f"{adata.n_obs:,} cells"
)


# ------------------------------------------------------------
# Mitochondrial filtering
# ------------------------------------------------------------

adata = adata[
    adata.obs["pct_counts_mt"] < 15,
    :
].copy()

print(
    f"After mitochondrial filter: "
    f"{adata.n_obs:,} cells"
)


# ------------------------------------------------------------
# Hemoglobin filtering
# ------------------------------------------------------------

adata = adata[
    adata.obs["pct_counts_hb"] < 5,
    :
].copy()

print(
    f"After hemoglobin filter: "
    f"{adata.n_obs:,} cells"
)


# In[60]:


n_cells_after_basic_qc = adata.n_obs
n_genes_after_basic_qc = adata.n_vars


# ## 1.5 Doublet detection using Scrublet
# 
# Potential doublets were detected using Scrublet. Doublet detection was
# performed separately for each sample because doublet rates and
# transcriptional distributions can differ between libraries.
# 
# Scrublet assigns each cell a doublet score and predicts whether the cell
# is a doublet based on simulated doublets generated from the observed
# expression data.

# In[61]:


# ------------------------------------------------------------
# Initialize Scrublet result columns
# ------------------------------------------------------------

adata.obs["doublet_score"] = 0.0
adata.obs["predicted_doublet"] = False


# ------------------------------------------------------------
# Run Scrublet independently for each sample
# ------------------------------------------------------------

for sample_id in adata.obs["sample"].unique():

    print(f"\nRunning Scrublet for {sample_id}...")

    sample_mask = (
        adata.obs["sample"] == sample_id
    )

    counts_matrix = adata[sample_mask].X

    scrub = scr.Scrublet(
        counts_matrix
    )

    doublet_scores, predicted_doublets = (
        scrub.scrub_doublets(
            verbose=False
        )
    )

    adata.obs.loc[
        sample_mask,
        "doublet_score"
    ] = doublet_scores

    adata.obs.loc[
        sample_mask,
        "predicted_doublet"
    ] = predicted_doublets


# In[62]:


# ------------------------------------------------------------
# Human-readable labels
# ------------------------------------------------------------

adata.obs["doublet_info"] = (
    adata.obs["predicted_doublet"]
    .map({
        False: "Singlet",
        True: "Doublet"
    })
)


# ------------------------------------------------------------
# Overall doublet counts
# ------------------------------------------------------------

doublet_counts = (
    adata.obs["predicted_doublet"]
    .value_counts()
)

print("Scrublet results:")
print(doublet_counts)


# ------------------------------------------------------------
# Overall doublet rate
# ------------------------------------------------------------

doublet_rate = (
    adata.obs["predicted_doublet"].mean()
    * 100
)

print(
    f"\nOverall predicted doublet rate: "
    f"{doublet_rate:.2f}%"
)


# ## 1.6 Inspect predicted doublets
# 
# The relationship between Scrublet predictions and common QC metrics was
# examined. Doublets are expected to show, in many cases, elevated numbers
# of detected genes and total UMI counts.

# In[63]:


sc.pl.violin(
    adata,
    "n_genes_by_counts",
    groupby="doublet_info",
    jitter=0.4,
    show=False
)

plt.savefig(
    figures_dir / "03_doublet_genes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# In[64]:


sc.pl.violin(
    adata,
    "total_counts",
    groupby="doublet_info",
    jitter=0.4,
    show=False
)

plt.savefig(
    figures_dir / "04_doublet_total_counts.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# In[65]:


sc.pl.violin(
    adata,
    "doublet_score",
    groupby="sample",
    jitter=0.4,
    rotation=45,
    show=False
)

plt.savefig(
    figures_dir / "05_scrublet_scores_by_sample.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# In[66]:


# ------------------------------------------------------------
# Per-sample doublet summary
# ------------------------------------------------------------

doublet_summary = (
    adata.obs
    .groupby("sample")
    .agg(
        cells=("sample", "size"),
        predicted_doublets=(
            "predicted_doublet",
            "sum"
        ),
        median_doublet_score=(
            "doublet_score",
            "median"
        )
    )
)

doublet_summary["doublet_rate_percent"] = (
    doublet_summary["predicted_doublets"]
    / doublet_summary["cells"]
    * 100
)

doublet_summary.to_csv(
    preprocessing_dir /
    "doublet_summary_by_sample.csv"
)

display(doublet_summary)


# ## 1.7 Remove Scrublet-predicted doublets
# 
# Cells classified as doublets by Scrublet were excluded from the
# downstream analysis. These cells are referred to as
# "Scrublet-predicted doublets" because doublet status is computationally
# inferred rather than experimentally confirmed.

# In[67]:


# ------------------------------------------------------------
# Record numbers before removal
# ------------------------------------------------------------

n_cells_before_doublet_removal = adata.n_obs

n_doublets = (
    adata.obs["predicted_doublet"].sum()
)

print(
    f"Cells before doublet removal: "
    f"{n_cells_before_doublet_removal:,}"
)

print(
    f"Scrublet-predicted doublets: "
    f"{n_doublets:,}"
)


# ------------------------------------------------------------
# Remove predicted doublets
# ------------------------------------------------------------

adata = adata[
    adata.obs["predicted_doublet"] == False,
    :
].copy()

print(
    f"Cells after doublet removal: "
    f"{adata.n_obs:,}"
)


# In[68]:


n_cells_after_scrublet = adata.n_obs
n_genes_after_scrublet = adata.n_vars


# In[28]:


# ------------------------------------------------------------
# Remove genes detected in fewer than 5 cells
# ------------------------------------------------------------

n_genes_before_filtering = adata.n_vars

sc.pp.filter_genes(
    adata,
    min_cells=5
)

n_genes_after_filtering = adata.n_vars

print(
    f"Genes before filtering: "
    f"{n_genes_before_filtering:,}"
)

print(
    f"Genes after filtering: "
    f"{n_genes_after_filtering:,}"
)


# ## 1.8 QC after filtering
# 
# QC metrics were visualized again after basic QC filtering, Scrublet
# doublet removal, and gene filtering.
# 
# The post-filtering distributions provide a direct assessment of whether
# the applied QC procedure successfully removed low-quality cells and
# extreme observations.

# In[69]:


# ------------------------------------------------------------
# Recalculate QC metrics
# ------------------------------------------------------------

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    inplace=True
)


# In[70]:


adata.obs[
    [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb"
    ]
].head()


# In[71]:


# ------------------------------------------------------------
# QC violin plots after filtering by sample
# ------------------------------------------------------------
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"],
    groupby="sample",
    jitter=0.4,
    multi_panel=True,
    show=False
)

plt.savefig(
    figures_dir / "06_qc_after_filtering_by_sample.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()


# In[72]:


# ------------------------------------------------------------
# QC violin plots after filtering by donor type
# ------------------------------------------------------------
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"],
    groupby="donor_type",
    jitter=0.4,
    multi_panel=True,
    show=False
)

plt.savefig(
    figures_dir / "07_qc_after_filtering_by_donor_type.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()


# In[73]:


# ------------------------------------------------------------
# Final cell counts by sample
# ------------------------------------------------------------

final_sample_counts = (
    adata.obs["sample"]
    .value_counts()
    .sort_index()
)

final_sample_counts.to_csv(
    preprocessing_dir /
    "final_cell_counts_by_sample.csv"
)

print("Final cell counts by sample:")
print(final_sample_counts)


# In[74]:


# ------------------------------------------------------------
# Final QC summary by sample
# ------------------------------------------------------------

final_qc_summary = (
    adata.obs
    .groupby("sample")
    [
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_hb"
        ]
    ]
    .median()
)

final_qc_summary.to_csv(
    preprocessing_dir /
    "qc_summary_after_filtering_by_sample.csv"
)

display(final_qc_summary)


# In[75]:


# ------------------------------------------------------------
# Final QC summary by donor type
# ------------------------------------------------------------

final_qc_summary = (
    adata.obs
    .groupby("donor_type")
    [
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_hb"
        ]
    ]
    .median()
)

final_qc_summary.to_csv(
    preprocessing_dir /
    "qc_summary_after_filtering_by_donor_type.csv"
)

display(final_qc_summary)


# In[76]:


# ------------------------------------------------------------
# Filtering summary: cells and genes
# ------------------------------------------------------------

filtering_summary = pd.DataFrame({
    "stage": [
        "Raw merged data",
        "After basic QC",
        "After Scrublet doublet removal",
        "After gene filtering"
    ],
    "cells": [
        n_cells_raw,
        n_cells_after_basic_qc,
        n_cells_after_scrublet,
        adata.n_obs
    ],
    "genes": [
        n_genes_raw,
        n_genes_after_basic_qc,
        n_genes_after_scrublet,
        n_genes_after_filtering
    ]
})

filtering_summary.to_csv(
    preprocessing_dir /
    "filtering_summary.csv",
    index=False
)

display(filtering_summary)


# ## 1.9 Save final QC-filtered dataset
# 
# The final QC-filtered AnnData object was saved and will serve as the
# input for all subsequent normalization, dimensionality-reduction,
# IRC-score, differential-expression, and pathway-enrichment analyses.

# In[77]:


# ------------------------------------------------------------
# Save final QC-filtered AnnData
# ------------------------------------------------------------

adata.write(
    preprocessing_dir /
    "merged_qc_filtered.h5ad"
)

print(
    "Final dataset saved to:"
)

print(
    preprocessing_dir /
    "merged_qc_filtered.h5ad"
)


# ------------------------------------------------------------
# Final dataset dimensions
# ------------------------------------------------------------

print("\n===================================")
print("FINAL DATASET")
print("===================================")

print(
    f"Cells: {adata.n_obs:,}"
)

print(
    f"Genes: {adata.n_vars:,}"
)

print(
    f"Samples: {adata.obs['sample'].nunique()}"
)

print("\nWorkflow complete.")


# In[78]:


adata

