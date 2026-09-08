#!/usr/bin/env python
# coding: utf-8

# # 3. Dimensionality Reduction, Clustering, and Immune-Lineage Marker Inspection
# 
# This section takes the processed HVG dataset generated in Section 2 and performs:
# 
# **PCA → neighborhood graph → UMAP → Leiden clustering → canonical immune-marker inspection**
# 
# ## Input
# 
# `results/preprocessing_tables/normalized_hvg_regressed_scaled.h5ad`
# 
# ## Main outputs
# 
# ### Figures
# 
# - `results/figures/09_pca_variance_ratio.png`
# - `results/figures/10_umap_clusters_donor_sample.png`
# - `results/figures/11_marker_dotplot_hvg.png`
# - `results/figures/12_marker_umap_hvg.png`
# - `results/figures/13_marker_dotplot_raw.png`
# - `results/figures/14_marker_umap_raw.png`
# 
# ### Tables
# 
# - `results/preprocessing_tables/pca_variance_ratio.csv`
# - `results/preprocessing_tables/leiden_cluster_counts.csv`
# - `results/preprocessing_tables/cluster_composition_by_sample.csv`
# - `results/preprocessing_tables/cluster_composition_by_donor_type.csv`
# - `results/preprocessing_tables/marker_availability_hvg.csv`
# - `results/preprocessing_tables/marker_availability_raw.csv`
# 
# ### Processed AnnData
# 
# - `results/preprocessing_tables/clustering_res1.0.h5ad`
# 
# > **Reproducibility:** No Google Colab, Google Drive, or user-specific absolute paths are used. The notebook detects the repository root automatically.
# 

# ## 3.1 Import libraries and locate the project
# 
# The notebook searches the current working directory and its parent directories for the repository structure. This makes the notebook portable after cloning.
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

input_file = preprocessing_dir / "normalized_hvg_regressed_scaled.h5ad"

print("Project root:", project_root)
print("Scanpy:", sc.__version__)
print("AnnData:", version("anndata"))
print("Pandas:", pd.__version__)
print("Input exists:", input_file.exists())


# ## 3.2 Load the Section 2 dataset
# 
# This input contains the 2,500 selected HVGs after normalization, log transformation, regression, and scaling.
# 
# The full normalized/log-transformed gene set should still be accessible through `adata.raw`, allowing canonical markers outside the HVG set to be examined.
# 

# In[2]:


if not input_file.exists():
    raise FileNotFoundError(
        f"Required input file was not found: {input_file}\n"
        "Run notebooks/02_normalization_hvg.ipynb first."
    )

adata = sc.read_h5ad(input_file)

print(adata)
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")
print("adata.raw available:", adata.raw is not None)


# ## 3.3 Principal Component Analysis
# 
# PCA reduces the HVG expression space into a smaller number of orthogonal dimensions while retaining major sources of variation.
# 
# The first 50 principal components are calculated using a fixed random seed for reproducibility.
# 

# In[3]:


sc.tl.pca(
    adata,
    n_comps=50,
    svd_solver="arpack",
    random_state=537
)

print("PCA shape:", adata.obsm["X_pca"].shape)


# ## 3.4 PCA variance explained
# 
# The PCA variance-ratio plot shows how much variation is explained by each principal component.
# 
# A CSV table is also saved so the variance contribution of every PC is recorded explicitly.
# 

# In[4]:


variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])

pca_variance_table = pd.DataFrame({
    "PC": np.arange(1, len(variance_ratio) + 1),
    "variance_ratio": variance_ratio,
    "cumulative_variance_ratio": np.cumsum(variance_ratio)
})

pca_variance_table.to_csv(
    preprocessing_dir / "pca_variance_ratio.csv",
    index=False
)

sc.pl.pca_variance_ratio(
    adata,
    n_pcs=50,
    log=True,
    show=False
)

plt.savefig(
    figures_dir / "09_pca_variance_ratio.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Saved:", preprocessing_dir / "pca_variance_ratio.csv")
print("Saved:", figures_dir / "09_pca_variance_ratio.png")

display(pca_variance_table.head(10))


# ## 3.5 Construct the neighborhood graph
# 
# Cell-to-cell similarity is calculated using the first 30 principal components.
# 
# Each cell is connected to its 20 nearest neighbors.
# 

# In[5]:


sc.pp.neighbors(
    adata,
    n_neighbors=20,
    n_pcs=30,
    random_state=159
)

print("Neighborhood graph constructed.")


# ## 3.6 Compute UMAP
# 
# UMAP provides a two-dimensional visualization of the neighborhood structure.
# 
# A fixed random seed is used so that the embedding is reproducible.
# 

# In[6]:


sc.tl.umap(
    adata,
    random_state=331
)

print("UMAP shape:", adata.obsm["X_umap"].shape)


# ## 3.7 Leiden clustering
# 
# Leiden clustering is performed on the neighborhood graph using resolution 1.0.
# 
# Cluster sizes and sample/donor composition are saved to tables for transparent reporting.
# 

# In[7]:


sc.tl.leiden(
    adata,
    resolution=1.0,
    key_added="leiden_res_1.0",
    random_state=537
)

n_clusters = adata.obs["leiden_res_1.0"].nunique()

print(f"Identified {n_clusters} clusters.")


# In[8]:


cluster_counts = (
    adata.obs["leiden_res_1.0"]
    .value_counts()
    .sort_index()
    .rename_axis("leiden_res_1.0")
    .reset_index(name="cells")
)

cluster_counts["fraction_of_cells"] = (
    cluster_counts["cells"] / adata.n_obs
)

cluster_counts.to_csv(
    preprocessing_dir / "leiden_cluster_counts.csv",
    index=False
)

cluster_by_sample = pd.crosstab(
    adata.obs["leiden_res_1.0"],
    adata.obs["sample"]
)

cluster_by_sample.to_csv(
    preprocessing_dir / "cluster_composition_by_sample.csv"
)

cluster_by_donor = pd.crosstab(
    adata.obs["leiden_res_1.0"],
    adata.obs["donor_type"]
)

cluster_by_donor.to_csv(
    preprocessing_dir / "cluster_composition_by_donor_type.csv"
)

print("Saved:", preprocessing_dir / "leiden_cluster_counts.csv")
print("Saved:", preprocessing_dir / "cluster_composition_by_sample.csv")
print("Saved:", preprocessing_dir / "cluster_composition_by_donor_type.csv")

display(cluster_counts)


# ## 3.8 Visualize clusters, donor type, and sample
# 
# These UMAPs help assess whether the major cell structures are shared across samples or dominated by individual samples/donor groups.
# 

# In[9]:


sc.pl.umap(
    adata,
    color=[
        "leiden_res_1.0",
        "donor_type",
        "sample"
    ],
    ncols=1,
    show=False
)

plt.savefig(
    figures_dir / "10_umap_clusters_donor_sample.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Saved:", figures_dir / "10_umap_clusters_donor_sample.png")


# ## 3.9 Define canonical immune-lineage markers
# 
# Canonical markers for major PBMC lineages are defined for preliminary biological interpretation of the Leiden clusters.
# 

# In[10]:


lineage_dict = {
    "CD4_T": ["CD3D", "CD4", "IL7R", "CCR7"],
    "CD8_T": ["CD3D", "CD8A", "CD8B", "GZMB"],
    "NK": ["NKG7", "GNLY", "NCAM1"],
    "B_cells": ["MS4A1", "CD19", "CD79A"],
    "Monocytes": ["CD14", "FCGR3A", "LYZ"],
    "DC": ["FCER1A", "CST3"]
}


# ## 3.10 Check marker availability in the HVG matrix
# 
# Some canonical markers may not be among the 2,500 selected HVGs.
# 
# Their availability is recorded explicitly before plotting.
# 

# In[11]:


marker_rows_hvg = []

for lineage, markers in lineage_dict.items():
    for gene in markers:
        marker_rows_hvg.append({
            "lineage": lineage,
            "gene": gene,
            "available_in_hvg_matrix": gene in adata.var_names
        })

marker_availability_hvg = pd.DataFrame(marker_rows_hvg)

marker_availability_hvg.to_csv(
    preprocessing_dir / "marker_availability_hvg.csv",
    index=False
)

available_markers = {
    lineage: [
        gene for gene in markers
        if gene in adata.var_names
    ]
    for lineage, markers in lineage_dict.items()
}

for lineage, markers in available_markers.items():
    print(f"{lineage}: {markers}")

print("Saved:", preprocessing_dir / "marker_availability_hvg.csv")


# ## 3.11 Marker dotplot using the HVG matrix
# 
# This dotplot displays canonical lineage markers that are present among the selected HVGs.
# 
# Dot size represents the fraction of cells expressing each gene, while color represents average expression.
# 

# In[12]:


available_markers_nonempty = {
    lineage: markers
    for lineage, markers in available_markers.items()
    if len(markers) > 0
}

if available_markers_nonempty:
    sc.pl.dotplot(
        adata,
        var_names=available_markers_nonempty,
        groupby="leiden_res_1.0",
        use_raw=False,
        show=False
    )

    plt.savefig(
        figures_dir / "11_marker_dotplot_hvg.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    print("Saved:", figures_dir / "11_marker_dotplot_hvg.png")
else:
    print("No canonical markers were available in the HVG matrix.")


# ## 3.12 Marker expression on UMAP using the HVG matrix
# 
# A compact set of common immune markers is visualized if those genes are present among the selected HVGs.
# 

# In[13]:


umap_markers_hvg = [
    gene for gene in [
        "CD3D",
        "CD4",
        "CD8A",
        "MS4A1",
        "CD14",
        "NKG7"
    ]
    if gene in adata.var_names
]

if umap_markers_hvg:
    sc.pl.umap(
        adata,
        color=umap_markers_hvg,
        use_raw=False,
        ncols=3,
        show=False
    )

    plt.savefig(
        figures_dir / "12_marker_umap_hvg.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    print("Saved:", figures_dir / "12_marker_umap_hvg.png")
else:
    print("None of the selected UMAP markers were available in the HVG matrix.")


# ## 3.13 Check marker availability in the full normalized expression matrix
# 
# `adata.raw` contains the complete normalized/log-transformed expression matrix saved before HVG subsetting in Section 2.
# 
# This allows canonical lineage markers that were not selected as HVGs to still be examined.
# 

# In[14]:


if adata.raw is None:
    raise ValueError(
        "adata.raw is missing. Section 2 should preserve the full normalized/log-transformed "
        "matrix in adata.raw before HVG subsetting."
    )

marker_rows_raw = []

for lineage, markers in lineage_dict.items():
    for gene in markers:
        marker_rows_raw.append({
            "lineage": lineage,
            "gene": gene,
            "available_in_raw": gene in adata.raw.var_names
        })

marker_availability_raw = pd.DataFrame(marker_rows_raw)

marker_availability_raw.to_csv(
    preprocessing_dir / "marker_availability_raw.csv",
    index=False
)

available_markers_raw = {
    lineage: [
        gene for gene in markers
        if gene in adata.raw.var_names
    ]
    for lineage, markers in lineage_dict.items()
}

for lineage, markers in available_markers_raw.items():
    print(f"{lineage}: {markers}")

print("Saved:", preprocessing_dir / "marker_availability_raw.csv")


# ## 3.14 Marker dotplot using the full normalized matrix
# 
# Using `use_raw=True` allows canonical markers outside the selected HVG set to be visualized across Leiden clusters.
# 

# In[15]:


available_markers_raw_nonempty = {
    lineage: markers
    for lineage, markers in available_markers_raw.items()
    if len(markers) > 0
}

if available_markers_raw_nonempty:
    sc.pl.dotplot(
        adata,
        var_names=available_markers_raw_nonempty,
        groupby="leiden_res_1.0",
        use_raw=True,
        show=False
    )

    plt.savefig(
        figures_dir / "13_marker_dotplot_raw.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    print("Saved:", figures_dir / "13_marker_dotplot_raw.png")
else:
    print("No canonical markers were available in adata.raw.")


# ## 3.15 Marker expression on UMAP using the full normalized matrix
# 
# Canonical immune markers are visualized using `adata.raw` so genes outside the 2,500 selected HVGs can still be inspected.
# 

# In[16]:


umap_markers_raw = [
    gene for gene in [
        "CD3D",
        "CD4",
        "CD8A",
        "MS4A1",
        "CD14",
        "NKG7"
    ]
    if gene in adata.raw.var_names
]

if umap_markers_raw:
    sc.pl.umap(
        adata,
        color=umap_markers_raw,
        use_raw=True,
        ncols=3,
        show=False
    )

    plt.savefig(
        figures_dir / "14_marker_umap_raw.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    print("Saved:", figures_dir / "14_marker_umap_raw.png")
else:
    print("None of the selected UMAP markers were available in adata.raw.")


# ## 3.16 Save the final clustered AnnData object
# 
# The final AnnData object contains:
# 
# - PCA coordinates
# - neighborhood graph
# - UMAP coordinates
# - Leiden cluster assignments
# - original sample metadata
# - `adata.raw` containing the full normalized/log-transformed expression matrix
# 
# This file becomes the handoff point for downstream annotation and IRC-related analyses.
# 

# In[17]:


final_file = preprocessing_dir / "clustering_res1.0.h5ad"

adata.write_h5ad(final_file)

print("Saved:", final_file)

print("\n==============================")
print("SECTION 3 FINAL DATASET")
print("==============================")
print(f"Cells: {adata.n_obs:,}")
print(f"HVG genes: {adata.n_vars:,}")
print(f"Leiden clusters: {adata.obs['leiden_res_1.0'].nunique()}")
print(f"PCA dimensions: {adata.obsm['X_pca'].shape[1]}")
print("UMAP available:", "X_umap" in adata.obsm)


# # Section 3 output summary
# 
# | Output | Meaning |
# |---|---|
# | `09_pca_variance_ratio.png` | Diagnostic plot showing variance explained by the first 50 principal components |
# | `pca_variance_ratio.csv` | Variance explained and cumulative variance explained by each principal component |
# | `leiden_cluster_counts.csv` | Number and fraction of cells assigned to each Leiden cluster |
# | `cluster_composition_by_sample.csv` | Cell counts for each sample within each Leiden cluster |
# | `cluster_composition_by_donor_type.csv` | Cell counts for healthy donors and patients within each Leiden cluster |
# | `10_umap_clusters_donor_sample.png` | UMAP visualization colored by Leiden cluster, donor type, and sample |
# | `marker_availability_hvg.csv` | Records whether canonical immune markers are present in the 2,500-HVG matrix |
# | `11_marker_dotplot_hvg.png` | Canonical marker expression across clusters using genes present in the HVG matrix |
# | `12_marker_umap_hvg.png` | Selected immune-marker expression across UMAP using the HVG matrix |
# | `marker_availability_raw.csv` | Records whether canonical immune markers are present in the complete normalized expression matrix |
# | `13_marker_dotplot_raw.png` | Canonical immune-lineage marker expression across clusters using `adata.raw` |
# | `14_marker_umap_raw.png` | Canonical immune-marker expression across UMAP using `adata.raw` |
# | `clustering_res1.0.h5ad` | Final Section 3 AnnData object containing PCA, neighbors, UMAP, and Leiden clustering results |
# 
# ## Interpretation
# 
# Section 2 prepared the quality-controlled expression matrix for multivariate analysis.
# 
# Section 3 reduces this expression space using PCA, constructs a cell-neighborhood graph, generates a UMAP representation, and identifies transcriptionally similar cell communities using Leiden clustering.
# 
# Canonical PBMC lineage markers are then inspected across these clusters. Marker expression from `adata.raw` is especially important because some biologically informative lineage markers may not have been selected among the 2,500 HVGs.
# 
# The final output is:
# 
# `clustering_res1.0.h5ad`
# 
# which can be used for downstream cell-type annotation and IRC-related analyses.
# 
