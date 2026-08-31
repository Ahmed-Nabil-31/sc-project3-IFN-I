#!/usr/bin/env python
# coding: utf-8

# In[2]:


import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from importlib.metadata import version

sc.settings.verbosity = 3
sc.settings.set_figure_params(
    dpi=120,
    facecolor="white",
    frameon=False
)

print("Scanpy version:", version("scanpy"))
print("AnnData version:", version("anndata"))


# In[3]:


data_dir = Path("../data")

print("Data directory:", data_dir.resolve())
print("Exists:", data_dir.exists())


# In[4]:


data_dir = Path("../data")

print("Data directory:", data_dir.resolve())
print("Exists:", data_dir.exists())

sample_dirs = sorted([
    p for p in data_dir.iterdir()
    if p.is_dir() and (p / "matrix.mtx.gz").exists()
])

print("Samples found:", [p.name for p in sample_dirs])


# In[6]:


adatas = {}

for sample_path in sample_dirs:
    sample_id = sample_path.name

    print(f"Loading {sample_id}...")

    adata_sample = sc.read_10x_mtx(
        path=sample_path,
        var_names="gene_symbols",
        cache=False
    )

    adata_sample.obs["sample"] = sample_id

    if sample_id.startswith("HD"):
        adata_sample.obs["donor_type"] = "Healthy_Donor"
    else:
        adata_sample.obs["donor_type"] = "Patient"

    adatas[sample_id] = adata_sample

print("Finished loading all samples.")


# In[7]:


for sample_id, adata_sample in adatas.items():
    print(
        f"{sample_id}: "
        f"{adata_sample.n_obs:,} cells × "
        f"{adata_sample.n_vars:,} genes"
    )


# In[8]:


adata = ad.concat(
    adatas,
    label="sample_batch",
    index_unique="-",
    join="outer"
)

adata.var_names_make_unique()

print(
    f"Merged dataset: "
    f"{adata.n_obs:,} cells × {adata.n_vars:,} genes"
)


# In[9]:


adata


# In[10]:


adata.obs.head()


# In[12]:


adata.obs.tail()


# In[12]:


from pathlib import Path

output_dir = Path("../results/preprocessing")
output_dir.mkdir(parents=True, exist_ok=True)


# In[13]:


adata.write_h5ad(output_dir / "merged_raw.h5ad")


# In[14]:


merged_file = output_dir / "merged_raw.h5ad"

print("Saved:", merged_file.resolve())
print(f"File size: {merged_file.stat().st_size / 1024**2:.2f} MB")


# In[23]:


adata.obs


# In[13]:


# Identify gene families used for QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")

adata.var["ribo"] = adata.var_names.str.startswith(
    ("RPS", "RPL")
)

adata.var["hb"] = adata.var_names.str.startswith(
    ("HBA", "HBB", "HBD", "HBE", "HBG", "HBM", "HBQ", "HBZ")
)

print("Mitochondrial genes:", adata.var["mt"].sum())
print("Ribosomal genes:", adata.var["ribo"].sum())
print("Hemoglobin genes:", adata.var["hb"].sum())


# In[14]:


sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=None,
    log1p=False,
    inplace=True
)


# In[15]:


sc.pl.violin(
    adata,
    [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_hb"
    ],
    groupby="donor_type",
    jitter=0.4,
    multi_panel=True
)


# In[16]:


sc.pl.violin(
    adata,
    [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_hb"
    ],
    groupby="sample",
    jitter=0.4,
    multi_panel=True
)


# In[26]:


adata.obs[
    [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb"
    ]
].describe()


# In[17]:


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

qc_summary


# In[ ]:


# Record cell count before filtering
n_before = adata.n_obs

# Filter cells with too few detected genes
sc.pp.filter_cells(adata, min_genes=300)

# Remove high mitochondrial cells
adata = adata[
    adata.obs["pct_counts_mt"] < 15.0,
    :
].copy()

# Remove cells with high hemoglobin contamination
adata = adata[
    adata.obs["pct_counts_hb"] < 5.0,
    :
].copy()

n_after = adata.n_obs

print(f"Cells before QC: {n_before:,}")
print(f"Cells after QC:  {n_after:,}")
print(f"Cells removed:   {n_before - n_after:,}")
print(f"Retention:       {n_after / n_before * 100:.2f}%")


# In[7]:


n_genes_before = adata.n_vars

sc.pp.filter_genes(
    adata,
    min_cells=5
)

n_genes_after = adata.n_vars

print(f"Genes before filtering: {n_genes_before:,}")
print(f"Genes after filtering:  {n_genes_after:,}")
print(f"Genes removed:          {n_genes_before - n_genes_after:,}")


# In[8]:


adata


# In[9]:


print(
    adata.obs["sample"].value_counts().sort_index()
)


# In[11]:


output_dir = Path("../results/preprocessing")
output_dir.mkdir(parents=True, exist_ok=True)

adata.write_h5ad(
    output_dir / "merged_qc_filtered.h5ad"
)

print(
    "Saved:",
    (output_dir / "merged_qc_filtered.h5ad").resolve()
)


# In[1]:


get_ipython().run_line_magic('history', '-f data_loading_qc.py')


# In[ ]:




