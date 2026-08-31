#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().system('pip install -q      "pandas==2.2.3"      "scanpy==1.11.5"      "anndata==0.12.19"      "gseapy==1.3.1"')


# In[ ]:





# In[ ]:


# Import Scanpy for single-cell RNA-seq analysis
import scanpy as sc

# Import AnnData for handling .h5ad files
import anndata as ad

# Import NumPy for numerical operations
import numpy as np

# Import pandas for handling metadata and tables
import pandas as pd

# Import Matplotlib for visualization
import matplotlib.pyplot as plt


# Configure Scanpy output
sc.settings.verbosity = 3

sc.settings.set_figure_params(
    dpi=120,
    facecolor="white",
    frameon=False
)


# Record the software versions used for this analysis
print("Scanpy:", sc.__version__)
print("AnnData:", ad.__version__)
print("Pandas:", pd.__version__)


# In[3]:


# Mount Google Drive so the processed .h5ad files
# can be accessed from persistent storage.
from google.colab import drive

drive.mount("/content/drive")


# In[4]:


# Load the dataset produced at the end of Step 8.4.
# This contains PCA, neighbors, UMAP, and Leiden clustering.

adata = sc.read_h5ad(
    "/content/drive/MyDrive/sc-project3-IFN-I/clustering_res1.0.h5ad"
)

print(adata)

print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")


# In[5]:


# Confirm that the Leiden cluster assignments from Step 8.4
# are present in the cell metadata.

print(adata.obs.columns.tolist())


# In[6]:


# Check whether the full normalized/log-transformed expression
# matrix was preserved in adata.raw during Step 8.3.

print(adata.raw)


# In[7]:


# Define the six core interferon-stimulated protein (ISP) genes
# identified by the publication as the IRC signature.

isg_signature = [
    "BST2",
    "EIF2AK2",
    "ISG15",
    "MX1",
    "IFIT3",
    "IRF7"
]


# In[8]:


# Check which IRC genes are present in the full normalized
# expression matrix stored in adata.raw.

valid_isgs = [
    gene for gene in isg_signature
    if gene in adata.raw.var_names
]

print(
    f"Calculating IRC score using "
    f"{len(valid_isgs)}/{len(isg_signature)} genes:"
)

print(valid_isgs)


# In[9]:


# Calculate the Interferon Response Capacity (IRC) module score.
#
# use_raw=True ensures that all six signature genes are taken
# from the full normalized/log-transformed expression matrix
# rather than only the 2,500 HVGs.

sc.tl.score_genes(
    adata,
    gene_list=valid_isgs,
    score_name="IRC_score",
    use_raw=True
)

print("IRC score calculated successfully.")


# In[10]:


adata.obs["IRC_score"]


# In[11]:


# Examine the distribution of IRC scores across all cells.

print(adata.obs["IRC_score"].describe())


# In[13]:


# Visualize the IRC score across the existing UMAP.
#
# The UMAP itself is unchanged; we are simply coloring
# cells according to their IRC score.

sc.pl.umap(
    adata,
    color=[
        "IRC_score",
        "donor_type",
        "leiden_res_1.0"
    ],
    cmap="viridis",
    ncols=1
)


# In[14]:


# Compare the distribution of IRC scores between samples.
# Each violin represents the distribution of IRC scores
# among cells belonging to that sample.

sc.pl.violin(
    adata,
    keys="IRC_score",
    groupby="sample",
    rotation=45
)


# In[15]:


# Save the AnnData object containing the IRC score.
#
# The IRC score is stored in:
# adata.obs["IRC_score"]

adata.write_h5ad(
    "/content/drive/MyDrive/sc-project3-IFN-I/ifn_irc_scored.h5ad"
)

print("8.6 IRC-scored dataset saved successfully.")


# In[16]:


# Select only cells from patients.
#
# This removes non-patient/donor cells from the analysis
# so that IRC-High vs IRC-Low stratification is performed
# specifically within the patient population.

patient_adata = adata[
    adata.obs["donor_type"] == "Patient"
].copy()


# Check how many patient cells remain
print(
    f"Patient cells: {patient_adata.n_obs:,}"
)


# In[17]:


# Calculate the median IRC score among patient cells.
#
# The median will be used as the cutoff to divide patient
# cells into two groups with relatively high or low
# interferon-response capacity.

median_score = patient_adata.obs["IRC_score"].median()

print(f"Median IRC score: {median_score:.4f}")


# In[18]:


# Assign each patient cell to an IRC group.
#
# Cells with IRC scores greater than or equal to the median
# are classified as IRC_High.
#
# Cells below the median are classified as IRC_Low.

patient_adata.obs["IRC_group"] = np.where(
    patient_adata.obs["IRC_score"] >= median_score,
    "IRC_High",
    "IRC_Low"
)


# Count cells in each IRC group
print(
    patient_adata.obs["IRC_group"].value_counts()
)


# In[19]:


# Perform differential expression analysis between
# IRC_High and IRC_Low patient cells.
#
# The Wilcoxon rank-sum test is used to identify genes
# whose expression differs between the two IRC groups.

sc.tl.rank_genes_groups(
    patient_adata,
    groupby="IRC_group",
    method="wilcoxon"
)


# In[20]:


# Extract the differential-expression results for
# the IRC_High group.
#
# Positive-ranked genes represent genes associated with
# higher expression in IRC_High cells.

de_df = sc.get.rank_genes_groups_df(
    patient_adata,
    group="IRC_High"
)


# Display the top 15 genes associated with IRC_High cells
print("Top upregulated genes in IRC-High patients:")
print(de_df.head(15))


# One thing I want you to be careful about
# 
# Your rulebook calls this:
# 
# Patient Stratification & Differential Expression
# 
# but the actual code stratifies patient cells, not necessarily entire patients.
# 
# That's an important distinction.
# 
# This:
# 
# patient_adata = adata[
#     adata.obs["donor_type"] == "Patient"
# ].copy()
# 
# means you're taking all cells belonging to patients and then assigning each cell to IRC-High or IRC-Low based on its individual IRC score.
# 
# So this:
# 
# IRC_High vs IRC_Low
# 
# is currently a cell-level comparison, not:
# 
# Patient 1 = IRC_High
# Patient 2 = IRC_Low
# Patient 3 = IRC_High
# ...
# 
# That's fine if it's what the publication/rulebook intends, but it's something we should explicitly recognize because it affects how the DE results should be interpreted.

# In[21]:


# Save the patient subset with IRC group assignments
# and differential-expression results.
#
# The DE results themselves are stored in adata.uns,
# while IRC_group is stored in patient_adata.obs.

patient_adata.write_h5ad(
    "/content/drive/MyDrive/sc-project3-IFN-I/patient_irc_DE.h5ad"
)

# Also save the differential-expression table as a CSV
# so that it can be easily opened and shared.

de_df.to_csv(
    "/content/drive/MyDrive/sc-project3-IFN-I/IRC_High_vs_Low_DE.csv",
    index=False
)

print("Step 8.7 results saved successfully.")

