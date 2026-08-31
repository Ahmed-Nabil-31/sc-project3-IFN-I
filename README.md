# Pre-Encoded Type I Interferon Responsiveness and Anti-PD-1 Therapy

## Group-3 Members

| Role | Name | Student ID | Email |
|---|---|---|---|
| **Group Leader** | **Suprokash Chakra Borty** | [Student ID] | suprokash.biochem@gmail.com |
| Member | Ahmed Nabil | [Student ID] | an7274003@gmail.com |
| Member | Md Osman Gani Bhuiyan | [Student ID] | mdosmanganibhuiyan@gmail.com |
| Member | Farah Ulfat | [Student ID] | farahulfattasnim@gmail.com |
| Member | Tahmeed Rezwan Shushmoy | [Student ID] | tahmeedrezwan555@gmail.com |
| Member | Shirajum Munira Oyshi | [Student ID] | oyshidu3062@gmail.com |
| Member | Tamanna Dilshad Phul | [Student ID] | tamannadilshad66@gmail.com |
| Member | Zahura Nasreen Akash | [Student ID] | zahuranasreen28@gmail.com |
| Member | Minhaz Abbasi | [Student ID] | mr.mahabub@gmail.com |




---

## Project Overview

### Scientific Question

Type I interferons (IFN-I) are important regulators of anti-tumor immunity, but persistent or excessive IFN-I signaling can also contribute to immune dysfunction and resistance to immunotherapy.

This project investigates **whether pre-existing IFN-I responsiveness in peripheral immune cells is associated with clinical outcome following anti-PD-1 therapy**, and explores the molecular programs associated with different IFN-I response states.

The analysis focuses on single-cell RNA-sequencing data from peripheral blood immune cells and follows a computational workflow covering quality control, normalization, highly variable gene selection, dimensionality reduction, immune-cell characterization, Interferon Response Capacity (IRC) scoring, differential expression, and pathway enrichment.

### Dataset

The dataset consists of single-cell gene-expression profiles from peripheral blood samples containing healthy-donor and cancer-patient immune cells. The analysis workflow identifies samples from healthy donors and patients and combines them into an AnnData object for downstream single-cell analysis.

The raw dataset is large (approximately **94,513 cells × 31,054 genes** before downstream feature selection). After normalization, the workflow identifies **2,500 highly variable genes (HVGs)** for dimensionality-reduction and clustering analyses.

The project uses the following six-gene IRC signature:

- `BST2`
- `EIF2AK2`
- `ISG15`
- `MX1`
- `IFIT3`
- `IRF7`

The IRC score is calculated from the normalized expression of these IFN-I-associated genes. Patient cells are subsequently divided into **IRC-High** and **IRC-Low** groups using the median IRC score as the cutoff.

---

## Analysis Workflow

The computational workflow is organized into five major stages:

```text
Raw 10x Genomics data
        │
        ▼
01. Data loading & quality control
        │
        ▼
02. Normalization & highly variable genes
        │
        ▼
03. PCA, neighborhood graph, UMAP & Leiden clustering
        │
        ▼
04. IFN-I / IRC score calculation
        │
        ▼
05. Differential expression & pathway enrichment
```

### Main analysis steps

1. **Data loading and QC**
   - Read 10x Genomics matrices.
   - Annotate samples as healthy donors or patients.
   - Merge samples.
   - Calculate mitochondrial, ribosomal, and hemoglobin metrics.
   - Remove low-quality cells.
   - Remove genes detected in fewer than five cells.

2. **Normalization and feature selection**
   - Store raw counts.
   - Normalize total counts to 10,000 per cell.
   - Apply log transformation.
   - Identify 2,500 highly variable genes.
   - Regress out total counts and mitochondrial percentage.
   - Scale the selected features.

3. **Dimensionality reduction and clustering**
   - Perform PCA using 50 components.
   - Construct a 20-nearest-neighbor graph using the first 30 PCs.
   - Generate a UMAP embedding.
   - Perform Leiden clustering at resolution 1.0.
   - Characterize clusters using canonical immune-cell markers such as `CD3D`, `CD4`, `CD8A`, `MS4A1`, `CD14`, and `NKG7`.

4. **IRC scoring**
   - Calculate the IFN-I response score using the six-gene IRC signature.
   - Visualize IRC scores across UMAP and samples.
   - Restrict the stratification analysis to patient cells.
   - Divide cells into IRC-High and IRC-Low groups using the median score.
   - Perform Wilcoxon-based differential expression analysis.

5. **Pathway enrichment**
   - Select significantly upregulated genes in IRC-High cells using:
     - `log2 fold change > 0.5`
     - `adjusted p-value < 0.05`
   - Use the top 100 significant genes for enrichment.
   - Analyze Reactome, KEGG, and Gene Ontology Biological Process gene sets.

---

# Quickstart Guide with `uv`

## 1. Prerequisites

Install:

- Git
- Python 3.10+
- `uv`

Check that they are available:

```bash
git --version
python --version
uv --version
```

If `uv` is not installed, follow the official installation instructions:

```bash
# Windows PowerShell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

For Linux/macOS:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Restart the terminal after installation if necessary.

---

## 2. Clone the repository

```bash
git clone https://github.com/Ahmed-Nabil-31/sc-project3-IFN-I.git
cd sc-project3-IFN-I
```

---

## 3. Create the `uv` environment

Create a project virtual environment:

```bash
uv venv
```

Activate it.

### Windows PowerShell

```powershell
.venv\Scripts\Activate.ps1
```

### Windows Command Prompt

```cmd
.venv\Scripts\activate
```

### Linux/macOS

```bash
source .venv/bin/activate
```

---

## 4. Install project dependencies

The repository contains a pinned `requirements.txt` file.

Install the dependencies with:

```bash
uv pip install -r requirements.txt
```

The analysis uses packages including:

- Scanpy
- AnnData
- NumPy
- Pandas
- Matplotlib
- GSEApy
- igraph
- Leidenalg
- Jupyter

---

## 5. Download the project data

The complete single-cell dataset is too large to be stored directly in this GitHub repository.

The project therefore provides the large files through the following Google Drive folder:

**Project data and large files:**

https://drive.google.com/drive/folders/1giT5FnZd7GkM8tmG28j1EH73148zT201?usp=sharing

Download the required files and place them according to the project structure.

For the first analysis stage, the raw 10x Genomics data should be available under:

```text
data/
├── HD.../
│   ├── matrix.mtx.gz
│   ├── barcodes.tsv.gz
│   └── features.tsv.gz
├── ...
└── patient_sample/
    ├── matrix.mtx.gz
    ├── barcodes.tsv.gz
    └── features.tsv.gz
```

Each sample directory must contain a `matrix.mtx.gz` file so that the first script can automatically detect the sample.

For later stages, the following processed files can be used to avoid repeating the most computationally expensive steps:

```text
merged_qc_filtered.h5ad
normalized_hvg_regressed_scaled.h5ad
clustering_res1.0.h5ad
ifn_irc_scored.h5ad
patient_irc_DE.h5ad
IRC_High_vs_Low_DE.csv
pathway_enrichment_results.csv
```

---

## 6. Run the analysis

The repository contains both notebooks and Python scripts:

```text
notebooks/
scripts/
results/
report/
```

The scripts correspond to:

```text
01_data_loading_qc.py
02_normalization_hvg.py
03_dimensionality_reduction.py
04_IRCscore.py
05_pathway_enrichment.py
```

### Start Jupyter Lab

```bash
uv run jupyter lab
```

Then open the notebooks in the following order:

```text
01_data_loading_qc.ipynb
02_normalization_hvg.ipynb
03_dimensionality_reduction.ipynb
04_IRCscore.ipynb
05_pathway_enrichment.ipynb
```

### Recommended execution strategy

The first stage can be run locally if the raw data are available:

```bash
uv run jupyter lab
```

For the later stages, the project workflow is designed to work well in **Google Colab**, particularly because the intermediate AnnData objects are large and the normalization, regression, PCA, clustering, and IRC analyses can require substantial RAM.

If using Colab, upload/copy the processed `.h5ad` and `.csv` files to the expected project directory and execute the notebooks sequentially.

> **Important:** The current notebooks contain some Google Colab/Google Drive-specific commands. Therefore, the README does not claim that every notebook can be executed directly with `uv run python` as a standalone script without modification.

---

# Key Findings

## 1. Peripheral IFN-I responsiveness is associated with anti-PD-1 outcome

The central biological finding is that **pre-treatment IFN-I responsiveness in peripheral effector T cells is associated with response to anti-PD-1 therapy**.

Patients with **lower IFN-I response capacity (IRC) in CD4 and CD8 effector T cells showed better long-term survival**, whereas high pre-existing IFN-I responsiveness was associated with poorer treatment outcome.

This supports the idea that IFN-I signaling is context-dependent: although acute IFN-I signaling can support anti-tumor immunity, a pre-existing highly IFN-I-responsive state may reflect a dysfunctional inflammatory state.

---

## 2. The IRC signature provides a compact representation of IFN-I responsiveness

The analysis focuses on six genes associated with IFN-I response:

```text
BST2
EIF2AK2
ISG15
MX1
IFIT3
IRF7
```

These genes are used to calculate an IRC score for individual cells.

The project therefore converts a high-dimensional single-cell expression profile into a biologically interpretable **IFN-I response score** that can be compared across cells, samples, and immune populations.

---

## 3. Single-cell analysis resolves distinct immune populations

PCA, neighborhood analysis, UMAP, and Leiden clustering reveal distinct immune-cell populations.

Canonical markers are used to identify major lineages, including:

- CD4 T cells
- CD8 T cells
- NK cells
- B cells
- Monocytes
- Dendritic cells

This cell-level resolution is important because the overall IFN-I response of total peripheral blood cells does not necessarily capture the biologically relevant response state of individual immune populations.

---

## 4. IRC-High and IRC-Low cells show different transcriptional programs

Differential expression analysis compares patient cells with high versus low IRC scores.

The pathway-enrichment stage examines genes that are significantly increased in IRC-High cells using:

```text
log2 fold change > 0.5
adjusted p-value < 0.05
```

The resulting genes are analyzed against:

- Reactome
- KEGG
- Gene Ontology Biological Process

The biological interpretation is consistent with the broader finding that highly IFN-I-responsive effector T cells can exhibit programs associated with activation, inflammatory signaling, and T-cell dysfunction.

---

## 5. High IFN-I responsiveness is linked to dysfunctional T-cell states

The underlying biological model suggests that heightened IFN-I responsiveness is not simply a marker of stronger anti-tumor immunity.

Instead, strong pre-existing IFN-I responsiveness can be associated with:

- inflammatory activation,
- suppressive signaling,
- T-cell exhaustion/dysfunction,
- altered metabolic programs, and
- reduced therapeutic benefit from PD-1 blockade.

Thus, **more IFN-I responsiveness is not necessarily better** in the context of chronic cancer-associated inflammation.

---

## 6. IFN-I responsiveness is associated with pre-existing cellular states

Multi-omic evidence from the underlying study indicates that differential IFN-I responsiveness can be associated with pre-existing transcriptional and chromatin states rather than simply reflecting differences in basal ISG expression.

This supports the concept that immune cells may be **epigenetically primed** for different responses to inflammatory stimulation.

---

# Important Interpretation Note

The current IRC stratification code operates at the **cell level**, not strictly at the whole-patient level.

Specifically, patient cells are assigned to `IRC_High` or `IRC_Low` according to each cell's IRC score:

```python
patient_adata.obs["IRC_group"] = np.where(
    patient_adata.obs["IRC_score"] >= median_score,
    "IRC_High",
    "IRC_Low"
)
```

Therefore, the differential-expression analysis should be described as a **cell-level IRC-High vs IRC-Low comparison** unless an explicit patient-level aggregation step is added.

This distinction is important when interpreting statistical results and making claims about patient prognosis.

---

# Repository Structure

```text
sc-project3-IFN-I/
│
├── notebooks/
│   ├── 01_data_loading_qc.ipynb
│   ├── 02_normalization_hvg.ipynb
│   ├── 03_dimensionality_reduction.ipynb
│   ├── 04_IRCscore.ipynb
│   └── 05_pathway_enrichment.ipynb
│
├── scripts/
│   ├── 01_data_loading_qc.py
│   ├── 02_normalization_hvg.py
│   ├── 03_dimensionality_reduction.py
│   ├── 04_IRCscore.py
│   └── 05_pathway_enrichment.py
│
├── results/
│   ├── figures/
│   └── preprocessing_tables/
│
├── report/
│   └── type-I-IFN-report-2.pdf
│
├── requirements.txt
├── Instructions for the reviewer.txt
├── LICENSE
└── README.md
```

---

# Reproducibility

The project records the main software dependencies in `requirements.txt` and uses fixed package versions for important analysis libraries.

For example:

```text
scanpy==1.11.5
anndata==0.12.19
pandas==2.2.3
gseapy==1.3.1
```

A fixed random seed is also used for PCA:

```python
random_state=537
```

This helps make the computational workflow more reproducible.

---

# Results and Figures

Generated figures are available under:

```text
results/figures/
```

The repository includes figures for:

- Data loading and quality control
- Normalization/HVG analysis
- PCA and dimensionality reduction
- UMAP and immune-cell characterization
- IRC scoring
- Differential expression/pathway enrichment

The complete project report is available at:

```text
report/type-I-IFN-report-2.pdf
```

---

# Reference

The biological framework for the project is based on:

**Boukhaled, G. M. et al. (2022). _Pre-encoded responsiveness to type I interferon in the peripheral immune system defines outcome of PD1 blockade therapy._ Nature Immunology, 23, 1273–1283.**

DOI: `10.1038/s41590-022-01262-7`

---

# License

This project is distributed under the license included in the repository.
