# Point-by-Point Response to Reviewers

**Manuscript Title**: Pre-Encoded Type I Interferon Responsiveness and Anti-PD-1 Therapy Outcomes  
**Authors**: Group 3 (Project 3) members
**Assigned Reviewers**: Group 2 (Drosophila Brain Project)  
**Date**: September 14, 2026  

Dear Reviewers and Mentors,

We sincerely thank Group 2 for their rigorous, constructive, and insightful evaluation of our manuscript and single-cell analysis workflow (`sc-project3-IFN-I`). The feedback regarding doublet detection, path portability, figure formatting, and statistical terminology has substantially improved the quality and reproducibility of our study.

Below, we provide our point-by-point responses detailing the analytical modifications, code updates, and text revisions implemented in the revised release (`v2.0-final`).

---

## Response to Major Comments

### Major Comment 1: Doublet Detection in Quality Control Workflow
> *Reviewer Comment*: "QC detects gene count, hemoglobin%, and mito% but does not include a doublet detection step. After filtering, the very high total_counts values increase the possibility that there are some doublets."

**Author Response**:  
We thank the reviewers for highlighting this important QC requirement. In our updated pipeline (scripts/01_data_loading_qc.py), we integrated Scrublet-based doublet detection independently for each sample. Scrublet doublet scores and predicted doublet labels were generated from the observed count matrices, and predicted doublets were removed prior to downstream normalization and clustering. This step removed 7,065 predicted doublets, reducing the dataset from 94,513 cells after basic QC to 87,448 cells. We additionally filtered cells with fewer than 300 detected genes, mitochondrial transcript content ≥15%, and hemoglobin transcript content ≥5%. No fixed upper UMI-count threshold was applied.

---

## Response to Minor Comments

### Minor Comment 1: Figure Organization and Labeling
> *Reviewer Comment*: "The figures could be neatly organized, and more labeling would be useful."

**Author Response**:  
We restructured the output directory (`results/figures/`) into sequential stage sub-folders (`01_qc`, `02_norm`, `03_dimred`, `04_irc`, `05_pathways`). Since the main labeling issue was the volcano plot that labelled the top genes, a table was provided just below it with the names of the top genes and there parameters. 

---

### Minor Comment 2: Figure Caption Conciseness
> *Reviewer Comment*: "The figure captions are quite lengthy and repetitive."

**Author Response**:  
We revised all figure captions in the manuscript to be concise and shortened them. Their description was explained in the discussion part rather then their captions, which were primarily in the Results section. 

---

### Minor Comment 3: Citation Style Consistency
> *Reviewer Comment*: "The reference section citation style is inconsistent."

**Author Response**:  
We standardized all references in the bibliography and in-text citations in APA 7th format.

---

### Minor Comment 4: Typo in Statistical Correction Method
> *Reviewer Comment*: "In section 2.5, method a-path is incorrectly written as 'Benjamini-Huang'. It should be 'Benjamini-Hochberg'."

**Author Response**:  
We corrected the typographical error in Section 2.5 of the report text from "Benjamini-Huang" to **"Benjamini-Hochberg"** multiple-testing correction.

---

## Response to Reproducibility & Execution Comments

### Execution Comment: Google Drive / Google Colab Hardcoded Paths
> *Reviewer Comment*: "There were some Google Drive/Google Colab paths provided, which caused issues during installation. Due to these errors, we were unable to run the notebook properly and generate figures and tables."

**Author Response**:  
We apologize for the execution failure caused by residual Colab/Drive path references. All notebooks and scripts in `scripts/` and `notebooks/` have been refactored to use dynamic, relative paths (`pathlib.Path`). We verified clean end-to-end execution within a fresh `uv` environment (`uv venv` and `uv pip install -r requirements.txt`) on a standard local machine. **(Just run the notebook folders sequentially after downloading the main data provided by our mentor and place it in a folder named data, within your local repository. As they contain huge files, they are excluded from the github repository.)**

---

## Summary of Code & Repository Updates

All code updates, structured figures, and manuscript revisions have been pushed to GitHub under release tag `v2.0-final`:
* `scripts/01_data_loading_qc.py`: Added Scrublet doublet detection and UMI upper bound filtering.
* `scripts/ & notebooks/`: Refactored file loading to relative project-root paths.
* `results/figures/`: Re-organized into stage-specific subdirectories with enhanced plot labels.
* `report/`: Corrected Section 2.5 typo ("Benjamini-Hochberg"), condensed figure captions, and unified citation styling.

Sincerely,  
**Group 3 Project Team**  
*(Ahmed Nabil, Suprokash Chakra Borty, Md Osman Gani Bhuiyan, Shirajum Munira Oyshi, Farah Ulfat, Tahmeed Rezwan Shushmoy, Tamanna Dilshad Phul, Zahura Nasreen Akash, Minhaz Abbasi)*
