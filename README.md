
---

Multi-Omics-Cancer-Subtype-Discovery Integration: A Hybrid R/Python
Pipeline
output: 
  github_document:
    toc: true
    toc_depth: 3

---


![Status](https://img.shields.io/badge/Status-In%20Progress-yellow)
![Language](https://img.shields.io/badge/Languages-R%20%7C%20Python-blue)
![Infrastructure](https://img.shields.io/badge/Infrastructure-HPC-orange)

## Project Overview

This repository hosts a **scalable hybrid workflow** for integrating
Mass Spectrometry Proteomics and Transcriptomics data to identify
molecular subtypes in cancer.

This pipeline utilises **R (Tidyverse)** for rigorous data engineering
and visualisation, and **Python (Scikit-Learn/PyTorch)** for
high-dimensional clustering and predictive modelling. The entire
workflow runs on a “headless” remote server (HPC) accessed via **SSH
Tunnelling**.

## Data Sources & Types

The analysis focuses on the **CPTAC (Clinical Proteomic Tumour Analysis
Consortium)** cohort.

| Data Type           | Technology                | Format                         | Dimensionality (Approx)         |
|:--------------------|:--------------------------|:-------------------------------|:--------------------------------|
| **Proteomics**      | TMT / Label-free LC-MS/MS | Abundance Ratio (Log)          | \~10,000 Proteins x 500 Samples |
| **Transcriptomics** | RNA-Seq (Illumina)        | FPKM / Counts                  | \~20,000 Genes x 500 Samples    |
| **Clinical**        | EMR / Pathology           | Mixed (Categorical/Continuous) | Survival, Grade, Stage          |

## Architecture

The project solves the “Infrastructure Gap” using a manual tunnelling
approach:

1.  **Remote Server (Backend):** High-Memory HPC Compute Node holding
    the 50GB+ datasets.
2.  **Local Client (Frontend):** SSH Tunnel forwarding port `8889`
    (JupyterLab) and `8787` (RStudio) to the local browser.
3.  **Environment:** A unified `conda` environment (`python_ai_env`)
    shared by both RStudio (`reticulate`) and JupyterLab.

## Analysis Workflow

### Phase 1: Data Ingestion & Engineering (R)

**Goal:** Cleanse and harmonise heterogeneous omics data.

* **Missing Value Imputation:** Handling missing Mass Spec values (MNAR/MAR) using `MSnbase` / `DEP`.
* **Normalisation:** TMM normalisation for RNA-seq and VSN/Median centring for Proteomics.
* **Batch Correction:**  Assessment and removal of technical batch effects using PCA/SVA.
* *Output:* Cleaned, matched `.csv` matrices.

### Phase 2: Dimensionality Reduction & Integration (Python)

**Goal:** Discover latent biological signals in high-dimensional space.
* **Feature Selection:** Variance filtering to remove static features.
* **Dimensionality Reduction:** UMAP and t-SNE projections to visualise patient manifold.
* **Integration:** Fusing RNA and Protein layers (Potential methods: MOFA+, Autoencoders).

### Phase 3: Subtype Discovery & Classification (Python)

**Goal:** Define and predict molecular subtypes. 
* **Clustering:** Unsupervised clustering (K-Means / Leiden) on integrated features to
define “Novel Subtypes.” 
* **Prediction:** Training Random Forest/ XGBoost classifiers to predict Clinical Survival from Omics features. 
* **Feature Importance:** Extracting the top 50 Biomarkers driving the separation.

### Phase 4: Visualisation & Reporting (R)

**Goal:** Publication-quality figures. 
* **Survival Analysis:** Kaplan-Meier curves (`survival`, `survminer`) comparing the new Subtypes.
* **Heatmaps:** ComplexHeatmap visualisation of the top differential proteins.
* **Pathway Analysis:** FGSEA / ORA to understand the biological context of the clusters.

## Directory Structure

``` bash
├── data/
│   ├── raw/                 # Original CPTAC downloads (Gitignored)
│   ├── processed/           # Cleaned matrices ready for ML
├── R/
│   ├── 01_data_cleaning.R   # Tidyverse scripts for QC
│   ├── 04_visualization.R   # ComplexHeatmap & Survival plots
├── python/
│   ├── 02_integration.ipynb # Jupyter Notebooks for UMAP/Clustering
│   ├── 03_classification.py # ML training scripts
├── envs/
│   └── environment.yml      # Conda environment specification
└── README.md                # Project Blueprint
```
