
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

This repository contains a **scalable hybrid workflow** that integrates:

- **Mass spectrometry proteomics** (protein abundance)
- **RNA-seq transcriptomics** (gene expression)
- **Clinical variables** (e.g., stage, grade, survival)

The goal is to **discover molecular subtypes in cancer**, interpret what separates them (biomarkers/pathways), and evaluate whether these subtypes relate to **clinical outcomes**.

The workflow uses:
- **R (tidyverse + bioconductor)** for data cleaning, QC, statistics, and publication-grade plots
- **Python (scikit-learn / PyTorch)** for dimensionality reduction, clustering, and predictive modelling

The pipeline is designed to run on a **headless HPC/remote server**, accessed locally via **SSH tunnelling** (RStudio + JupyterLab in the browser).

---

## Data Source

This analysis focuses on cohorts from **CPTAC (Clinical Proteomic Tumour Analysis Consortium)**.

| Data Type | Technology | Typical Input | Dimensionality (Approx.) |
|---|---|---|---|
| **Proteomics** | TMT / label-free LC–MS/MS | log-intensity / log-ratio matrix | ~5k–12k proteins × ~100–500 samples |
| **Transcriptomics** | RNA-seq (Illumina) | raw counts / normalised expression | ~15k–25k genes × ~100–500 samples |
| **Clinical** | patient metadata | categorical + continuous | survival, grade, stage, etc. |

> Note: exact dimensions vary by tumour type and CPTAC release.

---

## Architecture (HPC-Friendly Setup)

This project is built to work well with large datasets and remote compute:

1. **Remote server (backend):** HPC compute node with high memory + local storage for large CPTAC files  
2. **Local browser (frontend):** SSH tunnel forwards:
   - `8889` → JupyterLab
   - `8787` → RStudio Server
3. **Environment:** one shared `conda` env (e.g., `python_ai_env`)
   - R can call Python via `reticulate`
   - Python notebooks/scripts run in the same environment

---

## Analysis Workflow

### Phase 1 — Data Ingestion & Engineering (R)
**Goal:** produce clean, matched matrices ready for downstream modelling.

- **Proteomics QC + missingness review**
  - characterise missing values (common in mass spec)
  - impute using approaches supported by `DEP` / `MSnbase` when appropriate
- **Normalisation**
  - RNA-seq: start from counts and transform to ML-friendly values (e.g., log-scale normalised expression)
  - Proteomics: typical approaches include median/VSN-style scaling
- **Batch assessment**
  - visualise technical effects (PCA, metadata checks)
  - apply batch correction only if clearly required

**Output:** cleaned, aligned sample-by-feature matrices saved to `data/processed/`

---

### Phase 2 — Dimensionality Reduction & Integration (Python)
**Goal:** learn a compact representation of each patient/sample.

- **Feature filtering**
  - remove near-constant features and obvious noise
- **Dimensionality reduction**
  - PCA / UMAP / t-SNE for visualising sample structure
- **Integration**
  - start with a simple baseline (concatenate scaled omics layers)
  - optionally add a dedicated integration method later (e.g., factor models or neural embeddings)

**Output:** integrated latent matrix + embeddings used for clustering.

---

### Phase 3 — Subtype Discovery & Modelling (Python)
**Goal:** define subtypes and test whether they matter clinically.

- **Clustering**
  - k-means / Leiden (graph-based) on integrated features
  - compare results across settings to avoid “random” clusters
- **Subtype association with outcomes**
  - survival comparison across subtypes (later plotted in R)
- **Prediction (optional)**
  - train baseline models (Random Forest / XGBoost) for outcome prediction
  - extract feature importance to highlight candidate biomarkers

**Output:** subtype labels per patient + model outputs (metrics, feature importance).

---

### Phase 4 — Visualisation & Reporting (R)
**Goal:** publishable figures and interpretable biology.

- **Survival analysis**
  - Kaplan–Meier curves (`survival`, `survminer`) by subtype
- **Heatmaps**
  - top subtype-defining proteins/genes (`ComplexHeatmap`)
- **Pathway analysis**
  - FGSEA / ORA for biological interpretation

**Output:** figures + tables saved to `reports/` (or exported for manuscripts).

---

## Quality Checks

To keep results reliable, this repo will progressively add:

1. **RNA-seq inputs that are ML-friendly** (avoid relying on FPKM)
2. **Careful batch handling** (avoid removing real biology)
3. **A clear proteomics missingness strategy** (conservative vs exploratory)
4. **Cluster stability checks** (not just one run)
5. **Leakage-safe model evaluation** (proper cross-validation workflow)
6. **A small set of integration methods** (baseline + 1–2 stronger methods)

---

## Directory Structure

```bash
Multi-Omics-Cancer-Subtype-Discovery/
├── README.md
├── Multi-Omics-Cancer-Subtype-Discovery.Rproj
├── notebooks/                
├── R/                         
│   ├── 00_config.R
│   └── functions/
├── analysis/                  
│   └── report.Rmd            
├── results/                  
│   ├── figures/
│   └── tables/
└── data/                      
    ├── raw/                  
    └── processed/            


```

