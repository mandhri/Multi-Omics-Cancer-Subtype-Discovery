Multi-Omics Cancer Subtype Discovery (CPTAC UCEC)
================

- <a href="#project-overview" id="toc-project-overview">Project
  Overview</a>
- <a href="#project-overview-1" id="toc-project-overview-1">Project
  overview</a>
- <a href="#data-source" id="toc-data-source">Data Source</a>
- <a href="#architecture-hpc-remote-setup"
  id="toc-architecture-hpc-remote-setup">Architecture (HPC Remote
  Setup)</a>
- <a href="#analysis-workflow" id="toc-analysis-workflow">Analysis
  Workflow</a>
  - <a href="#1-data-ingestion--alignment-r"
    id="toc-1-data-ingestion--alignment-r">1 Data ingestion + alignment
    (R)</a>
  - <a href="#phase-2--dimensionality-reduction--integration-python"
    id="toc-phase-2--dimensionality-reduction--integration-python">Phase 2 —
    Dimensionality Reduction &amp; Integration (Python)</a>
  - <a href="#phase-3--subtype-discovery--modelling-python"
    id="toc-phase-3--subtype-discovery--modelling-python">Phase 3 — Subtype
    Discovery &amp; Modelling (Python)</a>
  - <a href="#phase-4--visualisation--reporting-r"
    id="toc-phase-4--visualisation--reporting-r">Phase 4 — Visualisation
    &amp; Reporting (R)</a>
- <a href="#quality-checks" id="toc-quality-checks">Quality Checks</a>
- <a href="#directory-structure" id="toc-directory-structure">Directory
  Structure</a>

![Status](https://img.shields.io/badge/Status-In%20Progress-yellow)
![Primary](https://img.shields.io/badge/Primary%20Pipeline-R-brightgreen)
![Supporting](https://img.shields.io/badge/Supporting-Python-blue)
![Infrastructure](https://img.shields.io/badge/Infrastructure-Remote%20%2F%20HPC-orange)

## Project Overview

## Project overview

This repository contains a **hybrid R/Python workflow** that integrates:

- **Mass spectrometry proteomics** (protein abundance)
- **RNA-seq transcriptomics** (gene expression)
- **Clinical variables** (grade, stage, survival)

The goal is to **discover molecular subtypes in cancer**, interpret what
separates them (candidate markers/pathways), and evaluate whether
subtypes relate to **clinical outcomes**.

**Tooling** - **R (tidyverse + Bioconductor)** is used for data
cleaning, QC, statistics, and publication-grade reporting (clinical
association + survival plots). - **Python (scikit-learn / PyTorch)** is
used for exploratory dimensionality reduction, clustering prototypes,
and baseline predictive modelling / feature ranking.

**Infrastructure** This project was executed across two environments: -
**Python (JupyterLab):** run on a headless remote server for exploratory
modelling and prototyping. - **R (Rmd/analysis):** run on a HPC
environment for data engineering, QC, clinical association testing, and
reporting.

------------------------------------------------------------------------

## Data Source

This analysis focuses on cohorts from **CPTAC (Clinical Proteomic Tumour
Analysis Consortium)**.

| Data Type           | Technology                | Typical Input                      | Dimensionality (Approx.)              |
|---------------------|---------------------------|------------------------------------|---------------------------------------|
| **Proteomics**      | TMT / label-free LC–MS/MS | log-intensity / log-ratio matrix   | \~5k–12k proteins × \~100–500 samples |
| **Transcriptomics** | RNA-seq (Illumina)        | raw counts / normalised expression | \~15k–25k genes × \~100–500 samples   |
| **Clinical**        | patient metadata          | categorical + continuous           | survival, grade, stage, etc.          |

------------------------------------------------------------------------

## Architecture (HPC Remote Setup)

This project was executed across two environments:

1.  **Headless remote server (Python / JupyterLab):**
    - Used for the Python notebooks (data download + quick modelling
      prototypes).
    - Accessed from a local browser via **SSH tunnelling** (for example,
      port `8889` for JupyterLab).
2.  **HPC environment (R / reporting + stats):**
    - Used for the R pipeline (data QC, clinical association testing,
      and Rmd reporting).
    - R scripts/Rmd were run directly in the HPC environment (separate
      from the Python runtime).
3.  **Outputs shared via files:**
    - Intermediate outputs are written to `data/processed/` (CSV/RDS)
      and reused across steps.

------------------------------------------------------------------------

## Analysis Workflow

### 1 Data ingestion + alignment (R)

**Goal:** produce clean, matched matrices ready for downstream
modelling.

- **Proteomics QC + missingness review**
  - characterise missing values (common in mass spec)
  - impute using approaches supported by `DEP` / `MSnbase` when
    appropriate
- **Normalisation**
  - RNA-seq: start from counts and transform to ML-friendly values
    (e.g., log-scale normalised expression)
  - Proteomics: typical approaches include median/VSN-style scaling
- **Batch assessment**
  - visualise technical effects (PCA, metadata checks)
  - apply batch correction only if clearly required

**Output:** cleaned, aligned sample-by-feature matrices saved to
`data/processed/`

------------------------------------------------------------------------

### Phase 2 — Dimensionality Reduction & Integration (Python)

**Goal:** learn a compact representation of each patient/sample.

- **Feature filtering**
  - remove near-constant features and obvious noise
- **Dimensionality reduction**
  - PCA / UMAP / t-SNE for visualising sample structure
- **Integration**
  - start with a simple baseline (concatenate scaled omics layers)
  - optionally add a dedicated integration method later (e.g., factor
    models or neural embeddings)

**Output:** integrated latent matrix + embeddings used for clustering.

------------------------------------------------------------------------

### Phase 3 — Subtype Discovery & Modelling (Python)

**Goal:** define subtypes and test whether they matter clinically.

- **Clustering**
  - k-means / Leiden (graph-based) on integrated features
  - compare results across settings to avoid “random” clusters
- **Subtype association with outcomes**
  - survival comparison across subtypes (later plotted in R)
- **Prediction (optional)**
  - train baseline models (Random Forest / XGBoost) for outcome
    prediction
  - extract feature importance to highlight candidate biomarkers

**Output:** subtype labels per patient + model outputs (metrics, feature
importance).

------------------------------------------------------------------------

### Phase 4 — Visualisation & Reporting (R)

**Goal:** publishable figures and interpretable biology.

- **Survival analysis**
  - Kaplan–Meier curves (`survival`, `survminer`) by subtype
- **Heatmaps**
  - top subtype-defining proteins/genes (`ComplexHeatmap`)
- **Pathway analysis**
  - FGSEA / ORA for biological interpretation

**Output:** figures + tables saved to `reports/` (or exported for
manuscripts).

------------------------------------------------------------------------

## Quality Checks

To keep results reliable, this repo will progressively add:

1.  **RNA-seq inputs that are ML-friendly** (avoid relying on FPKM)
2.  **Careful batch handling** (avoid removing real biology)
3.  **A clear proteomics missingness strategy** (conservative vs
    exploratory)
4.  **Cluster stability checks** (not just one run)
5.  **Leakage-safe model evaluation** (proper cross-validation workflow)
6.  **A small set of integration methods** (baseline + 1–2 stronger
    methods)

------------------------------------------------------------------------

## Directory Structure

``` bash
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
