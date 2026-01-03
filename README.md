# Multi-Omics Cancer Subtype Discovery (CPTAC UCEC)

![Status](https://img.shields.io/badge/Status-Complete-success)
![Data](https://img.shields.io/badge/Data-CPTAC%20(UCEC)-blueviolet)
![Pipeline](https://img.shields.io/badge/Pipeline-Python%20%2F%20R%20Hybrid-blue)
![Analysis](https://img.shields.io/badge/Analysis-Unsupervised%20Clustering-orange)

## Project Overview

This repository contains a comprehensive **hybrid R/Python workflow** for molecular subtype discovery in **Uterine Corpus Endometrial Carcinoma (UCEC)**. 

By integrating **Proteomics** (TMT mass spec) and **Transcriptomics** (RNA-seq) data, this project moves beyond standard clinical staging to identify molecularly distinct tumour subgroups. The workflow leverages the strengths of both ecosystems:
* **Python:** Used for data engineering (API ingestion), dimensionality reduction, and initial clustering prototypes.
* **R:** Used for rigorous statistical stability testing, "publication-ready" survival analysis, and biomarker discovery (`limma`).

### Key Findings
* **Prototype (Python):** Initial stability analysis identified a broad, stable 2-cluster split ($k=2$), separating tumours roughly by aggressiveness.
* **Refinement (R):** Advanced consensus clustering revealed a more granular **8-subtype structure ($k=8$)**.
* **Validation:** The $k=8$ subtypes showed a statistically significant association with **Histological Grade** (Fisher's Exact Test, $p < 0.001$), confirming biological relevance beyond technical batch effects.

---

## Architecture & Infrastructure

This project was executed in a high-performance computing (HPC) environment to handle large multi-omics matrices.

* **Remote Server (Python):** Data acquisition and heavy-lifting preprocessing were performed via JupyterLab, accessed via **SSH tunnelling** on a headless remote server.
* **HPC Node (R):** Statistical validation and reporting were executed in an RStudio Server environment within the HPC infrastructure.

---

## Data Source

This analysis focuses on cohorts from **CPTAC (Clinical Proteomic Tumour
Analysis Consortium)**.

| Data Type           | Technology                | Typical Input                      | Dimensionality (Approx.)              |
|---------------------|---------------------------|------------------------------------|---------------------------------------|
| **Proteomics**      | TMT / label-free LC–MS/MS | log-intensity / log-ratio matrix   | \~5k–12k proteins × \~100–500 samples |
| **Transcriptomics** | RNA-seq (Illumina)        | raw counts / normalised expression | \~15k–25k genes × \~100–500 samples   |
| **Clinical**        | patient metadata          | categorical + continuous           | survival, grade, stage, etc.          |

---

## Repository Structure & Workflow

The analysis is divided into two pipelines. To reproduce the findings, run the `python_pipeline` first (data acquisition), followed by the `r_analysis` (statistical validation).

### Phase 1: Python Pipeline (Data Engineering & Prototyping)
*Located in `/python_pipeline`*

| Script | Description | Original File |
| :--- | :--- | :--- |
| `01_download_data.ipynb` | Automated download of UCEC clinical & omics data via `cptac` API. | `1A_download_cptac.ipynb` |
| `02_align_samples.ipynb` | Aligns Clinical, Proteomics, and RNA indices (N=95 common samples). | `2A1_align_samples.ipynb` |
| `03_proteomics_qc.ipynb` | Missing value imputation (median) and dropout filtering (>40%). | `2A2_Proteomics_missingness_qc.ipynb` |
| `04_rna_filter.ipynb` | Variance filtering to remove near-constant genes. | `3A3_rna_sanity_filter.ipynb` |
| `05_baseline_model.ipynb` | PCA dimensionality reduction & baseline K-means modelling. | `4B_baseline_subtypes.ipynb` |
| `06_stability_check.ipynb` | Initial ARI (Adjusted Rand Index) stability testing. | `5B_cluster_stability.ipynb` |
| `07_parameter_tuning.ipynb` | Comparison of $k=2..6$. Identified stable $k=2$ prototype. | `6B_pick_k_by_stability.ipynb` |
| `08_prototype_labels.ipynb` | Generating labels for the $k=2$ prototype model. | `7B_final_subtypes_k2.ipynb` |
| `09_clinical_merge.ipynb` | Merging prototype subtypes with clinical metadata. | `8C_attach_subtypes_to_clinical.ipynb` |
| `10_surv_prep.ipynb` | Inspection of survival columns (OS/RFS). | `9C_survival_column_check.ipynb` |
| `11_km_plot_OS.ipynb` | Kaplan-Meier plotting for Overall Survival (Prototype). | `10C_kaplan_meier_OS.ipynb` |
| `12_km_plot_RFS.ipynb` | Kaplan-Meier plotting for Recurrence-Free Survival (Prototype). | `11C_kaplan_meier_RFS.ipynb` |

### Phase 2: R Pipeline (Statistical Refinement & Validation)
*Located in `/r_analysis`*

| Script | Description | Original File |
| :--- | :--- | :--- |
| `01_import_align.Rmd` | Re-importing raw data into the R ecosystem (`tidyverse`). | `2A1_align_samples` |
| `02_imputation_logic.Rmd` | Comparison of Median vs. MinProb imputation strategies. | `2A2_Proteomics_missingness_qc.Rmd` |
| `03_feature_selection.Rmd` | Variance filtering for zero-information genes. | `3A3_rna_sanity_filter.ipynb.Rmd` |
| `04_batch_detection.Rmd` | PCA visual inspection for technical batch effects. | `3A4_batch_detection_pca.Rmd` |
| `05_consensus_cluster.Rmd` | **Core Analysis:** Knee/Elbow plots and Silhouette analysis. | `4B_baseline_subtypes.Rmd` |
| `06_sensitivity_test.Rmd` | **Key Step:** Identified $k=8$ as optimal stable structure via ARI. | `5B_cluster_stability.Rmd` |
| `07_clinical_stats.Rmd` | **Validation:** Monte Carlo Chi-square testing against Histological Grade. | `6C_attach_subtypes_to_clinical.Rmd` |
| `08_survival_analysis.Rmd` | Final Kaplan-Meier curves and Log-rank testing for $k=8$. | `7C_Survival by subtype...Rmd` |

---

## Clustering Strategy: Prototype vs. Refinement

This project utilised a **two-stage clustering approach**, using Python for rapid prototyping and R for rigorous stability refinement.

### 1. Python Prototype (Standard K-Means)

The Python workflow utilised **single-fit K-Means clustering** and evaluated the optimal cluster number ($k$) using global variance metrics:

* **Elbow Method:** Identifies the point of diminishing returns in variance explanation.
* **Silhouette Score:** Measures cohesion (closeness to own cluster) vs. separation (distance to nearest neighbour).

**Result:**
Because these metrics are applied to a single fit of the data, they tend to emphasise the dominant global split. In this analysis, heuristics converged on a coarse **$k=2$ partition**.

> **Note:** The `k = 2` split roughly aligned with tumour grade (Low vs High).“Aggressiveness” was not used during clustering. This interpretation was assessed after clustering by comparing groups against clinical variables.

### 2. R Refinement (Consensus Clustering)

To detect stable substructure beyond the dominant global split, the R workflow used `ConsensusClusterPlus`, which performs resampling-based clustering.

To resolve granular phenotypes, the R workflow utilised `ConsensusClusterPlus`. This method employs **resampling-based clustering** to measure solution stability under perturbation.

**Example setup:**
* **Resampling:** Subsampled 80% of items across 1000 iterations.
* **Consensus:** Calculated how often pairs of samples clustered together across these iterations.
* **Matrix Generation:** Produced a consensus matrix to visualise membership consistency.


This approach prioritises **reproducibility** under sampling variation rather than simple variance explanation.

**Stability Metric:**
* **PAC (Proportion of Ambiguous Clustering):** A low PAC score indicates that clusters are well-resolved, with very few sample pairs falling into the "ambiguous" consensus range (i.e., pairs that sometimes cluster together and sometimes do not).

**Result:**
Stability diagnostics, specifically the minimisation of PAC, supported a more granular **$k=8$ solution**. This structure proved robust to resampling and provided a significantly higher-resolution view of biological heterogeneity compared to the prototype.

---

## Dependencies

**Python:**
* `cptac` (Data ingestion)
* `pandas`, `numpy` (Data manipulation)
* `scikit-learn` (PCA, K-Means, Metrics)
* `lifelines` (Survival Analysis)

**R:**
* `tidyverse` (Data manipulation)
* `ConsensusClusterPlus` (Robust Clustering)
* `limma` (Differential Expression)
* `survival`, `survminer` (Survival Stats & Viz)
* `ComplexHeatmap` (Visualisation)

---

## Directory structure (current)
```bash
Multi-Omics-Cancer-Subtype-Discovery/
├── README.md
├── Multi-Omics-Cancer-Subtype-Discovery.Rproj
├── notebooks/
│   ├── 1A_download_cptac.ipynb
│   ├── 2A1_align_samples.ipynb
│   ├── 2A2_Proteomics_missingness_qc.ipynb
│   ├── 3A3_rna_sanity_filter.ipynb
│   ├── 4B_baseline_subtypes.ipynb
│   ├── 5B_cluster_stability.ipynb
│   ├── 6B_pick_k_by_stability.ipynb
│   ├── 7B_final_subtypes_k2.ipynb
│   ├── 8C_attach_subtypes_to_clinical.ipynb
│   ├── 9C_survival_column_check.ipynb
│   ├── 10C_kaplan_meier_OS.ipynb
│   └── 11C_kaplan_meier_RFS.ipynb
├── R/
│   ├── 2A1_align_samples.Rmd
│   ├── 2A2_Proteomics_missingness_qc.Rmd
│   ├── 3A3_rna_sanity_filter.Rmd
│   ├── 3A4_batch_detection_pca.Rmd
│   ├── 4B_baseline_subtypes.Rmd
│   ├── 5B_cluster_stability.Rmd
│   ├── 6C_attach_subtypes_to_clinical.Rmd
│   └── 7C_survival_by_subtype_clinical_relevance.Rmd
├── results/
│   ├── figures/
│   └── tables/
└── data/
    ├── raw/
    └── processed/
```
---