# Single-Cell RNA-seq Analysis of Type 1 Diabetes Kidney Biopsies (GSE279086)

[![DOI](https://img.shields.io/badge/GEO-GSE279086-blue)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279086)
[![R](https://img.shields.io/badge/R-4.0%2B-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.8%2B-green)](https://www.python.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.0-orange)](https://satijalab.org/seurat/)

## 📋 Overview

This repository contains a  single-cell RNA-seq analysis pipeline for investigating kidney biopsies from Type 1 Diabetes (T1D) patients compared to Healthy Controls (HC). This analysis directly follows the methodology from Choi et al. (2024) "Attenuated kidney oxidative metabolism in young adults with type 1 diabetes" published in JCI.The analysis includes data preprocessing, quality control, batch correction using Harmony, clustering, and cell type annotation using CellTypist.

**Dataset:** GSE279086  
**Technology:** 10X Genomics scRNA-seq  
**Samples:** 40 samples (12 HC + 28 T1D)  
**Total Cells:** 24,813 cells  
**Genes:** 28,317 genes

---

## 🎯 Project Objectives

1. Process and quality control raw 10X Genomics data
2. Integrate multiple samples using Harmony batch correction
3. Identify cell populations through clustering
4. Annotate cell types using CellTypist
5. Compare cellular composition between HC and T1D conditions
6. Identify disease-associated gene expression patterns

---

## 📁 Repository Structure

```
GSE279086-T1D-Kidney-scRNAseq/
│
├── README.md                          # You are here
├── LICENSE                            # MIT License
├── .gitignore                         # Files not uploaded to GitHub
│
├── scripts/                           # 📜 ALL ANALYSIS CODE
│   ├── 01_download_and_prepare.Rmd    # Download 10X data, create Seurat objects
│   ├── 02_qc_filtering.Rmd            # QC metrics, filtering, visualization
│   ├── 03_Normalization.Rmd           # Normalization, PCA
│   ├── 04_Harmony Integration.Rmd     #Harmony,UMAP,Export files for CellTypist
│   ├── 05_celltypist_annotation.ipynb # PYTHON - Cell type annotation
│   └── 06_DEG and Pathways.Rmd        # T1D vs HC analysis, TCA cycle genes
│
├── data/                              # 💾 Small metadata files (UPLOADED)
│   └── metadata/                     
│       ├── GSE279086_metadata_full.csv    # Complete GEO metadata
│       └── GSE279086_metadata_parsed.csv  # GSM_ID → Condition (12 HC, 28 T1D)
│
├── outputs/                           # 📊 Processed results (SELECTIVE UPLOAD)
│   ├── tables/                       # ✅ UPLOAD - Summary statistics
│   ├── h5ad/                         # ❌ Not Uploaded - Python-compatible files
│   └── rds/                          # ❌ NOT UPLOADED - Too large (>10 GB)
│
├── plots/                            # 🎨 All figures (UPLOAD)
    ├── qc/                          # Quality control plots
    ├── integration/                 # Harmony before/after
    ├── annotation/                  # Cell type UMAPs
    └── differential_expression/     # TCA cycle, volcano plots



---

## 🔧 Installation & Setup

### R Environment

```r
# Install required R packages
install.packages("Seurat")
install.packages("harmony")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("patchwork")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
BiocManager::install("zellkonverter")
```

### Python Environment

In google colab

# Install packages
pip install scanpy
pip install celltypist
pip install anndata
pip install matplotlib
pip install seaborn
```

---

## 🚀 Quick Start

### 1. Data Preparation

```r
# Run the data preparation script
rmarkdown::render("scripts/01_data_preparation.Rmd")
```

This script:
- Loads raw 10X data from 40 samples
- Creates Seurat objects with QC metrics
- Merges samples
- Adds metadata (HC vs T1D conditions)

### 2. Harmony Integration

```r
# Run Harmony integration
rmarkdown::render("scripts/02_integration_harmony.Rmd")
```

This script:
- Normalizes and scales data
- Performs PCA
- Runs Harmony batch correction
- Performs clustering (33 clusters identified)
- Generates UMAP visualizations

### 3. Cell Type Annotation

```bash
# Run CellTypist annotation in Python
jupyter notebook scripts/03_celltypist_annotation.ipynb
```

This notebook:
- Loads Harmony-integrated h5ad file
- Runs CellTypist with immune cell models
- Assigns cell type labels
- Visualizes results

---

## 📊 Key Results

### Sample Distribution
- **Healthy Controls (HC):** 5,836 cells (23.5%)
- **Type 1 Diabetes (T1D):** 18,977 cells (76.5%)

### Quality Metrics
- **Mean genes/cell:** ~2,000
- **Mean UMIs/cell:** ~8,000
- **Mitochondrial %:** <10% (after filtering)

### Clustering
- **Pre-integration:** 33 clusters (PCA-based)
- **Post-Harmony:** Batch effects reduced
- **Cell types identified:** [To be added after CellTypist annotation]

---

## 📈 Workflow Overview


RAW DATA → QC → MERGE → HARMONY → UMAP → CELLTYPIST → T1D vs HC
 40 samples    ✓     ✓        ✓        ✓        ✓           ✓
               │      │        │        │        │           │
               ▼      ▼        ▼        ▼        ▼           ▼
            Remove   Add     Remove   See     Identify   Find TCA
            bad     condi-   batch    biol-   kidney     cycle
            cells   tions    effects  ogy     cell       genes
                                       types

```mermaid
graph TD
    A[Raw 10X Data<br/>40 samples] --> B[Quality Control<br/>Filter cells & genes]
    B --> C[Merge Samples<br/>24,813 cells]
    C --> D[Normalization<br/>LogNormalize]
    D --> E[Find Variable Features<br/>2,500 genes]
    E --> F[Scale Data<br/>PCA]
    F --> G[Harmony Integration<br/>Batch correction]
    G --> H[UMAP & Clustering<br/>33 clusters]
    H --> I[Export to h5ad<br/>For Python]
    I --> J[CellTypist Annotation<br/>Cell types]
    J --> K[Downstream Analysis<br/>DEGs, pathways]
```

---

## 📝 Methods Summary

### Quality Control
- **Filters applied:**
  - Genes detected in ≥3 cells
  - Cells with ≥200 features
  - Mitochondrial content < 20%
  - Ribosomal content calculated

### Normalization
- **Method:** LogNormalize (Seurat)
- **Scale factor:** 10,000
- **Variable features:** 2,500 genes (vst method)

### Batch Correction
- **Method:** Harmony
- **Batch variable:** Sample ID (orig.ident)
- **PCs used:** 1-50
- **Max iterations:** 20

### Clustering
- **Method:** Louvain algorithm (Seurat)
- **Resolution:** 0.8
- **Reduction:** Harmony
- **Dimensions:** 1-50

### Cell Type Annotation
- **Tool:** CellTypist
- **Model:** Immune_All_Low.pkl
- **Majority voting:** Enabled
 
### DEG & Pathway Analysis
- 
---

## 📂 Important Files to Keep

### Essential RDS Files
1. **`02_seurat_merged.rds`** (424.8 MB)
   - Merged object with all metadata
   - Use for re-running analyses from merged state

2. **`03_seurat_pca_umap_JOINED.rds`** (5 GB)
   - Normalized, PCA, UMAP, pre-Harmony
   - Use for comparing pre/post Harmony

3. **`04_harmony_integrated.rds`** (TBD)
   - Final Harmony-integrated object
   - **PRIMARY FILE for downstream analysis**

### H5AD File
1. **`04_harmony_integrated.h5ad`** (534.218 KB)
   - For Python/CellTypist
   - Compressed, no scale.data

### Metadata Files
1. **`GSE279086_metadata_full.csv`** (57.3 KB)
   - Complete sample metadata with conditions

### Plots Directory
- Keep all plots organized by analysis stage

---


## 🔬 Citation

If you use this analysis pipeline, please cite:

```
Dijith A R. (2026). Single-Cell RNA-seq Analysis of Type 1 Diabetes Kidney Biopsies (GSE279086). GitHub repository: 
https://github.com/ardijith/GSE279086-scRNAseq
```

**Original Dataset:**
```
Choi YJ, Richard G, Zhang G, et al. (2024). Attenuated kidney oxidative 
metabolism in young adults with type 1 diabetes. J Clin Invest. 
134(24):e183984. https://doi.org/10.1172/JCI183984
GEO Accession: GSE279086
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **Seurat** - Spatial analysis framework
- **Harmony** - Batch correction algorithm
- **CellTypist** - Automated cell type annotation
- **GEO** - Data repository (GSE279086)

---

## 📚 References

1. Hao et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*
2. Korsunsky et al. (2019). Fast, sensitive and accurate integration. *Nature Methods*
3. Domínguez Conde et al. (2022). Cross-tissue immune cell analysis. *Science*

---

**Last Updated:** February 11, 2026  
**Status:** ✅ Active Development
