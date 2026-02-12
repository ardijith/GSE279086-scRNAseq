# File Organization Guide for GSE279086 Project

## 📋 Complete File Inventory

### ✅ KEEP - Essential Files

#### **RDS Objects (Seurat)**

| File | Size | Purpose | Keep? |
|------|------|---------|-------|
| `02_seurat_merged.rds` | 424.8 MB | Merged object with metadata | ✅ YES |
| `03_GSE279086_seurat_pca_umap_JOINED.rds` | 5 GB | Pre-Harmony, normalized, PCA, UMAP | ✅ YES |
| `04_harmony_integrated.rds` | ~5 GB | **PRIMARY** - Harmony integrated | ✅ YES |

#### **H5AD Files (Python)**

| File | Size | Purpose | Keep? |
|------|------|---------|-------|
| `04_harmony_integrated.h5ad` | 534 KB | For CellTypist annotation | ✅ YES |

#### **Metadata & Summary Tables**

| File | Size | Purpose | Keep? |
|------|------|---------|-------|
| `GSE279086_metadata_full.csv` | 57.3 KB | Complete sample metadata | ✅ YES |
| `01_sample_summary.csv` | 3 KB | Sample statistics | ✅ YES |
| `pca_variance.csv` | 4.1 KB | PCA variance explained | ✅ YES |
| `GSE279086_QC_metrics.csv` | 1.1 KB | QC summary | ✅ YES |
| `qc_filtering_summary.csv` | 149 B | Filtering stats | ✅ YES |

#### **Scripts**

| File | Size | Purpose | Keep? |
|------|------|---------|-------|
| `03_10XtoSeurat_till_h5ad_conv.Rmd` | - | Main R analysis script | ✅ YES |
| `02_celltypist_GSE183276.ipynb` | - | Python CellTypist notebook | ✅ YES |

---

### ❌ DELETE - Redundant/Intermediate Files

#### **Duplicate RDS Files**

| File | Size | Reason to Delete |
|------|------|------------------|
| `GSE279086_merged.rds` | 101.9 MB | Old version with 81 layers issue |
| `GSE279086_seurat_processed.rds` | 5 GB | Redundant with 03_* file |
| `seurat_filtered.rds` | 103 MB | Intermediate QC step |
| `seurat_normalized.rds` | 5 GB | Intermediate normalization step |
| `seurat_with_conditions.rds` | 424.8 MB | Duplicate of 02_seurat_merged.rds |

#### **H5Seurat Files (Not Needed)**

| File | Size | Reason to Delete |
|------|------|------------------|
| `GSE279086_fixed.h5seurat` | 976 KB | Not needed (have RDS + h5ad) |
| `GSE279086_merged.h5seurat` | 100.9 KB | Not needed |
| `GSE279086_seurat.h5seurat` | 976.7 KB | Not needed |

#### **Other Files**

| File | Size | Reason to Delete |
|------|------|------------------|
| `.Rhistory` | 12 KB | R command history |
| `GSE279086_metadata_parsed.csv` | 1.5 KB | Redundant (keep full version) |

---

## 📁 Recommended Directory Structure

```
GSE279086-T1D-scRNAseq/
│
├── data/
│   ├── raw/                           # Original 10X data (gitignored)
│   │   ├── GSM8561110/
│   │   ├── GSM8561111/
│   │   └── ...
│   │
│   └── metadata/
│       └── GSE279086_metadata_full.csv
│
├── outputs/
│   ├── rds/
│   │   ├── 02_seurat_merged.rds
│   │   ├── 03_seurat_pca_umap_JOINED.rds
│   │   └── 04_harmony_integrated.rds
│   │
│   ├── h5ad/
│   │   └── 04_harmony_integrated.h5ad
│   │
│   └── tables/
│       ├── 01_sample_summary.csv
│       ├── pca_variance.csv
│       ├── GSE279086_QC_metrics.csv
│       └── qc_filtering_summary.csv
│
├── plots/
│   ├── qc/
│   │   └── cells_per_sample.png
│   │
│   ├── integration/
│   │   └── elbow_plot.png
│   │
│   └── annotation/
│       └── celltypist_umap.png
│
├── scripts/
│   ├── 01_data_preparation.Rmd
│   ├── 02_integration_harmony.Rmd
│   └── 03_celltypist_annotation.ipynb
│
├── environment/
│   ├── R_sessionInfo.txt
│   └── requirements.txt
│
├── README.md
├── LICENSE
└── .gitignore
```

---

## 🔧 Step-by-Step Organization

### Step 1: Create Directory Structure

```bash
# Create main directories
mkdir -p data/raw
mkdir -p data/metadata
mkdir -p outputs/rds
mkdir -p outputs/h5ad
mkdir -p outputs/tables
mkdir -p plots/qc
mkdir -p plots/integration
mkdir -p plots/annotation
mkdir -p scripts
mkdir -p environment
```

### Step 2: Move Files to Correct Locations

```bash
# Move RDS files
mv 02_seurat_merged.rds outputs/rds/
mv 03_GSE279086_seurat_pca_umap_JOINED.rds outputs/rds/
mv 04_harmony_integrated.rds outputs/rds/

# Move h5ad files
mv 04_harmony_integrated.h5ad outputs/h5ad/

# Move metadata
mv GSE279086_metadata_full.csv data/metadata/

# Move summary tables
mv 01_sample_summary.csv outputs/tables/
mv pca_variance.csv outputs/tables/
mv GSE279086_QC_metrics.csv outputs/tables/
mv qc_filtering_summary.csv outputs/tables/

# Move scripts
mv 03_10XtoSeurat_till_h5ad_conv.Rmd scripts/01_data_preparation.Rmd
mv 02_celltypist_GSE183276.ipynb scripts/03_celltypist_annotation.ipynb

# Move plots (you'll need to identify these)
# mv *.png plots/qc/  # or plots/integration/ or plots/annotation/
```

### Step 3: Delete Redundant Files

```bash
# Delete duplicate RDS files
rm GSE279086_merged.rds
rm GSE279086_seurat_processed.rds
rm seurat_filtered.rds
rm seurat_normalized.rds
rm seurat_with_conditions.rds

# Delete h5seurat files
rm GSE279086_fixed.h5seurat
rm GSE279086_merged.h5seurat
rm GSE279086_seurat.h5seurat

# Delete other files
rm .Rhistory
rm GSE279086_metadata_parsed.csv
```

### Step 4: Create Environment Files

```r
# In R, save session info
writeLines(capture.output(sessionInfo()), "environment/R_sessionInfo.txt")
```

```bash
# Create Python requirements.txt
cat > environment/requirements.txt << EOF
scanpy==1.9.3
celltypist==1.6.0
anndata==0.9.1
matplotlib==3.7.1
seaborn==0.12.2
pandas==2.0.1
numpy==1.24.3
scipy==1.10.1
EOF
```

---

## 📊 File Size Summary

### Essential Files (Keep)
- **Total RDS:** ~10.5 GB (3 files)
- **H5AD:** 534 KB (1 file)
- **Metadata/Tables:** ~66 KB (5 files)
- **Scripts:** ~50 KB (2-3 files)
- **Total to Keep:** ~10.6 GB

### Files to Delete
- **Total to Delete:** ~12 GB

### Space Savings
- **Before:** ~22.6 GB
- **After:** ~10.6 GB
- **Savings:** ~12 GB (53% reduction)

---

## 🎯 Priority Files for GitHub

### Upload to GitHub (With Git LFS for large files)

1. **All scripts** (small, essential)
2. **README.md** (documentation)
3. **Metadata files** (small CSV files)
4. **Summary tables** (small CSV files)
5. **Environment files** (reproducibility)
6. **Plots** (as images in plots/)

### Do NOT Upload (Too Large)

1. **RDS files** - Use Git LFS or host elsewhere (Zenodo, Figshare)
2. **H5AD files** - Use Git LFS or host elsewhere
3. **Raw 10X data** - Reference GEO accession instead

### Use Git LFS For

```bash
# Install Git LFS
git lfs install

# Track large files
git lfs track "*.rds"
git lfs track "*.h5ad"
git lfs track "*.h5seurat"

# Commit .gitattributes
git add .gitattributes
git commit -m "Configure Git LFS"
```

---

## 📝 Naming Conventions Used

### File Prefixes (Sequential)
- `01_*` - Initial data preparation
- `02_*` - Data merging
- `03_*` - Normalization, PCA, UMAP
- `04_*` - Harmony integration (final)

### File Types
- `.rds` - R Seurat objects
- `.h5ad` - Python AnnData objects
- `.csv` - Metadata and summary tables
- `.Rmd` - R Markdown analysis scripts
- `.ipynb` - Jupyter notebooks
- `.png` - Figures and plots

---

## ✅ Final Checklist

- [ ] Create directory structure
- [ ] Move essential files to correct locations
- [ ] Delete redundant files
- [ ] Rename scripts with sequential prefixes
- [ ] Create environment files
- [ ] Generate .gitignore
- [ ] Write comprehensive README.md
- [ ] Add LICENSE file
- [ ] Set up Git LFS for large files
- [ ] Make first commit
- [ ] Push to GitHub

---

**Space Saved:** ~12 GB  
**Files Kept:** 11 essential files  
**Files Deleted:** 8 redundant files
