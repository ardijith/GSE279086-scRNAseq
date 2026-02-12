# GitHub Setup Guide for GSE279086 Project

## 🚀 Quick Start - Push to GitHub

### Prerequisites
- Git installed on your computer
- GitHub account created
- Git configured with your credentials

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

---

## Step-by-Step Guide

### 1. Organize Your Files

```bash
cd ~/GSE279086

# Run the organization script
bash organize_project.sh

# Review files to delete
cat files_to_delete.txt

# Delete redundant files (optional but recommended)
bash delete_redundant_files.sh
```

### 2. Initialize Git Repository

```bash
# Initialize git
git init

# Add essential files
cp /path/to/README.md ./
cp /path/to/LICENSE ./
cp /path/to/.gitignore ./
```

### 3. Set Up Git LFS (for large files)

Git LFS allows you to version control large files without bloating your repository.

```bash
# Install Git LFS (if not already installed)
# On Mac: brew install git-lfs
# On Ubuntu: sudo apt-get install git-lfs
# On Windows: Download from https://git-lfs.github.com/

# Initialize Git LFS
git lfs install

# Track large files
git lfs track "*.rds"
git lfs track "*.h5ad"
git lfs track "*.h5seurat"
git lfs track "outputs/rds/*"
git lfs track "outputs/h5ad/*"

# Add .gitattributes (created by git lfs track)
git add .gitattributes
```

### 4. Create Repository on GitHub

1. Go to GitHub.com and log in
2. Click the **"+"** icon → **"New repository"**
3. Repository name: `GSE279086-T1D-scRNAseq`
4. Description: `Single-cell RNA-seq analysis of Type 1 Diabetes pancreatic islets`
5. Choose **Public** (for portfolio) or **Private**
6. **DO NOT** initialize with README (you already have one)
7. Click **"Create repository"**

### 5. Connect Local Repository to GitHub

```bash
# Add remote repository
git remote add origin https://github.com/YOUR_USERNAME/GSE279086-T1D-scRNAseq.git

# Verify remote
git remote -v
```

### 6. Make Initial Commit

```bash
# Check status
git status

# Add all files
git add .

# Check what will be committed
git status

# Commit with message
git commit -m "Initial commit: Complete scRNA-seq analysis pipeline"

# Push to GitHub
git branch -M main
git push -u origin main
```

---

## 📦 Alternative: Hosting Large Files Externally

If Git LFS is too complex or expensive, consider these alternatives:

### Option 1: Zenodo (Recommended for Data)

1. Go to [zenodo.org](https://zenodo.org)
2. Sign in with GitHub account
3. Upload your RDS and h5ad files
4. Get DOI (Digital Object Identifier)
5. Add DOI link to your README

**In README.md:**
```markdown
## Data Availability

Large processed data files are available at Zenodo:
- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

### Option 2: Figshare

Similar to Zenodo, but with different interface.

### Option 3: Google Drive (Quick & Easy)

1. Upload files to Google Drive
2. Make shareable
3. Add links to README

**In README.md:**
```markdown
## Data Files

Download processed data from Google Drive:
- [02_seurat_merged.rds](https://drive.google.com/file/d/XXXXX)
- [03_seurat_pca_umap_JOINED.rds](https://drive.google.com/file/d/XXXXX)
- [04_harmony_integrated.rds](https://drive.google.com/file/d/XXXXX)
- [04_harmony_integrated.h5ad](https://drive.google.com/file/d/XXXXX)
```

---

## 🔒 What to Include in GitHub

### ✅ DO Include

- [x] Scripts (`.Rmd`, `.ipynb`, `.R`, `.py`)
- [x] README.md
- [x] LICENSE
- [x] .gitignore
- [x] Metadata files (small CSVs)
- [x] Summary tables (small CSVs)
- [x] Plots (PNGs, JPGs)
- [x] Environment files (`requirements.txt`, etc.)
- [x] Documentation (methods, file descriptions)

### ❌ DON'T Include (Add to .gitignore)

- [ ] Large RDS files (use Git LFS or external hosting)
- [ ] Large h5ad files (use Git LFS or external hosting)
- [ ] Raw 10X data (reference GEO accession instead)
- [ ] Temporary files
- [ ] System files (.DS_Store, Thumbs.db)

---

## 📝 Update README with Data Links

After uploading large files to external hosting:

```bash
# Edit README.md to add download links
nano README.md

# Add section like:
## Data Availability

### Processed Data Files

Download from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Or individual files:
- `02_seurat_merged.rds` - [Download](https://zenodo.org/record/XXXXX/files/02_seurat_merged.rds)
- `03_seurat_pca_umap_JOINED.rds` - [Download](https://zenodo.org/record/XXXXX/files/03_seurat_pca_umap_JOINED.rds)
- `04_harmony_integrated.rds` - [Download](https://zenodo.org/record/XXXXX/files/04_harmony_integrated.rds)
- `04_harmony_integrated.h5ad` - [Download](https://zenodo.org/record/XXXXX/files/04_harmony_integrated.h5ad)

### Raw Data

Original 10X Genomics data available at GEO: [GSE279086](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279086)

# Commit changes
git add README.md
git commit -m "Add data availability links"
git push
```

---

## 🎨 Customize Your Repository

### Add Topics

On GitHub repository page:
1. Click **"Add topics"**
2. Add: `single-cell`, `rna-seq`, `type-1-diabetes`, `seurat`, `harmony`, `celltypist`, `bioinformatics`, `r`, `python`

### Add Description

Repository description: `Single-cell RNA-seq analysis pipeline for Type 1 Diabetes pancreatic islets using Seurat, Harmony, and CellTypist`

### Add Website

If you have a portfolio website, add it in the repository settings.

### Create Releases

After completing major milestones:

```bash
# Tag a release
git tag -a v1.0.0 -m "Complete analysis pipeline - initial release"
git push origin v1.0.0
```

Then create a formal release on GitHub with release notes.

---

## 🔧 Maintaining Your Repository

### Updating After Changes

```bash
# Check what changed
git status

# Add changed files
git add scripts/01_data_preparation.Rmd

# Or add all changes
git add .

# Commit with descriptive message
git commit -m "Update normalization parameters in data preparation"

# Push to GitHub
git push
```

### Creating Branches for New Features

```bash
# Create and switch to new branch
git checkout -b add-differential-expression

# Make changes, commit
git add .
git commit -m "Add differential expression analysis"

# Push branch
git push -u origin add-differential-expression

# On GitHub, create Pull Request to merge into main
```

---

## 📊 Adding a GitHub Pages Site (Optional)

Turn your README into a website:

1. Go to repository **Settings**
2. Scroll to **GitHub Pages**
3. Source: **main branch** → **/docs** (if you have docs folder)
4. Your site will be at: `https://YOUR_USERNAME.github.io/GSE279086-T1D-scRNAseq/`

---

## ✅ Final Checklist

Before making repository public:

- [ ] All scripts run without errors
- [ ] README is complete and informative
- [ ] Large files are handled (LFS or external)
- [ ] Sensitive information removed (passwords, API keys)
- [ ] License file included
- [ ] .gitignore properly configured
- [ ] Repository description and topics added
- [ ] Contact information updated in README
- [ ] Data availability clearly documented

---

## 🎯 Portfolio Tips

To make this repository stand out in your portfolio:

1. **Clear README:** Use the comprehensive README provided
2. **Professional plots:** Include high-quality figures
3. **Complete documentation:** Methods, file descriptions
4. **Reproducible:** Environment files, version info
5. **Clean code:** Well-commented scripts
6. **Results:** Include key findings in README
7. **Citation:** Proper attribution to tools and data
8. **Contact:** Make it easy for recruiters to reach you

---

## 🆘 Troubleshooting

### Push Rejected

```bash
# If remote has changes you don't have
git pull --rebase origin main
git push
```

### File Too Large Error

```bash
# Add to Git LFS
git lfs track "path/to/large/file.rds"
git add .gitattributes
git add path/to/large/file.rds
git commit --amend --no-edit
git push --force
```

### Remove Accidentally Committed File

```bash
# Remove from git but keep locally
git rm --cached path/to/file
echo "path/to/file" >> .gitignore
git add .gitignore
git commit -m "Remove large file from tracking"
git push
```

---

**Your repository is now ready for your portfolio!** 🎉
