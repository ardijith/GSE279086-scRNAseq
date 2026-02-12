# ============================================
# 04_cleanup_files.R
# Clean up GEO downloaded files for Seurat
# ============================================

# Set the base directory (UPDATE THIS PATH!)
base_dir <- "C:/ardij/Documents/GSE279086/raw_geo_data"

# Check if directory exists
if (!dir.exists(base_dir)) {
  stop("Directory not found: ", base_dir)
}

cat("Cleaning GEO files in:", base_dir, "\n\n")

# Get all GSM folders
gsm_folders <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
gsm_folders <- gsm_folders[grepl("GSM", basename(gsm_folders))]

cat("Found", length(gsm_folders), "GSM folders\n")

# Process each GSM folder
for (gsm_folder in gsm_folders) {
  cat("Processing", basename(gsm_folder), "...\n")
  
  # Set working directory to this GSM folder
  setwd(gsm_folder)
  
  # 1. Delete the processed files
  processed_files <- c(
    "*_barcodes_processed.tsv.gz",
    "*_features_processed.tsv.gz", 
    "*_matrix_processed.mtx.gz"
  )
  
  for (pattern in processed_files) {
    files_to_delete <- list.files(pattern = pattern, full.names = TRUE)
    if (length(files_to_delete) > 0) {
      file.remove(files_to_delete)
      cat("   Deleted", length(files_to_delete), "processed file(s)\n")
    }
  }
  
  # 2. Rename files to standard names
  # List of patterns and their new names
  rename_patterns <- list(
    c("*_barcodes.tsv.gz", "barcodes.tsv.gz"),
    c("*_features.tsv.gz", "features.tsv.gz"),
    c("*_matrix.mtx.gz", "matrix.mtx.gz")
  )
  
  for (pattern_new in rename_patterns) {
    pattern <- pattern_new[1]
    new_name <- pattern_new[2]
    
    matching_files <- list.files(pattern = pattern, full.names = TRUE)
    
    if (length(matching_files) == 1) {
      # Rename the file
      file.rename(matching_files, new_name)
      cat("   Renamed", basename(matching_files), "->", new_name, "\n")
    } else if (length(matching_files) > 1) {
      cat("   Warning: Multiple files match pattern", pattern, "\n")
    }
  }
  
  cat("   Cleaned up", basename(gsm_folder), "\n")
}

cat("\nCleanup complete! Each folder now contains:\n")
cat("barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz\n")






# ============================================
# 
#     Convert 10x into  Seurat
# ============================================



``{r Load the packages}
library(Seurat)
library(dplyr)
library(tidyr)
```


```{r Prepare the seurat object}
data_dir <- "~/GSE279086/raw_geo_data"

# List all GSM subdirectories
samples <- list.dirs(data_dir, recursive = FALSE)

seurat_list <- lapply(samples, function(sample_path) {
  sample_name <- basename(sample_path)
  message("Processing: ", sample_name)
  
  # Load data
  counts <- Read10X(data.dir = sample_path)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name,min.cells = 3,min.features = 200, assay = "RNA"
  )
  
  # Add metadata
  seurat_obj$sample <- sample_name
  
  # Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]|RPLP")
  
  # Handle any missing values (just in case)
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(
      percent.mt = replace_na(percent.mt, 0),
      percent.rb = replace_na(percent.rb, 0)
    )
  
  message(paste0("NaNs remaining in 'percent.mt': ", sum(is.nan(seurat_obj$percent.mt))))
  message(paste0("NaNs remaining in 'percent.rb': ", sum(is.nan(seurat_obj$percent.rb))))
  
  return(seurat_obj)
})

names(seurat_list) <- basename(samples)
names(seurat_list)
```


```{r Save individual samples into a seurat object}
output_dir <- file.path("~/GSE279086/input")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save each Seurat object
for (i in seq_along(seurat_list)) {
  sample_name <- seurat_list[[i]]@project.name
  save_path <- file.path(output_dir, paste0(sample_name, "_seurat.rds"))
  
  message("Saving: ", sample_name, " -> ", save_path)
  saveRDS(seurat_list[[i]], file = save_path)
}

message("All Seurat objects saved successfully in: ", output_dir)
```



