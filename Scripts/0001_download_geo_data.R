################################################################################
# DOWNLOAD scRNA-seq DATA FROM GEO (GSE279086)
# Customized for: C:/ardij/Documents
# Platform: Windows/Mac/Linux compatible
################################################################################

# === Configuration ===
cat("=== GEO Data Download Script ===\n\n")

# Mapping of GSM IDs to their sample codes
samples <- list(
    GSM8561110 = "S_1907_004567",
    GSM8561111 = "S_1907_004614",
    GSM8561112 = "S_1907_004802",
    GSM8561113 = "S_1907_004896",
    GSM8561114 = "S_1907_005225",
    GSM8561115 = "S_2006_004078",
    GSM8561116 = "S_2007_002809",
    GSM8561117 = "S_2007_002950",
    GSM8561118 = "S_2007_002997",
    GSM8561119 = "S_2007_003044",
    GSM8561120 = "S_2007_003091",
    GSM8561121 = "S_2007_003138",
    GSM8561122 = "S_2007_003793",
    GSM8561123 = "S_2007_003840",
    GSM8561124 = "S_2007_003934",
    GSM8561125 = "S_2007_004028",
    GSM8561126 = "S_2007_004216",
    GSM8561127 = "S_2007_004263",
    GSM8561128 = "S_2103_004019",
    GSM8561129 = "S_2103_004028",
    GSM8561130 = "S_2103_004037",
    GSM8561131 = "S_2103_004046",
    GSM8561132 = "S_2103_004064",
    GSM8561133 = "S_2103_004073",
    GSM8561134 = "S_2103_004091",
    GSM8561135 = "S_2103_004100",
    GSM8561136 = "S_2103_004109",
    GSM8561137 = "S_2103_004118",
    GSM8561138 = "S_2103_004127",
    GSM8561139 = "S_2103_004145",
    GSM8561140 = "S_2103_004154",
    GSM8561141 = "S_2103_004163",
    GSM8561142 = "S_2103_004181",
    GSM8561143 = "S_2107_023520",
    GSM8561144 = "S_2107_023538",
    GSM8561145 = "S_2107_023547",
    GSM8561146 = "S_2107_023556",
    GSM8561147 = "S_2107_023592",
    GSM8561148 = "S_2107_023610",
    GSM8561149 = "S_2107_023619"
)

# Base URL for GEO sample supplemental files
base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/samples"

# Files to download per GSM sample
files_to_download <- c(
    "barcodes.tsv.gz",
    "barcodes_processed.tsv.gz",
    "features.tsv.gz",
    "features_processed.tsv.gz",
    "matrix.mtx.gz",
    "matrix_processed.mtx.gz"
)

# === SET YOUR DOWNLOAD DIRECTORY HERE ===
# 🔽 YOUR CUSTOM DIRECTORY - Downloads will go to:
# C:/ardij/Documents/GSE279086/raw_geo_data

base_dir <- "C:/ardij/Documents/GSE279086/raw_geo_data"

# ⚠️ IMPORTANT: Make sure C:/ardij/Documents exists on your computer!
# If "ardij" folder doesn't exist, the script will try to create it

cat("Download directory:", base_dir, "\n\n")

# === Verify Directory Creation ===
# Try to create the directory
dir_created <- dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

# Check if directory exists now
if (!dir.exists(base_dir)) {
  cat("❌ ERROR: Cannot create directory!\n")
  cat("Please check:\n")
  cat("  1. Does C:/ardij exist?\n")
  cat("  2. Do you have write permission?\n")
  cat("  3. Is the path spelled correctly?\n\n")
  stop("Directory creation failed. Please fix the path and try again.")
} else {
  cat("✓ Directory verified and ready!\n\n")
}

# === Download Function ===
download_geo_file <- function(url, dest_file) {
  # Try to download, with error handling
  tryCatch({
    download.file(url, dest_file, mode = "wb", quiet = TRUE)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# === Download Process ===
cat("Starting download of", length(samples), "samples...\n")
cat("This may take 30-60 minutes depending on your internet speed.\n\n")

# Ask user to confirm before starting
cat("Press [Enter] to start downloading, or [Esc] to cancel: ")
readline()

total_files <- 0
downloaded_files <- 0
failed_files <- 0
skipped_files <- 0

# Record start time
start_time <- Sys.time()

for (gsm in names(samples)) {
  sample_code <- samples[[gsm]]
  
  # Create subdirectory for this GSM
  # e.g., GSM8561nnn
  subdir <- paste0(substr(gsm, 1, 7), "nnn")
  
  # Destination directory for this sample
  dest_dir <- file.path(base_dir, gsm)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("---------------------------------------------------\n")
  cat("Downloading:", gsm, "(", sample_code, ")\n")
  cat("Progress:", which(names(samples) == gsm), "of", length(samples), "\n\n")
  
  for (file in files_to_download) {
    total_files <- total_files + 1
    
    # Construct filename and URL
    filename <- paste0(gsm, "_", sample_code, "_", file)
    dest_file <- file.path(dest_dir, filename)
    url <- paste0(base_url, "/", subdir, "/", gsm, "/suppl/", filename)
    
    # Check if file already exists
    if (file.exists(dest_file)) {
      cat("  ✓ Already exists:", filename, "\n")
      skipped_files <- skipped_files + 1
      next
    }
    
    # Download file
    cat("  → Downloading:", filename, "... ")
    success <- download_geo_file(url, dest_file)
    
    if (success) {
      # Get file size
      size_mb <- round(file.size(dest_file) / 1024^2, 2)
      cat("✓ Done (", size_mb, "MB )\n", sep = "")
      downloaded_files <- downloaded_files + 1
    } else {
      cat("✗ FAILED\n")
      failed_files <- failed_files + 1
    }
  }
  
  cat("\n")
}

# Calculate elapsed time
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

# === Summary ===
cat("=======================================================\n")
cat("DOWNLOAD COMPLETE!\n")
cat("=======================================================\n\n")
cat("Summary:\n")
cat("  Total files:      ", total_files, "\n")
cat("  Downloaded:       ", downloaded_files, "\n")
cat("  Already existed:  ", skipped_files, "\n")
cat("  Failed:           ", failed_files, "\n")
cat("  Time taken:       ", round(elapsed_time, 1), "minutes\n\n")

if (failed_files > 0) {
  cat("⚠️ Some files failed to download. This could be due to:\n")
  cat("  - Internet connection issues\n")
  cat("  - Files not yet available on GEO server\n")
  cat("  - Incorrect URLs\n")
  cat("\n💡 TIP: You can re-run this script to retry failed downloads.\n")
  cat("         Already downloaded files will be skipped.\n\n")
}

cat("📁 Files saved to:", base_dir, "\n\n")

cat("🎉 Next steps:\n")
cat("1. Navigate to:", base_dir, "\n")
cat("2. Verify you have folders GSM8561110 through GSM8561149\n")
cat("3. Each folder should contain 6 .gz files\n")
cat("4. Use Seurat's Read10X() to load the data\n\n")

cat("Example code to load first sample:\n")
cat("─────────────────────────────────────────────────────────\n")
cat("library(Seurat)\n")
cat("data_path <- '", file.path(base_dir, "GSM8561110"), "'\n", sep = "")
cat("data <- Read10X(data.dir = data_path)\n")
cat("seurat_obj <- CreateSeuratObject(counts = data)\n")
cat("─────────────────────────────────────────────────────────\n\n")

################################################################################
# WHAT YOU DOWNLOADED
################################################################################

# This dataset is GSE279086 from GEO (Gene Expression Omnibus)
# It contains single-cell RNA-seq data from 40 samples
# 
# Each sample has 6 files:
# 1. barcodes.tsv.gz          - Cell barcodes (identifies each cell)
# 2. barcodes_processed.tsv.gz - Filtered/processed cell barcodes
# 3. features.tsv.gz           - Gene names and IDs
# 4. features_processed.tsv.gz - Filtered/processed gene list  
# 5. matrix.mtx.gz             - Count matrix (genes × cells)
# 6. matrix_processed.mtx.gz   - Filtered/processed count matrix
#
# Usually you'll want to use the "processed" versions for analysis

################################################################################
# NEXT SCRIPT TO RUN
################################################################################

# After downloading, you'll want to:
# 1. Load the data into Seurat
# 2. Perform quality control
# 3. Normalize and analyze
#
# Would you like me to create a script for that? 😊

################################################################################
