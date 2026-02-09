# =============================================================================
# Align All Covariate Rasters to Common Extent (Minimum NA Coverage)
# =============================================================================
# This script resamples rasters with extra NA pixels to match the valid extent
# of the reference raster (the one with minimum NAs).
# =============================================================================

library(terra)

# --- Configuration ---
HomeDir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"

input_folder <- file.path(HomeDir, "Data/data_in/soil_covariates_aligned")
output_folder <- file.path(HomeDir, "Data/data_in/soil_covariates_aligned_v2")

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# --- Load all rasters and check NA counts ---
cov_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)
cat("Found", length(cov_files), "raster files\n")

# Calculate NA count for each raster
na_counts <- sapply(cov_files, function(f) {
  r <- rast(f)
  global(is.na(r), sum)$sum
})

names(na_counts) <- basename(cov_files)
na_df <- data.frame(file = names(na_counts), na_count = na_counts, row.names = NULL)
na_df <- na_df[order(na_df$na_count), ]

cat("\n--- NA counts per raster (sorted) ---\n")
print(na_df)

# --- Find reference raster (minimum NAs) ---
min_na <- min(na_counts)
ref_file <- cov_files[which.min(na_counts)]
ref_name <- basename(ref_file)

cat("\n*** Reference raster (fewest NAs):", ref_name, "\n")
cat("    NA count:", format(min_na, big.mark = ","), "\n")

ref_rast <- rast(ref_file)

# Create a mask of valid pixels from reference
ref_mask <- !is.na(ref_rast)

# --- Identify rasters that need alignment (>5% more NAs than reference) ---
threshold <- min_na * 1.05
problem_files <- cov_files[na_counts > threshold]

cat("\n--- Rasters needing alignment (>5% more NAs than reference) ---\n")
for (f in problem_files) {
  extra_na <- na_counts[basename(f)] - min_na
  cat("  ", basename(f), ": +", format(extra_na, big.mark = ","), "extra NAs\n")
}

# --- Process all rasters ---
cat("\n=== Processing rasters ===\n")

for (f in cov_files) {
  fname <- basename(f)
  out_file <- file.path(output_folder, fname)
  
  r <- rast(f)
  
  if (f %in% problem_files) {
    cat("Aligning:", fname, "...")
    
    # Check if extents match
    if (!compareGeom(r, ref_rast, stopOnError = FALSE)) {
      # Resample to reference grid
      r_aligned <- resample(r, ref_rast, method = "bilinear")
    } else {
      r_aligned <- r
    }
    
    # Apply the reference mask (set extra NAs to NA in output)
    # This ensures all rasters have the same valid pixel extent
    r_aligned <- mask(r_aligned, ref_mask, maskvalues = 0)
    
    writeRaster(r_aligned, out_file, overwrite = TRUE)
    
    # Verify
    new_na <- global(is.na(r_aligned), sum)$sum
    cat(" done. NAs:", format(min_na, big.mark = ","), "->", format(new_na, big.mark = ","), "\n")
    
  } else {
    # Just copy the raster (already aligned)
    cat("Copying:", fname, "(already aligned)\n")
    writeRaster(r, out_file, overwrite = TRUE)
  }
}

# --- Final verification ---
cat("\n=== Verifying aligned rasters ===\n")

aligned_files <- list.files(output_folder, pattern = "\\.tif$", full.names = TRUE)
aligned_na <- sapply(aligned_files, function(f) {
  global(is.na(rast(f)), sum)$sum
})

names(aligned_na) <- basename(aligned_files)
unique_na <- unique(aligned_na)

if (length(unique_na) == 1) {
  cat("✓ SUCCESS: All rasters now have the same NA extent (", format(unique_na, big.mark = ","), " NAs)\n")
} else {
  cat("WARNING: NA counts still differ:\n")
  print(table(aligned_na))
}

cat("\nAligned rasters saved to:", output_folder, "\n")
cat("\nNext step: Update cov_folder in 4_1_prediction_map.R to use:\n")
cat('  cov_folder <- "', output_folder, '"\n', sep = "")
