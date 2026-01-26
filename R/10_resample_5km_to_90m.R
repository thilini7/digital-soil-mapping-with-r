# ============================================================================
# Script: 10_resample_5km_to_90m.R
# Purpose: Convert 5km resolution rasters to 90m and align to NSW/ACT extent
# Author: Auto-generated
# Date: 2026-01-23
# ============================================================================

# Load required libraries
library(terra)

# ============================================================================
# 1. Define file paths
# ============================================================================

# Set working directory to project root
project_root <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
setwd(project_root)

# Input directory with 5km rasters
input_dir <- file.path(project_root, "Data/data_in/soil_cov_5km")

# Reference raster for alignment (90m resolution, NSW/ACT extent)
reference_raster_path <- file.path(project_root, "Annual_90th_Percen_Temp_NSW_ACT_90m.tif")

# Output directory for resampled 90m rasters
output_dir <- file.path(project_root, "Data/data_out/soil_cov_90m")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 2. Load reference raster
# ============================================================================

cat("Loading reference raster to extract target properties...\n")
ref_raster <- rast(reference_raster_path)

# Extract properties from reference raster
ref_crs <- crs(ref_raster)
ref_extent <- ext(ref_raster)
ref_res <- res(ref_raster)

cat("Reference raster properties:\n")
cat("  CRS:", as.character(ref_crs), "\n")
cat("  Extent:", as.vector(ref_extent), "\n")
cat("  Resolution:", ref_res, "\n")
cat("  Dimensions:", dim(ref_raster), "\n\n")

# ============================================================================
# 3. Get list of input rasters
# ============================================================================

# Find all .tif files in the input directory
input_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

cat("Found", length(input_files), "raster files to process:\n")
for (f in input_files) {
  cat("  -", basename(f), "\n")
}
cat("\n")

# ============================================================================
# 4. Function to resample a single raster
# ============================================================================

resample_to_90m <- function(input_path, ref_raster, output_dir, method = "bilinear") {
  #' Resample a raster to match reference raster properties
  #' 
  #' @param input_path Path to input raster file
  #' @param ref_raster Reference raster (terra SpatRaster)
  #' @param output_dir Directory for output file

  #' @param method Resampling method ("bilinear" for continuous, "near" for categorical)
  #' @return Path to output raster
  
  file_name <- basename(input_path)
  output_name <- sub("\\.tif$", "_NSW_ACT_90m.tif", file_name)
  output_path <- file.path(output_dir, output_name)
  
  cat("Processing:", file_name, "\n")
  cat("  Loading raster...\n")
  
  # Load input raster
  input_raster <- rast(input_path)
  
  cat("  Original properties:\n")
  cat("    Dimensions:", dim(input_raster), "\n")
  cat("    Resolution:", res(input_raster), "\n")
  cat("    Extent:", as.vector(ext(input_raster)), "\n")
  cat("    CRS:", as.character(crs(input_raster)), "\n")
  
  # Get original statistics
  orig_stats <- global(input_raster, fun = c("min", "max", "mean"), na.rm = TRUE)
  cat("    Value range:", round(orig_stats$min, 4), "to", round(orig_stats$max, 4), "\n")
  cat("    Mean:", round(orig_stats$mean, 4), "\n")
  
  # Check if CRS matches reference, if not reproject
  if (!same.crs(input_raster, ref_raster)) {
    cat("  Reprojecting to match reference CRS...\n")
    input_raster <- project(input_raster, crs(ref_raster))
  }
  
  # Crop to reference extent (with some buffer to ensure full coverage)
  cat("  Cropping to NSW/ACT extent...\n")
  input_cropped <- crop(input_raster, ext(ref_raster), extend = TRUE)
  
  # Resample to match reference raster resolution and alignment
  cat("  Resampling from 5km to 90m resolution (method:", method, ")...\n")
  cat("    This may take a while for large rasters...\n")
  
  resampled_raster <- resample(
    input_cropped,
    ref_raster,
    method = method
  )
  
  cat("  Resampled properties:\n")
  cat("    Dimensions:", dim(resampled_raster), "\n")
  cat("    Resolution:", res(resampled_raster), "\n")
  cat("    Extent:", as.vector(ext(resampled_raster)), "\n")
  
  # Apply mask from reference raster to match exact coverage
  cat("  Applying mask from reference raster...\n")
  ref_mask <- !is.na(ref_raster)
  resampled_masked <- mask(resampled_raster, ref_mask, maskvalues = FALSE)
  
  # Set layer name (clean up the name)
  layer_name <- tools::file_path_sans_ext(output_name)
  names(resampled_masked) <- layer_name
  
  # Print final statistics
  final_stats <- global(resampled_masked, fun = c("min", "max", "mean", "sd"), na.rm = TRUE)
  cat("  Final statistics:\n")
  cat("    Min:", round(final_stats$min, 4), "\n")
  cat("    Max:", round(final_stats$max, 4), "\n")
  cat("    Mean:", round(final_stats$mean, 4), "\n")
  cat("    SD:", round(final_stats$sd, 4), "\n")
  
  # Count valid cells
  n_valid <- global(resampled_masked, fun = "notNA")
  n_total <- prod(dim(resampled_masked)[1:2])
  cat("    Valid cells:", n_valid$notNA, "of", n_total, "\n")
  
  # Write output raster
  cat("  Writing output to:", output_path, "\n")
  
  writeRaster(
    resampled_masked,
    filename = output_path,
    overwrite = TRUE,
    datatype = "FLT4S",  # 32-bit float
    gdal = c("COMPRESS=LZW", "TILED=YES")
  )
  
  cat("  Done!\n\n")
  
  return(output_path)
}

# ============================================================================
# 5. Process all rasters
# ============================================================================

cat("=" , rep("=", 70), "\n", sep = "")
cat("Starting batch resampling process...\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

# Track processing time
start_time <- Sys.time()

# Store output paths
output_paths <- character(length(input_files))

# Process each raster
for (i in seq_along(input_files)) {
  cat("Processing file", i, "of", length(input_files), "\n")
  cat("-", rep("-", 50), "\n", sep = "")
  
  output_paths[i] <- resample_to_90m(
    input_path = input_files[i],
    ref_raster = ref_raster,
    output_dir = output_dir,
    method = "bilinear"  # Use bilinear for continuous data
  )
}

end_time <- Sys.time()
processing_time <- difftime(end_time, start_time, units = "mins")

# ============================================================================
# 6. Summary and verification
# ============================================================================

cat("=" , rep("=", 70), "\n", sep = "")
cat("PROCESSING COMPLETE\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

cat("Total processing time:", round(processing_time, 2), "minutes\n\n")

cat("Output files created:\n")
for (op in output_paths) {
  cat("  -", basename(op), "\n")
}

cat("\nOutput directory:", output_dir, "\n\n")

# Verify all outputs
cat("Verifying output rasters...\n")
cat("-", rep("-", 50), "\n", sep = "")

for (op in output_paths) {
  output_check <- rast(op)
  cat("\n", basename(op), ":\n", sep = "")
  cat("  Dimensions:", dim(output_check), "\n")
  cat("  Resolution:", res(output_check), "\n")
  cat("  Matches reference dimensions:", all(dim(ref_raster) == dim(output_check)), "\n")
}

cat("\n=== Script completed successfully ===\n")
