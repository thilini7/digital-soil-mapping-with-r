# ============================================================================
# Script: 9_ascii_to_geotiff.R
# Purpose: Convert ASCII grid file to GeoTIFF using reference raster properties
# Author: Auto-generated
# Date: 2026-01-19
# ============================================================================

# Load required libraries
library(terra)

# ============================================================================
# 1. Define file paths
# ============================================================================

# Set working directory to project root
project_root <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
setwd(project_root)

# Input files
ascii_file <- file.path(project_root, "Data/data_in/asci_files/Temperature/Annual Max temp/mxtan.txt")
reference_raster_path <- file.path(project_root, "Annual_90th_Percen_Temp_NSW_ACT_90m.tif")

# Output file
output_raster_path <- file.path(project_root, "Data/data_out/Annual_Max_Temp_NSW_ACT_90m.tif")

# Create output directory if it doesn't exist
dir.create(dirname(output_raster_path), showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 2. Load reference raster to get target properties
# ============================================================================

cat("Loading reference raster to extract properties...\n")
ref_raster <- rast(reference_raster_path)

# Extract properties from reference raster
ref_crs <- crs(ref_raster)
ref_extent <- ext(ref_raster)
ref_res <- res(ref_raster)

cat("Reference raster properties:\n")
cat("  CRS:", as.character(ref_crs), "\n")
cat("  Extent:", as.vector(ref_extent), "\n")
cat("  Resolution:", ref_res, "\n")
cat("  Dimensions:", dim(ref_raster), "\n")

# ============================================================================
# 3. Load and process ASCII file
# ============================================================================

cat("\nLoading ASCII grid file...\n")

# Read ASCII grid - terra can read ESRI ASCII grid format directly
ascii_raster <- rast(ascii_file)

cat("Original ASCII raster properties:\n")
cat("  Dimensions:", dim(ascii_raster), "\n")
cat("  Resolution:", res(ascii_raster), "\n")
cat("  Extent:", as.vector(ext(ascii_raster)), "\n")

# Set CRS (ASCII files typically don't include CRS info, assume WGS84)
crs(ascii_raster) <- "EPSG:4326"

# Handle NoData values (99999.90 from the header)
# terra should automatically read nodata_value from header, but let's ensure
NAflag(ascii_raster) <- 99999.90

cat("\nOriginal raster value range (excluding NA):\n")
global_stats <- global(ascii_raster, fun = c("min", "max"), na.rm = TRUE)
cat("  Min:", global_stats$min, "\n")
cat("  Max:", global_stats$max, "\n")

# ============================================================================
# 4. Crop and resample to match reference raster
# ============================================================================

cat("\nCropping to NSW/ACT extent...\n")

# Crop to reference extent
ascii_cropped <- crop(ascii_raster, ref_extent)

cat("Cropped raster properties:\n")
cat("  Dimensions:", dim(ascii_cropped), "\n")
cat("  Extent:", as.vector(ext(ascii_cropped)), "\n")

cat("\nResampling to match reference raster resolution (90m)...\n")

# Resample to match reference raster resolution
# Using bilinear interpolation for continuous temperature data
ascii_resampled <- resample(
  ascii_cropped, 
  ref_raster, 
  method = "bilinear"
)

cat("Resampled raster properties:\n")
cat("  Dimensions:", dim(ascii_resampled), "\n")
cat("  Resolution:", res(ascii_resampled), "\n")
cat("  Extent:", as.vector(ext(ascii_resampled)), "\n")

# ============================================================================
# 5. Apply mask to match reference (optional - to match exact coverage)
# ============================================================================

cat("\nApplying mask from reference raster...\n")

# Create a mask from reference raster (where reference has valid values)
# This ensures the output has the same coverage as reference
ref_mask <- !is.na(ref_raster)
ascii_masked <- mask(ascii_resampled, ref_mask, maskvalues = FALSE)

# ============================================================================
# 6. Set layer name and metadata
# ============================================================================

names(ascii_masked) <- "Annual_Max_Temp_NSW_ACT_90m"

# Print final statistics
cat("\nFinal raster statistics:\n")
final_stats <- global(ascii_masked, fun = c("min", "max", "mean", "sd"), na.rm = TRUE)
cat("  Min:", round(final_stats$min, 2), "°C\n")
cat("  Max:", round(final_stats$max, 2), "°C\n")
cat("  Mean:", round(final_stats$mean, 2), "°C\n")
cat("  SD:", round(final_stats$sd, 2), "°C\n")

# Count cells
n_valid <- global(ascii_masked, fun = "notNA")
n_total <- prod(dim(ascii_masked)[1:2])
cat("  Valid cells:", n_valid$notNA, "of", n_total, "\n")

# ============================================================================
# 7. Write output raster
# ============================================================================

cat("\nWriting output raster to:", output_raster_path, "\n")

writeRaster(
  ascii_masked,
  filename = output_raster_path,
  overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float for temperature data
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

cat("\nRaster successfully saved!\n")

# ============================================================================
# 8. Verify output
# ============================================================================

cat("\nVerifying output raster...\n")
output_check <- rast(output_raster_path)
cat("Output raster info:\n")
print(output_check)

# Compare with reference
cat("\nComparison with reference raster:\n")
cat("  Reference dimensions:", dim(ref_raster), "\n")
cat("  Output dimensions:", dim(output_check), "\n")
cat("  Dimensions match:", all(dim(ref_raster) == dim(output_check)), "\n")

cat("\n=== Script completed successfully ===\n")
