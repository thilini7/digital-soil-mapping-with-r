# ============================================================================
# Script: 8_vector_to_raster_ldi.R
# Purpose: Convert vector shapefile (NSW_ACT_Landuse_merged.shp) to raster TIFF
#          using LDI values at 90m resolution
# Author: Auto-generated
# Date: 2026-01-19
# ============================================================================

# Load required libraries
library(sf)
library(terra)
library(dplyr)

# ============================================================================
# 1. Define file paths
# ============================================================================

# Set working directory to project root
project_root <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
setwd(project_root)

# Input files
shapefile_path <- file.path(project_root, "Data/data_in/shp_files/NSW_ACT_Landuse_merged.shp")
reference_raster_path <- file.path(project_root, "Annual_90th_Percen_Temp_NSW_ACT_90m.tif")

# Output file
output_raster_path <- file.path(project_root, "Data/data_out/LDI_NSW_ACT_90m.tif")

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
# 3. Load and prepare vector data
# ============================================================================

cat("\nLoading shapefile...\n")
land_use <- st_read(shapefile_path, quiet = TRUE)

cat("Shapefile loaded with", nrow(land_use), "features\n")
cat("Available columns:", paste(names(land_use), collapse = ", "), "\n")

# Check if LDI column exists
if (!"LDI" %in% names(land_use)) {
  stop("LDI column not found in shapefile. Please run 7_merge_nsw_act_shapefiles.R first.")
}

# Check LDI values
cat("\nLDI value summary:\n")
print(table(land_use$LDI, useNA = "always"))

# Transform to match reference CRS if needed
if (st_crs(land_use) != st_crs(4326)) {
  cat("\nTransforming shapefile to WGS 84 (EPSG:4326)...\n")
  land_use <- st_transform(land_use, crs = 4326)
}

# ============================================================================
# 4. Create empty raster template matching reference
# ============================================================================

cat("\nCreating raster template matching reference...\n")

# Create template raster with same properties as reference
template_raster <- rast(
  extent = ref_extent,
  resolution = ref_res,
  crs = ref_crs
)

cat("Template raster dimensions:", dim(template_raster), "\n")

# ============================================================================
# 5. Rasterize vector using LDI values
# ============================================================================

cat("\nRasterizing vector data (this may take several minutes)...\n")
start_time <- Sys.time()

# Convert sf to SpatVector for terra
land_use_vect <- vect(land_use)

# Rasterize using LDI field
ldi_raster <- rasterize(
  x = land_use_vect,
  y = template_raster,
  field = "LDI",
  fun = "min",  # If multiple polygons overlap, take minimum LDI (least disturbed)
  background = NA
)

end_time <- Sys.time()
cat("Rasterization completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# ============================================================================
# 6. Set raster metadata
# ============================================================================

# Set layer name
names(ldi_raster) <- "LDI_NSW_ACT_90m"

# Print raster statistics
cat("\nOutput raster summary:\n")
cat("  Dimensions:", dim(ldi_raster), "\n")
cat("  Resolution:", res(ldi_raster), "\n")
cat("  Extent:", as.vector(ext(ldi_raster)), "\n")
cat("  CRS:", as.character(crs(ldi_raster)), "\n")

# Value statistics
ldi_values <- values(ldi_raster, na.rm = TRUE)
cat("\nLDI raster value statistics:\n")
cat("  Min:", min(ldi_values, na.rm = TRUE), "\n")
cat("  Max:", max(ldi_values, na.rm = TRUE), "\n")
cat("  Mean:", round(mean(ldi_values, na.rm = TRUE), 2), "\n")
cat("  Non-NA cells:", sum(!is.na(ldi_values)), "\n")
cat("  NA cells:", sum(is.na(values(ldi_raster))), "\n")

print(table(ldi_values, useNA = "always"))

# ============================================================================
# 7. Write output raster
# ============================================================================

cat("\nWriting output raster to:", output_raster_path, "\n")

writeRaster(
  ldi_raster,
  filename = output_raster_path,
  overwrite = TRUE,
  datatype = "INT1U",  # Unsigned 8-bit integer (LDI values 1-6)
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

cat("\n=== Script completed successfully ===\n")
