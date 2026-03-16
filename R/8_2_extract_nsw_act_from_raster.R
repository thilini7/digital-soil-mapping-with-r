# ============================================================================
# Script: 8_2_extract_nsw_act_from_raster.R
# Purpose: Extract and align NSW + ACT raster from an Australia-wide raster
#          using the PM reference raster grid
# ============================================================================

# Load required libraries
library(terra)

# ============================================================================
# 1. Define file paths
# ============================================================================

input_raster_path <- file.path(
  "Data", "data_in", "raster_australia",
  "NLUM_v7_250_ALUMV8_2015_16_alb.tif"
)
reference_raster_path <- file.path(
  "Data", "data_in", "soil_covariates_for_geno_pheno_soil_aligned",
  "PM_radmap_v4_2019_filtered_dose_GAPFilled.tif"
)
output_raster_path <- file.path(
  "Data", "data_out", "soil_cov_90m",
  "NLUM_v7_250_ALUMV8_2015_16_alb_NSW_ACT_aligned.tif"
)

if (!file.exists(input_raster_path)) {
  stop("Input raster not found: ", input_raster_path)
}
if (!file.exists(reference_raster_path)) {
  stop("Reference raster not found: ", reference_raster_path)
}

dir.create(dirname(output_raster_path), showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 2. Load rasters
# ============================================================================

cat("Loading Australia-wide raster...\n")
input_raster <- rast(input_raster_path)

cat("Loading PM reference raster...\n")
ref_raster <- rast(reference_raster_path)

cat("\nInput raster info:\n")
cat("  CRS        : ", as.character(crs(input_raster)), "\n", sep = "")
cat("  Resolution : ", paste(res(input_raster), collapse = " x "), "\n", sep = "")
cat("  Extent     : ", paste(as.vector(ext(input_raster)), collapse = ", "), "\n", sep = "")

cat("\nReference raster info:\n")
cat("  CRS        : ", as.character(crs(ref_raster)), "\n", sep = "")
cat("  Resolution : ", paste(res(ref_raster), collapse = " x "), "\n", sep = "")
cat("  Extent     : ", paste(as.vector(ext(ref_raster)), collapse = ", "), "\n", sep = "")

# ============================================================================
# 3. Align raster to PM reference grid
# ============================================================================

cat("\nProjecting and aligning raster to PM reference grid...\n")
# Use nearest-neighbor because this is a categorical raster.
aligned_raster <- project(input_raster, ref_raster, method = "near")

cat("Masking aligned raster to the reference raster footprint...\n")
nsw_act_raster <- mask(aligned_raster, ref_raster)

names(nsw_act_raster) <- "NSW_ACT_Raster"

# ============================================================================
# 4. Write output
# ============================================================================

cat("\nWriting output raster to: ", output_raster_path, "\n", sep = "")
writeRaster(
  nsw_act_raster,
  filename = output_raster_path,
  overwrite = TRUE,
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

# ============================================================================
# 5. Verify output
# ============================================================================

cat("\nOutput raster summary:\n")
cat("  Name       : ", names(nsw_act_raster), "\n", sep = "")
cat("  Resolution : ", paste(res(nsw_act_raster), collapse = " x "), "\n", sep = "")
cat("  Extent     : ", paste(as.vector(ext(nsw_act_raster)), collapse = ", "), "\n", sep = "")
cat("  Dimensions : ", paste(dim(nsw_act_raster), collapse = " x "), "\n", sep = "")

non_na_cells <- global(!is.na(nsw_act_raster), "sum", na.rm = TRUE)[1, 1]
cat("  Non-NA cells: ", non_na_cells, "\n", sep = "")

cat("\n=== Completed successfully ===\n")
