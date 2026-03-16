# ============================================================================
# Script: 8_1_batch_shp_to_tiff_with_reference.R
# Purpose: Convert one .shp file to GeoTIFF using a reference raster grid
# ============================================================================

# Load required libraries
library(sf)
library(terra)

# ============================================================================
# 1. Define paths and options
# ============================================================================

shapefile_path <- file.path("Data", "data_in", "land_nswlanduse2017v1p5", "NSWLanduse2017_Ver1_5_20230921.shp")
reference_raster_path <- file.path(
  "Data", "data_in", "soil_covariates_for_geno_pheno_soil_aligned",
  "PM_radmap_v4_2019_filtered_dose_GAPFilled.tif"
)
output_dir <- file.path("Data", "data_out", "shp_to_tiff")

# Preferred field order. If none found, the script picks the first non-geometry field.
preferred_fields <- c("LDI", "landuse", "Landuse", "LANDUSE", "CLASS", "Class")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(shapefile_path)) {
  stop("Shapefile not found: ", shapefile_path)
}
if (!file.exists(reference_raster_path)) {
  stop("Reference raster not found: ", reference_raster_path)
}

# ============================================================================
# 2. Load reference raster (target grid/alignment)
# ============================================================================

cat("Loading reference raster...\n")
ref_raster <- rast(reference_raster_path)

cat("Reference raster info:\n")
cat("  CRS:", as.character(crs(ref_raster)), "\n")
cat("  Resolution:", paste(res(ref_raster), collapse = " x "), "\n")
cat("  Extent:", paste(as.vector(ext(ref_raster)), collapse = ", "), "\n")
cat("  Dimensions:", paste(dim(ref_raster), collapse = " x "), "\n")

template_raster <- ref_raster
template_raster[] <- NA

# ============================================================================
# 3. Load shapefile and pick raster field
# ============================================================================

cat("\nLoading shapefile...\n")
# Allow GDAL to recreate missing .shx index files when possible.
Sys.setenv(SHAPE_RESTORE_SHX = "YES")

shx_path <- paste0(tools::file_path_sans_ext(shapefile_path), ".shx")
if (!file.exists(shx_path)) {
  cat("Warning: .shx file not found. GDAL will attempt to restore it: ", shx_path, "\n", sep = "")
}

sf_obj <- st_read(shapefile_path, quiet = TRUE)
cat("Shapefile loaded with ", nrow(sf_obj), " features\n", sep = "")

# Select a rasterization field from available columns.
pick_field <- function(sf_obj, preferred) {
  cols <- names(sf_obj)
  candidates <- cols[cols != attr(sf_obj, "sf_column")]

  hit <- preferred[preferred %in% candidates]
  if (length(hit) > 0) return(hit[1])

  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

field_name <- pick_field(sf_obj, preferred_fields)
if (is.na(field_name)) {
  stop("No non-geometry attribute fields found.")
}
cat("Raster field: ", field_name, "\n", sep = "")

# ============================================================================
# 4. Convert shapefile to TIFF
# ============================================================================

if (st_crs(sf_obj) != st_crs(as.character(crs(ref_raster)))) {
  cat("Reprojecting to reference CRS...\n")
  sf_obj <- st_transform(sf_obj, crs = as.character(crs(ref_raster)))
}

vec_obj <- vect(sf_obj)

cat("Rasterizing...\n")
out_raster <- rasterize(
  x = vec_obj,
  y = template_raster,
  field = field_name,
  background = NA
)

out_name <- paste0(tools::file_path_sans_ext(basename(shapefile_path)), "_", field_name, ".tif")
out_path <- file.path(output_dir, out_name)

writeRaster(
  out_raster,
  filename = out_path,
  overwrite = TRUE,
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

non_na <- global(!is.na(out_raster), "sum", na.rm = TRUE)[1, 1]
cat("Output: ", out_path, "\n", sep = "")
cat("Non-NA cells: ", non_na, "\n", sep = "")

# ============================================================================
# 5. Summary
# ============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("Single shapefile conversion complete\n")
cat("  Input shapefile  : ", shapefile_path, "\n", sep = "")
cat("  Output TIFF      : ", out_path, "\n", sep = "")
cat("  Output directory : ", output_dir, "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
