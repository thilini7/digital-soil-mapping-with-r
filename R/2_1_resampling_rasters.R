# --- Install packages if needed ---
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")

library(terra)

# Set terra options - create a proper temp directory
temp_dir <- file.path(getwd(), "temp_terra")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
terraOptions(memfrac = 0.8, tempdir = temp_dir)

# ----------- USER CONFIG -----------
cov_path <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_in/soil_covariates" # source folder with raw rasters
out_path <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_in/soil_covariates_aligned_v2"# output folder for aligned rasters
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

ref_raster_path <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_in/soil_covariates/Relief_TPI_3s.tif"
ref_rast <- rast(ref_raster_path)
cat("Reference raster:\n", ref_raster_path, "\n")
cat("  Extent:", paste(as.vector(ext(ref_rast)), collapse = " "), "\n")
cat("  Resolution:", paste(res(ref_rast), collapse = " "), "\n")
cat("  CRS:", crs(ref_rast), "\n\n")

cov_files <- list.files(cov_path, pattern = ".tif$", full.names = TRUE)
cat("Found", length(cov_files), "covariate rasters in", cov_path, "\n\n")

# ----------- ALIGN & SAVE LOOP -----------
for (i in seq_along(cov_files)) {
  f <- cov_files[i]
  out_file <- file.path(out_path, basename(f))
  
  cat(sprintf("[%d/%d] Processing: %s\n", i, length(cov_files), f))
  
  tryCatch({
    r <- rast(f)
    
    # Align to reference: use a single-step approach
    cat("  - Aligning to reference grid...\n")
    
    # Create output raster matching reference, then fill with resampled values
    # This writes directly to output file avoiding temp file issues
    if (same.crs(r, ref_rast)) {
      # Same CRS: crop then resample directly to file
      r_crop <- crop(r, ref_rast, snap = "out")
      resample(r_crop, ref_rast, method = "bilinear", filename = out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    } else {
      # Different CRS: project directly to file
      project(r, ref_rast, method = "bilinear", align = TRUE, filename = out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    }
    
    cat("  ✓ Saved aligned raster to:", out_file, "\n\n")
    
    # Clean up
    rm(r)
    if (exists("r_crop")) rm(r_crop)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("  ✗ ERROR:", conditionMessage(e), "\n\n")
  })
}

cat("🎉 All rasters aligned and saved in:\n", out_path, "\n")

# Clean up temp directory
unlink(temp_dir, recursive = TRUE)
