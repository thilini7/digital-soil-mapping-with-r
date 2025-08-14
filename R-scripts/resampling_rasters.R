# --- Install packages if needed ---
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")

library(terra)

# ----------- USER CONFIG -----------
cov_path <- "D:/InputData/Environmental_Covariate/Covariates" # source folder with raw rasters
out_path <- "D:/InputData/Environmental_Covariate/Covariates_Aligned"# output folder for aligned rasters
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

ref_raster_path <- "D:/InputData/Environmental_Covariate/Covariates/All/Annual_10th_Percen_RF_NSW_ACT_90m.tif"
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
  r <- rast(f)
  
  # Project if CRS is different
  if (!crs(r) == crs(ref_rast)) {
    cat("  - Projecting to reference CRS...\n")
    r <- project(r, ref_rast)
  }
  # Crop and resample to reference raster
  cat("  - Cropping to reference extent...\n")
  r_crop <- crop(r, ref_rast)
  cat("  - Resampling to reference grid...\n")
  # Use method = "bilinear" for continuous, "near" for categorical
  r_resamp <- resample(r_crop, ref_rast, method = "bilinear")
  
  # Save aligned raster
  writeRaster(r_resamp, out_file, overwrite = TRUE)
  cat("  ✓ Saved aligned raster to:", out_file, "\n\n")
}

cat("🎉 All rasters aligned and saved in:\n", out_path, "\n")
