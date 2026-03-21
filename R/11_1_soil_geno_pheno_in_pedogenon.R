# =============================================================================
# Calculate mean soil property per pedogenon class from Geno/Pheno soil maps
# =============================================================================
# For each soil property (pH, Bulk_Density, AWC) and soil type (Geno, Pheno):
#   1. Load pedogenon raster (1380 classes)
#   2. Load genosoil / phenosoil property prediction raster
#   3. Compute zonal mean of the property within each pedogenon class
#   4. Create a new raster where every pixel in a pedogenon class receives
#      that class's mean property value
#   5. Save results as GeoTIFF and export summary CSV
# =============================================================================

# --- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
})

set.seed(42)

# --- Terra options -----------------------------------------------------------
temp_dir <- file.path(getwd(), "temp_terra")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
terraOptions(memfrac = 0.8, tempdir = temp_dir)

# =============================================================================
# CONFIGURATION
# =============================================================================
depth_layer <- "X0.30cm"

soil_properties <- c("pH", "Bulk_Density", "AWC")


soil_types <- c("Geno", "Pheno")

pedogenon_path <- file.path("Data", "data_in_v2", "genosoil_phenosoil",
                            "pedogenon_nsw.tif")

# Build a named list of input prediction raster paths
# Geno: Models_Geno/mod.cubist.<property>/preds/<property>_pred_X0.30cm.tif
# Pheno: Models_Pheno/mod.cubist.<property>/preds/<property>_pred_X0.30cm.tif
get_pred_path <- function(soil_type, soil_property, depth) {
  model_dir <- ifelse(soil_type == "Geno", "Models_Geno", "Models_Pheno")
  file.path(model_dir, paste0("mod.cubist.", soil_property), "preds",
            paste0(soil_property, "_pred_", depth, ".tif"))
}

# Output directory
output_base <- file.path("Data", "data_out", "pedogenon_zonal_means")
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Print configuration
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PEDOGENON ZONAL MEAN — SOIL PROPERTY ASSIGNMENT\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("Depth layer      : ", depth_layer, "\n", sep = "")
cat("Soil properties  : ", paste(soil_properties, collapse = ", "), "\n", sep = "")
cat("Soil types       : ", paste(soil_types, collapse = ", "), "\n", sep = "")
cat("Pedogenon raster : ", pedogenon_path, "\n", sep = "")
cat("Output directory : ", output_base, "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

# =============================================================================
# 1. Load pedogenon raster
# =============================================================================
cat("Loading pedogenon raster...\n")
pedogenon_rast <- rast(pedogenon_path)
cat("  CRS        : ", crs(pedogenon_rast, describe = TRUE)$name, "\n", sep = "")
cat("  Resolution : ", paste(res(pedogenon_rast), collapse = " x "), "\n", sep = "")
cat("  Extent     : ", paste(round(as.vector(ext(pedogenon_rast)), 4), collapse = ", "), "\n", sep = "")

# Count unique pedogenon classes (excluding NA)
ped_vals <- unique(values(pedogenon_rast, mat = FALSE))
ped_vals <- ped_vals[!is.na(ped_vals)]
cat("  Unique pedogenon classes: ", length(ped_vals), "\n\n", sep = "")

# =============================================================================
# 2. Process each soil type × soil property combination
# =============================================================================
all_summaries <- list()

for (soil_type in soil_types) {
  for (soil_property in soil_properties) {

    combo_label <- paste0(soil_type, " — ", soil_property)
    cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")
    cat("Processing: ", combo_label, " (", depth_layer, ")\n", sep = "")
    cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")

    # --- 2a. Load prediction raster ---
    pred_path <- get_pred_path(soil_type, soil_property, depth_layer)
    if (!file.exists(pred_path)) {
      cat("  ** SKIPPED — file not found: ", pred_path, "\n\n", sep = "")
      next
    }
    cat("  Prediction raster: ", pred_path, "\n", sep = "")
    pred_rast <- rast(pred_path)

    cat("  Pred CRS        : ", crs(pred_rast, describe = TRUE)$name, "\n", sep = "")
    cat("  Pred Resolution : ", paste(res(pred_rast), collapse = " x "), "\n", sep = "")

    # --- 2b. Align rasters if needed ---
    # The pedogenon and prediction rasters must share extent, resolution, and CRS.
    needs_align <- !same.crs(pedogenon_rast, pred_rast) ||
                   !all(res(pedogenon_rast) == res(pred_rast)) ||
                   ext(pedogenon_rast) != ext(pred_rast)

    if (needs_align) {
      cat("  Aligning prediction raster to pedogenon grid...\n")
      pred_rast <- resample(pred_rast, pedogenon_rast, method = "bilinear")
      cat("  Alignment complete.\n")
    } else {
      cat("  Rasters already aligned.\n")
    }

    # --- 2c. Zonal mean ---
    cat("  Computing zonal mean per pedogenon class...\n")
    zonal_df <- zonal(pred_rast, pedogenon_rast, fun = "mean", na.rm = TRUE)
    colnames(zonal_df) <- c("pedogenon_class", "mean_value")

    # Remove rows where mean is NA (classes with no valid data)
    zonal_df <- zonal_df[!is.na(zonal_df$mean_value), ]

    cat("  Classes with valid mean: ", nrow(zonal_df), " / ", length(ped_vals), "\n", sep = "")
    cat("  Mean value range       : [", round(min(zonal_df$mean_value), 4),
        ", ", round(max(zonal_df$mean_value), 4), "]\n", sep = "")

    # --- 2d. Reclassify pedogenon raster with mean values ---
    cat("  Reclassifying pedogenon raster with zonal means...\n")
    rcl_matrix <- as.matrix(zonal_df[, c("pedogenon_class", "mean_value")])
    result_rast <- classify(pedogenon_rast, rcl = rcl_matrix, others = NA)

    # --- 2e. Save output raster ---
    out_dir <- file.path(output_base, soil_type)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    out_tif <- file.path(out_dir,
                         paste0(soil_property, "_", tolower(soil_type),
                                "_pedogenon_mean_", depth_layer, ".tif"))

    writeRaster(result_rast, out_tif, overwrite = TRUE,
                gdal = c("COMPRESS=LZW"))
    cat("  Raster saved : ", out_tif, "\n", sep = "")

    # --- 2f. Save zonal summary CSV ---
    out_csv <- file.path(out_dir,
                         paste0(soil_property, "_", tolower(soil_type),
                                "_pedogenon_zonal_summary_", depth_layer, ".csv"))
    write.csv(zonal_df, out_csv, row.names = FALSE)
    cat("  CSV saved    : ", out_csv, "\n", sep = "")

    # --- 2g. Store summary for final report ---
    all_summaries[[combo_label]] <- list(
      soil_type       = soil_type,
      soil_property   = soil_property,
      n_classes       = nrow(zonal_df),
      mean_min        = min(zonal_df$mean_value),
      mean_max        = max(zonal_df$mean_value),
      overall_mean    = mean(zonal_df$mean_value),
      output_raster   = out_tif,
      output_csv      = out_csv
    )

    cat("\n")
  }
}

# =============================================================================
# 3. Final summary report
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

summary_table <- do.call(rbind, lapply(all_summaries, function(s) {
  data.frame(
    Soil_Type     = s$soil_type,
    Property      = s$soil_property,
    N_Classes     = s$n_classes,
    Mean_Min      = round(s$mean_min, 4),
    Mean_Max      = round(s$mean_max, 4),
    Overall_Mean  = round(s$overall_mean, 4),
    stringsAsFactors = FALSE
  )
}))

print(summary_table, row.names = FALSE)

# Save combined summary
summary_csv <- file.path(output_base, paste0("pedogenon_zonal_summary_all_", depth_layer, ".csv"))
write.csv(summary_table, summary_csv, row.names = FALSE)
cat("\nCombined summary saved: ", summary_csv, "\n", sep = "")

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("DONE\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
