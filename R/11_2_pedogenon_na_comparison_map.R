# =============================================================================
# Compare Geno vs Pheno: Highlight pedogenon classes with NA (removed) values
# =============================================================================
# For each soil property, identify pedogenon classes that have no valid zonal
# mean in Geno, Pheno, or both. Produce a categorical map showing:
#   1 = Valid in both Geno and Pheno
#   2 = NA in Geno only
#   3 = NA in Pheno only
#   4 = NA in both
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

pedogenon_path <- file.path("Data", "data_in_v2", "genosoil_phenosoil",
                            "pedogenon_nsw.tif")

get_pred_path <- function(soil_type, soil_property, depth) {
  model_dir <- ifelse(soil_type == "Geno", "Models_Geno", "Models_Pheno")
  file.path(model_dir, paste0("mod.cubist.", soil_property), "preds",
            paste0(soil_property, "_pred_", depth, ".tif"))
}

output_base <- file.path("Data", "data_out", "pedogenon_na_comparison")
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load pedogenon raster
# =============================================================================
cat("Loading pedogenon raster...\n")
pedogenon_rast <- rast(pedogenon_path)

ped_vals <- unique(values(pedogenon_rast, mat = FALSE))
ped_vals <- ped_vals[!is.na(ped_vals)]
cat("  Unique pedogenon classes: ", length(ped_vals), "\n\n", sep = "")

# =============================================================================
# 2. For each soil property, compare Geno vs Pheno NA patterns
# =============================================================================
for (soil_property in soil_properties) {

  cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat("Property: ", soil_property, " (", depth_layer, ")\n", sep = "")
  cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

  # --- Load and align both prediction rasters --------------------------------
  zonal_list <- list()

  for (soil_type in c("Geno", "Pheno")) {
    pred_path <- get_pred_path(soil_type, soil_property, depth_layer)
    if (!file.exists(pred_path)) {
      cat("  ** SKIPPED — file not found: ", pred_path, "\n")
      next
    }

    pred_rast <- rast(pred_path)

    needs_align <- !same.crs(pedogenon_rast, pred_rast) ||
                   !all(res(pedogenon_rast) == res(pred_rast)) ||
                   ext(pedogenon_rast) != ext(pred_rast)

    if (needs_align) {
      cat("  Aligning ", soil_type, " raster to pedogenon grid...\n", sep = "")
      pred_rast <- resample(pred_rast, pedogenon_rast, method = "bilinear")
    }

    zdf <- zonal(pred_rast, pedogenon_rast, fun = "mean", na.rm = TRUE)
    colnames(zdf) <- c("pedogenon_class", "mean_value")
    zdf <- zdf[!is.na(zdf$mean_value), ]

    zonal_list[[soil_type]] <- zdf
    cat("  ", soil_type, " valid classes: ", nrow(zdf), "\n", sep = "")
  }

  if (length(zonal_list) < 2) {
    cat("  ** Cannot compare — missing one or both rasters. Skipping.\n\n")
    next
  }

  # --- Identify NA classes per type ------------------------------------------
  geno_classes  <- zonal_list[["Geno"]]$pedogenon_class
  pheno_classes <- zonal_list[["Pheno"]]$pedogenon_class

  all_classes <- data.frame(pedogenon_class = ped_vals) %>%
    mutate(
      geno_valid  = pedogenon_class %in% geno_classes,
      pheno_valid = pedogenon_class %in% pheno_classes,
      category = case_when(
         geno_valid &  pheno_valid ~ 1L,  # Valid in both
        !geno_valid &  pheno_valid ~ 2L,  # NA in Geno only
         geno_valid & !pheno_valid ~ 3L,  # NA in Pheno only
        !geno_valid & !pheno_valid ~ 4L   # NA in both
      )
    )

  # --- Print summary ---------------------------------------------------------
  cat("\n  Category summary:\n")
  cat("    1 - Valid in both  : ", sum(all_classes$category == 1), " classes\n", sep = "")
  cat("    2 - NA in Geno only: ", sum(all_classes$category == 2), " classes\n", sep = "")
  cat("    3 - NA in Pheno only:", sum(all_classes$category == 3), " classes\n", sep = "")
  cat("    4 - NA in both     : ", sum(all_classes$category == 4), " classes\n", sep = "")

  # List the specific classes that differ
  geno_only_na  <- all_classes$pedogenon_class[all_classes$category == 2]
  pheno_only_na <- all_classes$pedogenon_class[all_classes$category == 3]
  both_na       <- all_classes$pedogenon_class[all_classes$category == 4]

  if (length(geno_only_na) > 0) {
    cat("\n  Pedogenon classes NA in Geno only: ",
        paste(geno_only_na, collapse = ", "), "\n", sep = "")
  }
  if (length(pheno_only_na) > 0) {
    cat("  Pedogenon classes NA in Pheno only: ",
        paste(pheno_only_na, collapse = ", "), "\n", sep = "")
  }
  if (length(both_na) > 0) {
    cat("  Pedogenon classes NA in both: ",
        paste(both_na, collapse = ", "), "\n", sep = "")
  }

  # --- Build categorical raster ----------------------------------------------
  rcl_matrix <- as.matrix(all_classes[, c("pedogenon_class", "category")])
  na_map <- classify(pedogenon_rast, rcl = rcl_matrix, others = NA)

  # Assign category labels
  levels(na_map) <- data.frame(
    id    = 1:4,
    label = c("Valid in both",
              "NA in Geno only",
              "NA in Pheno only",
              "NA in both")
  )

  # --- Save raster -----------------------------------------------------------
  out_tif <- file.path(output_base,
                       paste0(soil_property, "_na_comparison_", depth_layer, ".tif"))
  writeRaster(na_map, out_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  cat("\n  Raster saved: ", out_tif, "\n", sep = "")

  # --- Save CSV summary ------------------------------------------------------
  out_csv <- file.path(output_base,
                       paste0(soil_property, "_na_comparison_", depth_layer, ".csv"))
  write.csv(all_classes, out_csv, row.names = FALSE)
  cat("  CSV saved   : ", out_csv, "\n", sep = "")

  # --- Plot map --------------------------------------------------------------
  out_png <- file.path(output_base,
                       paste0(soil_property, "_na_comparison_", depth_layer, ".png"))

  png(out_png, width = 1200, height = 900, res = 150)

  plot(na_map,
       col  = c("grey80",      # 1: Valid in both
                "red",          # 2: NA in Geno only
                "blue",         # 3: NA in Pheno only
                "purple"),      # 4: NA in both
       main = paste0(soil_property, " — Pedogenon NA Comparison (", depth_layer, ")"),
       mar  = c(3, 3, 3, 8),
       plg  = list(title = "Category"))

  dev.off()
  cat("  Map saved   : ", out_png, "\n\n", sep = "")
}

cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("DONE\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
