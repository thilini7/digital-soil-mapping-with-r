# ==============================================================================
# Script: 12_calculate_awc.R
# Purpose: Calculate Available Water Capacity (AWC) using pedotransfer functions
# Reference: Francos et al. (2024) "Soil production capital – Relating available
#            water capacity to farmland price" Soil Security 16, 100157
#            PTFs from Searle et al. (2023)
# Input:  Sand, Clay, Silt, Bulk Density, Organic Carbon rasters (0-30 cm)
# Output: DUL, LL, AWC (%), AWC (mm) rasters
# ==============================================================================

library(terra)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

# Input directory containing soil property rasters
input_dir <- file.path("Data", "data_in_v2", "Sand_and_clay_and_silt_with_oc_and_bd")

# Input raster paths
sand_path <- file.path(input_dir, "NSW_Sand0_30_mean_220901reduced.tif")
clay_path <- file.path(input_dir, "NSW_Clay0_30_mean_220901reduced.tif")
silt_path <- file.path(input_dir, "NSW_Silt0_30_mean_220901reduced.tif")
bd_path   <- file.path(input_dir, "Bulk_Density_pred_X0.30cm.tif")
oc_path   <- file.path(input_dir, "Organic_Carbon_pred_X0.30cm.tif")

# Output directory
output_dir <- file.path("Data", "data_out", "AWC")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Depth of the soil layer in mm (0-30 cm = 300 mm)
depth_mm <- 300

# ==============================================================================
# 2. LOAD RASTERS
# ==============================================================================

cat("Loading input rasters...\n")

sand <- rast(sand_path)
clay <- rast(clay_path)
silt <- rast(silt_path)
bd   <- rast(bd_path)
oc   <- rast(oc_path)

cat("  Sand: ", sand_path, "\n")
cat("  Clay: ", clay_path, "\n")
cat("  Silt: ", silt_path, "\n")
cat("  BD:   ", bd_path, "\n")
cat("  OC:   ", oc_path, "\n\n")

# ==============================================================================
# 3. ALIGN RASTERS TO COMMON EXTENT AND RESOLUTION
# ==============================================================================

cat("Aligning rasters to common extent and resolution...\n")

# Use the smallest common extent
ref <- sand
clay <- resample(clay, ref, method = "bilinear")
silt <- resample(silt, ref, method = "bilinear")
bd   <- resample(bd,   ref, method = "bilinear")
oc   <- resample(oc,   ref, method = "bilinear")

cat("  Aligned to extent: ", as.vector(ext(ref)), "\n")
cat("  Resolution:        ", res(ref), "\n")
cat("  CRS:               ", crs(ref, describe = TRUE)$name, "\n\n")

# ==============================================================================
# 4. CALCULATE AWC USING PEDOTRANSFER FUNCTIONS (Searle et al., 2023)
# ==============================================================================

# Eq. (1): DUL (%) = 48.98 - 21.86*BD + 0.36*Clay - 0.06*Sand
#                     - 0.19*sqrt(Silt) + 2.90*sqrt(OC)
#
# Eq. (2): LL (%)  = 17.40 - 10.05*BD + 0.34*Clay - 0.02*Sand + 0.18*Silt
#
# AWC (%) = DUL - LL
# AWC (mm) = AWC(%) / 100 * depth_mm

cat("Calculating DUL (Drained Upper Limit)...\n")
DUL <- 48.98 - 21.86 * bd + 0.36 * clay - 0.06 * sand -
  0.19 * sqrt(silt) + 2.90 * sqrt(oc)

cat("Calculating LL (Crop Lower Limit)...\n")
LL <- 17.40 - 10.05 * bd + 0.34 * clay - 0.02 * sand + 0.18 * silt

cat("Calculating AWC...\n")
AWC_pct <- DUL - LL

# Clamp AWC to physically reasonable range (0-100%)
AWC_pct <- clamp(AWC_pct, lower = 0, upper = 100)

# Convert AWC percentage to mm for the 0-30 cm layer
AWC_mm <- AWC_pct / 100 * depth_mm

cat("  Depth layer: 0-30 cm (", depth_mm, " mm)\n")
cat("  AWC (%) range:  ", round(global(AWC_pct, "min", na.rm = TRUE)[[1]], 2), " - ",
    round(global(AWC_pct, "max", na.rm = TRUE)[[1]], 2), "\n")
cat("  AWC (mm) range: ", round(global(AWC_mm, "min", na.rm = TRUE)[[1]], 2), " - ",
    round(global(AWC_mm, "max", na.rm = TRUE)[[1]], 2), "\n\n")

# ==============================================================================
# 5. SAVE OUTPUT RASTERS
# ==============================================================================

cat("Saving output rasters...\n")

write_opts <- list(gdal = c("COMPRESS=DEFLATE"))

writeRaster(DUL, file.path(output_dir, "DUL_pct_X0.30cm.tif"),
            overwrite = TRUE, wopt = write_opts)
writeRaster(LL, file.path(output_dir, "LL_pct_X0.30cm.tif"),
            overwrite = TRUE, wopt = write_opts)
writeRaster(AWC_pct, file.path(output_dir, "AWC_pct_X0.30cm.tif"),
            overwrite = TRUE, wopt = write_opts)
writeRaster(AWC_mm, file.path(output_dir, "AWC_mm_X0.30cm.tif"),
            overwrite = TRUE, wopt = write_opts)

cat("  Saved 4 rasters to: ", output_dir, "\n\n")

# ==============================================================================
# 6. SUMMARY STATISTICS
# ==============================================================================

cat("===== SUMMARY STATISTICS =====\n\n")

summary_stats <- data.frame(
  Layer = c("DUL (%)", "LL (%)", "AWC (%)", "AWC (mm)"),
  Min = round(c(
    global(DUL, "min", na.rm = TRUE)[[1]],
    global(LL, "min", na.rm = TRUE)[[1]],
    global(AWC_pct, "min", na.rm = TRUE)[[1]],
    global(AWC_mm, "min", na.rm = TRUE)[[1]]
  ), 2),
  Mean = round(c(
    global(DUL, "mean", na.rm = TRUE)[[1]],
    global(LL, "mean", na.rm = TRUE)[[1]],
    global(AWC_pct, "mean", na.rm = TRUE)[[1]],
    global(AWC_mm, "mean", na.rm = TRUE)[[1]]
  ), 2),
  Max = round(c(
    global(DUL, "max", na.rm = TRUE)[[1]],
    global(LL, "max", na.rm = TRUE)[[1]],
    global(AWC_pct, "max", na.rm = TRUE)[[1]],
    global(AWC_mm, "max", na.rm = TRUE)[[1]]
  ), 2)
)
print(summary_stats)

# ==============================================================================
# 7. VISUALISATION
# ==============================================================================

cat("\nGenerating plots...\n")

# Save plots to PDF
pdf(file.path(output_dir, "AWC_Maps.pdf"), width = 12, height = 8)

par(mfrow = c(2, 2), mar = c(3, 3, 3, 5))

plot(DUL, main = "DUL (%) - Drained Upper Limit\n0-30 cm",
     col = hcl.colors(50, "Blues 3"))
plot(LL, main = "LL (%) - Crop Lower Limit\n0-30 cm",
     col = hcl.colors(50, "Reds 3"))
plot(AWC_pct, main = "AWC (%) - Available Water Capacity\n0-30 cm",
     col = hcl.colors(50, "Viridis"))
plot(AWC_mm, main = "AWC (mm)\n0-30 cm",
     col = hcl.colors(50, "Viridis"))

dev.off()

cat("  Plots saved to: ", file.path(output_dir, "AWC_Maps.pdf"), "\n")

# Also display in R
par(mfrow = c(2, 2), mar = c(3, 3, 3, 5))
plot(DUL, main = "DUL (%)", col = hcl.colors(50, "Blues 3"))
plot(LL, main = "LL (%)", col = hcl.colors(50, "Reds 3"))
plot(AWC_pct, main = "AWC (%)", col = hcl.colors(50, "Viridis"))
plot(AWC_mm, main = "AWC (mm)", col = hcl.colors(50, "Viridis"))

cat("\n===== COMPLETE =====\n")
cat("Reference: Francos et al. (2024) Soil Security 16, 100157\n")
cat("PTFs:      Searle et al. (2023)\n")
cat("Data:      0-30 cm depth for NSW\n")
