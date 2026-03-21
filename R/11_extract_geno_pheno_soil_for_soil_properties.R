# =============================================================================
# Extract random points from Geno/Pheno Soil raster, then extract predicted values
# =============================================================================
# Steps:
#   1. Load geno_pheno_soil raster
#   2. Generate random points within the raster extent
#   3. Extract geno_pheno_soil values at those points
#   4. Filter out NA and NoData values (-9999)
#   5. Export CSV with id, Longitude, Latitude, Geno_Pheno_Soil
#   6. Extract soil property predictions at the filtered points
#   7. Export final CSV with id, Longitude, Latitude, Geno_Pheno_Soil, prediction
# =============================================================================

# --- Install required packages (only if not already installed) ---
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(terra)
library(dplyr)

set.seed(42)  # For reproducibility

# INPUT: Change soil_property and depth_layer to match your data
soil_property <- "AWC"       # Options: Organic_Carbon, Nitrogen, Phosphorus, pH, Bulk_Density, CEC, EC, Sum_of_Bases
depth_layer   <- "X0.30cm"  # Options: X0.30cm, X30.60cm, X60.100cm, X100.200cm
geno_pheno_soil_type <- "genosoil"

# ---- Input paths ----
geno_pheno_soil_path <- file.path("Data", "data_in_v2", "genosoil_phenosoil", "genosoil_nsw.tif")
pred_path     <- file.path("Models_v2", paste0("mod.cubist.", soil_property), "preds",
                           paste0(soil_property, "_pred_", depth_layer, ".tif"))
#This path for different raster files
# pred_path     <- file.path("Data", "data_out", "AWC",
#                                   paste0(soil_property, "_pct_", depth_layer, ".tif"))

# ---- Output paths ----
output_dir <- file.path("Data", "data_out", "geno_pheno_soil_extraction", soil_property, depth_layer)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

geno_pheno_soil_csv <- file.path(output_dir, paste0(geno_pheno_soil_type, "_random_", soil_property, "_", depth_layer, ".csv"))
final_csv    <- file.path(output_dir, paste0(geno_pheno_soil_type, "_", soil_property, "_", depth_layer, ".csv"))

# ---- Print configuration ----
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("GENO and PHENO SOIL / SOIL PROPERTY EXTRACTION - CONFIGURATION\n")
cat("Soil Property    : ", soil_property, "\n", sep = "")
cat("Depth Layer      : ", depth_layer, "\n", sep = "")
cat("Geno_Pheno soil Raster  : ", geno_pheno_soil_path, "\n", sep = "")
cat("Prediction Raster: ", pred_path, "\n", sep = "")
cat("Output (geno_pheno): ", geno_pheno_soil_csv, "\n", sep = "")
cat("Output (final)   : ", final_csv, "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

# =============================================================================
# 1. Load geno_pheno_soil raster
# =============================================================================
cat("Loading geno_pheno_soil raster...\n")
geno_pheno_soil_rast <- rast(geno_pheno_soil_path)
cat("  CRS : ", crs(geno_pheno_soil_rast, describe = TRUE)$name, "\n", sep = "")
cat("  Res : ", paste(res(geno_pheno_soil_rast), collapse = " x "), "\n", sep = "")

# Quick look at value range to verify NoData values
cat("\nGeno_Pheno_Soil raster summary:\n")
print(summary(geno_pheno_soil_rast))

# Check for common NoData sentinel values
vals_sample <- values(geno_pheno_soil_rast, mat = FALSE)
vals_sample <- vals_sample[!is.na(vals_sample)]
cat("\nValue range: [", min(vals_sample), ", ", max(vals_sample), "]\n", sep = "")

nodata_candidates <- c(-9999, -9999.0, -3.4e+38)
for (nd in nodata_candidates) {
  n <- sum(vals_sample == nd, na.rm = TRUE)
  if (n > 0) cat("  Found", n, "pixels with value", nd, "(likely NoData)\n")
}

# =============================================================================
# 2. Iteratively sample random points until we have >= min_valid valid points
# =============================================================================
min_valid    <- 8000   # Minimum number of valid (non-NA, non-NoData) points required
batch_size   <- 16000  # Points to sample per iteration
max_iter     <- 20     # Safety limit to avoid infinite loops

cat("\nTarget valid points: ", min_valid, "\n", sep = "")

geno_pheno_soil_filtered <- data.frame()
total_sampled <- 0

for (iter in seq_len(max_iter)) {
  cat("  Iteration ", iter, ": sampling ", batch_size, " points ... ", sep = "")
  
  pts <- spatSample(geno_pheno_soil_rast, size = batch_size, method = "random",
                    as.points = TRUE, na.rm = FALSE)
  total_sampled <- total_sampled + nrow(pts)
  
  coords <- crds(pts)
  vals   <- values(pts)
  
  batch_df <- data.frame(
    Longitude        = coords[, "x"],
    Latitude         = coords[, "y"],
    Geno_Pheno_Soil  = vals[, 1]
  )
  
  # Filter out NA and NoData (-9999) values
  batch_valid <- batch_df %>%
    filter(!is.na(Geno_Pheno_Soil), Geno_Pheno_Soil != -9999)
  
  geno_pheno_soil_filtered <- rbind(geno_pheno_soil_filtered, batch_valid)
  
  cat(nrow(batch_valid), " valid (cumulative: ", nrow(geno_pheno_soil_filtered), ")\n", sep = "")
  
  if (nrow(geno_pheno_soil_filtered) >= min_valid) break
}

# Trim to exactly min_valid if we overshot (keeps it clean)
if (nrow(geno_pheno_soil_filtered) > min_valid) {
  geno_pheno_soil_filtered <- geno_pheno_soil_filtered[seq_len(min_valid), ]
}

# Assign sequential IDs
geno_pheno_soil_filtered$id <- seq_len(nrow(geno_pheno_soil_filtered))
geno_pheno_soil_filtered <- geno_pheno_soil_filtered[, c("id", "Longitude", "Latitude", "Geno_Pheno_Soil")]

cat("\nSampling summary:\n")
cat("  Total sampled      : ", total_sampled, "\n", sep = "")
cat("  Valid points kept  : ", nrow(geno_pheno_soil_filtered), "\n", sep = "")

if (nrow(geno_pheno_soil_filtered) < min_valid) {
  warning("Could not reach ", min_valid, " valid points after ", max_iter,
          " iterations. Got ", nrow(geno_pheno_soil_filtered), ".")
}

# =============================================================================
# 5. Export genosoil CSV
# =============================================================================
write.csv(geno_pheno_soil_filtered, geno_pheno_soil_csv, row.names = FALSE)
cat("\nGeno_Pheno_Soil CSV saved to: ", geno_pheno_soil_csv, "\n", sep = "")

# =============================================================================
# 6. Load pH prediction raster and extract values at filtered points
# =============================================================================
cat("\nLoading", soil_property, "prediction raster...\n")
pred_rast <- rast(pred_path)
cat("  CRS : ", crs(pred_rast, describe = TRUE)$name, "\n", sep = "")
cat("  Res : ", paste(res(pred_rast), collapse = " x "), "\n", sep = "")

# Convert filtered points to SpatVector for extraction
pts_vect <- vect(geno_pheno_soil_filtered, geom = c("Longitude", "Latitude"),
                 crs = crs(geno_pheno_soil_rast))

# Reproject points if CRS differs between the two rasters
if (!same.crs(geno_pheno_soil_rast, pred_rast)) {
  cat("  Reprojecting points to prediction raster CRS...\n")
  pts_vect <- project(pts_vect, crs(pred_rast))
}

cat("Extracting", soil_property, "values at", nrow(geno_pheno_soil_filtered), "points...\n")
pred_vals <- terra::extract(pred_rast, pts_vect)

# =============================================================================
# 7. Build final dataframe and export
# =============================================================================
# Dynamic column name: e.g. "X0.30cm"
pred_col_name <- paste0(depth_layer)

final_df <- geno_pheno_soil_filtered
final_df[[pred_col_name]] <- pred_vals[, 2]  # Column 2 = raster values (column 1 = ID)

cat("\n", soil_property, " extraction summary:\n", sep = "")
cat("  NA values          : ", sum(is.na(final_df[[pred_col_name]])), "\n", sep = "")
cat("  Valid values       : ", sum(!is.na(final_df[[pred_col_name]])), "\n", sep = "")

# Remove rows with NA in the prediction column
final_df <- final_df[!is.na(final_df[[pred_col_name]]), ]
final_df$id <- seq_len(nrow(final_df))
cat("  After removing NA  : ", nrow(final_df), " rows\n", sep = "")

write.csv(final_df, final_csv, row.names = FALSE)
cat("\nFinal CSV saved to: ", final_csv, "\n", sep = "")

# =============================================================================
# 8. Plot sampled points on NSW/ACT map
# =============================================================================
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("nswgeo", quietly = TRUE)) install.packages("nswgeo")

library(sf)
library(ggplot2)
library(nswgeo)

cat("\nGenerating map of sampled points on NSW/ACT...\n")

# Load NSW boundary
nsw_boundary <- nswgeo::nsw
if (!inherits(nsw_boundary, "sf")) nsw_boundary <- st_as_sf(nsw_boundary)

# Build the map
map_plot <- ggplot() +
  geom_sf(data = nsw_boundary, fill = "lightblue", color = "black", alpha = 0.3) +
  geom_point(data = final_df, aes(x = Longitude, y = Latitude, color = .data[[pred_col_name]]),
             size = 0.8, alpha = 0.7) +
  scale_color_viridis_c(name = paste0(soil_property, "\n(", depth_layer, ")"),
                        na.value = "grey50") +
  labs(
    title = paste("Sampled Points -", soil_property, depth_layer),
    subtitle = paste(nrow(final_df), "points from", geno_pheno_soil_type, "NSW"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Save the map
map_file <- file.path(output_dir, paste0("map_", geno_pheno_soil_type, "_", soil_property, "_", depth_layer, ".png"))
ggsave(map_file, plot = map_plot, width = 10, height = 8, dpi = 300)
cat("Map saved to: ", map_file, "\n", sep = "")

# Display the map
print(map_plot)

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("DONE\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
