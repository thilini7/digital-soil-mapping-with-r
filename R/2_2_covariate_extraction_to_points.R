# --- Install required packages (only if not already installed) ---
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(terra)
library(dplyr)

# =============================================================================
# FILE PATHS - CENTRALIZED CONFIGURATION
# =============================================================================
HomeDir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
setwd(HomeDir)

# INPUT: Change soil_property to match your data
soil_property <- "Phosphorus"  # Options: pH, OC, BD, CEC, EC, Clay, etc.

# Input paths
cov_path <- file.path(HomeDir, "Data", "data_in", "soil_covariates_aligned_v2")
soil_points_csv <- file.path(HomeDir, "Data", "data_out", "splined_data", soil_property,
                             paste0(soil_property, "_splined_NSW_1991_2020_Data.csv"))

# Output path
output_dir <- file.path(HomeDir, "Data", "data_out", "Soil_data_with_covariates")
output_file <- file.path(output_dir, paste0(soil_property, "_with_covariates_new.csv"))

# Print configuration
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("COVARIATE EXTRACTION - CONFIGURATION\n")
cat("Soil Property: ", soil_property, "\n", sep = "")
cat("Input Points: ", basename(soil_points_csv), "\n", sep = "")
cat("Covariates Directory: ", cov_path, "\n", sep = "")
cat("Output File: ", basename(output_file), "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

# =============================================================================
# LOAD AND EXTRACT COVARIATE DATA
# =============================================================================

# 1. Load covariates
cat("Loading covariate rasters...\n")
cov_files <- list.files(cov_path, pattern = ".tif$", full.names = TRUE)
if (length(cov_files) == 0) {
  stop("No .tif files found in: ", cov_path)
}

cov_stack <- rast(cov_files)

# Use filenames as layer names to avoid duplicates from internal layer names
names(cov_stack) <- tools::file_path_sans_ext(basename(cov_files))

cat("Loaded", nlyr(cov_stack), "covariate layers.\n")

# 2. Load soil point data (must have columns: Longitude, Latitude, property)
cat("Loading soil point data...\n")
soil_points <- read.csv(soil_points_csv)
cat("Loaded", nrow(soil_points), "soil points.\n")

# 3. Convert soil points to SpatVector for terra extraction
pts_vect <- vect(soil_points, geom = c("Longitude", "Latitude"), crs = crs(cov_stack))

# 4. Extract all raster values at point locations (one call, very fast)
cat("Extracting covariate values...\n")
cov_values <- terra::extract(cov_stack, pts_vect)
# Remove the first column (ID), as extract adds a point ID column
cov_values <- cov_values[, -1, drop = FALSE]

# 5. Combine soil property + extracted covariates
soil_with_covs <- cbind(soil_points, cov_values)

# 6. Save dataset for modeling
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(soil_with_covs, output_file, row.names = FALSE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("✅ EXTRACTION COMPLETE\n")
cat("Output file: ", output_file, "\n", sep = "")
cat("Rows: ", nrow(soil_with_covs), " | Columns: ", ncol(soil_with_covs), "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")
