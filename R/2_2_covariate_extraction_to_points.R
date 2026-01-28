# --- Install required packages (only if not already installed) ---
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(terra)
library(dplyr)

# 1. Define path to aligned covariates folder
cov_path <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_in/soil_covariates_aligned"
cov_files <- list.files(cov_path, pattern = ".tif$", full.names = TRUE)

# 2. Load all aligned raster covariates into a SpatRaster stack (super fast in terra)
cov_stack <- rast(cov_files)

# Use filenames as layer names to avoid duplicates from internal layer names
names(cov_stack) <- tools::file_path_sans_ext(basename(cov_files))

cat("Loaded", nlyr(cov_stack), "covariate layers.\n")

# 3. Load soil point data (must have columns: Longitude, Latitude, property)
soil_points <- read.csv("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/splined_data/Phosphorus/Phosphorus_splined_NSW_1991_2020_Data.csv")
cat("Loaded", nrow(soil_points), "soil points.\n")

# 4. Convert soil points to SpatVector for terra extraction
pts_vect <- vect(soil_points, geom = c("Longitude", "Latitude"), crs = crs(cov_stack))

# 5. Extract all raster values at point locations (one call, very fast)
cov_values <- terra::extract(cov_stack, pts_vect)
# Remove the first column (ID), as extract adds a point ID column
cov_values <- cov_values[, -1, drop = FALSE]

# 6. Combine soil property + extracted covariates
soil_with_covs <- cbind(soil_points, cov_values)

# 7. Save dataset for modeling in a separate folder
output_dir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/Soil_data_with_covariates"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(output_dir, "Phosphorus_with_covariates_new.csv")
write.csv(soil_with_covs, output_file, row.names = FALSE)

cat("✅ File saved to:", output_file, "\n")
