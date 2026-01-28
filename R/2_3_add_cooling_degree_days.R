# --- Add Cooling Degree Days covariate to existing dataset ---
# This script adds the missing climate covariate that has future projections available

library(terra)

# 1. Load existing soil data with covariates
soil_data <- read.csv("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/Soil_data_with_covariates/OC_with_covariates_new.csv")
cat("Loaded", nrow(soil_data), "soil points\n")

# 2. Load the Cooling Degree Days raster (from original covariates folder)
cdd_path <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_in/soil_covariates/Cooling_Degree_Days_Annual_90m.tif"

if (!file.exists(cdd_path)) {
  stop("Cooling Degree Days raster not found at: ", cdd_path)
}

cdd_rast <- rast(cdd_path)
cat("Loaded Cooling Degree Days raster\n")

# 3. Extract CDD values at soil point locations
pts_vect <- vect(soil_data, geom = c("Longitude", "Latitude"), crs = crs(cdd_rast))
cdd_values <- terra::extract(cdd_rast, pts_vect)

# 4. Add to dataframe
soil_data$Cooling_Degree_Days_Annual_90m <- cdd_values[, 2]  # First column is ID

# Check for NAs
n_na <- sum(is.na(soil_data$Cooling_Degree_Days_Annual_90m))
cat("Points with NA for Cooling Degree Days:", n_na, "of", nrow(soil_data), "\n")

# 5. Save updated dataset
output_file <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/Soil_data_with_covariates/OC_with_covariates_new.csv"
write.csv(soil_data, output_file, row.names = FALSE)
cat("✅ Updated dataset saved with Cooling Degree Days covariate\n")

# Show column names
cat("\nCurrent covariates:\n")
cov_cols <- setdiff(names(soil_data), c("id", "Longitude", "Latitude", 
                                         "X0.5cm", "X5.15cm", "X15.30cm", 
                                         "X30.60cm", "X60.100cm", "X100.200cm"))
cat(paste(cov_cols, collapse = "\n"))
cat("\n\nTotal covariates:", length(cov_cols), "\n")
