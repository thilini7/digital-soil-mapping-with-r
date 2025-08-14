library(terra)

# --- Paths ---
cov_folder <- "/Users/neo/Development/Data-Science/data-science/Data/data_in/soil_covariates_aligned" 
model_file <- "/Users/neo/Development/Data-Science/data-science/mod.regression.oc/models/X0.5cm.rds"
output_tif <- "/Users/neo/Development/Data-Science/data-science/maps/oc_pred_X0.5cm_v1.tif"

# --- Load your covariate names (should match your model!) ---
# e.g. from your previous script or saved .RData
load("/Users/neo/Development/Data-Science/data-science/Data/data_out/RData/oc_covs_regression.RData")
# Now cov_names is available

# --- List and match raster files to covariate names ---
cov_files <- list.files(cov_folder, pattern = "\\.tif$", full.names = TRUE)
base_names <- tools::file_path_sans_ext(basename(cov_files))
ordered_files <- cov_files[match(cov_names, base_names)]

print(setdiff(cov_names, base_names))

# --- Error check: missing files? ---
if (any(is.na(ordered_files))) {
  cat("Missing raster files for these covariates:\n")
  print(cov_names[is.na(ordered_files)])
  stop("Some covariate rasters are missing! Aborting.")
}

# --- Stack rasters in model order and assign names ---
cov_stack <- rast(ordered_files)
names(cov_stack) <- cov_names

# --- Check: Layer names match? ---
if (!all(names(cov_stack) == cov_names)) {
  stop("Raster layer names do not match covariate names! Check file names and order.")
}
cat("Covariate stack ready!\n")

# --- Load the trained random forest model for a given depth ---
model <- readRDS(model_file)

# --- Define a wrapper for ranger prediction ---
predict_fun <- function(model, data, ...) {
  predict(model, data = data)$predictions
}

# --- Predict SOC for all pixels in NSW ---
cat("Predicting SOC map...\n")
soc_pred <- terra::predict(cov_stack, model, fun = predict_fun, na.rm = TRUE, progress = "text")

# --- Save and plot the result ---
writeRaster(soc_pred, output_tif, overwrite = TRUE)
plot(soc_pred, main = "Predicted OC 0–5cm")

cat("Map saved to:", output_tif, "\n")
