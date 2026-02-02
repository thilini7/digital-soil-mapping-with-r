library(terra)
library(Cubist)

# --- Configuration ---
HomeDir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
soil_property <- "Organic_Carbon"
depth_name <- "X0.5cm"  # Change to desired depth: X0.5cm, X5.15cm, X15.30cm, X30.60cm, X60.100cm, X100.200cm

# --- Paths ---
cov_folder <- file.path(HomeDir, "Data/data_in/soil_covariates_aligned")
model_file <- file.path(HomeDir, paste0("mod.cubist.", soil_property), "models", paste0(depth_name, ".rds"))
output_tif <- file.path(HomeDir, paste0("mod.cubist.", soil_property), "preds", paste0(soil_property, "_pred_", depth_name, ".tif"))

# --- Load your covariate names (should match your model!) ---
load(file.path(HomeDir, "Data/data_out/RData", paste0(soil_property, "_covs_regression.RData")))
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

# --- Load the trained Cubist model for a given depth ---
model <- readRDS(model_file)
cat("Loaded model:", model_file, "\n")

# --- Define a wrapper for Cubist prediction ---
predict_fun <- function(model, data, ...) {
  predict(model, newdata = data)
}

# --- Predict SOC for all pixels ---
cat("Predicting", soil_property, "map for depth:", depth_name, "...\n")
soc_pred <- terra::predict(cov_stack, model, fun = predict_fun, na.rm = TRUE, progress = "text")

# --- Back-transform if model was trained on log scale ---
# Uncomment the next line if use_log_transform = TRUE in training
# soc_pred <- exp(soc_pred)

# --- Save and plot the result ---
writeRaster(soc_pred, output_tif, overwrite = TRUE)
plot(soc_pred, main = paste("Predicted", soil_property, depth_name))

cat("Map saved to:", output_tif, "\n")
