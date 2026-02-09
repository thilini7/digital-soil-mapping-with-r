library(terra)
library(Cubist)

# --- Configuration ---
HomeDir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
ModelsDir <- file.path(HomeDir, "Models")
soil_property <- "Organic_Carbon"
depth_name <- "X0.5cm"  # Change to desired depth: X0.5cm, X5.15cm, X15.30cm, X30.60cm, X60.100cm, X100.200cm

# --- Paths ---
cov_folder <- file.path(HomeDir, "Data/data_in/soil_covariates_aligned_v2")
model_file <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "models", paste0(depth_name, ".rds"))
output_tif <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "preds", paste0(soil_property, "_pred_", depth_name, ".tif"))

# --- Load your covariate names (for reference) ---
load(file.path(HomeDir, "Data/data_out/RData", paste0(soil_property, "_covs_regression.RData")))
# Note: The model may use fewer covariates due to NZV/correlation filtering

# --- Load the trained Cubist model for a given depth ---
model <- readRDS(model_file)
cat("Loaded model:", model_file, "\n")

# --- Get the covariates actually used by the model ---
# Cubist stores variable names in the model object
model_cov_names <- model$vars$all
if (is.null(model_cov_names)) {
  # Fallback: try to get from coefficients or names
  model_cov_names <- names(model$coefficients)
}
if (is.null(model_cov_names)) {
  # Last resort: use cov_names from RData (may cause issues)
  model_cov_names <- cov_names
  warning("Could not extract covariate names from model - using full cov_names")
}

cat("Model uses", length(model_cov_names), "covariates:\n")
print(model_cov_names)

# --- List and match raster files to MODEL's covariate names ---
cov_files <- list.files(cov_folder, pattern = "\\.tif$", full.names = TRUE)
base_names <- tools::file_path_sans_ext(basename(cov_files))

# Check for missing rasters
missing_covs <- setdiff(model_cov_names, base_names)
if (length(missing_covs) > 0) {
  cat("WARNING: Missing raster files for these model covariates:\n")
  print(missing_covs)
  stop("Cannot predict - missing covariate rasters!")
}

ordered_files <- cov_files[match(model_cov_names, base_names)]

# --- Stack rasters in model order and assign names ---
cov_stack <- rast(ordered_files)
names(cov_stack) <- model_cov_names

cat("Covariate stack ready with", nlyr(cov_stack), "layers\n")

# --- Check for NA pixels in covariates ---
na_count <- global(is.na(cov_stack), sum)
cat("NA pixels per covariate:\n")
print(na_count[na_count$sum > 0, , drop = FALSE])

# --- Identify which covariate has the most NAs (causing extra white spots) ---
min_na <- min(na_count$sum)
max_na <- max(na_count$sum)
extra_na_covs <- rownames(na_count)[na_count$sum > min_na * 1.1]  # >10% more NAs than minimum

if (length(extra_na_covs) > 0) {
  cat("\n*** WARNING: These covariates have more NAs than others (may cause white spots):\n")
  for (cov in extra_na_covs) {
    extra <- na_count[cov, "sum"] - min_na
    cat("  ", cov, ": +", format(extra, big.mark = ","), "extra NA pixels\n")
  }
  cat("\nConsider checking alignment of these rasters or excluding them.\n\n")
}

# --- Define a wrapper for Cubist prediction ---
predict_fun <- function(model, data, ...) {
  predict(model, newdata = data)
}

# --- Predict SOC for all pixels (OPTIMIZED) ---
cat("Predicting", soil_property, "map for depth:", depth_name, "...\n")
cat("Raster dimensions:", nrow(cov_stack), "x", ncol(cov_stack), "=", 
    format(ncell(cov_stack), big.mark = ","), "pixels\n")

t0 <- Sys.time()

# Note: Parallel prediction with Cubist can be problematic
# Using single-core for reliability (still faster than base R loops)
cat("Running prediction (single-core for Cubist compatibility)...\n")

soc_pred <- terra::predict(
  cov_stack, 
  model, 
  fun = predict_fun, 
  na.rm = TRUE, 
  progress = "text"
)

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat("Prediction completed in", elapsed, "minutes\n")

# --- Back-transform if model was trained on log scale ---
# Uncomment the next line if use_log_transform = TRUE in training
# soc_pred <- exp(soc_pred)

# --- Save and plot the result ---
writeRaster(soc_pred, output_tif, overwrite = TRUE)
plot(soc_pred, main = paste("Predicted", soil_property, depth_name))

cat("Map saved to:", output_tif, "\n")
