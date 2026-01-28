# --- Libraries and setup ---
suppressPackageStartupMessages({
  library(ranger)
  library(caret)
  library(dplyr)
  library(doParallel)
})

# Reproducibility
set.seed(42)

# =============================================================================
# FILE PATHS - CENTRALIZED CONFIGURATION
# =============================================================================
HomeDir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
setwd(HomeDir)

# INPUT: Change soil_property to match your data file
soil_property <- "EC"  # Options: pH, OC, BD, CEC, EC, Clay, etc.

# Input data paths
data_input_path <- file.path(HomeDir, "Data/data_out/RData", 
                             paste0(soil_property, "_covs_regression.RData"))
soil_covariates_csv <- file.path(HomeDir, "Data/data_out/Soil_data_with_covariates",
                                 paste0(soil_property, "_with_covariates_new.csv"))

# Output model directory name
mod_type <- paste0("mod.regression.", soil_property)

# Output subdirectories
dir_models <- file.path(HomeDir, mod_type, "models")
dir_cv <- file.path(HomeDir, mod_type, "cv")
dir_metrics <- file.path(HomeDir, mod_type, "metrics")
dir_preds <- file.path(HomeDir, mod_type, "preds")
dir_importance <- file.path(HomeDir, mod_type, "importance")

# Print configuration
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("RANDOM FOREST MODEL TRAINING - CONFIGURATION\n")
cat("Soil Property: ", soil_property, "\n", sep = "")
cat("Input RData: ", basename(data_input_path), "\n", sep = "")
cat("Output Directory: ", file.path(HomeDir, mod_type), "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

# Optional helpers
if (file.exists("R/3_3_func.R")) {
  source("R/3_3_func.R")
}

# Create output directories
for (dir in list(dir_models, dir_cv, dir_metrics, dir_preds, dir_importance)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# --- Load your regression-ready data (must provide df_conc and cov_names) ---
load(data_input_path)

# Sanity checks
stopifnot(exists("df_conc"), exists("cov_names"))
depth_cols <- c("X0.5cm", "X5.15cm", "X15.30cm", "X30.60cm", "X60.100cm", "X100.200cm")
missing_depths <- setdiff(depth_cols, names(df_conc))
if (length(missing_depths)) {
  stop("These depth columns are missing in df_conc: ", paste(missing_depths, collapse = ", "))
}

start_col <- which(names(df_conc) == depth_cols[1])
end_col   <- which(names(df_conc) == depth_cols[length(depth_cols)])
cat("start_col =", start_col, "end_col =", end_col, "\n")

# --- Parallel setup ---
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
on.exit({
  try(stopCluster(cl), silent = TRUE)
  registerDoSEQ()
}, add = TRUE)
cat("Running on", n_cores, "cores\n")

# --- Controls / knobs ---
retrain_on_all <- TRUE     # set FALSE to keep the train-only model
num_trees      <- 500      # increased for better accuracy (was 200)
use_log_transform <- TRUE  # set TRUE to use log-transformed OC values
use_spatial_cv <- TRUE     # set TRUE for spatial cross-validation (recommended for soil mapping)
add_interactions <- TRUE   # set TRUE to add feature interactions

# Install blockCV if needed for spatial CV
if (use_spatial_cv) {
  if (!requireNamespace("blockCV", quietly = TRUE)) {
    install.packages("blockCV")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    install.packages("sf")
  }
  library(blockCV)
  library(sf)
}

# Standard CV control (used if spatial CV is disabled or as fallback)
ctrl <- trainControl(
  method = "cv",
  number = 10,              # increased to 10-fold CV for more stable estimates
  allowParallel = TRUE,
  savePredictions = "final" # save predictions for analysis
)

# caret's ranger: tune mtry and min.node.size (splitrule fixed to "variance" for regression)
# Expanded tuning grid for better optimization
n_covariates <- length(cov_names)
# Test mtry from sqrt(p) to p/2 (common range for RF)
mtry_values <- unique(c(
  floor(sqrt(n_covariates)),      # sqrt(p) - default
  floor(n_covariates / 3),        # p/3
  floor(n_covariates / 2),        # p/2
  min(10, n_covariates)           # fixed value
))
mtry_values <- sort(mtry_values[mtry_values > 0])

rg.tuneGrid <- expand.grid(
  mtry = mtry_values,
  splitrule = "variance",
  min.node.size = c(3, 5, 10)     # added smaller node size for more complex trees
)

cat("Number of covariates:", n_covariates, "\n")
cat("mtry values to tune:", paste(mtry_values, collapse = ", "), "\n")
cat("Using log-transformed data:", use_log_transform, "\n")

# Select data based on log-transform option
if (use_log_transform && exists("df_conc_log")) {
  df_model <- df_conc_log
  cat(">>> Using LOG-TRANSFORMED OC values\n")
} else {
  df_model <- df_conc
  cat(">>> Using ORIGINAL OC values\n")
}

# --- Feature Engineering: Add interaction terms ---
if (add_interactions) {
  cat("\n>>> Adding feature interactions...\n")
  
  # Check which covariates exist in the data
  has_twi <- "Relief_TWI_3s" %in% cov_names
  has_rf <- "Annual_Mean_RF_90m" %in% cov_names
  has_temp <- "Annual_Max_Temp_90m" %in% cov_names
  has_dem <- "Relief_Dems_3s_Mosaic" %in% cov_names
  has_slope <- "Relief_Slope_Perc" %in% cov_names
  
  # TWI x Rainfall interaction (moisture accumulation potential)
  if (has_twi && has_rf) {
    df_model$TWI_x_Rainfall <- df_model$Relief_TWI_3s * df_model$Annual_Mean_RF_90m / 1000
    cov_names <- c(cov_names, "TWI_x_Rainfall")
    cat("  + Added TWI_x_Rainfall\n")
  }
  
  # Elevation x Temperature interaction (climate gradient)
  if (has_dem && has_temp) {
    df_model$DEM_x_Temp <- df_model$Relief_Dems_3s_Mosaic * df_model$Annual_Max_Temp_90m / 1000
    cov_names <- c(cov_names, "DEM_x_Temp")
    cat("  + Added DEM_x_Temp\n")
  }
  
  # Slope x Rainfall interaction (erosion potential)
  if (has_slope && has_rf) {
    df_model$Slope_x_Rainfall <- df_model$Relief_Slope_Perc * df_model$Annual_Mean_RF_90m / 1000
    cov_names <- c(cov_names, "Slope_x_Rainfall")
    cat("  + Added Slope_x_Rainfall\n")
  }
  
  cat("Total covariates after feature engineering:", length(cov_names), "\n")
}

# --- Training loop over depth columns ---
for (cc in seq(start_col, end_col)) {
  depth_name <- names(df_model)[cc]
  cat("\n==============================\n")
  cat(">>> Working on", depth_name, "\n")
  
  # --- Prepare data ---
  df <- df_model %>%
    dplyr::select(all_of(c(depth_name, cov_names))) %>%
    filter(complete.cases(.)) %>%
    filter(is.finite(.data[[depth_name]]))  # handle log-transformed data
  
  n_all <- nrow(df)
  cat("Records after cleaning:", n_all, "\n")
  if (n_all < 50) {
    cat("Skipping", depth_name, "- too few records for split/CV.\n")
    next
  }
  
  # --- 80/20 split (stratified on response) ---
  idx_train <- caret::createDataPartition(y = df[[depth_name]], p = 0.80, list = FALSE)
  df_train  <- df[idx_train, , drop = FALSE]
  df_test   <- df[-idx_train, , drop = FALSE]
  cat("Train:", nrow(df_train), " Test:", nrow(df_test), "\n")
  
  # --- Set up cross-validation (spatial or random) ---
  t0 <- Sys.time()
  formulaString <- as.formula(paste0(depth_name, " ~ ", paste(cov_names, collapse = " + ")))
  
  if (use_spatial_cv) {
    cat("  Setting up spatial block CV...\n")
    tryCatch({
      # Need coordinates - load from original data
      orig_df <- read.csv(soil_covariates_csv)
      train_coords <- orig_df[idx_train, c("Longitude", "Latitude")]
      train_coords <- train_coords[complete.cases(train_coords), ]
      
      # Match rows with df_train
      if (nrow(train_coords) >= nrow(df_train)) {
        train_coords <- train_coords[1:nrow(df_train), ]
      }
      
      # Create sf object
      train_sf <- st_as_sf(cbind(df_train, train_coords), 
                           coords = c("Longitude", "Latitude"), 
                           crs = 4326)
      
      # Create spatial blocks (~50km blocks for NSW scale)
      sb <- cv_spatial(
        x = train_sf,
        k = 10,                    # 10-fold spatial CV
        size = 50000,              # 50km block size
        selection = "random",
        iteration = 50,
        progress = FALSE
      )
      
      # Convert to caret-compatible format
      ctrl_spatial <- trainControl(
        method = "cv",
        number = 10,
        index = sb$folds_ids,
        allowParallel = TRUE,
        savePredictions = "final"
      )
      ctrl <- ctrl_spatial
      cat("  ✓ Using spatial block CV (50km blocks)\n")
    }, error = function(e) {
      cat("  ✗ Spatial CV failed:", conditionMessage(e), "\n")
      cat("  → Falling back to random CV\n")
    })
  }
  
  # --- Train model with tuning ---
  cat("  Training Random Forest model...\n")
  t.mrg <- caret::train(
    formulaString,
    data = df_train,
    method = "ranger",
    trControl = ctrl,
    tuneGrid = rg.tuneGrid,
    num.trees = num_trees,          # passed via ... to ranger::ranger
    importance = "impurity",
    quantreg = TRUE,                # quantile RF enabled (point preds still fine)
    keep.inbag = TRUE
  )
  cat("CV done in", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")
  
  # Save CV grid + best tune
  write.csv(t.mrg$results,
            file = file.path(dir_cv, paste0("cv_grid_", depth_name, ".csv")),
            row.names = FALSE)
  saveRDS(t.mrg$bestTune,
          file = file.path(dir_cv, paste0("bestTune_", depth_name, ".rds")))
  
  # --- Evaluate on 20% hold-out ---
  pred_test <- predict(t.mrg, newdata = df_test)
  obs_test  <- df_test[[depth_name]]
  
  RMSE <- sqrt(mean((pred_test - obs_test)^2))
  MAE  <- mean(abs(pred_test - obs_test))
  R2   <- {
    cxy <- cov(pred_test, obs_test)
    if (is.na(cxy) || sd(pred_test) == 0 || sd(obs_test) == 0) NA_real_ else cor(pred_test, obs_test)^2
  }
  
  metrics <- data.frame(
    depth   = depth_name,
    n_train = nrow(df_train),
    n_test  = nrow(df_test),
    RMSE    = RMSE,
    MAE     = MAE,
    R2      = R2,
    stringsAsFactors = FALSE
  )
  write.csv(metrics,
            file = file.path(dir_metrics, paste0("metrics_", depth_name, ".csv")),
            row.names = FALSE)
  
  # Save test predictions (paired)
  preds_out <- data.frame(
    depth = depth_name,
    obs   = obs_test,
    pred  = pred_test
  )
  write.csv(preds_out,
            file = file.path(dir_preds, paste0("pred_test_", depth_name, ".csv")),
            row.names = FALSE)
  
  # --- Final model (either refit on all data or keep train-only) ---
  if (retrain_on_all) {
    bt <- t.mrg$bestTune
    final_model <- ranger::ranger(
      dependent.variable.name = depth_name,
      data = df,                            # ALL data for highest capacity
      mtry = bt$mtry,
      splitrule = bt$splitrule,
      min.node.size = bt$min.node.size,
      num.trees = num_trees,
      importance = "impurity",
      quantreg = TRUE,
      keep.inbag = TRUE
    )
  } else {
    final_model <- t.mrg$finalModel        # trained on train only
  }
  
  # Save model
  saveRDS(final_model,
          file = file.path(dir_models, paste0(depth_name, ".rds")))
  cat("Saved model for", depth_name, "\n")
  
  # --- Variable importance export (from final model if available) ---
  vi <- tryCatch({
    as.data.frame(final_model$variable.importance) %>%
      tibble::rownames_to_column(var = "covariate") %>%
      dplyr::rename(importance = `final_model$variable.importance`) %>%
      arrange(desc(importance))
  }, error = function(e) NULL)
  
  if (!is.null(vi)) {
    write.csv(vi,
              file = file.path(dir_importance, paste0("importance_", depth_name, ".csv")),
              row.names = FALSE)
  }
}

cat("\nAll done!\n")
