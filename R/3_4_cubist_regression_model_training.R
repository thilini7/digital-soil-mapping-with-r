# =============================================================================
# Cubist Regression Model Training for Digital Soil Mapping
# =============================================================================
# This script trains Cubist regression models for soil property prediction
# Cubist is a rule-based regression model that combines decision trees with 
# linear regression models at the terminal nodes
# =============================================================================

# --- Libraries and setup ---
suppressPackageStartupMessages({
  library(Cubist)      # Main Cubist package
  library(caret)       # For model training framework
  library(dplyr)       # Data manipulation
  library(doParallel)  # Parallel processing
  library(tibble)      # For rownames_to_column
})

# Reproducibility
set.seed(42)

# --- Paths ---
setwd("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/")
HomeDir <- getwd()
mod.type <- "mod.cubist.Phosphorus"  # Change to your target property

# Optional helpers
if (file.exists("R/3_3_func.R")) {
  source("R/3_3_func.R")
}

# Create output directories
dir.create(file.path(mod.type), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "cv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "metrics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "preds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "importance"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "rules"), recursive = TRUE, showWarnings = FALSE)

# --- Load your regression-ready data ---
# Must provide df_conc (data frame) and cov_names (character vector of covariate names)
load("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/RData/Phosphorus_covs_regression.RData")

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
fast_mode <- TRUE           # Set TRUE for faster training (reduced tuning grid)
use_log_transform <- TRUE   # Set TRUE for log-transformed values (common for OC)
use_spatial_cv <- FALSE     # Set FALSE for faster training (spatial CV is slow)
add_interactions <- FALSE   # Set FALSE for faster training

# Cubist-specific parameters
# committees: number of boosting iterations (like ensemble of rule-based models)
# neighbors: number of nearest neighbors for prediction adjustment (0 = no adjustment)
if (fast_mode) {
  committees_values <- c(1, 10, 20)    # Reduced grid for speed
  neighbors_values <- c(0, 5)           # Reduced grid for speed
  cv_folds <- 5                         # Fewer folds
} else {
  committees_values <- c(1, 5, 10, 20, 50)
  neighbors_values <- c(0, 3, 5, 9)
  cv_folds <- 10
}

# Install spatial CV packages if needed
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
  number = cv_folds,        # CV folds (5 for fast mode, 10 otherwise)
  allowParallel = TRUE,
  savePredictions = "final"
)

# Cubist tuning grid
# committees: like boosting, each subsequent model focuses on residuals
# neighbors: applies instance-based correction using training data neighbors
cubist_tuneGrid <- expand.grid(
  committees = committees_values,
  neighbors = neighbors_values
)

cat("Fast mode:", fast_mode, "\n")
cat("Number of covariates:", length(cov_names), "\n")
cat("Committees to tune:", paste(committees_values, collapse = ", "), "\n")
cat("Neighbors to tune:", paste(neighbors_values, collapse = ", "), "\n")
cat("CV folds:", cv_folds, "\n")
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

# --- Function to extract Cubist rules ---
extract_cubist_rules <- function(model, depth_name, output_dir) {
  tryCatch({
    # Get the rules summary
    rules_summary <- summary(model)
    
    # Capture the printed output
    rules_text <- capture.output(print(rules_summary))
    
    # Write rules to text file
    writeLines(rules_text, 
               file.path(output_dir, paste0("rules_", depth_name, ".txt")))
    
    cat("  ✓ Saved Cubist rules to text file\n")
    return(TRUE)
  }, error = function(e) {
    cat("  ✗ Failed to extract rules:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# --- Function to get variable importance from Cubist ---
get_cubist_importance <- function(model) {
  tryCatch({
    # Cubist provides usage statistics
    usage <- varImp(model)
    
    if (!is.null(usage)) {
      vi_df <- data.frame(
        covariate = rownames(usage$importance),
        importance = usage$importance$Overall
      ) %>%
        arrange(desc(importance))
      return(vi_df)
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

# =============================================================================
# Training loop over depth columns
# =============================================================================
all_metrics <- list()

for (cc in seq(start_col, end_col)) {
  depth_name <- names(df_model)[cc]
  cat("\n==============================\n")
  cat(">>> Working on", depth_name, "\n")
  
  # --- Prepare data ---
  df <- df_model %>%
    dplyr::select(all_of(c(depth_name, cov_names))) %>%
    filter(complete.cases(.)) %>%
    filter(is.finite(.data[[depth_name]]))
  
  n_all <- nrow(df)
  cat("Records after cleaning:", n_all, "\n")
  if (n_all < 50) {
    cat("Skipping", depth_name, "- too few records for split/CV.\n")
    next
  }
  
  # --- 80/20 split ---
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
      # Check if coordinates exist in df_model
      if (all(c("Longitude", "Latitude") %in% names(df_model))) {
        # Use coordinates from df_model directly (preferred method)
        df_with_coords <- df_model %>%
          dplyr::select(all_of(c(depth_name, cov_names, "Longitude", "Latitude"))) %>%
          filter(complete.cases(.)) %>%
          filter(is.finite(.data[[depth_name]]))
        
        train_coords <- df_with_coords[idx_train, c("Longitude", "Latitude")]
      } else {
        # Load from original CSV and match by row indices BEFORE filtering
        orig_df <- read.csv("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/Soil_data_with_covariates/OC_with_covariates_new.csv")
        
        # Add coordinates to df_model first, then filter together
        df_with_coords <- df_model
        if (nrow(orig_df) == nrow(df_model)) {
          df_with_coords$Longitude <- orig_df$Longitude
          df_with_coords$Latitude <- orig_df$Latitude
        } else {
          stop("Row count mismatch between df_model and original CSV")
        }
        
        # Filter with coordinates
        df_with_coords <- df_with_coords %>%
          dplyr::select(all_of(c(depth_name, cov_names, "Longitude", "Latitude"))) %>%
          filter(complete.cases(.)) %>%
          filter(is.finite(.data[[depth_name]]))
        
        train_coords <- df_with_coords[idx_train, c("Longitude", "Latitude")]
      }
      
      # Verify row counts match
      if (nrow(train_coords) != nrow(df_train)) {
        stop(paste("Row mismatch: train_coords has", nrow(train_coords), 
                   "rows but df_train has", nrow(df_train), "rows"))
      }
      
      # Create sf object
      train_sf <- st_as_sf(cbind(df_train, train_coords), 
                           coords = c("Longitude", "Latitude"), 
                           crs = 4326)
      
      # Create spatial blocks (~50km blocks for NSW scale)
      sb <- cv_spatial(
        x = train_sf,
        k = cv_folds,
        size = 50000,
        selection = "random",
        iteration = if (fast_mode) 10 else 50,
        progress = FALSE
      )
      
      # Convert to caret-compatible format
      ctrl_spatial <- trainControl(
        method = "cv",
        number = cv_folds,
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
  
  # --- Train Cubist model with tuning ---
  cat("  Training Cubist regression model...\n")
  t.cubist <- caret::train(
    formulaString,
    data = df_train,
    method = "cubist",
    trControl = ctrl,
    tuneGrid = cubist_tuneGrid
  )
  cat("CV done in", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")
  
  # Print best parameters
  cat("  Best parameters: committees =", t.cubist$bestTune$committees, 
      ", neighbors =", t.cubist$bestTune$neighbors, "\n")
  
  # Save CV grid + best tune
  write.csv(t.cubist$results,
            file = file.path(HomeDir, mod.type, "cv", paste0("cv_grid_", depth_name, ".csv")),
            row.names = FALSE)
  saveRDS(t.cubist$bestTune,
          file = file.path(HomeDir, mod.type, "cv", paste0("bestTune_", depth_name, ".rds")))
  
  # --- Evaluate on 20% hold-out ---
  pred_test <- predict(t.cubist, newdata = df_test)
  obs_test  <- df_test[[depth_name]]
  
  # Calculate metrics
  RMSE <- sqrt(mean((pred_test - obs_test)^2))
  MAE  <- mean(abs(pred_test - obs_test))
  R2   <- {
    cxy <- cov(pred_test, obs_test)
    if (is.na(cxy) || sd(pred_test) == 0 || sd(obs_test) == 0) NA_real_ else cor(pred_test, obs_test)^2
  }
  
  # Bias
  bias <- mean(pred_test - obs_test)
  
  # Lin's Concordance Correlation Coefficient (CCC)
  CCC <- tryCatch({
    mean_obs <- mean(obs_test)
    mean_pred <- mean(pred_test)
    var_obs <- var(obs_test)
    var_pred <- var(pred_test)
    cov_op <- cov(obs_test, pred_test)
    2 * cov_op / (var_obs + var_pred + (mean_obs - mean_pred)^2)
  }, error = function(e) NA_real_)
  
  metrics <- data.frame(
    depth      = depth_name,
    n_train    = nrow(df_train),
    n_test     = nrow(df_test),
    committees = t.cubist$bestTune$committees,
    neighbors  = t.cubist$bestTune$neighbors,
    RMSE       = round(RMSE, 4),
    MAE        = round(MAE, 4),
    R2         = round(R2, 4),
    Bias       = round(bias, 4),
    CCC        = round(CCC, 4),
    stringsAsFactors = FALSE
  )
  
  all_metrics[[depth_name]] <- metrics
  
  write.csv(metrics,
            file = file.path(HomeDir, mod.type, "metrics", paste0("metrics_", depth_name, ".csv")),
            row.names = FALSE)
  
  cat("  Metrics: RMSE =", round(RMSE, 4), 
      ", MAE =", round(MAE, 4), 
      ", R² =", round(R2, 4), "\n")
  
  # Save test predictions (paired)
  preds_out <- data.frame(
    depth    = depth_name,
    observed = obs_test,
    predicted = pred_test,
    residual = obs_test - pred_test
  )
  write.csv(preds_out,
            file = file.path(HomeDir, mod.type, "preds", paste0("pred_test_", depth_name, ".csv")),
            row.names = FALSE)
  
  # --- Final model trained on ALL data ---
  cat("  Training final model on all data...\n")
  bt <- t.cubist$bestTune
  
  # Prepare x and y for cubist() function
  x_all <- df[, cov_names, drop = FALSE]
  y_all <- df[[depth_name]]
  
  final_model <- cubist(
    x = x_all,
    y = y_all,
    committees = bt$committees
  )
  
  # Apply neighbors correction if specified
  if (bt$neighbors > 0) {
    final_model$neighbors <- bt$neighbors
    attr(final_model, "neighbors") <- bt$neighbors
  }
  
  # Save model
  saveRDS(final_model,
          file = file.path(HomeDir, mod.type, "models", paste0(depth_name, ".rds")))
  cat("  ✓ Saved model for", depth_name, "\n")
  
  # --- Extract and save Cubist rules ---
  extract_cubist_rules(final_model, depth_name, 
                       file.path(HomeDir, mod.type, "rules"))
  
  # --- Variable importance export ---
  vi <- get_cubist_importance(t.cubist)
  
  if (!is.null(vi)) {
    write.csv(vi,
              file = file.path(HomeDir, mod.type, "importance", paste0("importance_", depth_name, ".csv")),
              row.names = FALSE)
    
    # Print top 10 important variables
    cat("  Top 10 important covariates:\n")
    top10 <- head(vi, 10)
    for (i in 1:nrow(top10)) {
      cat("    ", i, ". ", top10$covariate[i], " (", round(top10$importance[i], 2), ")\n", sep = "")
    }
  }
}

# =============================================================================
# Summary of all models
# =============================================================================
cat("\n==============================\n")
cat(">>> SUMMARY OF ALL CUBIST MODELS\n")
cat("==============================\n")

if (length(all_metrics) > 0) {
  summary_df <- do.call(rbind, all_metrics)
  
  cat("\nPerformance metrics across depths:\n")
  print(summary_df)
  
  # Save combined metrics
  write.csv(summary_df,
            file = file.path(HomeDir, mod.type, "metrics", "all_depths_metrics.csv"),
            row.names = FALSE)
  
  cat("\nMean performance:\n")
  cat("  Mean RMSE:", round(mean(summary_df$RMSE, na.rm = TRUE), 4), "\n")
  cat("  Mean MAE:", round(mean(summary_df$MAE, na.rm = TRUE), 4), "\n")
  cat("  Mean R²:", round(mean(summary_df$R2, na.rm = TRUE), 4), "\n")
  cat("  Mean CCC:", round(mean(summary_df$CCC, na.rm = TRUE), 4), "\n")
}

cat("\n==============================\n")
cat("All Cubist models trained successfully!\n")
cat("Output directory:", file.path(HomeDir, mod.type), "\n")
cat("==============================\n")

# =============================================================================
# Cubist Model Prediction Function (for use in mapping)
# =============================================================================
# Usage example for prediction:
# 
# # Load the trained model
# model <- readRDS("mod.cubist.OC/models/X0.5cm.rds")
# 
# # Load covariate raster stack
# covariates <- terra::rast("path/to/covariates.tif")
# 
# # Predict using the model
# # Note: Cubist with neighbors requires training data for prediction
# prediction <- predict(covariates, model, 
#                       neighbors = model$neighbors,
#                       na.rm = TRUE)
# 
# # If using log-transformed model, back-transform
# if (use_log_transform) {
#   prediction <- exp(prediction)
# }
# =============================================================================
