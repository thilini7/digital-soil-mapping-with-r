# --- Libraries and setup ---
suppressPackageStartupMessages({
  library(ranger)
  library(caret)
  library(dplyr)
  library(doParallel)
})

# Reproducibility
set.seed(42)

# --- Paths ---
setwd("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/")
HomeDir <- getwd()
mod.type <- "mod.regression.OC"

# Optional helpers
if (file.exists("R/3_3_func.R")) {
  source("R/3_3_func.R")
}

# Create output dirs
dir.create(file.path(mod.type), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "cv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "metrics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "preds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mod.type, "importance"), recursive = TRUE, showWarnings = FALSE)

# --- Load your regression-ready data (must provide df_conc and cov_names) ---
load("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/RData/OC_covs_regression.RData")

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
num_trees      <- 200      # used for both tuning fit and final refit

ctrl <- trainControl(
  method = "cv",
  number = 5,
  allowParallel = TRUE
)

# caret's ranger: tune mtry and min.node.size (splitrule fixed to "variance" for regression)
# Ensure mtry doesn't exceed number of covariates
n_covariates <- length(cov_names)
max_mtry <- min(8, n_covariates)
mtry_values <- unique(pmin(c(2, 5, 8), max_mtry))
mtry_values <- mtry_values[mtry_values > 0]  # Keep only positive values

rg.tuneGrid <- expand.grid(
  mtry = mtry_values,
  splitrule = "variance",
  min.node.size = c(5, 10)
)

cat("Number of covariates:", n_covariates, "\n")
cat("mtry values to tune:", paste(mtry_values, collapse = ", "), "\n")

# --- Training loop over depth columns ---
for (cc in seq(start_col, end_col)) {
  depth_name <- names(df_conc)[cc]
  cat("\n==============================\n")
  cat(">>> Working on", depth_name, "\n")
  
  # --- Prepare data ---
  df <- df_conc %>%
    dplyr::select(all_of(c(depth_name, cov_names))) %>%
    filter(complete.cases(.)) %>%
    filter(.data[[depth_name]] > 0, .data[[depth_name]] < 200)
  
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
  
  # --- Tuning on training set only ---
  t0 <- Sys.time()
  formulaString <- as.formula(paste0(depth_name, " ~ ", paste(cov_names, collapse = " + ")))
  
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
            file = file.path(HomeDir, mod.type, "cv", paste0("cv_grid_", depth_name, ".csv")),
            row.names = FALSE)
  saveRDS(t.mrg$bestTune,
          file = file.path(HomeDir, mod.type, "cv", paste0("bestTune_", depth_name, ".rds")))
  
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
            file = file.path(HomeDir, mod.type, "metrics", paste0("metrics_", depth_name, ".csv")),
            row.names = FALSE)
  
  # Save test predictions (paired)
  preds_out <- data.frame(
    depth = depth_name,
    obs   = obs_test,
    pred  = pred_test
  )
  write.csv(preds_out,
            file = file.path(HomeDir, mod.type, "preds", paste0("pred_test_", depth_name, ".csv")),
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
          file = file.path(HomeDir, mod.type, "models", paste0(depth_name, ".rds")))
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
              file = file.path(HomeDir, mod.type, "importance", paste0("importance_", depth_name, ".csv")),
              row.names = FALSE)
  }
}

cat("\nAll done!\n")
