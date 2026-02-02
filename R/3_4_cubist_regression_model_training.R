# =============================================================================
# Cubist Regression Model Training for Digital Soil Mapping (WITH FULL LOGGING)
# =============================================================================

suppressPackageStartupMessages({
  library(Cubist)
  library(caret)
  library(dplyr)
  library(doParallel)
  library(tibble)
})

set.seed(42)

# =============================================================================
# PATHS
# =============================================================================
HomeDir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
setwd(HomeDir)

soil_property <- "Organic_Carbon"

data_input_path <- file.path(
  HomeDir, "Data/data_out/RData",
  paste0(soil_property, "_covs_regression.RData")
)

soil_covariates_csv <- file.path(
  HomeDir, "Data/data_out/Soil_data_with_covariates",
  paste0(soil_property, "_with_covariates_new.csv")
)

mod_type <- paste0("mod.cubist.", soil_property)

dirs <- list(
  models     = file.path(HomeDir, mod_type, "models"),
  cv         = file.path(HomeDir, mod_type, "cv"),
  metrics    = file.path(HomeDir, mod_type, "metrics"),
  preds      = file.path(HomeDir, mod_type, "preds"),
  importance = file.path(HomeDir, mod_type, "importance"),
  rules      = file.path(HomeDir, mod_type, "rules"),
  logs       = file.path(HomeDir, mod_type, "logs")
)

lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================
load(data_input_path)
stopifnot(exists("df_conc"), exists("cov_names"))

depth_cols <- c("X0.5cm", "X5.15cm", "X15.30cm",
                "X30.60cm", "X60.100cm", "X100.200cm")

start_col <- which(names(df_conc) == depth_cols[1])
end_col   <- which(names(df_conc) == depth_cols[length(depth_cols)])

# =============================================================================
# PARALLEL
# =============================================================================
cl <- makePSOCKcluster(max(1, parallel::detectCores() - 1))
registerDoParallel(cl)
on.exit({ try(stopCluster(cl), silent = TRUE); registerDoSEQ() }, add = TRUE)

# =============================================================================
# SETTINGS
# =============================================================================
use_log_transform <- TRUE
cv_folds <- 10
cv_repeats <- 3  # Repeated CV for more stable estimates

# Expanded tuning grid for better optimization
committees_values <- c(1, 5, 10, 20, 50, 100)
neighbors_values  <- c(0, 1, 3, 5, 7, 9)  # Added 1 and 7

# Relaxed correlation cutoff (0.95 instead of 0.9) to keep more predictors
correlation_cutoff <- 0.95

cubist_grid <- expand.grid(
  committees = committees_values,
  neighbors  = neighbors_values
)

# =============================================================================
# HELPERS
# =============================================================================

log_line <- function(logfile, ...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = logfile, append = TRUE)
}

# Cubist: variable usage (%) from a caret Cubist model (robust)
get_cubist_usage_pct <- function(caret_model) {
  u <- caret_model$finalModel$usage
  if (is.null(u)) return(NULL)
  
  # Handle data.frame format (most common for Cubist)
  if (is.data.frame(u)) {
    # Get the first numeric column for usage values
    num_cols <- names(u)[sapply(u, is.numeric)]
    if (length(num_cols) == 0) return(NULL)
    
    # Prefer "Conditions" or "Model" column if available
    use_col <- if ("Conditions" %in% num_cols) "Conditions" 
               else if ("Model" %in% num_cols) "Model" 
               else num_cols[1]
    
    usage <- u[[use_col]]
    
    # Get covariate names - check for "Variable" column first, then rownames
    if ("Variable" %in% names(u)) {
      cov_names <- as.character(u$Variable)
    } else {
      cov_names <- rownames(u)
    }
    
    if (is.null(cov_names) || length(cov_names) != length(usage)) {
      cov_names <- paste0("Var_", seq_along(usage))
    }
    names(usage) <- cov_names
    
  } else if (is.matrix(u)) {
    usage <- u[, 1]
    cov_names <- rownames(u)
    if (is.null(cov_names)) cov_names <- paste0("Var_", seq_along(usage))
    names(usage) <- cov_names
    usage <- as.numeric(usage)
    
  } else {
    usage <- u
    if (is.null(names(usage))) names(usage) <- paste0("Var_", seq_along(usage))
  }
  
  # Filter to finite values

  usage <- usage[is.finite(usage)]
  if (length(usage) == 0) return(NULL)
  
  # Sort by usage

  usage <- sort(usage, decreasing = TRUE)
  
  # Calculate percentage (handle zero sum case)
  total <- sum(usage)
  if (total == 0 || !is.finite(total)) {
    pct <- rep(0, length(usage))
  } else {
    pct <- 100 * usage / total
  }
  
  data.frame(
    covariate = names(usage),
    usage = as.numeric(usage),
    pct = as.numeric(pct),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# METRICS STORAGE
# =============================================================================
all_metrics <- list()
top_cov_summary <- list()

# =============================================================================
# LOOP OVER DEPTHS
# =============================================================================
for (cc in seq(start_col, end_col)) {
  
  depth_name <- names(df_conc)[cc]
  logfile <- file.path(dirs$logs, paste0("log_", depth_name, ".txt"))
  writeLines(character(0), con = logfile)
  
  log_line(logfile, "==============================")
  log_line(logfile, ">>> Depth: ", depth_name)
  log_line(logfile, "Timestamp: ", as.character(Sys.time()))
  log_line(logfile, "use_log_transform: ", use_log_transform)
  log_line(logfile, "cv_folds: ", cv_folds)
  log_line(logfile, "Tune committees: ", paste(committees_values, collapse = ", "))
  log_line(logfile, "Tune neighbors:  ", paste(neighbors_values, collapse = ", "))
  
  df_raw <- if (use_log_transform && exists("df_conc_log")) df_conc_log else df_conc
  
  df <- df_raw %>%
    dplyr::select(all_of(c(depth_name, cov_names))) %>%
    filter(complete.cases(.), is.finite(.data[[depth_name]]))
  
  log_line(logfile, "Records after cleaning: ", nrow(df))
  
  # Log target variable statistics for diagnostics
  y_vals <- df[[depth_name]]
  log_line(logfile, "Target stats: min=", round(min(y_vals), 4),
           " | max=", round(max(y_vals), 4),
           " | mean=", round(mean(y_vals), 4),
           " | sd=", round(sd(y_vals), 4),
           " | CV%=", round(100 * sd(y_vals) / abs(mean(y_vals)), 2))
  
  if (nrow(df) < 100) {
    log_line(logfile, "Skipping — too few samples")
    next
  }
  
  # --- Train / Test split
  idx <- createDataPartition(df[[depth_name]], p = 0.8, list = FALSE)
  train <- df[idx, , drop = FALSE]
  test  <- df[-idx, , drop = FALSE]
  log_line(logfile, "Train: ", nrow(train), " | Test: ", nrow(test))
  
  # --- Adaptive CV: use fewer folds and repeats for small samples
  n_train <- nrow(train)
  if (n_train < 500) {
    actual_folds <- 5
    actual_repeats <- 5  # More repeats to compensate for fewer folds
    log_line(logfile, "Small sample size - using 5-fold CV with 5 repeats")
  } else if (n_train < 1000) {
    actual_folds <- 5
    actual_repeats <- 3
    log_line(logfile, "Medium sample size - using 5-fold CV with 3 repeats")
  } else {
    actual_folds <- cv_folds
    actual_repeats <- cv_repeats
  }
  
  ctrl <- trainControl(
    method = "repeatedcv",
    number = actual_folds,
    repeats = actual_repeats,
    allowParallel = TRUE,
    savePredictions = "final"
  )
  
  # --- Remove NZV predictors (per depth) + LOG NAMES
  nzv <- nearZeroVar(train[, cov_names, drop = FALSE], saveMetrics = TRUE)
  
  cov_nzv_keep <- rownames(nzv)[!nzv$nzv]
  cov_nzv_drop <- rownames(nzv)[nzv$nzv]
  
  log_line(logfile, "NZV kept (", length(cov_nzv_keep), "):")
  log_line(logfile, paste(cov_nzv_keep, collapse = ", "))
  
  if (length(cov_nzv_drop) > 0) {
    log_line(logfile, "NZV dropped (", length(cov_nzv_drop), "):")
    log_line(logfile, paste(cov_nzv_drop, collapse = ", "))
  } else {
    log_line(logfile, "NZV dropped (0): none")
  }
  
  # --- Split numeric vs categorical (LOG NAMES)
  is_num <- sapply(train[, cov_nzv_keep, drop = FALSE], is.numeric)
  
  cov_num <- cov_nzv_keep[is_num]
  cov_cat <- cov_nzv_keep[!is_num]
  
  log_line(logfile, "Numeric covariates (", length(cov_num), "):")
  log_line(logfile, if (length(cov_num) > 0) paste(cov_num, collapse = ", ") else "none")
  
  log_line(logfile, "Categorical covariates (", length(cov_cat), "):")
  log_line(logfile, if (length(cov_cat) > 0) paste(cov_cat, collapse = ", ") else "none")
  
  # --- Remove highly correlated predictors (NUMERIC ONLY) + LOG NAMES
  # Using relaxed cutoff (0.95) to keep more predictors
  cov_corr_drop <- character(0)
  if (length(cov_num) >= 2) {
    corr <- cor(train[, cov_num, drop = FALSE], use = "pairwise.complete.obs")
    drop_idx <- findCorrelation(corr, cutoff = correlation_cutoff)
    if (length(drop_idx) > 0) {
      cov_corr_drop <- cov_num[drop_idx]
      cov_num <- cov_num[-drop_idx]
    }
  }
  
  if (length(cov_corr_drop) > 0) {
    log_line(logfile, "Correlation dropped (|r| > ", correlation_cutoff, ") (", length(cov_corr_drop), "):")
    log_line(logfile, paste(cov_corr_drop, collapse = ", "))
  } else {
    log_line(logfile, "Correlation dropped (0): none")
  }
  
  # --- Final covariates used in model
  cov_use <- c(cov_num, cov_cat)
  log_line(logfile, "FINAL covariates used (", length(cov_use), "):")
  log_line(logfile, paste(cov_use, collapse = ", "))
  
  # --- Model formula
  formula <- as.formula(
    paste(depth_name, "~", paste(cov_use, collapse = " + "))
  )
  
  # --- Train Cubist
  log_line(logfile, "Training Cubist...")
  t0 <- Sys.time()
  
  t.cubist <- train(
    formula,
    data = train,
    method = "cubist",
    trControl = ctrl,
    tuneGrid = cubist_grid,
    control = cubistControl(
      extrapolation = 0.5,
      unbiased = TRUE
    )
  )
  
  log_line(logfile, "CV training time (sec): ",
           round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  
  bt <- t.cubist$bestTune
  log_line(logfile, "BestTune: committees=", bt$committees, " | neighbors=", bt$neighbors)
  
  # Save CV grid + best tune
  write.csv(
    t.cubist$results,
    file.path(dirs$cv, paste0("cv_grid_", depth_name, ".csv")),
    row.names = FALSE
  )
  saveRDS(
    bt,
    file.path(dirs$cv, paste0("bestTune_", depth_name, ".rds"))
  )
  
  # --- Variable usage (%) + log top covariate
  usage_df <- get_cubist_usage_pct(t.cubist)
  if (!is.null(usage_df) && nrow(usage_df) > 0) {
    write.csv(
      usage_df,
      file.path(dirs$importance, paste0("usage_pct_", depth_name, ".csv")),
      row.names = FALSE
    )
    
    top1 <- usage_df[1, ]
    log_line(logfile, "Top covariate (Cubist usage): ", top1$covariate,
             " | ", round(top1$pct, 2), "%")
    
    top_cov_summary[[depth_name]] <- data.frame(
      depth = depth_name,
      top_covariate = top1$covariate,
      pct = round(top1$pct, 2),
      stringsAsFactors = FALSE
    )
    
    log_line(logfile, "Top 10 covariates by usage (%):")
    top10 <- utils::head(usage_df[, c("covariate", "pct")], 10)
    cat(paste(capture.output(print(top10)), collapse = "\n"),
        "\n", file = logfile, append = TRUE)
    print(top10)
  } else {
    log_line(logfile, "No usage statistics available from Cubist model.")
  }
  
  # --- Test prediction
  pred <- predict(t.cubist, newdata = test)
  obs  <- test[[depth_name]]
  
  # --- Metrics (R²: SSE/TSS)
  SSE <- sum((obs - pred)^2)
  TSS <- sum((obs - mean(obs))^2)
  
  R2  <- 1 - SSE / TSS
  RMSE <- sqrt(mean((obs - pred)^2))
  MAE  <- mean(abs(obs - pred))
  Bias <- mean(pred - obs)
  
  CCC <- 2 * cov(obs, pred) /
    (var(obs) + var(pred) + (mean(obs) - mean(pred))^2)
  
  # --- Back-transform metrics (if trained in log scale)
  if (use_log_transform) {
    obs_bt  <- exp(obs)
    pred_bt <- exp(pred)
    R2_bt <- 1 - sum((obs_bt - pred_bt)^2) /
      sum((obs_bt - mean(obs_bt))^2)
  } else {
    R2_bt <- NA_real_
  }
  
  metrics <- data.frame(
    depth = depth_name,
    n_train = nrow(train),
    n_test  = nrow(test),
    committees = bt$committees,
    neighbors  = bt$neighbors,
    RMSE = round(RMSE, 4),
    MAE  = round(MAE, 4),
    R2   = round(R2, 4),
    R2_backtrans = round(R2_bt, 4),
    Bias = round(Bias, 4),
    CCC  = round(CCC, 4),
    stringsAsFactors = FALSE
  )
  
  all_metrics[[depth_name]] <- metrics
  
  write.csv(
    metrics,
    file.path(dirs$metrics, paste0("metrics_", depth_name, ".csv")),
    row.names = FALSE
  )
  
  # Save test predictions
  preds_out <- data.frame(
    depth = depth_name,
    observed = obs,
    predicted = pred,
    residual = obs - pred,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    preds_out,
    file.path(dirs$preds, paste0("pred_test_", depth_name, ".csv")),
    row.names = FALSE
  )
  
  log_line(logfile, "Metrics: RMSE=", round(RMSE, 4),
           " | MAE=", round(MAE, 4),
           " | R2=", round(R2, 4),
           " | R2_backtrans=", round(R2_bt, 4),
           " | Bias=", round(Bias, 4),
           " | CCC=", round(CCC, 4))
  
  # --- Final model (neighbors properly applied)
  final_model <- cubist(
    x = df[, cov_use, drop = FALSE],
    y = df[[depth_name]],
    committees = bt$committees,
    neighbors  = bt$neighbors
  )
  
  saveRDS(
    final_model,
    file.path(dirs$models, paste0(depth_name, ".rds"))
  )
  
  # Save caret model too (optional but useful)
  saveRDS(
    t.cubist,
    file.path(dirs$models, paste0("caret_", depth_name, ".rds"))
  )
  
  # Variable importance via caret (secondary)
  vi <- varImp(t.cubist)$importance %>%
    rownames_to_column("covariate") %>%
    arrange(desc(Overall))
  
  write.csv(
    vi,
    file.path(dirs$importance, paste0("importance_", depth_name, ".csv")),
    row.names = FALSE
  )
  
  log_line(logfile, "Saved model: ", file.path(dirs$models, paste0(depth_name, ".rds")))
  log_line(logfile, "Saved usage%: ", file.path(dirs$importance, paste0("usage_pct_", depth_name, ".csv")))
  log_line(logfile, "Saved metrics: ", file.path(dirs$metrics, paste0("metrics_", depth_name, ".csv")))
  log_line(logfile, "Done depth: ", depth_name)
}

# =============================================================================
# SUMMARY
# =============================================================================
if (length(all_metrics) > 0) {
  summary_df <- do.call(rbind, all_metrics)
  
  write.csv(
    summary_df,
    file.path(dirs$metrics, "all_depths_metrics.csv"),
    row.names = FALSE
  )
  
  cat("\n==============================\n")
  cat(">>> PERFORMANCE SUMMARY\n")
  cat("==============================\n")
  print(summary_df)
} else {
  cat("\nNo depth models were trained (too few records?)\n")
}

# Top covariate summary across depths
if (length(top_cov_summary) > 0) {
  top_cov_df <- do.call(rbind, top_cov_summary)
  
  write.csv(
    top_cov_df,
    file.path(dirs$importance, "top_covariate_by_depth.csv"),
    row.names = FALSE
  )
  
  cat("\n==============================\n")
  cat(">>> TOP COVARIATE BY DEPTH (Cubist usage %)\n")
  cat("==============================\n")
  print(top_cov_df)
}

cat("\n✓ Cubist modelling completed successfully\n")
