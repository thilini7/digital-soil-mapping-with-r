# Read CSV
df <- read.csv("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/Soil_data_with_covariates/Phosphorus_with_covariates_new.csv")

# Identify columns
depth_cols <- c("X0.5cm", "X5.15cm", "X15.30cm", "X30.60cm", "X60.100cm", "X100.200cm")

# Get all potential covariates
all_covariates <- setdiff(names(df), c("id", "Longitude", "Latitude", depth_cols))

# Exclude covariates that change in the future (for climate projection studies)
# - FPAR (fractional vegetation) - changes with climate
# - Clay content - if present, changes are uncertain
excluded_covariates <- c(
  "FPAR_pct_max_NSW_ACT_90m",
  "FPAR_pct_mean_NSW_ACT_90m", 
  "FPAR_pct_median_NSW_ACT_90m",
  "FPAR_pct_min_NSW_ACT_90m"
  # Add any clay-related covariates here if present
)

# Keep only static (terrain, parent material) and climate covariates (have future projections)
covariate_cols <- setdiff(all_covariates, excluded_covariates)
cat("Excluded covariates (change in future):", paste(excluded_covariates, collapse = ", "), "\n")
cat("Using", length(covariate_cols), "covariates for modeling:\n")
cat(paste(covariate_cols, collapse = ", "), "\n\n")

# Regression data
df_conc <- df[, c(depth_cols, covariate_cols)]

# Remove rows where all depth columns are NA
df_conc <- df_conc[rowSums(!is.na(df_conc[, depth_cols])) > 0, ]

# Remove rows with NA in any covariate column
df_conc <- df_conc[complete.cases(df_conc[, covariate_cols]), ]

# 4. Create log-transformed version (often improves model accuracy for soil OC)
df_conc_log <- df_conc
for (col in depth_cols) {
  # log1p handles zeros: log1p(x) = log(1 + x)
  df_conc_log[[col]] <- log1p(df_conc_log[[col]])
}
cat("Created log-transformed data (df_conc_log) for improved modeling\n")

# Classification data (presence/absence)
df_preab <- df_conc
for (col in depth_cols) {
  df_preab[[col]] <- ifelse(!is.na(df_conc[[col]]) & df_conc[[col]] > 0, 1, 0)
  df_preab[[col]] <- factor(df_preab[[col]], levels = c(0, 1))
}

# Save
cov_names <- covariate_cols
output_dir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/RData"
if (!dir.exists(output_dir)) {
  # If it doesn't exist, create it
  dir.create(output_dir, recursive = TRUE)
  cat(paste("Directory", output_dir, "created.\n"))
}

save(df_preab, df_conc, df_conc_log, cov_names, 
     file = file.path(output_dir, "Phosphorus_covs_regression.RData"))

cat("\nSaved: df_conc (original), df_conc_log (log-transformed), df_preab (classification), cov_names\n")
