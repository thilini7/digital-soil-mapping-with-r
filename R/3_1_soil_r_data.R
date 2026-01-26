# Read CSV
df <- read.csv("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_out/Soil_data_with_covariates/OC_with_covariates_new.csv")

# Identify columns
depth_cols <- c("X0.5cm", "X5.15cm", "X15.30cm", "X30.60cm", "X60.100cm", "X100.200cm")
covariate_cols <- setdiff(names(df), c("id", "Longitude", "Latitude", depth_cols))

# Define thresholds for filtering depth values
depth_threshold_lower <- 0.01  # lower limit - change as needed
depth_threshold_upper <- 10  # upper limit - change as needed

# Regression datap
df_conc <- df[, c(depth_cols, covariate_cols)]

# 1. Set depth values outside threshold range to NA
for (col in depth_cols) {
  df_conc[[col]][df_conc[[col]] < depth_threshold_lower | df_conc[[col]] > depth_threshold_upper] <- NA
}

# 2. Remove rows where all depth columns are NA
df_conc <- df_conc[rowSums(!is.na(df_conc[, depth_cols])) > 0, ]

# 3. Remove rows with NA in any covariate column
df_conc <- df_conc[complete.cases(df_conc[, covariate_cols]), ]

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

save(df_preab, df_conc, cov_names, 
     file = file.path(output_dir, "OC_covs_regression.RData"))
