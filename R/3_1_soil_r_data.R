# Read CSV
df <- read.csv("/Users/neo/Development/Data-Science/data-science/Data/data_in/Soil_with_raster/OC_with_covariates.csv")

# Identify columns
depth_cols <- c("X0.5cm", "X5.15cm", "X15.30cm", "X30.60cm", "X60.100cm", "X100.200cm")
covariate_cols <- setdiff(names(df), c("id", "Longitude", "Latitude", depth_cols))

# Define threshold for filtering depth values
depth_threshold <- 10  # change this as needed

# Regression data
df_conc <- df[, c(depth_cols, covariate_cols)]

# 1. Set depth values > threshold to NA
for (col in depth_cols) {
  df_conc[[col]][df_conc[[col]] > depth_threshold] <- NA
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
save(df_preab, df_conc, cov_names, 
     file = "/Users/neo/Development/Data-Science/data-science/Data/data_out/RData/oc_covs_regression.RData")
