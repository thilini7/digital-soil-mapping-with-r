# =============================================================================
# SOIL PROPERTY RANGE ANALYSIS & STATISTICS
# Calculates percentile-based ranges and summary statistics by depth
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(e1071)   # skewness()
  library(stringr)
})

# =============================================================================
# FILE PATHS - CENTRALIZED CONFIGURATION
# =============================================================================
HomeDir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
setwd(HomeDir)

# INPUT: Soil property data file
soil_property <- "EC"  #Options: Organic_Carbon, Nitrogen, Phosphorus, pH, Bulk_Density, CEC, EC, Clay, Sum_of_Bases etc.

# =============================================================================
# PROPERTY-SPECIFIC PERCENTILE SETTINGS
# =============================================================================
# Different properties have different distributions - adjust percentiles accordingly
percentile_settings <- list(
  "Organic_Carbon" = c(lower = 0.01, upper = 0.99),  # P1-P99
  "pH"             = c(lower = 0.01, upper = 0.99),  # P1-P99
  "Bulk_Density"   = c(lower = 0.01, upper = 0.99),  # P1-P99
  "Clay"           = c(lower = 0.01, upper = 0.99),  # P1-P99
  "CEC"            = c(lower = 0.01, upper = 0.99),  # P1-P99
  "EC"             = c(lower = 0.01, upper = 0.95),  # P1-P95 (very skewed)
  "Nitrogen"       = c(lower = 0.01, upper = 0.99),  # P1-P99
  "Phosphorus"     = c(lower = 0.01, upper = 0.95),  # P1-P95 (fertilizer hotspots)
  "Sum_of_Bases"   = c(lower = 0.01, upper = 0.99)   # P1-P99
)

# Get percentiles for current property (default P1-P99 if not specified)
use_raw_range <- FALSE
if (soil_property == "EC") {
  # EC uses fixed raw value range (dS/m) instead of percentiles
  # because EC is extremely skewed with a long tail
  lower_raw <- 0.1    # Minimum valid EC (dS/m)
  upper_raw <- 15    # Maximum valid EC (dS/m)
  use_raw_range <- TRUE
  cat("Note: EC uses fixed range [", lower_raw, "-", upper_raw, "] dS/m instead of percentiles\n")
} else if (soil_property %in% names(percentile_settings)) {
  lower_pct <- percentile_settings[[soil_property]]["lower"]
  upper_pct <- percentile_settings[[soil_property]]["upper"]
} else {
  lower_pct <- 0.01
  upper_pct <- 0.99
  cat("Note: Using default P1-P99 for", soil_property, "\n")
}

# Input file path - adjust to match your data location
file_path <- file.path(HomeDir, "Data/data_out/Soil_data_with_covariates_v2",
                       paste0(soil_property, "_with_covariates_new.csv"))

# Output directory for stats
output_dir <- file.path(HomeDir, "Data/data_out/Soil_property_ranges")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("SOIL PROPERTY RANGE ANALYSIS\n")
cat("Property: ", soil_property, "\n", sep = "")
cat("Input file: ", basename(file_path), "\n", sep = "")
if (use_raw_range) {
  cat("Fixed range: ", lower_raw, " - ", upper_raw, " (raw values)\n", sep = "")
} else {
  cat("Percentile range: P", lower_pct*100, "-P", upper_pct*100, "\n", sep = "")
}
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

# =============================================================================
# LOAD DATA
# =============================================================================
df <- read_csv(file_path, show_col_types = FALSE)
cat("Loaded", nrow(df), "records\n\n")

# =============================================================================
# IDENTIFY SOIL PROPERTY COLUMNS FOR DEPTH INTERVALS
# =============================================================================
# Your column names: X0.30cm, X30.60cm, X60.100cm, X100.200cm
depth_cols <- c(
  "0-30"    = "X0.30cm",
  "30-60"   = "X30.60cm", 
  "60-100"  = "X60.100cm",
  "100-200" = "X100.200cm"
)

# Verify columns exist
missing_cols <- depth_cols[!depth_cols %in% names(df)]
if (length(missing_cols) > 0) {
  cat("WARNING: Missing columns:\n")
  print(missing_cols)
  cat("\nAvailable columns:\n")
  print(names(df))
  stop("Please update depth_cols to match your column names.")
}

cat("Using depth columns:\n")
print(depth_cols)
cat("\n")

# =============================================================================
# CREATE LONG FORMAT TABLE
# =============================================================================
prop_long <- bind_rows(lapply(names(depth_cols), function(d) {
  col <- depth_cols[[d]]
  tibble(depth = d, value = df[[col]])
})) %>%
  mutate(value = as.numeric(value))

# =============================================================================
# CALCULATE SUMMARY STATISTICS + SKEWNESS + PERCENTILES
# =============================================================================
prop_stats <- prop_long %>%
  group_by(depth) %>%
  summarise(
    n_total   = n(),
    n_missing = sum(is.na(value)),
    n_used    = sum(!is.na(value)),
    min       = suppressWarnings(min(value, na.rm = TRUE)),
    q_lower   = if (use_raw_range) lower_raw else suppressWarnings(quantile(value, lower_pct, na.rm = TRUE, names = FALSE)),
    q05       = suppressWarnings(quantile(value, 0.05, na.rm = TRUE, names = FALSE)),
    median    = suppressWarnings(median(value, na.rm = TRUE)),
    mean      = suppressWarnings(mean(value, na.rm = TRUE)),
    q95       = suppressWarnings(quantile(value, 0.95, na.rm = TRUE, names = FALSE)),
    q_upper   = if (use_raw_range) upper_raw else suppressWarnings(quantile(value, upper_pct, na.rm = TRUE, names = FALSE)),
    max       = suppressWarnings(max(value, na.rm = TRUE)),
    skewness  = suppressWarnings(e1071::skewness(value, na.rm = TRUE, type = 2)),
    .groups = "drop"
  ) %>%
  # Calculate percent that would be capped
  left_join(
    prop_long %>%
      group_by(depth) %>%
      summarise(
        q_lower = if (use_raw_range) lower_raw else quantile(value, lower_pct, na.rm = TRUE, names = FALSE),
        q_upper = if (use_raw_range) upper_raw else quantile(value, upper_pct, na.rm = TRUE, names = FALSE),
        n_used = sum(!is.na(value)),
        n_capped = sum(!is.na(value) & (value < q_lower | value > q_upper)),
        pct_capped = 100 * n_capped / n_used,
        .groups = "drop"
      ) %>% 
      select(depth, n_capped, pct_capped),
    by = "depth"
  ) %>%
  arrange(factor(depth, levels = c("0-30", "30-60", "60-100", "100-200")))

# =============================================================================
# DISPLAY RECOMMENDED RANGE
# =============================================================================
best_range <- prop_stats %>%
  transmute(
    depth,
    P1  = round(q_lower, 3),
    P99 = round(q_upper, 3),
    capped_n = n_capped,
    capped_pct = round(pct_capped, 2)
  )

# If P1 is 0, set it to 0.01 (avoid zero values for log transforms)
best_range$P1[best_range$P1 == 0] <- 0.01

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
if (use_raw_range) {
  cat("RECOMMENDED", soil_property, "TRAINING RANGE (fixed:", lower_raw, "-", upper_raw, "dS/m)\n")
} else {
  cat("RECOMMENDED", soil_property, "TRAINING RANGE (P", lower_pct*100, "-P", upper_pct*100, ")\n", sep = " ")
}
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")
print(as.data.frame(best_range))

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("FULL SUMMARY STATISTICS BY DEPTH\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")
print(as.data.frame(prop_stats))

# =============================================================================
# GENERATE CODE SNIPPET FOR 3_1_soil_r_data.R
# =============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("CODE SNIPPET FOR 3_1_soil_r_data.R (copy-paste ready)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

if (use_raw_range) {
  cat(sprintf('# Apply valid range filters for %s (fixed range: %g - %g dS/m)\n', soil_property, lower_raw, upper_raw))
  cat(sprintf('if (soil_property == "%s") {\n', soil_property))
  cat(sprintf('  cat("Applying %s fixed range filters [%g - %g dS/m]...\\n")\n', soil_property, lower_raw, upper_raw))
} else {
  cat(sprintf('# Apply valid range filters for %s (P%g-P%g percentiles)\n', soil_property, lower_pct*100, upper_pct*100))
  cat(sprintf('if (soil_property == "%s") {\n', soil_property))
  cat(sprintf('  cat("Applying %s range filters (P%g-P%g)...\\n")\n', soil_property, lower_pct*100, upper_pct*100))
}
cat('  n_before <- nrow(df_conc)\n')
cat('  \n')

for (i in seq_along(depth_cols)) {
  d <- names(depth_cols)[i]
  col <- depth_cols[[i]]
  p_low <- round(prop_stats$q_lower[i], 2)
  p_high <- round(prop_stats$q_upper[i], 2)
  cat(sprintf('  df_conc$%s[df_conc$%s < %.2f | df_conc$%s > %.2f] <- NA\n', 
              col, col, p_low, col, p_high))
}

cat('  \n')
cat('  # Remove rows where all depth columns are now NA\n')
cat('  df_conc <- df_conc[rowSums(!is.na(df_conc[, depth_cols])) > 0, ]\n')
cat('  cat("Removed", n_before - nrow(df_conc), "rows with out-of-range values\\n")\n')
cat('}\n')

# =============================================================================
# OPTIONAL: CREATE CAPPED (WINSORIZED) COLUMNS
# =============================================================================
cap_to_range <- function(x) {
  if (use_raw_range) {
    p_low  <- lower_raw
    p_high <- upper_raw
  } else {
    p_low  <- quantile(x, lower_pct, na.rm = TRUE, names = FALSE)
    p_high <- quantile(x, upper_pct, na.rm = TRUE, names = FALSE)
  }
  x <- ifelse(is.na(x), NA, pmin(pmax(x, p_low), p_high))
  x
}

for (d in names(depth_cols)) {
  col <- depth_cols[[d]]
  newcol <- paste0(soil_property, "_", gsub("-", "_", d), "_capped")
  df[[newcol]] <- cap_to_range(as.numeric(df[[col]]))
  df[[paste0(newcol, "_log")]] <- log1p(df[[newcol]])  # log-transformed target
}

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
# Save statistics to CSV
stats_file <- file.path(output_dir, paste0(soil_property, "_range_stats.csv"))
write_csv(prop_stats, stats_file)
cat("\n\nStatistics saved to:", stats_file, "\n")

# Save best range to CSV
range_file <- file.path(output_dir, paste0(soil_property, "_best_range.csv"))
write_csv(best_range, range_file)
cat("Best range saved to:", range_file, "\n")

# Optional: Save capped data
# capped_file <- file.path(output_dir, paste0(soil_property, "_capped_data.csv"))
# write_csv(df, capped_file)
# cat("Capped data saved to:", capped_file, "\n")

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
