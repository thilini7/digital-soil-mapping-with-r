# =============================================================================
# Script: Extract Specific Lab Method Data from ANSIS Combined Dataset
# =============================================================================
# This script extracts soil data for specific laboratory methods from the
# ANSIS_combined_with_LDI.csv file and saves them as separate CSV files.
#
# Common Lab Methods (ObservedProperty codes):
# pH:
#   - 4A1: pH in water (1:5 soil:water)
#
# Organic Carbon:
#   - 6A1: Organic carbon - Walkley-Black
#
# Electrical Conductivity:
#   - 3A1: EC (1:5 soil:water)
#
# Cation Exchange Capacity:
#   - 15A1: CEC - Ammonium acetate at pH 7
#   - 15C1: CEC - Alcoholic ammonium chloride
#   - 15D3: CEC - Compulsive exchange
#   - 15E1: ECEC
#
# Exchangeable Cations:
#   - 15F1: Exchangeable Ca
#   - 15F2: Exchangeable Mg
#   - 15F3: Exchangeable Na
#   - 15G1: Exchangeable K
#   - 15M1: Exchangeable bases sum
#   - 15N1: Exchangeable acidity
#
# Particle Size:
#   - 2B1: Particle size - clay content
#
# Bulk Density:
#   - 503.01, 503.03: Bulk density methods
#
# Nitrogen:
#   - 7C2b_NH4_N: Ammonium N - 2 M KCl
#
# Phosphorus:
#   - 9E1: Available P - Olsen
# 
# Moisture content:
#   - 2B1: As received moisture content
# =============================================================================

# --- Install required packages (only if not already installed) ---
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(dplyr)
library(readr)
library(tidyr)

# =============================================================================
# 1. Define file paths and parameters
# =============================================================================
base_dir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"

# Input file (with or without LDI)
input_file <- file.path(base_dir, "Data/data_out/ansis_lab_measurements/All Sites/ANSIS_combined_with_LDI.csv")

# If LDI file doesn't exist, use original
if (!file.exists(input_file)) {
  input_file <- file.path(base_dir, "Data/data_out/ansis_lab_measurements/All Sites/ANSIS_combined.csv")
  cat("Note: Using original file without LDI\n")
}

# Output directory
output_dir <- file.path(base_dir, "Data/data_out/lab_method_extracts")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. Define lab method groups for extraction
# =============================================================================

# Define which lab methods to extract (modify as needed)
lab_method_groups <- list(
  
  # pH measurements
  pH = list(
    methods = c("4A1"),
    description = "pH measurements (water)",
    output_name = "pH"
  ),
  
  # Organic Carbon
  OC = list(
    methods = c("6A1"),
    description = "Organic Carbon (Walkley-Black)",
    output_name = "Organic_Carbon"
  ),
  
  # Electrical Conductivity
  EC = list(
    methods = c("3A1"),
    description = "Electrical Conductivity",
    output_name = "EC"
  ),
  
  # Cation Exchange Capacity
  CEC = list(
    methods = c("15A1", "15C1", "15D3", "15E1"),
    description = "Cation Exchange Capacity",
    output_name = "CEC"
  ),
  
  # Exchangeable Cations
  ExCa = list(
    methods = c("15F1"),
    description = "Exchangeable Calcium",
    output_name = "Exchangeable_Ca"
  ),
  
  ExMg = list(
    methods = c("15F2"),
    description = "Exchangeable Magnesium",
    output_name = "Exchangeable_Mg"
  ),
  
  ExNa = list(
    methods = c("15F3"),
    description = "Exchangeable Sodium",
    output_name = "Exchangeable_Na"
  ),
  
  ExK = list(
    methods = c("15G1"),
    description = "Exchangeable Potassium",
    output_name = "Exchangeable_K"
  ),
  
  # Bulk Density
  BD = list(
    methods = c("503.01", "503.03"),
    description = "Bulk Density",
    output_name = "Bulk_Density"
  ),
  
  # Nitrogen
  TN = list(
    methods = c("7C2b_NH4_N"),
    description = "Ammonium N (2 M KCl)",
    output_name = "Nitrogen"
  ),
  
  # Available Phosphorus
  P = list(
    methods = c("9E1"),
    description = "Available Phosphorus (Olsen)",
    output_name = "Phosphorus"
  ),
  
  # Moisture Content (As received)
  Moisture = list(
    methods = c("2B1"),
    description = "As Received Moisture Content",
    output_name = "Moisture"
  )
)

# =============================================================================
# 3. Function to extract and save lab method data
# =============================================================================

extract_lab_method <- function(data, method_info, output_dir) {
  
  methods <- method_info$methods
  description <- method_info$description
  output_name <- method_info$output_name
  
  # Filter data for specified methods
  extracted <- data %>%
    filter(ObservedProperty %in% methods)
  
  if (nrow(extracted) == 0) {
    cat("   No data found for", description, "\n")
    return(NULL)
  }
  
  # Save to CSV
  output_file <- file.path(output_dir, paste0(output_name, "_data.csv"))
  write_csv(extracted, output_file)
  
  cat("   ", description, ":\n", sep = "")
  cat("      Records:", nrow(extracted), "\n")
  cat("      Unique locations:", n_distinct(extracted$Location_ID), "\n")
  cat("      Methods:", paste(unique(extracted$ObservedProperty), collapse = ", "), "\n")
  cat("      Saved to:", basename(output_file), "\n")
  
  return(extracted)
}

# =============================================================================
# 4. Function to create wide format (pivot) data
# =============================================================================

create_wide_format <- function(data, output_dir, output_name) {
  
  # Create a simplified dataset with key columns
  wide_data <- data %>%
    select(Location_ID, Longitude, Latitude, UpperDepth, LowerDepth, 
           ObservedProperty, Value, Units, SampleDate,
           any_of(c("LDI", "LandUse_Secondary", "LandUse_Tertiary"))) %>%
    # Create unique row identifier
    mutate(row_id = row_number()) %>%
    # Pivot to wide format
    pivot_wider(
      id_cols = c(Location_ID, Longitude, Latitude, UpperDepth, LowerDepth, 
                  SampleDate, any_of(c("LDI", "LandUse_Secondary", "LandUse_Tertiary"))),
      names_from = ObservedProperty,
      values_from = Value,
      values_fn = list(Value = mean)  # Average if multiple values
    )
  
  # Save wide format
  output_file <- file.path(output_dir, paste0(output_name, "_wide.csv"))
  write_csv(wide_data, output_file)
  
  cat("   Wide format saved:", basename(output_file), "\n")
  cat("      Records:", nrow(wide_data), "\n")
  
  return(wide_data)
}

# =============================================================================
# 5. Load data
# =============================================================================
cat("=" , rep("=", 59), "\n", sep = "")
cat("Extracting Lab Method Data from ANSIS Dataset\n")
cat("=" , rep("=", 59), "\n\n", sep = "")

cat("1. Loading data...\n")
cat("   File:", input_file, "\n")

soil_data <- read_csv(input_file, show_col_types = FALSE)
cat("   Loaded", nrow(soil_data), "records\n")

# Show available methods
cat("\n2. Available lab methods in dataset:\n")
method_counts <- soil_data %>%
  count(ObservedProperty, sort = TRUE) %>%
  head(20)
print(method_counts)

# =============================================================================
# 6. Extract data for each lab method group
# =============================================================================
cat("\n3. Extracting data by lab method group...\n\n")

extracted_list <- list()

for (group_name in names(lab_method_groups)) {
  method_info <- lab_method_groups[[group_name]]
  extracted_list[[group_name]] <- extract_lab_method(
    soil_data, 
    method_info, 
    output_dir
  )
  cat("\n")
}

# =============================================================================
# 7. Create combined datasets for common analyses
# =============================================================================
cat("4. Creating combined datasets...\n\n")

# Combine all exchangeable cations
cat("   Combining exchangeable cations data...\n")
ex_cations <- soil_data %>%
  filter(ObservedProperty %in% c("15F1", "15F2", "15F3", "15G1", "15M1", "15N1"))

if (nrow(ex_cations) > 0) {
  write_csv(ex_cations, file.path(output_dir, "Exchangeable_Cations_all.csv"))
  cat("      Saved: Exchangeable_Cations_all.csv (", nrow(ex_cations), " records)\n", sep = "")
  
  # Create wide format
  create_wide_format(ex_cations, output_dir, "Exchangeable_Cations")
}

# Create all-properties wide format for unique locations
cat("\n   Creating site-level summary with all properties...\n")

# Get one record per location-depth combination with key properties
site_summary <- soil_data %>%
  filter(ObservedProperty %in% c("4A1", "6A1", "3A1", 
                                  "15D3", "15E1", "7C2b_NH4_N", "9E1")) %>%
  select(Location_ID, Longitude, Latitude, UpperDepth, LowerDepth,
         ObservedProperty, Value, 
         any_of(c("LDI", "LandUse_Secondary"))) %>%
  pivot_wider(
    id_cols = c(Location_ID, Longitude, Latitude, UpperDepth, LowerDepth,
                any_of(c("LDI", "LandUse_Secondary"))),
    names_from = ObservedProperty,
    values_from = Value,
    values_fn = list(Value = mean)
  )

write_csv(site_summary, file.path(output_dir, "Site_Summary_wide.csv"))
cat("      Saved: Site_Summary_wide.csv (", nrow(site_summary), " records)\n", sep = "")

# =============================================================================
# 8. Summary
# =============================================================================
cat("\n", "=" , rep("=", 59), "\n", sep = "")
cat("SUMMARY\n")
cat("=" , rep("=", 59), "\n", sep = "")
cat("Output directory:", output_dir, "\n\n")

cat("Files created:\n")
output_files <- list.files(output_dir, pattern = "\\.csv$")
for (f in output_files) {
  file_path <- file.path(output_dir, f)
  file_size <- file.info(file_path)$size / 1024  # KB
  cat("  ", f, " (", round(file_size, 1), " KB)\n", sep = "")
}

cat("\n✅ Extraction complete!\n")
