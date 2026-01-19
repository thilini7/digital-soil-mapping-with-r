# =============================================================================
# Script: Add LDI (Land Disturbance Index) to NSW Land Use Shapefile
# =============================================================================
# This script assigns LDI values to each polygon in the NSW Land Use shapefile
# based on the Secondary land use classification.
#
# LDI Classification (Gray et al., 2015):
# - LDI 1: No effective disturbance (e.g., national park, nature reserve)
# - LDI 2: Limited disturbance, minor native vegetation clearing 
#          (e.g., selective logging, production forestry)
# - LDI 3: Moderate disturbance, moderate native vegetation clearing, 
#          light grazing in woodland, hardwood plantation
# - LDI 4: High disturbance, complete native vegetation clearing 
#          (e.g., native and improved pasture, softwood plantation)
# - LDI 5: Very high disturbance (e.g., improved pasture with moderate 
#          cropping, orchards, viticulture)
# - LDI 6: Extreme disturbance, predominant cropping (rain-fed or irrigated)
# =============================================================================

# --- Load required packages ---
library(sf)
library(dplyr)

# =============================================================================
# 1. Define file paths
# =============================================================================
base_dir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"

shapefile_input <- file.path(base_dir, "Data/data_in/shp_files/NSW/NSWLanduse2017_Ver1_5_20230921.shp")
shapefile_output <- file.path(base_dir, "Data/data_in/shp_files/NSW/NSWLanduse2017_with_LDI.shp")

# =============================================================================
# 2. Define LDI mapping based on NSW Land Use Secondary Classification
# =============================================================================

create_ldi_mapping <- function() {
  # Create a lookup table for Secondary land use codes to LDI values
  ldi_map <- data.frame(
    Secondary = c(
      # Conservation and minimal use (LDI 1)
      "1.1.0 Nature conservation",
      "1.2.0 Managed resource protection", 
      "1.3.0 Other minimal use",
      
      # Production forestry (LDI 2)
      "2.2.0 Production native forestry",
      
      # Light grazing and hardwood plantation (LDI 3)
      "2.1.0 Grazing native vegetation",
      "3.1.0 Plantation forests",
      
      # Modified pastures and softwood (LDI 4)
      "3.2.0 Grazing modified pastures",
      "3.6.0 Land in transition",
      "4.6.0 Irrigated land in transition",
      
      # High intensity - horticulture, orchards, vineyards (LDI 5)
      "3.4.0 Perennial horticulture",
      "3.5.0 Seasonal horticulture",
      "4.2.0 Grazing irrigated modified pastures",
      "4.4.0 Irrigated perennial horticulture",
      "4.5.0 Irrigated seasonal horticulture",
      "5.1.0 Intensive horticulture",
      "5.2.0 Intensive animal production",
      
      # Extreme disturbance - cropping (LDI 6)
      "3.3.0 Cropping",
      "4.1.0 Irrigated plantation forests",
      "4.3.0 Irrigated cropping",
      
      # Urban/Industrial - assign LDI 6 (extreme disturbance)
      "5.3.0 Manufacturing and industrial",
      "5.4.0 Residential and farm infrastructure",
      "5.5.0 Services",
      "5.6.0 Utilities",
      "5.7.0 Transport and communication",
      "5.8.0 Mining",
      "5.9.0 Waste treatment and disposal",
      
      # Water bodies - assign NA
      "6.1.0 Lake",
      "6.2.0 Reservoir/dam",
      "6.3.0 River",
      "6.4.0 Channel/aqueduct",
      "6.5.0 Marsh/wetland",
      "6.6.0 Estuary/coastal waters"
    ),
    LDI_value = c(
      # Conservation (LDI 1)
      1, 1, 1,
      # Production forestry (LDI 2)
      2,
      # Light grazing (LDI 3)
      3, 3,
      # Modified pastures (LDI 4)
      4, 4, 4,
      # Horticulture (LDI 5)
      5, 5, 5, 5, 5, 5, 5,
      # Cropping (LDI 6)
      6, 6, 6,
      # Urban/Industrial (LDI 6)
      6, 6, 6, 6, 6, 6, 6,
      # Water bodies (NA)
      NA, NA, NA, NA, NA, NA
    ),
    stringsAsFactors = FALSE
  )
  
  return(ldi_map)
}

# =============================================================================
# 3. Load NSW Land Use shapefile
# =============================================================================
cat("=" , rep("=", 69), "\n", sep = "")
cat("Adding LDI (Land Disturbance Index) to NSW Land Use Shapefile\n")
cat("=" , rep("=", 69), "\n\n", sep = "")

cat("1. Loading NSW Land Use shapefile...\n")
cat("   File:", shapefile_input, "\n")

landuse_sf <- st_read(shapefile_input, quiet = TRUE)
cat("   Loaded", nrow(landuse_sf), "polygons\n")
cat("   Columns:", paste(names(landuse_sf), collapse = ", "), "\n")

# =============================================================================
# 4. Check existing Secondary column
# =============================================================================
cat("\n2. Checking Secondary land use classification column...\n")

if (!"Secondary" %in% names(landuse_sf)) {
  stop("ERROR: 'Secondary' column not found in shapefile!")
}

# Show unique Secondary values
unique_secondary <- unique(landuse_sf$Secondary)
cat("   Found", length(unique_secondary), "unique Secondary land use categories\n")
cat("   Sample categories:\n")
print(head(unique_secondary, 10))

# =============================================================================
# 5. Create and apply LDI mapping
# =============================================================================
cat("\n3. Creating LDI mapping...\n")
ldi_mapping <- create_ldi_mapping()
cat("   Defined", nrow(ldi_mapping), "land use to LDI mappings\n")

# Remove existing LDI column if present (with NA values)
if ("LDI" %in% names(landuse_sf)) {
  cat("   Removing existing empty LDI column...\n")
  landuse_sf <- landuse_sf[, setdiff(names(landuse_sf), "LDI")]
}

# Join LDI values to land use data
cat("\n4. Joining LDI values to polygons...\n")
landuse_sf <- landuse_sf %>%
  left_join(ldi_mapping, by = "Secondary") %>%
  rename(LDI = LDI_value)

# =============================================================================
# 6. Report LDI distribution
# =============================================================================
cat("\n5. LDI value distribution:\n")
ldi_summary <- landuse_sf %>%
  st_drop_geometry() %>%
  group_by(LDI) %>%
  summarise(
    Count = n(),
    Percentage = round(n() / nrow(landuse_sf) * 100, 2)
  ) %>%
  arrange(LDI)

print(as.data.frame(ldi_summary))

# Count unmapped categories
unmapped <- landuse_sf %>%
  st_drop_geometry() %>%
  filter(is.na(LDI)) %>%
  pull(Secondary) %>%
  unique()

if (length(unmapped) > 0) {
  cat("\n   WARNING: Some Secondary categories were not mapped to LDI:\n")
  print(unmapped)
}

# =============================================================================
# 7. Save updated shapefile
# =============================================================================
cat("\n6. Saving updated shapefile...\n")
cat("   Output:", shapefile_output, "\n")

st_write(landuse_sf, shapefile_output, delete_dsn = TRUE)

cat("\n" , rep("=", 69), "\n", sep = "")
cat("✅ SUCCESS! LDI column added to shapefile\n")
cat("   Total polygons:", nrow(landuse_sf), "\n")
cat("   Polygons with LDI values:", sum(!is.na(landuse_sf$LDI)), "\n")
cat("   Polygons without LDI (water bodies):", sum(is.na(landuse_sf$LDI)), "\n")
cat(rep("=", 69), "\n", sep = "")
