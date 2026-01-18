# =============================================================================
# Script: Add LDI (Land Disturbance Index) to Soil Point Data
# =============================================================================
# This script extracts LDI values from NSW Land Use shapefile and adds them
# to the ANSIS soil point data.
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

# --- Install required packages (only if not already installed) ---
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

library(sf)
library(dplyr)
library(readr)

# =============================================================================
# 1. Define file paths
# =============================================================================
base_dir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"

shapefile_path <- file.path(base_dir, "Data/data_in/shp_files/NSWLanduse2017_Ver1_5_20230921.shp")
csv_input_path <- file.path(base_dir, "Data/data_out/ansis_lab_measurements/All Sites/ANSIS_combined.csv")
csv_output_path <- file.path(base_dir, "Data/data_out/ansis_lab_measurements/All Sites/ANSIS_combined_with_LDI.csv")

# =============================================================================
# 2. Define LDI mapping based on NSW Land Use Secondary Classification
# =============================================================================
# Mapping based on Gray et al. (2015) and Australian Soil and Land Survey Field Handbook

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
      "3.1.0 Plantation forests",  # Note: Could be 3 or 4 depending on type
      
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
      
      # Water bodies - assign NA (or could use LDI 1)
      "6.1.0 Lake",
      "6.2.0 Reservoir/dam",
      "6.3.0 River",
      "6.4.0 Channel/aqueduct",
      "6.5.0 Marsh/wetland",
      "6.6.0 Estuary/coastal waters"
    ),
    LDI = c(
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
# 3. Load soil point data
# =============================================================================
cat("=" , rep("=", 59), "\n", sep = "")
cat("Adding LDI (Land Disturbance Index) to soil point data\n")
cat("=" , rep("=", 59), "\n\n", sep = "")

cat("1. Loading soil point data...\n")
soil_df <- read_csv(csv_input_path, show_col_types = FALSE)
cat("   Loaded", nrow(soil_df), "records\n")
cat("   Columns:", paste(names(soil_df)[1:10], collapse = ", "), "...\n")

# =============================================================================
# 4. Get unique locations to reduce processing
# =============================================================================
cat("\n2. Extracting unique locations...\n")
unique_locs <- soil_df %>%
  select(Longitude, Latitude) %>%
  distinct()
cat("   Found", nrow(unique_locs), "unique locations\n")

# =============================================================================
# 5. Convert soil points to sf object
# =============================================================================
cat("\n3. Converting to spatial points...\n")
points_sf <- st_as_sf(unique_locs, 
                      coords = c("Longitude", "Latitude"),
                      crs = 4283)  # GDA94
cat("   Created sf object with", nrow(points_sf), "points\n")

# Keep original coordinates for later merge
points_sf$Longitude <- unique_locs$Longitude
points_sf$Latitude <- unique_locs$Latitude

# =============================================================================
# 6. Load NSW Land Use shapefile
# =============================================================================
cat("\n4. Loading NSW Land Use shapefile (this may take several minutes)...\n")
cat("   File:", shapefile_path, "\n")

landuse_sf <- st_read(shapefile_path, quiet = TRUE)
cat("   Loaded", nrow(landuse_sf), "land use polygons\n")
cat("   CRS:", st_crs(landuse_sf)$input, "\n")

# =============================================================================
# 7. Add LDI values to land use shapefile based on Secondary classification
# =============================================================================
cat("\n5. Mapping land use to LDI values...\n")

ldi_mapping <- create_ldi_mapping()

# Remove original empty LDI column if it exists
if ("LDI" %in% names(landuse_sf)) {
  landuse_sf <- landuse_sf %>% select(-LDI)
}

# Join LDI values to land use data
landuse_sf <- landuse_sf %>%
  left_join(ldi_mapping, by = "Secondary")

# Check LDI distribution
cat("   LDI value distribution in land use data:\n")
ldi_counts <- table(landuse_sf$LDI, useNA = "ifany")
print(ldi_counts)

# =============================================================================
# 8. Ensure CRS match and fix invalid geometries
# =============================================================================
cat("\n6. Checking coordinate reference systems and fixing geometries...\n")
if (st_crs(points_sf) != st_crs(landuse_sf)) {
  cat("   Transforming points to match land use CRS...\n")
  points_sf <- st_transform(points_sf, st_crs(landuse_sf))
}
cat("   CRS matched\n")

# Disable S2 spherical geometry to avoid edge issues
sf_use_s2(FALSE)
cat("   Disabled S2 spherical geometry processing\n")

# Fix invalid geometries
cat("   Repairing invalid geometries in land use data...\n")
landuse_sf <- st_make_valid(landuse_sf)
cat("   Geometries repaired\n")

# =============================================================================
# 9. Perform spatial join
# =============================================================================
cat("\n7. Performing spatial join (this may take a while)...\n")

# Spatial join - find which polygon each point falls within
joined_sf <- st_join(points_sf, 
                     landuse_sf[, c("Secondary", "Tertiary", "LDI")],
                     join = st_within,
                     left = TRUE)

# Handle duplicates (point in multiple polygons - take first)
joined_df <- joined_sf %>%
  st_drop_geometry() %>%
  distinct(Longitude, Latitude, .keep_all = TRUE)

cat("   Joined", nrow(joined_df), "unique locations\n")
cat("   Points with LDI:", sum(!is.na(joined_df$LDI)), "\n")
cat("   Points without LDI:", sum(is.na(joined_df$LDI)), "\n")

# =============================================================================
# 9b. Apply nearest neighbor for points without LDI (fallback)
# =============================================================================
na_count <- sum(is.na(joined_df$LDI))

if (na_count > 0) {
  cat("\n7b. Applying nearest neighbor for", na_count, "points without LDI...\n")
  
  # Get points that need nearest neighbor assignment
  na_points <- joined_sf %>%
    filter(is.na(LDI)) %>%
    select(Longitude, Latitude) %>%
    distinct()
  
  if (nrow(na_points) > 0) {
    cat("   Finding nearest polygons for", nrow(na_points), "unique locations...\n")
    
    # Find nearest polygon for each NA point
    # Use st_nearest_feature to find index of nearest polygon
    nearest_idx <- st_nearest_feature(na_points, landuse_sf)
    
    # Extract LDI values from nearest polygons
    nearest_ldi <- landuse_sf$LDI[nearest_idx]
    nearest_secondary <- landuse_sf$Secondary[nearest_idx]
    nearest_tertiary <- landuse_sf$Tertiary[nearest_idx]
    
    # Create lookup table for NA points
    na_lookup <- data.frame(
      Longitude = na_points$Longitude,
      Latitude = na_points$Latitude,
      LDI_nearest = nearest_ldi,
      Secondary_nearest = nearest_secondary,
      Tertiary_nearest = nearest_tertiary
    )
    
    # Update joined_df with nearest neighbor values
    joined_df <- joined_df %>%
      left_join(na_lookup, by = c("Longitude", "Latitude")) %>%
      mutate(
        # Use nearest neighbor value if original LDI is NA
        LDI = ifelse(is.na(LDI), LDI_nearest, LDI),
        Secondary = ifelse(is.na(Secondary), Secondary_nearest, Secondary),
        Tertiary = ifelse(is.na(Tertiary), Tertiary_nearest, Tertiary)
      ) %>%
      select(-LDI_nearest, -Secondary_nearest, -Tertiary_nearest)
    
    cat("   Nearest neighbor assignment complete\n")
    cat("   Points with LDI after nearest neighbor:", sum(!is.na(joined_df$LDI)), "\n")
    cat("   Points still without LDI:", sum(is.na(joined_df$LDI)), "\n")
  }
}

# =============================================================================
# 10. Merge LDI values back to original soil data
# =============================================================================
cat("\n8. Merging LDI values back to original data...\n")

# Rename columns to avoid confusion
joined_df <- joined_df %>%
  rename(LandUse_Secondary = Secondary,
         LandUse_Tertiary = Tertiary)

# Merge with original data
result_df <- soil_df %>%
  left_join(joined_df[, c("Longitude", "Latitude", "LDI", 
                          "LandUse_Secondary", "LandUse_Tertiary")],
            by = c("Longitude", "Latitude"))

# =============================================================================
# 11. Save output
# =============================================================================
cat("\n9. Saving results...\n")
write_csv(result_df, csv_output_path)
cat("   Saved to:", csv_output_path, "\n")

# =============================================================================
# 12. Summary
# =============================================================================
cat("\n", "=" , rep("=", 59), "\n", sep = "")
cat("SUMMARY\n")
cat("=" , rep("=", 59), "\n", sep = "")
cat("Total records:", nrow(result_df), "\n")
cat("Records with LDI:", sum(!is.na(result_df$LDI)), "\n")
cat("Records without LDI:", sum(is.na(result_df$LDI)), "\n")
cat("\nLDI distribution in output:\n")
print(table(result_df$LDI, useNA = "ifany"))

cat("\nLDI Classification Reference:\n")
cat("  LDI 1: No effective disturbance (national parks, reserves)\n")
cat("  LDI 2: Limited disturbance (production forestry)\n")
cat("  LDI 3: Moderate disturbance (grazing native vegetation)\n")
cat("  LDI 4: High disturbance (modified pastures)\n")
cat("  LDI 5: Very high disturbance (horticulture, orchards)\n")
cat("  LDI 6: Extreme disturbance (cropping, urban)\n")

cat("\n✅ Complete!\n")
