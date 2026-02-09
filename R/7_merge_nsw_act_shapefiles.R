# =============================================================================
# Script: Merge NSW and ACT Land Use Shapefiles with LDI Calculation
# =============================================================================
# This script merges the NSW Land Use shapefile with the ACT Land Use shapefile,
# harmonizing column names and adding LDI values to both before merging.
# Can also work with a single shapefile if only one is available.
# =============================================================================

# --- Load required packages ---
library(sf)
library(dplyr)

# =============================================================================
# 1. Define file paths
# =============================================================================
base_dir <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"

nsw_shapefile <- file.path(base_dir, "Data/data_in/shp_files/NSW/NSWLanduse2017_Ver1_5_20230921.shp")
act_shapefile <- file.path(base_dir, "Data/data_in/shp_files/ACT/ACT_landuse_final_2012.shp")
merged_output <- file.path(base_dir, "Data/data_in/shp_files/NSW_ACT_Landuse_merged_new.shp")

# Check which files exist
nsw_exists <- file.exists(nsw_shapefile)
act_exists <- file.exists(act_shapefile)

cat("=" , rep("=", 69), "\n", sep = "")
cat("Processing Land Use Shapefiles with LDI Calculation\n")
cat("=" , rep("=", 69), "\n\n", sep = "")
cat("NSW shapefile exists:", nsw_exists, "\n")
cat("ACT shapefile exists:", act_exists, "\n")

if (!nsw_exists && !act_exists) {
  stop("ERROR: Neither NSW nor ACT shapefile found!")
}

# =============================================================================
# 2. Define LDI mapping for NSW Secondary classification
# =============================================================================
# LDI Classification (Gray et al., 2015):
# - LDI 1: No effective disturbance (e.g., nature reserve)
# - LDI 2: Limited disturbance (e.g., light native forestry)
# - LDI 3: Moderate disturbance (e.g., woodland grazing) - grazing with retained native woodland
# - LDI 4: High disturbance (e.g., grazing on pasture) - cleared native vegetation for pasture
# - LDI 5: Very high disturbance (e.g., cropping/grazing, horticulture)
# - LDI 6: Extreme disturbance (e.g., cropping)

create_ldi_mapping_nsw <- function() {
  ldi_map <- data.frame(
    Secondary = c(
      # LDI 1: No effective disturbance - Conservation and minimal use
      "1.1.0 Nature conservation",
      "1.2.0 Managed resource protection", 
      "1.3.0 Other minimal use",
      
      # LDI 2: Limited disturbance - Production forestry (selective logging)
      "2.2.0 Production native forestry",
      
      # LDI 3: Moderate disturbance - Woodland grazing (native vegetation retained)
      "2.1.0 Grazing native vegetation",
      
      # LDI 4: High disturbance - Grazing on pasture (native vegetation cleared), plantations
      "3.1.0 Plantation forests",
      "3.2.0 Grazing modified pastures",
      "3.6.0 Land in transition",
      "4.6.0 Irrigated land in transition",
      
      # LDI 5: Very high disturbance - Orchards, viticulture, horticulture
      "3.4.0 Perennial horticulture",
      "3.5.0 Seasonal horticulture",
      "4.2.0 Grazing irrigated modified pastures",
      "4.4.0 Irrigated perennial horticulture",
      "4.5.0 Irrigated seasonal horticulture",
      "5.1.0 Intensive horticulture",
      "5.2.0 Intensive animal production",
      
      # LDI 6: Extreme disturbance - Predominant cropping, urban/industrial
      "3.3.0 Cropping",
      "4.1.0 Irrigated plantation forests",
      "4.3.0 Irrigated cropping",
      "5.3.0 Manufacturing and industrial",
      "5.4.0 Residential and farm infrastructure",
      "5.5.0 Services",
      "5.6.0 Utilities",
      "5.7.0 Transport and communication",
      "5.8.0 Mining",
      "5.9.0 Waste treatment and disposal",
      
      # Water bodies (NA - excluded from analysis)
      "6.1.0 Lake",
      "6.2.0 Reservoir/dam",
      "6.3.0 River",
      "6.4.0 Channel/aqueduct",
      "6.5.0 Marsh/wetland",
      "6.6.0 Estuary/coastal waters"
    ),
    LDI = c(
      # LDI 1: Conservation
      1, 1, 1,
      # LDI 2: Production forestry
      2,
      # LDI 3: Woodland grazing
      3,
      # LDI 4: Grazing on pasture, plantations
      4, 4, 4, 4,
      # LDI 5: Horticulture, orchards
      5, 5, 5, 5, 5, 5, 5,
      # LDI 6: Cropping, urban/industrial
      6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
      # Water bodies
      NA, NA, NA, NA, NA, NA
    ),
    stringsAsFactors = FALSE
  )
  return(ldi_map)
}

# =============================================================================
# 3. Define LDI mapping for ACT Secondary classification
# =============================================================================

create_ldi_mapping_act <- function() {
  ldi_map <- data.frame(
    SECONDARY_ = c(
      # LDI 1: No effective disturbance
      "1.1 Nature conservation",
      "1.2 Managed resource protection",
      "1.3 Other minimal use",
      
      # LDI 2: Limited disturbance - Production forestry
      "2.2 Production native forestry",
      
      # LDI 3: Moderate disturbance - Woodland grazing (native vegetation retained)
      "2.1 Grazing native vegetation",
      
      # LDI 4: High disturbance - Grazing on pasture (native vegetation cleared), plantations
      "3.1 Plantation forestry",
      "3.2 Grazing modified pastures",
      "3.6 Land in transition",
      "4.6 Irrigated land in transition",
      
      # LDI 5: Very high disturbance - Horticulture
      "3.4 Perennial horticulture",
      "3.5 Seasonal horticulture",
      "4.2 Grazing irrigated modified pastures",
      "4.4 Irrigated perennial horticulture",
      "4.5 Irrigated seasonal horticulture",
      "5.1 Intensive horticulture",
      "5.2 Intensive animal husbandry",
      
      # LDI 6: Extreme disturbance - Cropping, urban/industrial
      "3.3 Cropping",
      "4.1 Irrigated plantation forestry",
      "4.3 Irrigated cropping",
      "5.3 Manufacturing and industrial",
      "5.4 Residential and farm infrastructure",
      "5.5 Services",
      "5.6 Utilities",
      "5.7 Transport and communication",
      "5.8 Mining",
      "5.9 Waste treatment and disposal",
      
      # Water bodies (NA)
      "6.1 Lake",
      "6.2 Reservoir/dam",
      "6.3 River",
      "6.4 Channel/aqueduct",
      "6.5 Marsh/Wetland",
      "6.6 Estuary/coastal waters"
    ),
    LDI = c(
      # LDI 1
      1, 1, 1,
      # LDI 2
      2,
      # LDI 3: Woodland grazing
      3,
      # LDI 4: Grazing on pasture, plantations
      4, 4, 4, 4,
      # LDI 5
      5, 5, 5, 5, 5, 5, 5,
      # LDI 6
      6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
      # Water
      NA, NA, NA, NA, NA, NA
    ),
    stringsAsFactors = FALSE
  )
  return(ldi_map)
}

# =============================================================================
# 4. Function to calculate LDI from code prefix (fallback for unmapped values)
# =============================================================================
# Fallback LDI based on land use code pattern (Gray et al., 2015):

calculate_ldi_from_code <- function(secondary_code) {
  if (is.na(secondary_code) || secondary_code == "") return(NA)
  
  code <- trimws(secondary_code)
  
 # Extract the numeric code (e.g., "1.1", "3.3", "5.4")
  code_pattern <- regmatches(code, regexpr("^[0-9]+\\.[0-9]+", code))
  if (length(code_pattern) == 0) {
    # Try just first digit
    first_char <- substr(code, 1, 1)
    ldi <- case_when(
      first_char == "1" ~ 1L,   # Conservation -> LDI 1
      first_char == "2" ~ 2L,   # Native forestry -> LDI 2
      first_char == "3" ~ 4L,   # Dryland agriculture -> LDI 4 (default)
      first_char == "4" ~ 5L,   # Irrigated -> LDI 5
      first_char == "5" ~ 6L,   # Intensive -> LDI 6
      first_char == "6" ~ NA_integer_,  # Water -> NA
      TRUE ~ NA_integer_
    )
    return(ldi)
  }
  
  code_num <- code_pattern[1]
  
  # Detailed mapping based on secondary code
  ldi <- case_when(
    # 1.x: Conservation -> LDI 1
    grepl("^1\\.", code_num) ~ 1L,
    
    # 2.1: Grazing native vegetation (woodland grazing) -> LDI 3
    grepl("^2\\.1", code_num) ~ 3L,
    # 2.2: Production forestry -> LDI 2
    grepl("^2\\.2", code_num) ~ 2L,
    # 2.x other -> LDI 2
    grepl("^2\\.", code_num) ~ 2L,
    
    # 3.1: Plantation -> LDI 4 (native vegetation cleared)
    grepl("^3\\.1", code_num) ~ 4L,
    # 3.2: Modified pastures (grazing on pasture) -> LDI 4
    grepl("^3\\.2", code_num) ~ 4L,
    # 3.3: Cropping -> LDI 6
    grepl("^3\\.3", code_num) ~ 6L,
    # 3.4, 3.5: Horticulture -> LDI 5
    grepl("^3\\.[45]", code_num) ~ 5L,
    # 3.6: Land in transition -> LDI 4
    grepl("^3\\.6", code_num) ~ 4L,
    # 3.x other -> LDI 4
    grepl("^3\\.", code_num) ~ 4L,
    
    # 4.1: Irrigated plantation -> LDI 6
    grepl("^4\\.1", code_num) ~ 6L,
    # 4.2: Irrigated pastures -> LDI 5
    grepl("^4\\.2", code_num) ~ 5L,
    # 4.3: Irrigated cropping -> LDI 6
    grepl("^4\\.3", code_num) ~ 6L,
    # 4.4, 4.5: Irrigated horticulture -> LDI 5
    grepl("^4\\.[45]", code_num) ~ 5L,
    # 4.6: Land in transition -> LDI 4
    grepl("^4\\.6", code_num) ~ 4L,
    # 4.x other -> LDI 5
    grepl("^4\\.", code_num) ~ 5L,
    
    # 5.x: Intensive uses -> LDI 6
    grepl("^5\\.", code_num) ~ 6L,
    
    # 6.x: Water -> NA
    grepl("^6\\.", code_num) ~ NA_integer_,
    
    TRUE ~ NA_integer_
  )
  return(ldi)
}

# Vectorized version
calculate_ldi_from_code_vec <- Vectorize(calculate_ldi_from_code)

# =============================================================================
# 5. Load shapefiles (only those that exist)
# =============================================================================

nsw_cols <- NULL
act_cols <- NULL

# Load NSW if exists
if (nsw_exists) {
  cat("\n1. Loading NSW Land Use shapefile...\n")
  nsw_sf <- st_read(nsw_shapefile, quiet = TRUE)
  cat("   Loaded", nrow(nsw_sf), "NSW polygons\n")
  cat("   CRS:", st_crs(nsw_sf)$input, "\n")
  
  # Add LDI values
  cat("   Adding LDI values to NSW...\n")
  ldi_nsw <- create_ldi_mapping_nsw()
  if ("LDI" %in% names(nsw_sf)) {
    nsw_sf <- nsw_sf[, setdiff(names(nsw_sf), "LDI")]
  }
  nsw_sf <- nsw_sf %>% left_join(ldi_nsw, by = "Secondary")
  
  # Handle unmapped Secondary categories using fallback logic
  unmapped_nsw <- sum(is.na(nsw_sf$LDI) & !grepl("^6\\.", nsw_sf$Secondary))
  if (unmapped_nsw > 0) {
    cat("   ⚠️ Found", unmapped_nsw, "polygons with unmapped Secondary codes\n")
    cat("   Applying fallback LDI calculation based on primary category...\n")
    
    # Get unmapped categories
    unmapped_cats <- nsw_sf %>%
      st_drop_geometry() %>%
      filter(is.na(LDI) & !grepl("^6\\.", Secondary)) %>%
      pull(Secondary) %>%
      unique()
    cat("   Unmapped categories:\n")
    print(unmapped_cats)
    
    # Apply fallback calculation for unmapped values (excluding water bodies starting with 6.)
    nsw_sf <- nsw_sf %>%
      mutate(LDI = ifelse(is.na(LDI) & !grepl("^6\\.", Secondary),
                          calculate_ldi_from_code_vec(Secondary),
                          LDI))
  }
  
  # Harmonize columns
  nsw_cols <- nsw_sf %>%
    dplyr::transmute(
      Secondary = Secondary,
      Tertiary = Tertiary,
      Area_Ha = Area_Ha,
      LDI = LDI,
      Source = "NSW",
      geometry = geometry
    )
  cat("   NSW LDI assigned:", sum(!is.na(nsw_cols$LDI)), "/", nrow(nsw_cols), "polygons\n")
}

# Load ACT if exists
if (act_exists) {
  cat("\n2. Loading ACT Land Use shapefile...\n")
  act_sf <- st_read(act_shapefile, quiet = TRUE)
  cat("   Loaded", nrow(act_sf), "ACT polygons\n")
  cat("   CRS:", st_crs(act_sf)$input, "\n")
  
  # Add LDI values
  cat("   Adding LDI values to ACT...\n")
  ldi_act <- create_ldi_mapping_act()
  act_sf <- act_sf %>% left_join(ldi_act, by = "SECONDARY_")
  
  # Handle unmapped Secondary categories using fallback logic
  unmapped_act <- sum(is.na(act_sf$LDI) & !grepl("^6\\.", act_sf$SECONDARY_))
  if (unmapped_act > 0) {
    cat("   ⚠️ Found", unmapped_act, "polygons with unmapped Secondary codes\n")
    cat("   Applying fallback LDI calculation based on primary category...\n")
    
    # Get unmapped categories
    unmapped_cats <- act_sf %>%
      st_drop_geometry() %>%
      filter(is.na(LDI) & !grepl("^6\\.", SECONDARY_)) %>%
      pull(SECONDARY_) %>%
      unique()
    cat("   Unmapped categories:\n")
    print(unmapped_cats)
    
    # Apply fallback calculation
    act_sf <- act_sf %>%
      mutate(LDI = ifelse(is.na(LDI) & !grepl("^6\\.", SECONDARY_),
                          calculate_ldi_from_code_vec(SECONDARY_),
                          LDI))
  }
  
  # Harmonize columns
  act_cols <- act_sf %>%
    dplyr::transmute(
      Secondary = SECONDARY_,
      Tertiary = TERTIARY_V,
      Area_Ha = Area_ha,
      LDI = LDI,
      Source = "ACT",
      geometry = geometry
    )
  cat("   ACT LDI assigned:", sum(!is.na(act_cols$LDI)), "/", nrow(act_cols), "polygons\n")
}

# =============================================================================
# 6. Merge or use single shapefile
# =============================================================================
cat("\n3. Preparing final output...\n")

if (!is.null(nsw_cols) && !is.null(act_cols)) {
  # Both exist - transform ACT to match NSW CRS and merge
  cat("   Transforming ACT CRS to match NSW...\n")
  act_cols <- st_transform(act_cols, st_crs(nsw_cols))
  
  cat("   Merging NSW and ACT shapefiles...\n")
  merged_sf <- rbind(nsw_cols, act_cols)
  cat("   Total merged polygons:", nrow(merged_sf), "\n")
  cat("   - NSW polygons:", sum(merged_sf$Source == "NSW"), "\n")
  cat("   - ACT polygons:", sum(merged_sf$Source == "ACT"), "\n")
  
} else if (!is.null(nsw_cols)) {
  # Only NSW exists
  cat("   Using NSW shapefile only (ACT not found)\n")
  merged_sf <- nsw_cols
  cat("   Total polygons:", nrow(merged_sf), "\n")
  
} else {
  # Only ACT exists
  cat("   Using ACT shapefile only (NSW not found)\n")
  merged_sf <- act_cols
  cat("   Total polygons:", nrow(merged_sf), "\n")
}

# =============================================================================
# 7. Report LDI distribution
# =============================================================================
cat("\n4. LDI value distribution in final data:\n")
ldi_summary <- merged_sf %>%
  st_drop_geometry() %>%
  group_by(LDI, Source) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Source, values_from = Count, values_fill = 0)

# Calculate total based on which columns exist
if ("NSW" %in% names(ldi_summary) && "ACT" %in% names(ldi_summary)) {
  ldi_summary <- ldi_summary %>% mutate(Total = NSW + ACT)
} else if ("NSW" %in% names(ldi_summary)) {
  ldi_summary <- ldi_summary %>% mutate(Total = NSW)
} else {
  ldi_summary <- ldi_summary %>% mutate(Total = ACT)
}
ldi_summary <- ldi_summary %>% arrange(LDI)

print(as.data.frame(ldi_summary))

# Report any remaining unmapped categories (non-water with NA LDI)
remaining_na <- merged_sf %>%
  st_drop_geometry() %>%
  filter(is.na(LDI) & !grepl("^6\\.", Secondary))

if (nrow(remaining_na) > 0) {
  cat("\n   ⚠️ WARNING: Still have", nrow(remaining_na), "non-water polygons without LDI:\n")
  print(unique(remaining_na$Secondary))
}

# =============================================================================
# 8. Save output shapefile
# =============================================================================
cat("\n5. Saving output shapefile...\n")
cat("   Output:", merged_output, "\n")

st_write(merged_sf, merged_output, delete_dsn = TRUE)

# Final summary
sources_used <- unique(merged_sf$Source)
cat("\n" , rep("=", 69), "\n", sep = "")
cat("✅ SUCCESS! Land use shapefile processed\n")
cat("   Sources:", paste(sources_used, collapse = " + "), "\n")
cat("   Total polygons:", nrow(merged_sf), "\n")
cat("   Polygons with LDI values:", sum(!is.na(merged_sf$LDI)), "\n")
cat("   Water body polygons (LDI=NA):", sum(is.na(merged_sf$LDI)), "\n")
cat("   Output file:", merged_output, "\n")
cat(rep("=", 69), "\n", sep = "")
