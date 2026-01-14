# ==============================================================================
# Script: Convert ANSIS JSON to TERN Soil Data Federator CSV Format
# Description: Converts ANSIS JSON soil data files to a CSV format similar to 
#              TERNSoilDataFederator output. Processes all soil properties.
# Author: Digital Soil Mapping Project
# Date: 2026-01-12
# Updated: 2026-01-14 - Fixed path handling, JSON parsing, performance issues
# ==============================================================================

# Load required libraries
library(jsonlite)
library(stringr)
library(dplyr)

# ==============================================================================
# Utility: Safe extraction from nested lists
# ==============================================================================

#' Safely extract a value from a nested list structure
#' @param x List to extract from
#' @param ... Keys to traverse
#' @param default Default value if extraction fails
#' @return Extracted value or default
safe_extract <- function(x, ..., default = NULL) {
  keys <- list(...)
  tryCatch({
    for (key in keys) {
      if (is.null(x)) return(default)
      if (is.numeric(key)) {
        if (length(x) < key) return(default)
        x <- x[[key]]
      } else {
        if (!key %in% names(x)) return(default)
        x <- x[[key]]
      }
    }
    if (is.null(x)) default else x
  }, error = function(e) default)
}

# ==============================================================================
# Helper Functions (adapted from ANSISUtils-main)
# ==============================================================================

#' Get Site ID from ANSIS site data
#' @param siteAsList List containing site data
#' @return Character string with site ID
getSiteID <- function(siteAsList) {
  scopeID <- safe_extract(siteAsList, "scopedIdentifier", 1, "value")
  if (is.null(scopeID)) {
    sid <- safe_extract(siteAsList, "id", default = paste0("unknown_", Sys.time()))
  } else {
    auth <- safe_extract(siteAsList, "scopedIdentifier", 1, "authority", default = "unknown")
    sid <- paste0(auth, "+", scopeID)
  }
  return(sid)
}

#' Get Site Location (coordinates) from ANSIS site data
#' @param siteAsList List containing site data
#' @return List with X (Longitude) and Y (Latitude)
getSiteLocation <- function(siteAsList) {
  geom <- safe_extract(siteAsList, "geometry")
  if (is.null(geom)) {
    return(list(X = NA_real_, Y = NA_real_))
  }
  
  # Handle different geometry formats
  srid <- NULL
  if (is.character(geom)) {
    srid <- geom
  } else if (is.list(geom) && !is.null(geom$result)) {
    srid <- geom$result
  } else if (is.list(geom) && length(geom) > 0) {
    srid <- geom[[1]]
  }
  
  if (is.null(srid) || !is.character(srid)) {
    return(list(X = NA_real_, Y = NA_real_))
  }
  
  # Parse POINT coordinates from WKT format

  # Format: "SRID=4283;POINT(143.58154632872 -35.0968873016132)"
  # Also handles: "POINT(143.58 -35.09)" without SRID prefix

  tryCatch({
    # Check if it looks like WKT POINT
    if (!grepl("POINT", srid, ignore.case = TRUE)) {
      return(list(X = NA_real_, Y = NA_real_))
    }
    
    # Extract coordinates between parentheses
    coords_match <- regmatches(srid, regexec("POINT\\s*\\(\\s*([^)]+)\\s*\\)", srid, ignore.case = TRUE))
    if (length(coords_match[[1]]) < 2) {
      return(list(X = NA_real_, Y = NA_real_))
    }
    
    coord_str <- trimws(coords_match[[1]][2])
    coords <- as.numeric(strsplit(coord_str, "\\s+")[[1]])
    
    if (length(coords) >= 2 && !any(is.na(coords[1:2]))) {
      return(list(X = coords[1], Y = coords[2]))
    }
    return(list(X = NA_real_, Y = NA_real_))
  }, error = function(e) {
    return(list(X = NA_real_, Y = NA_real_))
  })
}

#' Get sample date from site visit
#' @param siteAsList List containing site data
#' @param visit_index Which site visit to use (default 1)
#' @return Character string with date
getSampleDate <- function(siteAsList, visit_index = 1) {
  tryCatch({
    # Check if siteVisit exists and has entries
    site_visits <- safe_extract(siteAsList, "siteVisit")
    if (is.null(site_visits) || length(site_visits) == 0) {
      return(NA_character_)
    }
    
    # Clamp visit_index to valid range
    visit_index <- min(visit_index, length(site_visits))
    dt <- safe_extract(site_visits, visit_index, "startedAtTime")
    
    if (!is.null(dt) && is.character(dt) && nzchar(dt)) {
      # Convert from ISO format to dd-mm-yyyy
      date_str <- strsplit(dt, "T")[[1]][1]
      date_obj <- as.Date(date_str)
      return(format(date_obj, "%d-%m-%Y"))
    }
    return(NA_character_)
  }, error = function(e) {
    return(NA_character_)
  })
}

#' Check if a layer property is a lab measurement
#' @param layer List containing layer data
#' @param prop Property name to check
#' @return Logical TRUE if it's a lab property with valid usedProcedure
isLabProperty <- function(layer, prop) {
  tryCatch({
    prop_data <- safe_extract(layer, prop)
    if (is.null(prop_data) || length(prop_data) == 0) {
      return(FALSE)
    }
    
    # Check first element for usedProcedure
    used_proc <- safe_extract(prop_data, 1, "usedProcedure")
    
    # Must be non-null and non-empty string
    if (!is.null(used_proc) && is.character(used_proc) && nzchar(used_proc)) {
      return(TRUE)
    }
    return(FALSE)
  }, error = function(e) {
    return(FALSE)
  })
}

#' Extract method code from usedProcedure URI/string
#' @param procedure_str The raw usedProcedure string (may be URI, prefixed code, etc.)
#' @return Clean method code (e.g., "4B2")
extractMethodCode <- function(procedure_str) {
  if (is.null(procedure_str) || !is.character(procedure_str) || !nzchar(procedure_str)) {
    return(NA_character_)
  }
  
  # Remove common prefixes
  code <- procedure_str
  code <- gsub("^scm:", "", code)
  code <- gsub("^spmile:", "", code)
  
  # If it's a URI, extract the last path segment
  if (grepl("^https?://", code) || grepl("/", code)) {
    parts <- strsplit(code, "/")[[1]]
    code <- tail(parts[parts != ""], 1)
    if (length(code) == 0) code <- procedure_str
  }
  
  # Remove any remaining prefixes like "method:" or "procedure:"
  code <- gsub("^(method|procedure|code):", "", code, ignore.case = TRUE)
  
  # Extract just the method code pattern (e.g., "4B2", "6A1", "15A1")
  # Pattern: digits followed by letter(s) optionally followed by more digits
  method_match <- regmatches(code, regexec("([0-9]+[A-Za-z]+[0-9]*)", code))
  if (length(method_match[[1]]) >= 2) {
    return(toupper(method_match[[1]][2]))
  }
  
  # If no standard pattern found, return cleaned string
  return(trimws(code))
}

#' Extract lab values from a soil layer
#' @param layer List containing layer data
#' @param prop Property name
#' @param ud Upper depth (in metres from ANSIS)
#' @param ld Lower depth (in metres from ANSIS)
#' @param depth_unit Target depth unit ("cm" or "m")
#' @return Data frame with extracted values
extractLabValues <- function(layer, prop, ud, ld, depth_unit = "m") {
  results <- list()
  
  tryCatch({
    prop_data <- safe_extract(layer, prop)
    
    if (is.null(prop_data) || length(prop_data) == 0) return(data.frame())
    
    # Convert depths from metres to cm if needed
    if (depth_unit == "cm") {
      ud_out <- if (!is.na(ud)) ud * 100 else NA_real_
      ld_out <- if (!is.na(ld)) ld * 100 else NA_real_
    } else {
      ud_out <- ud
      ld_out <- ld
    }
    
    # Handle multiple measurements for the same property
    for (i in seq_along(prop_data)) {
      measurement <- prop_data[[i]]
      if (is.null(measurement)) next
      
      used_proc <- safe_extract(measurement, "usedProcedure")
      if (is.null(used_proc)) next
      
      procedure <- extractMethodCode(used_proc)
      
      value <- safe_extract(measurement, "result", "value")
      unit <- safe_extract(measurement, "result", "unit")
      
      if (!is.null(value)) {
        results[[length(results) + 1]] <- data.frame(
          UpperDepth = ud_out,
          LowerDepth = ld_out,
          DepthUnit = depth_unit,
          PropertyType = "LaboratoryMeasurement",
          ObservedProperty = procedure,
          RawProcedure = used_proc,  # Keep original for debugging
          Value = as.numeric(value),
          Units = if (is.null(unit)) NA_character_ else gsub("^unit:", "", unit),
          stringsAsFactors = FALSE
        )
      }
    }
  }, error = function(e) {
    # Silent error handling
  })
  
  if (length(results) == 0) return(data.frame())
  bind_rows(results)
}

#' Parse all soil layers from a site
#' @param siteAsList List containing site data
#' @param visit_index Which site visit to process (default 1, use NULL for all)
#' @param depth_unit Target depth unit ("cm" or "m")
#' @return Data frame with all layer data
parseSoilLayers <- function(siteAsList, visit_index = 1, depth_unit = "m") {
  all_results <- list()
  
  tryCatch({
    site_visits <- safe_extract(siteAsList, "siteVisit")
    
    if (is.null(site_visits) || length(site_visits) == 0) {
      return(data.frame())
    }
    
    # Determine which visits to process
    if (is.null(visit_index)) {
      visit_indices <- seq_along(site_visits)
    } else {
      visit_indices <- min(visit_index, length(site_visits))
    }
    
    for (v_idx in visit_indices) {
      soil_profile <- safe_extract(site_visits, v_idx, "soilProfile")
      if (is.null(soil_profile) || length(soil_profile) == 0) next
      
      # Handle multiple soil profiles if present
      for (p_idx in seq_along(soil_profile)) {
        slsl <- safe_extract(soil_profile, p_idx, "soilLayer")
        if (is.null(slsl) || length(slsl) == 0) next
        
        for (i in seq_along(slsl)) {
          layer <- slsl[[i]]
          if (is.null(layer)) next
          
          # Get depth values (ANSIS stores in metres)
          ud <- safe_extract(layer, "depthUpper", "result", "value")
          ld <- safe_extract(layer, "depthLower", "result", "value")
          
          if (is.null(ud)) ud <- NA_real_
          if (is.null(ld)) ld <- NA_real_
          
          # Convert to numeric
          ud <- as.numeric(ud)
          ld <- as.numeric(ld)
          
          # Get all property names in this layer
          propNames <- names(layer)
          if (is.null(propNames)) next
          
          # Known lab properties (case-sensitive as per ANSIS schema)
          labProps <- c("ph", "electricalConductivity", "organicCarbon", 
                        "totalOrganicCarbon", "cationExchangeCapacity",
                        "exchangeableCalcium", "exchangeableMagnesium",
                        "exchangeablePotassium", "exchangeableSodium",
                        "exchangeableAluminium", "extractablePhosphorus",
                        "moistureContent", "bulkDensity", "clay", "silt", "sand",
                        "totalNitrogen", "dryAggregates", "linearShrinkage",
                        "soilWaterCharacteristic", "texture", "coarseFragments",
                        "fieldTexture", "pH", "EC")
          
          for (prop in propNames) {
            if (prop %in% labProps || isLabProperty(layer, prop)) {
              labVals <- extractLabValues(layer, prop, ud, ld, depth_unit)
              if (nrow(labVals) > 0) {
                labVals$VisitIndex <- v_idx
                labVals$ProfileIndex <- p_idx
                labVals$LayerIndex <- i
                all_results[[length(all_results) + 1]] <- labVals
              }
            }
          }
        }
      }
    }
  }, error = function(e) {
    message(paste("Error parsing layers:", e$message))
  })
  
  if (length(all_results) == 0) return(data.frame())
  bind_rows(all_results)
}

# ==============================================================================
# Main Conversion Function
# ==============================================================================

#' Convert ANSIS JSON to TERN Soil Data Federator CSV format
#' @param json_file Path to the input JSON file
#' @param output_file Path for the output CSV file
#' @param observed_property Optional: filter to specific property code (e.g., "4B2" for pH)
#' @param depth_unit Target depth unit ("cm" or "m", default "m")
#' @return Data frame with converted data
convertANSIStoTERNFormat <- function(json_file, output_file = NULL, 
                                      observed_property = NULL,
                                      depth_unit = "m") {
  
  message("Reading ANSIS JSON file: ", json_file)
  
  # Validate file exists
  if (!file.exists(json_file)) {
    stop("JSON file not found: ", json_file)
  }
  
  # Read JSON file
  json_data <- fromJSON(json_file, simplifyDataFrame = FALSE)
  
  # Get number of sites - handle different JSON structures
  sites_data <- safe_extract(json_data, "data")
  if (is.null(sites_data) || length(sites_data) == 0) {
    # Try alternative structure (array at root)
    if (is.list(json_data) && length(json_data) > 0 && is.null(names(json_data))) {
      sites_data <- json_data
    } else {
      sites_data <- list()
    }
  }
  nsites <- length(sites_data)
  
  if (nsites == 0) {
    message("Warning: No sites found in JSON file")
    return(data.frame())
  }
  
  message(paste("Processing", nsites, "sites..."))
  
  # Collect all rows in a list (much faster than rbind in loop)
  all_rows <- list()
  row_counter <- 0
  
  # Process each site
  for (k in seq_len(nsites)) {
    site <- sites_data[[k]]
    if (is.null(site)) next
    
    # Get site info
    siteID <- getSiteID(site)
    location <- getSiteLocation(site)
    sampleDate <- getSampleDate(site)
    
    # Parse soil layers for lab values
    layerData <- parseSoilLayers(site, depth_unit = depth_unit)
    
    if (nrow(layerData) > 0) {
      # Create unique sample ID per site
      sample_id <- paste0("S", k)
      
      for (i in seq_len(nrow(layerData))) {
        row <- layerData[i, ]
        row_counter <- row_counter + 1
        
        # Create unique observation ID
        obs_id <- paste0(siteID, "_", 
                         row$VisitIndex, "_",
                         row$ProfileIndex, "_",
                         row$LayerIndex, "_",
                         row$ObservedProperty, "_",
                         row_counter)
        
        # Create unique layer ID
        layer_id <- paste0("V", row$VisitIndex, 
                           "_P", row$ProfileIndex,
                           "_L", row$LayerIndex)
        
        all_rows[[row_counter]] <- data.frame(
          DataStore = "ANSIS",
          Dataset = "ANSIS_Federator",
          Provider = "ANSIS",
          Observation_ID = obs_id,
          SampleID = sample_id,
          SampleDate = sampleDate,
          Longitude = location$X,
          Latitude = location$Y,
          UpperDepth = row$UpperDepth,
          LowerDepth = row$LowerDepth,
          DepthUnit = row$DepthUnit,
          PropertyType = row$PropertyType,
          ObservedProperty = row$ObservedProperty,
          RawProcedure = row$RawProcedure,
          Value = row$Value,
          Units = row$Units,
          QualCollection = NA_character_,
          QualSpatialAggregation = NA_character_,
          QualManagement = NA_character_,
          QualSpatialAccuracy = NA_character_,
          Location_ID = siteID,
          Layer_ID = layer_id,
          ExtractTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Progress indicator
    if (k %% 50 == 0 || k == nsites) {
      message(paste("Processed", k, "of", nsites, "sites"))
    }
  }
  
  # Combine all rows at once (much faster than iterative rbind)
  if (length(all_rows) == 0) {
    message("Warning: No lab measurements found in any sites")
    return(data.frame())
  }
  
  output_df <- bind_rows(all_rows)
  
  # Filter by observed property if specified
  if (!is.null(observed_property)) {
    # Support partial matching for flexibility
    matches <- grepl(observed_property, output_df$ObservedProperty, ignore.case = TRUE)
    output_df <- output_df[matches, ]
    message(paste("Filtered to property:", observed_property, 
                  "- Rows:", nrow(output_df)))
  }
  
  # Save to CSV if output file specified (no row names)
  if (!is.null(output_file)) {
    write.csv(output_df, output_file, row.names = FALSE)
    message(paste("Data saved to:", output_file))
  }
  
  message(paste("Conversion complete! Total rows:", nrow(output_df)))
  return(output_df)
}

# ==============================================================================
# Specific Property Extraction Functions
# ==============================================================================

#' Extract pH data (4B1, 4B2 methods) from ANSIS JSON
#' @param json_file Path to the input JSON file
#' @param output_file Path for the output CSV file
#' @param method pH method code (e.g., "4B1", "4B2", or NULL for all)
#' @return Data frame with pH data
extractPhData <- function(json_file, output_file = NULL, method = NULL) {
  
  all_data <- convertANSIStoTERNFormat(json_file)
  
  # Filter for pH measurements (use partial matching for flexibility)
  ph_methods <- c("4A1", "4B1", "4B2", "4B3")
  ph_pattern <- paste(ph_methods, collapse = "|")
  ph_data <- all_data[grepl(ph_pattern, all_data$ObservedProperty, ignore.case = TRUE), ]
  
  if (!is.null(method)) {
    ph_data <- ph_data[grepl(method, ph_data$ObservedProperty, ignore.case = TRUE), ]
  }
  
  if (!is.null(output_file)) {
    write.csv(ph_data, output_file, row.names = FALSE)
    message(paste("pH data saved to:", output_file))
  }
  
  return(ph_data)
}

# ==============================================================================
# Main Execution
# ==============================================================================

# ==============================================================================
# Main Execution - Process Soil Property Data
# ==============================================================================

#' Determine project root directory robustly
#' @return Absolute path to project root
getProjectRoot <- function() {
  # Helper to normalize paths consistently across platforms
  normPath <- function(p) {
    if (is.null(p) || !nzchar(p)) return(NULL)
    tryCatch({
      # Convert backslashes to forward slashes first (Windows compatibility)
      p <- gsub("\\\\", "/", p)
      # Remove trailing slashes
      p <- sub("/$", "", p)
      # Use normalizePath if path exists, otherwise return cleaned path
      if (file.exists(p)) {
        normalizePath(p, winslash = "/", mustWork = FALSE)
      } else {
        p
      }
    }, error = function(e) p)
  }
  
  # Helper to check if directory looks like project root
  isProjectRoot <- function(dir) {
    if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(FALSE)
    markers <- c("Data/data_in/ansis_data", "R/0_convert_ansis_json_to_csv.R")
    all(sapply(markers, function(m) {
      test_path <- file.path(dir, m)
      # Also try with backslash conversion for Windows
      file.exists(test_path) || file.exists(gsub("/", "\\\\", test_path))
    }))
  }
  
  # Method 1: Try rstudioapi if available and running in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    tryCatch({
      if (rstudioapi::isAvailable()) {
        script_path <- rstudioapi::getSourceEditorContext()$path
        if (!is.null(script_path) && nzchar(script_path)) {
          script_path <- normPath(script_path)
          script_dir <- dirname(script_path)
          
          # Check if we're in R/ subdirectory (handle both / and \)
          if (grepl("[/\\\\]R$", script_dir, ignore.case = FALSE)) {
            candidate <- dirname(script_dir)
            candidate <- normPath(candidate)
            if (isProjectRoot(candidate)) {
              return(candidate)
            }
          }
          
          # Maybe script is directly in project root
          if (isProjectRoot(script_dir)) {
            return(normPath(script_dir))
          }
        }
      }
    }, error = function(e) NULL)
  }
  
  # Method 2: Check current working directory and parents
  cwd <- normPath(getwd())
  
  # Check current directory

  if (isProjectRoot(cwd)) {
    return(cwd)
  }
  
  # Check parent directory (if running from R/ subdirectory)
  parent <- dirname(cwd)
  if (isProjectRoot(parent)) {
    return(normPath(parent))
  }
  
  # Check grandparent (in case of deeper nesting)
  grandparent <- dirname(parent)
  if (isProjectRoot(grandparent)) {
    return(normPath(grandparent))
  }
  
  # Method 3: Use here package if available
  if (requireNamespace("here", quietly = TRUE)) {
    tryCatch({
      here_root <- normPath(here::here())
      if (isProjectRoot(here_root)) {
        return(here_root)
      }
    }, error = function(e) NULL)
  }
  
  # Method 4: Try to find project root by walking up directory tree
  search_dir <- cwd
  for (i in 1:10) {  # Limit search depth
    if (isProjectRoot(search_dir)) {
      return(normPath(search_dir))
    }
    parent_dir <- dirname(search_dir)
    # Stop if we've reached filesystem root
    if (parent_dir == search_dir) break
    search_dir <- parent_dir
  }
  
  # Fallback: return current directory with warning
  warning("Could not reliably determine project root. Using current directory: ", cwd)
  return(cwd)
}

# Only run if this script is being executed directly (not sourced)
if (sys.nframe() == 0 || !exists("SOURCED_ONLY")) {
  
  # Helper function to check if path is absolute
  isAbsolutePath <- function(path) {
    if (is.null(path) || !nzchar(path)) return(FALSE)
    # Windows: starts with drive letter (C:/, D:\, etc.)
    # Unix/Mac: starts with /
    grepl("^[A-Za-z]:[/\\\\]", path) || grepl("^/", path)
  }
  
  # Helper function to resolve path (handle both absolute and relative)
  resolvePath <- function(path, base_dir) {
    if (is.null(path) || !nzchar(path)) return(NULL)
    # Normalize slashes
    path <- gsub("\\\\", "/", path)
    # Remove leading slash if it's not a Unix absolute path and not a drive letter
    # This handles cases like "/Data/..." which should be relative
    if (grepl("^/[^/]", path) && !grepl("^/[A-Za-z]:", path)) {
      # Check if it looks like a relative path with accidental leading slash
      if (!dir.exists(path)) {
        path <- sub("^/", "", path)
      }
    }
    # If already absolute, use as-is
    if (isAbsolutePath(path)) {
      return(normalizePath(path, winslash = "/", mustWork = FALSE))
    }
    # Otherwise, combine with base directory
    return(normalizePath(file.path(base_dir, path), winslash = "/", mustWork = FALSE))
  }
  
  # Determine and set project root
  project_root <- getProjectRoot()
  setwd(project_root)
  message("Working directory: ", getwd())
  
  # ============================================================================
  # CONFIGURATION - Change these paths for each property
  # Paths can be:
  #   - Absolute: "D:/ANSIS_DATA/..." or "/home/user/data/..."
  #   - Relative to project root: "Data/data_in/..." or "./Data/..."
  # ============================================================================
  
  # Input: Directory containing JSON files for the property
  json_dir <- "Data/data_in/ansis_data/New Baseline (1990-2020)/Moisture Content"
  
  # Output: Directory where CSV files will be saved
  output_dir <- "Data/data_out/ansis_lab_measurements/Moisture Content"
  
  # ============================================================================
  
  # Convert to absolute path and validate
  json_dir_abs <- resolvePath(json_dir, project_root)
  
  message("\n========================================")
  message("ANSIS JSON to CSV Conversion")
  message("========================================")
  message("Input directory (configured): ", json_dir)
  message("Input directory (resolved): ", json_dir_abs)
  
  # Check if directory exists
  if (!dir.exists(json_dir_abs)) {
    stop("Input directory does not exist: ", json_dir_abs, 
         "\nPlease check the path and ensure files are in place.")
  }
  
  # Find JSON files (case-insensitive for cross-platform support)
  json_files <- list.files(json_dir_abs, 
                           pattern = "\\.[jJ][sS][oO][nN]$", 
                           full.names = TRUE,
                           recursive = TRUE)  # Also check subdirectories
  
  if (length(json_files) == 0) {
    # List what IS in the directory for debugging
    all_files <- list.files(json_dir_abs, full.names = FALSE, recursive = TRUE)
    message("Files found in directory: ", 
            if (length(all_files) > 0) paste(head(all_files, 10), collapse = ", ") else "(empty)")
    stop("No JSON files found in: ", json_dir_abs)
  }
  
  message("Output directory: ", output_dir)
  message("JSON files found: ", length(json_files))
  for (f in json_files) {
    message("  - ", basename(f))
  }
  
  # Process all JSON files - collect in list for performance
  all_data_list <- list()
  
  for (input_json in json_files) {
    message("\nProcessing: ", basename(input_json))
    
    # Convert all properties from JSON
    json_data <- convertANSIStoTERNFormat(json_file = input_json)
    
    if (nrow(json_data) > 0) {
      # Include all LaboratoryMeasurement records
      lab_data <- json_data[json_data$PropertyType == "LaboratoryMeasurement", ]
      if (nrow(lab_data) > 0) {
        all_data_list[[length(all_data_list) + 1]] <- lab_data
      }
    }
  }
  
  # Combine all data at once (faster than iterative rbind)
  all_data <- if (length(all_data_list) > 0) bind_rows(all_data_list) else data.frame()
  
  if (nrow(all_data) == 0) {
    message("\nWarning: No laboratory measurements extracted from any files.")
    message("This may indicate:")
    message("  - JSON structure doesn't match expected ANSIS format")
    message("  - No lab measurements in the data")
    message("  - Property names don't match expected values")
  } else {
    # Create output directory if not exists
    output_dir_abs <- resolvePath(output_dir, project_root)
    if (!dir.exists(output_dir_abs)) {
      dir.create(output_dir_abs, recursive = TRUE)
    }
    
    # Get unique lab methods
    lab_methods <- unique(all_data$ObservedProperty)
    lab_methods <- lab_methods[!is.na(lab_methods)]
    
    message("\n=== Generating CSV files by method ===")
    
    # Save separate CSV for each lab method (no row names)
    for (method in lab_methods) {
      method_data <- all_data[all_data$ObservedProperty == method, ]
      # Clean method name for filename
      safe_method <- gsub("[^A-Za-z0-9_-]", "_", method)
      output_file <- file.path(output_dir_abs, paste0("ANSIS-", safe_method, ".csv"))
      write.csv(method_data, output_file, row.names = FALSE)
      message("Saved: ", basename(output_file), " (", nrow(method_data), " records)")
    }
    
    # Also save combined file
    output_csv_all <- file.path(output_dir_abs, "ANSIS_combined.csv")
    write.csv(all_data, output_csv_all, row.names = FALSE)
    
    # Display summary
    message("\n========================================")
    message("=== SUMMARY ===")
    message("========================================")
    message("JSON files processed: ", length(json_files))
    message("Total records extracted: ", nrow(all_data))
    message("CSV files generated: ", length(lab_methods) + 1)
    message("Output directory: ", output_dir_abs)
    message("\nLab methods found:")
    print(table(all_data$ObservedProperty, useNA = "ifany"))
  }
}
