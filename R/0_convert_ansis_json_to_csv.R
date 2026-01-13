# ==============================================================================
# Script: Convert ANSIS JSON to TERN Soil Data Federator CSV Format
# Description: Converts ANSIS JSON soil data files to a CSV format similar to 
#              TERNSoilDataFederator output. Processes all soil properties.
# Author: Digital Soil Mapping Project
# Date: 2026-01-12
# Updated: 2026-01-13 - Generic script for all soil property data
# ==============================================================================

# Load required libraries
library(jsonlite)
library(stringr)
library(dplyr)

# ==============================================================================
# Helper Functions (adapted from ANSISUtils-main)
# ==============================================================================

#' Get Site ID from ANSIS site data
#' @param siteAsList List containing site data
#' @return Character string with site ID
getSiteID <- function(siteAsList) {
  scopeID <- siteAsList$scopedIdentifier[[1]]$value
  if (is.null(scopeID)) {
    sid <- siteAsList$id
  } else {
    auth <- siteAsList$scopedIdentifier[[1]]$authority
    sid <- paste0(auth, "+", scopeID)
  }
  return(sid)
}

#' Get Site Location (coordinates) from ANSIS site data
#' @param siteAsList List containing site data
#' @return List with X (Longitude) and Y (Latitude)
getSiteLocation <- function(siteAsList) {
  geom <- siteAsList$geometry
  if (is.null(geom)) {
    return(list(X = NA, Y = NA))
  }
  
  # Handle different geometry formats
  if (is.character(geom)) {
    srid <- geom
  } else if (is.list(geom) && !is.null(geom$result)) {
    srid <- geom$result
  } else if (is.list(geom) && length(geom) > 0) {
    srid <- geom[[1]]
  } else {
    return(list(X = NA, Y = NA))
  }
  
  # Parse POINT coordinates from WKT format
  # Format: "SRID=4283;POINT(143.58154632872 -35.0968873016132)"
  tryCatch({
    bits <- str_split(srid, "[(]")
    bits2 <- str_remove(bits[[1]][2], "[)]")
    bits3 <- str_split(bits2, " ")[[1]]
    ol <- list()
    ol$X <- as.numeric(bits3[1])
    ol$Y <- as.numeric(bits3[2])
    return(ol)
  }, error = function(e) {
    return(list(X = NA, Y = NA))
  })
}

#' Get sample date from site visit
#' @param siteAsList List containing site data
#' @return Character string with date
getSampleDate <- function(siteAsList) {
  tryCatch({
    dt <- siteAsList$siteVisit[[1]]$startedAtTime
    if (!is.null(dt)) {
      # Convert from ISO format to dd-mm-yyyy
      date_obj <- as.Date(str_split(dt, "T")[[1]][1])
      return(format(date_obj, "%d-%m-%Y"))
    }
    return(NA)
  }, error = function(e) {
    return(NA)
  })
}

#' Check if a layer property is a lab measurement
#' @param layer List containing layer data
#' @param prop Property name to check
#' @return Logical TRUE if it's a lab property
isLabProperty <- function(layer, prop) {
  tryCatch({
    if (length(layer[[prop]][[1]]$usedProcedure) > -1 & 
        !is.null(layer[[prop]][[1]]$usedProcedure)) {
      return(TRUE)
    }
  }, error = function(e) {
    return(FALSE)
  })
  return(FALSE)
}

#' Extract lab values from a soil layer
#' @param layer List containing layer data
#' @param prop Property name
#' @param ud Upper depth
#' @param ld Lower depth
#' @return Data frame with extracted values
extractLabValues <- function(layer, prop, ud, ld) {
  results <- data.frame()
  
  tryCatch({
    prop_data <- layer[[prop]]
    
    if (is.null(prop_data)) return(results)
    
    # Handle multiple measurements for the same property
    for (i in 1:length(prop_data)) {
      measurement <- prop_data[[i]]
      
      if (!is.null(measurement$usedProcedure)) {
        procedure <- str_remove(measurement$usedProcedure, "scm:")
        procedure <- str_remove(procedure, "spmile:")
        
        value <- measurement$result$value
        unit <- measurement$result$unit
        
        if (!is.null(value)) {
          r <- data.frame(
            UpperDepth = ud,
            LowerDepth = ld,
            PropertyType = "LaboratoryMeasurement",
            ObservedProperty = procedure,
            Value = value,
            Units = ifelse(is.null(unit), NA, str_remove(unit, "unit:")),
            stringsAsFactors = FALSE
          )
          results <- rbind(results, r)
        }
      }
    }
  }, error = function(e) {
    # Silent error handling
  })
  
  return(results)
}

#' Parse all soil layers from a site
#' @param siteAsList List containing site data
#' @return Data frame with all layer data
parseSoilLayers <- function(siteAsList) {
  allResults <- data.frame()
  
  tryCatch({
    slsl <- siteAsList[["siteVisit"]][[1]]$soilProfile[[1]]$soilLayer
    
    if (is.null(slsl)) return(allResults)
    
    for (i in 1:length(slsl)) {
      layer <- slsl[[i]]
      
      # Get depth values (convert from meters to appropriate units)
      ud <- layer$depthUpper$result$value
      ld <- layer$depthLower$result$value
      
      if (is.null(ud)) ud <- NA
      if (is.null(ld)) ld <- NA
      
      # Get all property names
      propNames <- names(layer)
      
      # Filter to lab properties
      labProps <- c("ph", "electricalConductivity", "organicCarbon", 
                    "totalOrganicCarbon", "cationExchangeCapacity",
                    "exchangeableCalcium", "exchangeableMagnesium",
                    "exchangeablePotassium", "exchangeableSodium",
                    "exchangeableAluminium", "extractablePhosphorus",
                    "moistureContent", "bulkDensity", "clay", "silt", "sand",
                    "totalNitrogen", "dryAggregates", "linearShrinkage",
                    "soilWaterCharacteristic", "texture")
      
      for (prop in propNames) {
        if (prop %in% labProps || isLabProperty(layer, prop)) {
          labVals <- extractLabValues(layer, prop, ud, ld)
          if (nrow(labVals) > 0) {
            allResults <- rbind(allResults, labVals)
          }
        }
      }
    }
  }, error = function(e) {
    message(paste("Error parsing layers:", e$message))
  })
  
  return(allResults)
}

# ==============================================================================
# Main Conversion Function
# ==============================================================================

#' Convert ANSIS JSON to TERN Soil Data Federator CSV format
#' @param json_file Path to the input JSON file
#' @param output_file Path for the output CSV file
#' @param observed_property Optional: filter to specific property code (e.g., "4B2" for pH)
#' @return Data frame with converted data
convertANSIStoTERNFormat <- function(json_file, output_file = NULL, 
                                      observed_property = NULL) {
  
  message("Reading ANSIS JSON file...")
  
  # Read JSON file
  json_data <- fromJSON(json_file, simplifyDataFrame = FALSE)
  
  # Initialize output data frame matching TERN format
  output_df <- data.frame(
    DataStore = character(),
    Dataset = character(),
    Provider = character(),
    Observation_ID = character(),
    SampleID = character(),
    SampleDate = character(),
    Longitude = numeric(),
    Latitude = numeric(),
    UpperDepth = numeric(),
    LowerDepth = numeric(),
    PropertyType = character(),
    ObservedProperty = character(),
    Value = numeric(),
    Units = character(),
    QualCollection = character(),
    QualSpatialAggregation = character(),
    QualManagement = character(),
    QualSpatialAccuracy = character(),
    Location_ID = character(),
    Layer_ID = character(),
    ExtractTime = character(),
    stringsAsFactors = FALSE
  )
  
  # Get number of sites
  nsites <- length(json_data$data)
  message(paste("Processing", nsites, "sites..."))
  
  # Process each site
  for (k in 1:nsites) {
    site <- json_data$data[[k]]
    
    # Get site info
    siteID <- getSiteID(site)
    location <- getSiteLocation(site)
    sampleDate <- getSampleDate(site)
    
    # Parse soil layers for lab values
    layerData <- parseSoilLayers(site)
    
    if (nrow(layerData) > 0) {
      # Create layer IDs
      unique_depths <- unique(layerData[, c("UpperDepth", "LowerDepth")])
      unique_depths <- unique_depths[order(unique_depths$UpperDepth), ]
      
      for (i in 1:nrow(layerData)) {
        row <- layerData[i, ]
        
        # Find layer ID based on depth order
        layer_idx <- which(unique_depths$UpperDepth == row$UpperDepth & 
                             unique_depths$LowerDepth == row$LowerDepth)
        if (length(layer_idx) == 0) layer_idx <- i
        
        new_row <- data.frame(
          DataStore = "ANSIS",
          Dataset = "ANSIS_Federator",
          Provider = "ANSIS",
          Observation_ID = NA,
          SampleID = "1",
          SampleDate = sampleDate,
          Longitude = location$X,
          Latitude = location$Y,
          UpperDepth = row$UpperDepth,
          LowerDepth = row$LowerDepth,
          PropertyType = row$PropertyType,
          ObservedProperty = row$ObservedProperty,
          Value = row$Value,
          Units = row$Units,
          QualCollection = NA,
          QualSpatialAggregation = NA,
          QualManagement = NA,
          QualSpatialAccuracy = NA,
          Location_ID = siteID,
          Layer_ID = as.character(layer_idx),
          ExtractTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
          stringsAsFactors = FALSE
        )
        
        output_df <- rbind(output_df, new_row)
      }
    }
    
    # Progress indicator
    if (k %% 10 == 0 || k == nsites) {
      message(paste("Processed", k, "of", nsites, "sites"))
    }
  }
  
  # Filter by observed property if specified
  if (!is.null(observed_property)) {
    output_df <- output_df[output_df$ObservedProperty == observed_property, ]
    message(paste("Filtered to property:", observed_property, 
                  "- Rows:", nrow(output_df)))
  }
  
  # Save to CSV if output file specified
  if (!is.null(output_file)) {
    write.csv(output_df, output_file, row.names = TRUE)
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
  
  # Filter for pH measurements
  ph_methods <- c("4A1", "4B1", "4B2", "4B3")
  ph_data <- all_data[all_data$ObservedProperty %in% ph_methods, ]
  
  if (!is.null(method)) {
    ph_data <- ph_data[ph_data$ObservedProperty == method, ]
  }
  
  if (!is.null(output_file)) {
    write.csv(ph_data, output_file, row.names = TRUE)
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

# Only run if this script is being executed directly (not sourced)
if (sys.nframe() == 0 || !exists("SOURCED_ONLY")) {
  
  # Get the script directory and set working directory to project root
  script_dir <- tryCatch({
    dirname(rstudioapi::getSourceEditorContext()$path)
  }, error = function(e) {
    "."
  })
  
  # Set working directory to project root
  project_root <- tryCatch({
    if (grepl("/R$", script_dir)) {
      dirname(script_dir)
    } else if (file.exists("Data/data_in/ansis_data")) {
      getwd()
    } else if (file.exists("../Data/data_in/ansis_data")) {
      normalizePath("..")
    } else {
      getwd()
    }
  }, error = function(e) {
    getwd()
  })
  
  setwd(project_root)
  message("Working directory: ", getwd())
  
  # ============================================================================
  # CONFIGURATION - Change these paths for each property
  # ============================================================================
  
  # Input: Directory containing JSON files for the property
  json_dir <- "Data/data_in/ansis_data/ansis_oc_sites"
  
  # Output: Directory where CSV files will be saved
  output_dir <- "Data/data_out/ansis_lab_measurements/oc"
  
  # ============================================================================
  
  # Find JSON files
  json_files <- list.files(json_dir, pattern = "\\.json$", full.names = TRUE)
  
  if (length(json_files) == 0) {
    stop("No JSON files found in: ", json_dir)
  }
  
  message("\n========================================")
  message("ANSIS JSON to CSV Conversion")
  message("========================================")
  message("Input directory: ", json_dir)
  message("Output directory: ", output_dir)
  message("JSON files found: ", length(json_files))
  
  # Process all JSON files
  all_data <- data.frame()
  
  for (input_json in json_files) {
    message("\nProcessing: ", basename(input_json))
    
    # Convert all properties from JSON
    json_data <- convertANSIStoTERNFormat(json_file = input_json)
    
    # Include all LaboratoryMeasurement records
    lab_data <- json_data[json_data$PropertyType == "LaboratoryMeasurement", ]
    
    all_data <- rbind(all_data, lab_data)
  }
  
  # Create output directory if not exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get unique lab methods
  lab_methods <- unique(all_data$ObservedProperty)
  
  message("\n=== Generating CSV files by method ===")
  
  # Save separate CSV for each lab method
  for (method in lab_methods) {
    method_data <- all_data[all_data$ObservedProperty == method, ]
    output_file <- file.path(output_dir, paste0("ANSIS-", method, ".csv"))
    write.csv(method_data, output_file, row.names = TRUE)
    message("Saved: ", basename(output_file), " (", nrow(method_data), " records)")
  }
  
  # Also save combined file
  output_csv_all <- file.path(output_dir, "ANSIS_combined.csv")
  write.csv(all_data, output_csv_all, row.names = TRUE)
  
  # Display summary
  message("\n========================================")
  message("=== SUMMARY ===")
  message("========================================")
  message("JSON files processed: ", length(json_files))
  message("Total records extracted: ", nrow(all_data))
  message("CSV files generated: ", length(lab_methods) + 1)
  message("\nLab methods found:")
  print(table(all_data$ObservedProperty))
}
