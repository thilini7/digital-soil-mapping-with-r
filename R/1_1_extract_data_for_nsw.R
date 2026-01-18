# ---- Load required libraries ----
library(ranger)
library(sf)
library(raster)
library(mapview)
library(tidyr)
library(dplyr)
library(aqp)
library(ggplot2)
library(lubridate)
library(nswgeo)

# ---- Set working directory ----
setwd("/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r")

# Platform-independent paths
data_in_dir <- file.path("Data", "data_out", "lab_method_extracts")
data_output_filter_prefix <- file.path("Data", "data_out", "splined_data", "carbon", "Carbon_filter_NSW")
data_output_splined_prefix <- file.path("Data", "data_out", "splined_data", "carbon", "Carbon_splined_NSW")

# Ensure output directories exist
dir.create(dirname(data_output_filter_prefix), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(data_output_splined_prefix), recursive = TRUE, showWarnings = FALSE)

# ---- Load helper functions ----
source("R/1_2_extract_data_for_func.R")

# ---- Load and combine data ----
files <- list.files(path = data_in_dir, pattern = "Organic_Carbon.*\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No Organic_Carbon CSV files found in: ", data_in_dir)
}
dfs <- lapply(files, read.csv)

# Required columns for analysis
required_cols <- c("Latitude", "Longitude", "UpperDepth", "LowerDepth", "Value", "LDI")

dfs <- lapply(dfs, function(x) {
  # Create Location_ID if it doesn't exist (use alternative identifier columns)
  if (!"Location_ID" %in% names(x)) {
    if ("Soilprofileid" %in% names(x)) {
      x$Location_ID <- as.character(x$Soilprofileid)
    } else if ("Observation_ID" %in% names(x)) {
      x$Location_ID <- as.character(x$Observation_ID)
    } else {
      x$Location_ID <- as.character(seq_len(nrow(x)))
    }
  }
  
  # Check if all required columns exist
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    warning(paste("Skipping file - missing columns:", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # Standardize column types
  x$Location_ID <- as.character(x$Location_ID)
  if ("SampleDate" %in% names(x)) x$SampleDate <- as.character(x$SampleDate)
  
  num_cols <- c("UpperDepth", "LowerDepth", "Value", "Latitude", "Longitude", "LDI")
  for (col in num_cols) {
    if (col %in% names(x)) {
      x[[col]] <- suppressWarnings(as.numeric(as.character(x[[col]])))
    }
  }
  
  # Keep only needed columns
  keep_cols <- intersect(names(x), c("Location_ID", "SampleDate", required_cols))
  x[, keep_cols, drop = FALSE]
})

# Remove NULL entries (files that were skipped)
dfs <- Filter(Negate(is.null), dfs)

if (length(dfs) == 0) {
  stop("No valid files found with required columns: ", paste(required_cols, collapse = ", "))
}

df <- dplyr::bind_rows(dfs)
if ("X" %in% names(df)) df$X <- NULL

# ---- Basic cleaning ----
cat("🔍 Checking for NAs...\n")
rows_before <- nrow(df)
df <- df[complete.cases(df[, c("Latitude", "Longitude", "UpperDepth", "LowerDepth", "Value", "LDI")]), ]
cat("Removed", rows_before - nrow(df), "rows with NA in key columns.\n")

# ---- Convert depth units (m -> cm) ----
cat("Converting depths from meters to centimeters...\n")
df$UpperDepth <- df$UpperDepth * 100
df$LowerDepth <- df$LowerDepth * 100
cat("After conversion, UpperDepth min/max:", min(df$UpperDepth), max(df$UpperDepth), "\n")
cat("After conversion, LowerDepth min/max:", min(df$LowerDepth), max(df$LowerDepth), "\n")

# Remove impossible/outlier depths and values (after conversion)
df <- df[df$UpperDepth >= 0 & df$LowerDepth > 0 & df$LowerDepth <= 200, ]
df <- df[df$Value >= 0 & df$Value < 100, ]

# ---- Fix reversed/zero-thickness horizons ----
reversed <- df$LowerDepth <= df$UpperDepth
cat("⚠ Dropping", sum(reversed), "rows where LowerDepth <= UpperDepth (reversed or zero-thickness).\n")
df <- df[!reversed, ]

# ---- NSW boundary filter ----
filter_in_nsw <- function(data) {
  nsw_boundary <- nswgeo::nsw
  if (!inherits(nsw_boundary, "sf")) nsw_boundary <- st_as_sf(nsw_boundary)
  data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
  if (st_crs(nsw_boundary) != st_crs(data_sf)) {
    nsw_boundary <- st_transform(nsw_boundary, st_crs(data_sf))
  }
  filtered_data <- st_filter(data_sf, nsw_boundary)
  st_drop_geometry(filtered_data)
}
df <- filter_in_nsw(df)

# ---- SampleDate cleaning ----
if ("SampleDate" %in% names(df)) {
  df <- df[!grepl("NA-NA-NA", df$SampleDate), ]
  df$SampleDate <- trimws(df$SampleDate)
  df$SampleDate <- gsub("\\s+", " ", df$SampleDate)
  df$SampleDate <- gsub(" ([AP]M)$", "\\1", df$SampleDate, ignore.case = TRUE)
  df$SampleDate <- parse_date_time(
    df$SampleDate,
    orders = c("ymd", "mdy HMS p", "mdy HM p", "mdy", "dmy HMS", "dmy HM", "dmy"),
    exact = FALSE
  )
  df$SampleDate <- as.Date(df$SampleDate)
  
  start_date <- as.Date("1970-01-01")
  end_date   <- as.Date("2000-01-01")
  df <- filter(df, SampleDate >= start_date & SampleDate <= end_date)
}

# ---- Save filtered data ----
if (nrow(df) > 0) {
  data_output_filter_file <- paste0(data_output_filter_prefix, "_1991_2020_Data.csv")
  write.csv(df, data_output_filter_file, row.names = FALSE)
  cat("✅ Filtering complete! Saved as:", data_output_filter_file, "\n")
} else {
  stop("⚠ No records found inside NSW or in the date range! Check your filters and input data.\n")
}

# ---- Run equal-area spline ----
gsm.depths <- c(0, 5, 15, 30, 60, 100, 200)
eaFit <- ea_spline(
  df[, c("Location_ID", "UpperDepth", "LowerDepth", "Value")],
  var.name = "Value",
  d = t(gsm.depths),
  lam = 0.01,
  vlow = 0,
  show.progress = FALSE
)

data_splined <- eaFit$harmonised %>%
  mutate(across(where(is.numeric), ~na_if(., -9999))) %>%
  dplyr::select(-`soil depth`)

# ---- Merge coords and LDI ----
df_coords <- df %>%
  group_by(Location_ID) %>%
  summarise(across(c(Longitude, Latitude, LDI), mean, na.rm = TRUE), .groups = "drop")

data_splined$id <- as.character(data_splined$id)
df_coords$Location_ID <- as.character(df_coords$Location_ID)
data_splined <- merge(data_splined, df_coords, by.x = "id", by.y = "Location_ID")

# ---- Save splined data ----
data_output_splined_file <- paste0(data_output_splined_prefix, "_1991_2020_Data.csv")
write.csv(data_splined, data_output_splined_file, row.names = FALSE)
cat("✅ Splined data saved as:", data_output_splined_file, "\n")

# ---- Optional clean NA ----
cols_to_clean <- names(data_splined)[grepl("^X", names(data_splined))]
if (length(cols_to_clean) > 0) {
  cleaned_file <- paste0(data_output_splined_prefix, "_cleaned_1991_2020_Data.csv")
  data_splined <- na_fill(data_splined, cutoffs = 100, colnames = cols_to_clean)
  write.csv(data_splined, cleaned_file, row.names = FALSE)
  cat("✅ Cleaned splined data saved as:", cleaned_file, "\n")
}

# ---- Quick map ----
if (nrow(data_splined) > 0) {
  nsw_boundary <- nswgeo::nsw
  if (!inherits(nsw_boundary, "sf")) nsw_boundary <- st_as_sf(nsw_boundary)
  ggplot() +
    geom_sf(data = nsw_boundary, fill = "lightblue", color = "black", alpha = 0.5) +
    geom_point(data = data_splined, aes(x = Longitude, y = Latitude), color = "red", size = 1) +
    labs(title = "Filtered Data Points in NSW", x = "Longitude", y = "Latitude") +
    theme_minimal()
}
