library(terra)

# --- Configuration ---
HomeDir      <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
ModelsDir    <- file.path(HomeDir, "Models")
soil_property <- "CEC"
depth_name   <- "X0.5cm"
TilesDir     <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "preds", paste0("tiles_", depth_name))
output_tif   <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "preds", paste0(soil_property, "_pred_", depth_name, ".tif"))

# Clamping: Set min/max bounds to prevent unrealistic extrapolation
# Adjust these values based on the soil property and expected realistic range
clamp_min <- 0       # Minimum expected value (e.g., CEC can't be negative)
clamp_max <- 100     # Maximum expected value (adjust based on training data range)

# Ensure output directory exists
output_dir <- dirname(output_tif)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 1. Find the Tile Files Directly ---
# We look for the tiles again to ensure we have the correct, absolute paths
tile_files <- list.files(TilesDir, pattern = "\\.tif$", full.names = TRUE)

# Safety Check: Did we find them?
if (length(tile_files) == 0) {
  stop("ERROR: No .tif files found in: ", TilesDir, "\nCheck your folder path!")
} else {
  cat("Found", length(tile_files), "tiles. Assembling map...\n")
}

# --- 2. Create a fresh Virtual Map in Memory ---
# This creates a VRT object 'live' so links are guaranteed to be correct
v_live <- vrt(tile_files)

# --- 3. Apply Clamping ---
cat("Applying clamping (min:", clamp_min, ", max:", clamp_max, ")...\n")
v_clamped <- clamp(v_live, lower = clamp_min, upper = clamp_max)

# --- 4. Write to Final TIF ---
cat("Converting to GeoTIFF (This may take a few minutes)...\n")
writeRaster(
  v_clamped, 
  output_tif, 
  overwrite = TRUE,
  datatype = "FLT4S",  
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
)

cat("Success! Map saved to:", output_tif, "\n")