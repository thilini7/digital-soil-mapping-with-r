library(terra)

# --- Configuration ---
HomeDir      <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
ModelsDir    <- file.path(HomeDir, "Models")
soil_property <- "CEC"
depth_name   <- "X0.5cm"
TilesDir     <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "preds", paste0("tiles_", depth_name))
output_tif   <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "preds", paste0(soil_property, "_pred_", depth_name, ".tif"))

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

# --- 3. Write to Final TIF ---
cat("Converting to GeoTIFF (This may take a few minutes)...\n")
writeRaster(
  v_live, 
  output_tif, 
  overwrite = TRUE,
  datatype = "FLT4S",  
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
)

cat("Success! Map saved to:", output_tif, "\n")