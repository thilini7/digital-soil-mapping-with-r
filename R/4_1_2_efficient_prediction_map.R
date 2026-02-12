# ==============================================================================
# SAFE MODE: TILE-BASED PREDICTION
# Best for low RAM or massive areas. Saves progress automatically.
# ==============================================================================

library(terra)
library(Cubist)
library(caret)

# --- 1. CONFIGURATION ---
HomeDir       <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
ModelsDir     <- file.path(HomeDir, "Models")
CovariatesDir <- file.path(HomeDir, "Data/data_in/soil_covariates_aligned_v2")
soil_property <- "CEC"
depth_name    <- "X0.5cm"
model_file    <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "models", paste0(depth_name, ".rds"))

# Back-transformation setting: Set to TRUE if model was trained on log1p-transformed data
# (check 3_4_cubist_regression_model_training.R -> use_log_transform setting)
use_log_transform <- TRUE  # Set to FALSE for properties like pH where log may not be appropriate

# Output folder for tiles
TilesDir      <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "preds", paste0("tiles_", depth_name))
if(!dir.exists(TilesDir)) dir.create(TilesDir, recursive=TRUE)

# Final merged output
final_output  <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "preds", paste0(soil_property, "_pred_", depth_name, ".tif"))

# Number of splits (e.g., 4 means 4x4 = 16 total tiles)
# Increase this number if you still get crashes.
SPLIT_FACTOR <- 10 

# --- 2. LOAD MODEL & DATA ---
cat("Loading Model and Covariates...\n")
model <- readRDS(model_file)

# Extract names based on model type
if (inherits(model, "train")) {
  model_cov_names <- model$coefnames
} else {
  model_cov_names <- model$vars$all
  if (is.null(model_cov_names)) model_cov_names <- names(model$coefficients)
}
# Fallback
if (is.null(model_cov_names)) {
  load(file.path(HomeDir, "Data/data_out/RData", paste0(soil_property, "_covs_regression.RData")))
  model_cov_names <- cov_names
}

# Load Raster Stack (Lazy load - uses no RAM yet)
cov_files <- list.files(CovariatesDir, pattern = "\\.tif$", full.names = TRUE)
base_names <- tools::file_path_sans_ext(basename(cov_files))
cov_stack <- rast(cov_files[match(model_cov_names, base_names)])
names(cov_stack) <- model_cov_names

# --- 3. CREATE TILES ---
# We define a grid of extents covering the whole area
total_ext <- ext(cov_stack)
x_edges <- seq(total_ext$xmin, total_ext$xmax, length.out = SPLIT_FACTOR + 1)
y_edges <- seq(total_ext$ymin, total_ext$ymax, length.out = SPLIT_FACTOR + 1)

cat("Map split into", SPLIT_FACTOR * SPLIT_FACTOR, "tiles.\n")

# --- 4. TILE LOOP ---
# This loop processes one tile at a time
count <- 1
total_tiles <- SPLIT_FACTOR * SPLIT_FACTOR

for (i in 1:SPLIT_FACTOR) {
  for (j in 1:SPLIT_FACTOR) {
    
    # A. Define Tile Filename
    tile_filename <- file.path(TilesDir, paste0("tile_", i, "_", j, ".tif"))
    
    # B. Skip if already exists (Resume feature)
    if (file.exists(tile_filename)) {
      cat(sprintf("[%d/%d] Tile %d-%d exists. Skipping.\n", count, total_tiles, i, j))
      count <- count + 1
      next
    }
    
    cat(sprintf("[%d/%d] Processing Tile %d-%d... ", count, total_tiles, i, j))
    
    # C. Crop Stack to Tile Extent
    # We add a small buffer? No, exact crop is fine for pixel-based prediction
    e <- ext(x_edges[i], x_edges[i+1], y_edges[j], y_edges[j+1])
    tile_stack <- crop(cov_stack, e)
    
    # D. Check if tile has data (it might be empty ocean/border)
    # Using global is fast check
    if (all(is.na(values(tile_stack[[1]], mat=FALSE)[1:100]))) { 
      # Deeper check if edge is NA
      # If you want to be safe, just predict. terra handles NA fast.
    }
    
    tryCatch({
      # E. Predict Tile
      # Notes: CORES=1 to prevent RAM explosion.
      pred_tile <- terra::predict(
        object = tile_stack,
        model = model,
        fun = function(model, data) {
          # Force neighbors=0 for speed/safety. 
          # Remove neighbors=0 if you absolutely need high accuracy and have time.
          predict(model, newdata = as.data.frame(data), neighbors = 0) 
        },
        na.rm = TRUE
      )
      
      # F. Back-transform if model was trained on log1p scale
      if (use_log_transform) {
        pred_tile <- expm1(pred_tile)  # expm1(x) = exp(x) - 1, reverses log1p()
      }
      
      # G. Write tile to file
      writeRaster(pred_tile, tile_filename, overwrite = TRUE,
                  wopt = list(gdal=c("COMPRESS=DEFLATE")))
      cat("Done.\n")
      
    }, error = function(e) {
      cat("Error on tile (might be empty/ocean):", e$message, "\n")
    })
    
    # H. CLEANUP MEMORY (Crucial)
    rm(tile_stack)
    if (exists("pred_tile")) rm(pred_tile)
    gc() # Force R to release RAM
    
    count <- count + 1
  }
}

# --- 5. MERGE TILES ---
cat("All tiles processed. Merging into final VRT...\n")

# Get list of created tiles
tile_files <- list.files(TilesDir, pattern = "\\.tif$", full.names = TRUE)

if (length(tile_files) > 0) {
  # Create a Virtual Raster (VRT) - instantaneous, takes no space
  vrt_file <- file.path(ModelsDir, paste0("mod.cubist.", soil_property), "preds", paste0(soil_property, "_pred_", depth_name, ".vrt"))
  v <- vrt(tile_files, vrt_file, overwrite=TRUE)
  
  # Optional: Save as single TIF (might take time/space)
  # writeRaster(v, final_output, overwrite=TRUE)
  
  cat("Success! Map available at:", vrt_file, "\n")
  plot(v, main="Final Tiled Prediction")
} else {
  stop("No tiles were generated.")
}