library(terra)

output_tif <- "/Users/neo/Development/Data-Science/data-science/maps/OC_pred_X0.5cm.tif"

# Load raster
soc_map <- rast(output_tif)

# Stronger brown/yellow/red at low SOC, then green → blue
soc_cols <- colorRampPalette(c(
  "#8c2d04", # dark brown
  "#cc4c02", # reddish brown
  "#ec7014", # orange-brown
  "#fe9929", # orange
  "#fed98e", # light orange/yellow
  "#ffffbf", # pale yellow
  "#d9f0a3", # pale green
  "#91cf60", # green
  "#41ab5d", # dark green
  "#1a9850", # deep green
  "#0571b0", # blue
  "#08306b"  # dark blue
))(200)

# Plot
plot(soc_map, col = soc_cols, main = "Predicted SOC 0–5 cm (Warm Low, Cool High)")
