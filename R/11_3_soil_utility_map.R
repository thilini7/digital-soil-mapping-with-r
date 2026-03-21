# =============================================================================
# soil Utility Map & Graph — Gaussian Bell-Shape Function
# =============================================================================
# Publication-quality utility maps and graphs following the style of
# Hunakunti et al. (2025) — Int. Soil and Water Conservation Research.
#
# Applies a Gaussian utility function to the soil pedogenon mean raster
# (Genosoil) to produce:
#   1. A utility graph (soil vs Utility) with pedogenon class overlay
#   2. A spatial utility map with NSW boundary, scale bar, north arrow
#   3. A utility distribution histogram
#   4. Per-class utility CSV
#
# Gaussian utility:  U(soil) = exp( -0.5 * ((soil - μ) / σ)^2 )
#   μ = 6.5  (optimum soil for plant nutrient availability)
#   σ = 1    (spread)
#
# Reference:
#   - documents/Estimating-surrogates-utility.pdf
#   - documents/Application-of-fuzzy(soil shape).pdf
#   - documents/1-s2.0-S2095633925000711-main.pdf
# =============================================================================

# --- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(tidyterra)
  library(sf)
  library(ggspatial)
  library(cowplot)
})

set.seed(42)

# --- Terra options -----------------------------------------------------------
temp_dir <- file.path(getwd(), "temp_terra")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
terraOptions(memfrac = 0.8, tempdir = temp_dir)

# =============================================================================
# CONFIGURATION
# =============================================================================
depth_layer <- "X0.30cm"

# Gaussian parameters
mu    <- 6.5
sigma <- 1

# Input pH pedogenon mean rasters (output from 11_1 script)
input_base <- file.path("Data", "data_out", "pedogenon_zonal_means")

ph_inputs <- list(
  Geno  = file.path(input_base, "Geno",
                    paste0("pH_geno_pedogenon_mean_", depth_layer, ".tif"))
)

# NSW state boundary for map overlay
# Use rnaturalearth or ozmaps if available; fall back to raster extent
nsw_boundary <- tryCatch({
  if (requireNamespace("ozmaps", quietly = TRUE)) {
    oz <- ozmaps::ozmap_states
    nsw <- oz[oz$NAME == "New South Wales", ]
    st_transform(nsw, 4326)
  } else {
    NULL
  }
}, error = function(e) NULL)

# Output directory
output_base <- file.path("Data", "data_out", "pH_utility")
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

# Publication theme
theme_pub <- theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle    = element_text(size = 10, hjust = 0.5, colour = "grey30"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9),
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

# =============================================================================
# Gaussian utility function
# =============================================================================
gaussian_utility <- function(x, mu, sigma) {
  exp(-0.5 * ((x - mu) / sigma)^2)
}

# =============================================================================
# 1. Apply utility function to raster
# =============================================================================
for (soil_type in names(ph_inputs)) {

  cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat("Processing: ", soil_type, " pH utility\n", sep = "")
  cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

  in_path <- ph_inputs[[soil_type]]
  if (!file.exists(in_path)) {
    cat("  ** SKIPPED — file not found: ", in_path, "\n\n", sep = "")
    next
  }

  cat("  Input raster : ", in_path, "\n", sep = "")
  ph_rast <- rast(in_path)

  # --- Apply Gaussian utility ------------------------------------------------
  cat("  Applying Gaussian utility (mu=", mu, ", sigma=", sigma, ")...\n", sep = "")
  utility_rast <- app(ph_rast, fun = function(x) {
    gaussian_utility(x, mu = mu, sigma = sigma)
  })
  names(utility_rast) <- "Utility"

  # --- Summary statistics ----------------------------------------------------
  util_vals <- values(utility_rast, mat = FALSE)
  util_vals <- util_vals[!is.na(util_vals)]
  cat("  Utility range  : [", round(min(util_vals), 4), ", ",
      round(max(util_vals), 4), "]\n", sep = "")
  cat("  Utility mean   : ", round(mean(util_vals), 4), "\n", sep = "")
  cat("  Utility median : ", round(median(util_vals), 4), "\n", sep = "")

  # --- Save utility raster (GeoTIFF) ----------------------------------------
  out_tif <- file.path(output_base,
                       paste0("pH_utility_", tolower(soil_type), "_",
                              depth_layer, ".tif"))
  writeRaster(utility_rast, out_tif, overwrite = TRUE,
              gdal = c("COMPRESS=LZW"))
  cat("  Raster saved   : ", out_tif, "\n", sep = "")

  # --- Per-class utility CSV -------------------------------------------------
  zonal_csv <- file.path(input_base, soil_type,
                         paste0("pH_", tolower(soil_type),
                                "_pedogenon_zonal_summary_", depth_layer, ".csv"))
  zonal_df <- NULL
  if (file.exists(zonal_csv)) {
    zonal_df <- read.csv(zonal_csv)
    zonal_df$utility <- gaussian_utility(zonal_df$mean_value, mu, sigma)
    zonal_df <- zonal_df[order(zonal_df$pedogenon_class), ]

    out_util_csv <- file.path(output_base,
                              paste0("pH_utility_per_class_", tolower(soil_type),
                                     "_", depth_layer, ".csv"))
    write.csv(zonal_df, out_util_csv, row.names = FALSE)
    cat("  Class CSV      : ", out_util_csv, "\n", sep = "")
    cat("    (", nrow(zonal_df), " classes, utility range: [",
        round(min(zonal_df$utility), 4), ", ",
        round(max(zonal_df$utility), 4), "])\n", sep = "")
  } else {
    cat("  ** Zonal summary CSV not found, skipping class-level CSV\n")
  }

  # =========================================================================
  # 2. Utility Graph (publication style — like Fig. 12 in Hunakunti et al.)
  # =========================================================================
  cat("\n  Generating utility graph...\n")

  curve_df <- data.frame(
    pH      = seq(2, 12, by = 0.01),
    Utility = gaussian_utility(seq(2, 12, by = 0.01), mu, sigma)
  )

  p_curve <- ggplot() +
    # Shaded area under curve
    geom_area(data = curve_df, aes(x = pH, y = Utility),
              fill = "#2166AC", alpha = 0.08) +
    # Main curve
    geom_line(data = curve_df, aes(x = pH, y = Utility),
              colour = "#2166AC", linewidth = 1.2) +
    # Overlay pedogenon class points if available
    {if (!is.null(zonal_df))
      geom_point(data = zonal_df,
                 aes(x = mean_value, y = utility),
                 colour = "#B2182B", size = 1.2, alpha = 0.5)
    } +
    # Optimum line
    geom_vline(xintercept = mu, linetype = "dashed", colour = "grey40",
               linewidth = 0.5) +
    annotate("text", x = mu + 0.4, y = 0.97,
             label = paste0("mu == ", mu), parse = TRUE,
             size = 3.5, colour = "grey30") +
    # ±1σ markers
    geom_segment(aes(x = mu - sigma, xend = mu - sigma, y = 0, yend = exp(-0.5)),
                 linetype = "dotted", colour = "grey60") +
    geom_segment(aes(x = mu + sigma, xend = mu + sigma, y = 0, yend = exp(-0.5)),
                 linetype = "dotted", colour = "grey60") +
    annotate("text", x = mu + sigma + 0.5, y = exp(-0.5) + 0.03,
             label = "1*sigma", parse = TRUE, size = 3, colour = "grey50") +
    scale_x_continuous(breaks = seq(2, 12, by = 1), limits = c(2, 12),
                       expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.05),
                       expand = c(0, 0)) +
    labs(
      x        = "pH",
      y        = "Utility",
      title    = paste0("pH Capacity Utility — ", soil_type, "soil"),
      subtitle = bquote("Gaussian:" ~ mu == .(mu) * "," ~ sigma == .(sigma))
    ) +
    theme_pub

  # Add legend annotation for class points
  if (!is.null(zonal_df)) {
    p_curve <- p_curve +
      annotate("point", x = 10.5, y = 0.92,
               colour = "#B2182B", size = 2) +
      annotate("text", x = 10.8, y = 0.92,
               label = paste0("Pedogenon classes (n=", nrow(zonal_df), ")"),
               size = 3, hjust = 0, colour = "grey30")
  }

  curve_png <- file.path(output_base,
                         paste0("pH_utility_graph_", tolower(soil_type), "_",
                                depth_layer, ".png"))
  ggsave(curve_png, p_curve, width = 8, height = 5, dpi = 300)
  cat("  Graph saved    : ", curve_png, "\n", sep = "")

  # =========================================================================
  # 3. Utility Map (publication style — like Fig. 11 in Hunakunti et al.)
  # =========================================================================
  cat("  Generating utility map...\n")

  # Downsample raster for faster ggplot rendering
  target_ncell <- 2e6
  agg_factor <- max(1, round(sqrt(ncell(utility_rast) / target_ncell)))
  if (agg_factor > 1) {
    cat("  Aggregating raster by factor ", agg_factor, " for plotting...\n", sep = "")
    plot_rast <- aggregate(utility_rast, fact = agg_factor, fun = "mean",
                           na.rm = TRUE)
  } else {
    plot_rast <- utility_rast
  }

  p_map <- ggplot() +
    geom_spatraster(data = plot_rast) +
    scale_fill_whitebox_c(
      palette = "muted",
      limits  = c(0, 1),
      breaks  = seq(0, 1, by = 0.2),
      labels  = sprintf("%.1f", seq(0, 1, by = 0.2)),
      na.value = "transparent",
      name    = "Utility"
    ) +
    # NSW boundary overlay
    {if (!is.null(nsw_boundary))
      geom_sf(data = nsw_boundary, fill = NA, colour = "black",
              linewidth = 0.4)
    } +
    # Scale bar and north arrow
    annotation_scale(
      location = "bl", width_hint = 0.25,
      pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
      style = "ticks", text_cex = 0.7
    ) +
    annotation_north_arrow(
      location = "tr", which_north = "true",
      height = unit(1.2, "cm"), width = unit(1.0, "cm"),
      pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
      style = north_arrow_fancy_orienteering(text_size = 7)
    ) +
    labs(
      title    = paste0("pH Capacity Utility Map — ", soil_type, "soil"),
      subtitle = bquote("Gaussian:" ~ mu == .(mu) * "," ~ sigma == .(sigma) ~
                         " | Depth: " * .(gsub("X", "", depth_layer))),
      x = "Longitude", y = "Latitude"
    ) +
    coord_sf(expand = FALSE) +
    theme_pub +
    theme(
      legend.position  = "right",
      legend.key.height = unit(1.8, "cm"),
      legend.key.width  = unit(0.4, "cm"),
      axis.text         = element_text(size = 8)
    )

  map_png <- file.path(output_base,
                       paste0("pH_utility_map_", tolower(soil_type), "_",
                              depth_layer, ".png"))
  ggsave(map_png, p_map, width = 10, height = 8, dpi = 300)
  cat("  Map saved      : ", map_png, "\n", sep = "")

  # =========================================================================
  # 4. Utility Distribution Histogram
  # =========================================================================
  cat("  Generating histogram...\n")

  hist_df <- data.frame(utility = util_vals)
  med_val <- median(util_vals)

  p_hist <- ggplot(hist_df, aes(x = utility)) +
    geom_histogram(bins = 60, fill = "#2166AC", colour = "white",
                   linewidth = 0.2, alpha = 0.85) +
    geom_vline(xintercept = med_val, linetype = "dashed",
               colour = "#B2182B", linewidth = 0.8) +
    annotate("text", x = med_val + 0.04, y = Inf, vjust = 2,
             label = paste0("Median = ", round(med_val, 3)),
             colour = "#B2182B", size = 3.5, hjust = 0) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    labs(
      x        = "Utility (0–1)",
      y        = "Pixel count",
      title    = paste0("pH Utility Distribution — ", soil_type, "soil"),
      subtitle = paste0("Depth: ", gsub("X", "", depth_layer))
    ) +
    theme_pub

  hist_png <- file.path(output_base,
                        paste0("pH_utility_hist_", tolower(soil_type), "_",
                               depth_layer, ".png"))
  ggsave(hist_png, p_hist, width = 8, height = 5, dpi = 300)
  cat("  Histogram saved: ", hist_png, "\n", sep = "")

  # =========================================================================
  # 5. Combined panel (graph + map side by side)
  # =========================================================================
  cat("  Generating combined panel...\n")

  p_combined <- plot_grid(
    p_curve, p_map,
    ncol = 2, rel_widths = c(0.45, 0.55),
    labels = c("(a)", "(b)"), label_size = 12
  )

  combined_png <- file.path(output_base,
                            paste0("pH_utility_combined_", tolower(soil_type),
                                   "_", depth_layer, ".png"))
  ggsave(combined_png, p_combined, width = 16, height = 7, dpi = 300)
  cat("  Combined saved : ", combined_png, "\n\n", sep = "")
}

# =============================================================================
# Summary
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("DONE — pH Utility outputs in: ", output_base, "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
