# =============================================================================
# Soil Utility Map & Graph — Capacity and Condition
# =============================================================================
# Publication-quality utility maps and graphs following the style of
# Hunakunti et al. (2025) — Int. Soil and Water Conservation Research.
#
# Supports two modes:
#   CAPACITY  — Francos et al. (2025) Eq. 1 for pH:
#               U(pH) = exp( -(pH - 6.5)^2 / 2 )
#   CONDITION — two methods:
#     pH: Francos et al. (2025) Eqs. 1–3
#         U_geno/pheno = exp(-(pH - 6.5)^2 / 2)       [Gaussian capacity]
#         U_condition  = exp(-k*(U_geno - U_pheno))    [exponential decay]
#         capped at 1 when phenosoil ≥ genosoil utility
#
# Reference:
#   - Ng et al. (2024) Soil Research — pH Gaussian μ=6.5, σ=1
#   - McBratney & Odeh (1997) Geoderma — fuzzy membership functions
#   - Evangelista et al. (2025) Soil Security — relative difference
#   - Francos et al. (2025) Catena 259:109372 — pH capacity & condition Eqs. 1–3
#   - Hunakunti et al. (2025) ISWCR — capacity-condition framework
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
soil_property <- "AWC"              # <-- "pH", "Bulk_Density", or "AWC"
depth_layer   <- "X0.30cm"
utility_mode  <- "condition"       # <-- "capacity" or "condition"

# =============================================================================
# Utility function definitions per soil property
# =============================================================================
# capacity_cfg:  Francos et al. (2025) Eq. 1 for pH
#   U(pH) = exp( -(pH - 6.5)^2 / 2 )
#
# condition_cfg: Francos et al. (2025) for pH (method = "utility_diff"),
#   or Evangelista et al. (2025) linear ratio for other properties
# =============================================================================
utility_config <- list(

  pH = list(
    capacity_cfg = list(
      fn     = function(x, params) exp(-(x - params$mu)^2 / 2),
      params = list(mu = 6.5),
      x_lab  = "pH_capacity",
      x_range = c(2, 12),
      subtitle_expr = function(p) bquote("Gaussian:" ~ mu == .(p$mu))
    ),
    condition_cfg = list(
      method = "utility_diff",
      fn     = function(x, params) ifelse(x < 0, 1, exp(-params$k * x)),
      params = list(k = 10, mu = 6.5),
      x_lab  = expression(Delta*"U  ("*U[geno] - U[pheno]*")"),
      x_range = c(-0.5, 1),
      subtitle_expr = function(p) bquote("Exponential Decay:  k =" ~ .(p$k))
    )
  ),

  Bulk_Density = list(
    capacity_cfg = list(
      fn     = function(x, params) 1 / (1 + exp(params$a * (x - params$b))),
      params = list(a = 10, b = 1.5),
      x_lab  = expression("Bulk Density (g/cm"^3*")"),
      x_range = c(0.5, 2.2),
      subtitle_expr = function(p) bquote("Less is better: b =" ~ .(p$b))
    ),
    condition_cfg = list(
      method = "utility_diff",
      fn     = function(x, params) ifelse(x < 0, 1, exp(-params$k * x)),
      params = list(k = 10),
      x_lab  = expression(Delta*"U  ("*U[geno] - U[pheno]*")"),
      x_range = c(-0.5, 1),
      subtitle_expr = function(p) bquote("Exponential Decay:  k =" ~ .(p$k))
    )
  ),

  AWC = list(
    capacity_cfg = list(
      fn     = function(x, params) 1 / (1 + exp(-params$c * (x - params$d))),
      params = list(c = 1.0, d = 12),
      x_lab  = "AWC (%)",
      x_range = c(0, 40),
      subtitle_expr = function(p) bquote("More is better: d =" ~ .(p$d))
    ),
    condition_cfg = list(
      method = "utility_diff",
      fn     = function(x, params) ifelse(x < 0, 1, exp(-params$k * x)),
      params = list(k = 10),
      x_lab  = expression(Delta*"U  ("*U[geno] - U[pheno]*")"),
      x_range = c(-0.5, 1),
      subtitle_expr = function(p) bquote("Exponential Decay:  k =" ~ .(p$k))
    )
  )
)

# --- Resolve config ----------------------------------------------------------
if (!soil_property %in% names(utility_config)) {
  stop("No utility config for '", soil_property,
       "'. Available: ", paste(names(utility_config), collapse = ", "))
}
if (!utility_mode %in% c("capacity", "condition")) {
  stop("utility_mode must be 'capacity' or 'condition'")
}

cfg_key     <- paste0(utility_mode, "_cfg")
prop_cfg    <- utility_config[[soil_property]][[cfg_key]]
utility_fn  <- prop_cfg$fn
util_params <- prop_cfg$params
is_utility_diff <- isTRUE(prop_cfg$method == "utility_diff")

# When utility_diff condition, retrieve the capacity function for U_geno / U_pheno
if (is_utility_diff) {
  cap_cfg     <- utility_config[[soil_property]][["capacity_cfg"]]
  capacity_fn <- cap_cfg$fn
  cap_params  <- cap_cfg$params
}

# --- Resolve input raster path -----------------------------------------------
input_base <- file.path("Data", "data_out", "pedogenon_zonal_means")

if (utility_mode == "capacity") {
  soil_type    <- "Geno"
  input_raster <- file.path(input_base, "Geno",
                            paste0(soil_property, "_geno_pedogenon_mean_",
                                   depth_layer, ".tif"))
  zonal_csv    <- file.path(input_base, "Geno",
                            paste0(soil_property, "_geno_pedogenon_zonal_summary_",
                                   depth_layer, ".csv"))
  mode_label   <- "Capacity"
} else {
  soil_type    <- "Condition"
  mode_label   <- "Condition"
  zonal_csv    <- file.path("Data", "data_out", "condition_maps", soil_property,
                            paste0(soil_property, "_condition_per_class_",
                                   depth_layer, ".csv"))
  if (is_utility_diff) {
    # Francos et al. (2025): need both genosoil and phenosoil rasters
    geno_raster  <- file.path(input_base, "Geno",
                              paste0(soil_property, "_geno_pedogenon_mean_",
                                     depth_layer, ".tif"))
    pheno_raster <- file.path(input_base, "Pheno",
                              paste0(soil_property, "_pheno_pedogenon_mean_",
                                     depth_layer, ".tif"))
    input_raster <- NULL
  } else {
    input_raster <- file.path("Data", "data_out", "condition_maps", soil_property,
                              paste0(soil_property, "_condition_ratio_",
                                     depth_layer, ".tif"))
  }
}

# NSW state boundary
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
output_base <- file.path("Data", "data_out",
                         paste0(soil_property, "_utility"))
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
# Utility function wrapper
# =============================================================================
compute_utility <- function(x) utility_fn(x, util_params)

# =============================================================================
# 1. Load and process
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("Processing: ", soil_property, " ", mode_label, " Utility\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

cat("  Mode         : ", mode_label, "\n", sep = "")

if (is_utility_diff) {
  # --- Utility-difference condition approach ----------------------------------
  if (!file.exists(geno_raster))  stop("Genosoil raster not found: ",  geno_raster)
  if (!file.exists(pheno_raster)) stop("Phenosoil raster not found: ", pheno_raster)

  cat("  Method       : Utility-difference condition\n", sep = "")
  cat("  Geno raster  : ", geno_raster,  "\n", sep = "")
  cat("  Pheno raster : ", pheno_raster, "\n", sep = "")
  cat("  Condition k  : ", util_params$k, "\n", sep = "")

  geno_rast  <- rast(geno_raster)
  pheno_rast <- rast(pheno_raster)
  k <- util_params$k

  cat("  Computing U_geno  (capacity function)...\n")
  u_geno  <- app(geno_rast,  fun = function(x) capacity_fn(x, cap_params))
  cat("  Computing U_pheno (capacity function)...\n")
  u_pheno <- app(pheno_rast, fun = function(x) capacity_fn(x, cap_params))

  cat("  Computing U_condition (exponential decay, k=", k, ")...\n", sep = "")
  u_diff <- u_geno - u_pheno
  utility_rast <- app(u_diff, fun = function(d) ifelse(d < 0, 1, exp(-k * d)))
  names(utility_rast) <- "Utility"

} else {
  # --- Standard single-raster approach ---------------------------------------
  if (!file.exists(input_raster)) stop("File not found: ", input_raster)
  cat("  Input raster : ", input_raster, "\n", sep = "")
  if (length(util_params) > 0) {
    cat("  Parameters   : ", paste(names(util_params), util_params, sep = "=", collapse = ", "), "\n", sep = "")
  } else {
    cat("  Parameters   : Linear (ratio = score, Evangelista et al. 2025)\n", sep = "")
  }

  prop_rast <- rast(input_raster)

  # --- Apply utility function ------------------------------------------------
  cat("  Applying utility function...\n")
  utility_rast <- app(prop_rast, fun = function(x) compute_utility(x))
  names(utility_rast) <- "Utility"
}

# --- Summary statistics ------------------------------------------------------
util_vals <- values(utility_rast, mat = FALSE)
util_vals <- util_vals[!is.na(util_vals)]
cat("  Utility range  : [", round(min(util_vals), 4), ", ",
    round(max(util_vals), 4), "]\n", sep = "")
cat("  Utility mean   : ", round(mean(util_vals), 4), "\n", sep = "")
cat("  Utility median : ", round(median(util_vals), 4), "\n", sep = "")

# --- Save utility raster (GeoTIFF) ------------------------------------------
out_tif <- file.path(output_base,
                     paste0(soil_property, "_", tolower(mode_label),
                            "_utility_", depth_layer, ".tif"))
writeRaster(utility_rast, out_tif, overwrite = TRUE,
            gdal = c("COMPRESS=LZW"))
cat("  Raster saved   : ", out_tif, "\n", sep = "")

# --- Per-class utility CSV ---------------------------------------------------
zonal_df <- NULL
if (file.exists(zonal_csv)) {
  zonal_df <- read.csv(zonal_csv)

  if (is_utility_diff) {
    # Compute per-class condition utilities via capacity function
    zonal_df$u_geno  <- capacity_fn(zonal_df$geno_mean, cap_params)
    zonal_df$u_pheno <- capacity_fn(zonal_df$pheno_mean, cap_params)
    zonal_df$u_diff  <- zonal_df$u_geno - zonal_df$u_pheno
    zonal_df$utility <- ifelse(zonal_df$u_diff < 0, 1,
                               exp(-util_params$k * zonal_df$u_diff))
    point_x_col <- "u_diff"
  } else if (utility_mode == "capacity") {
    zonal_df$utility <- compute_utility(zonal_df$mean_value)
    point_x_col <- "mean_value"
  } else {
    zonal_df$utility <- compute_utility(zonal_df$condition_ratio)
    point_x_col <- "condition_ratio"
  }
  zonal_df <- zonal_df[order(zonal_df$pedogenon_class), ]

  out_util_csv <- file.path(output_base,
                            paste0(soil_property, "_", tolower(mode_label),
                                   "_utility_per_class_", depth_layer, ".csv"))
  write.csv(zonal_df, out_util_csv, row.names = FALSE)
  cat("  Class CSV      : ", out_util_csv, "\n", sep = "")
  cat("    (", nrow(zonal_df), " classes, utility range: [",
      round(min(zonal_df$utility, na.rm = TRUE), 4), ", ",
      round(max(zonal_df$utility, na.rm = TRUE), 4), "])\n", sep = "")
} else {
  cat("  ** CSV not found: ", zonal_csv, ", skipping class-level CSV\n", sep = "")
  point_x_col <- NULL
}

# =============================================================================
# 2. Utility Graph
# =============================================================================
cat("\n  Generating utility graph...\n")

x_seq <- seq(prop_cfg$x_range[1], prop_cfg$x_range[2], length.out = 1000)
curve_df <- data.frame(value = x_seq, Utility = compute_utility(x_seq))

p_curve <- ggplot() +
  geom_line(data = curve_df, aes(x = value, y = Utility),
            colour = "grey60", linewidth = 0.6, linetype = "dashed") +
  {if (utility_mode == "condition")
    geom_vline(xintercept = if (is_utility_diff) 0 else 1,
               linetype = "dotted", colour = "forestgreen",
               linewidth = 0.7)
  } +
  {if (!is.null(zonal_df) && !is.null(point_x_col))
    geom_point(data = zonal_df,
               aes(x = .data[[point_x_col]], y = utility),
               colour = "#B2182B", size = 1.2, alpha = 0.5)
  } +
  scale_x_continuous(limits = prop_cfg$x_range, expand = c(0.01, 0)) +
  {if (utility_mode == "condition" && !is_utility_diff)
    scale_y_continuous(expand = c(0, 0))
  else
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.05),
                       expand = c(0, 0))
  } +
  labs(
    x        = prop_cfg$x_lab,
    y        = "Utility",
    title    = paste0(soil_property, " ", mode_label, " Utility"),
    subtitle = prop_cfg$subtitle_expr(util_params)
  ) +
  theme_pub

curve_png <- file.path(output_base,
                       paste0(soil_property, "_", tolower(mode_label),
                              "_utility_graph_", depth_layer, ".png"))
ggsave(curve_png, p_curve, width = 8, height = 5, dpi = 300)
cat("  Graph saved    : ", curve_png, "\n", sep = "")

# =============================================================================
# 3. Utility Map
# =============================================================================
cat("  Generating utility map...\n")

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
  {if (utility_mode == "condition" && !is_utility_diff)
    scale_fill_gradient2(
      low      = "#B2182B",
      mid      = "white",
      high     = "#2166AC",
      midpoint = 1,
      na.value = "transparent",
      name     = "Utility"
    )
  else
    scale_fill_whitebox_c(
      palette = "muted",
      limits  = c(0, 1),
      breaks  = seq(0, 1, by = 0.2),
      labels  = sprintf("%.1f", seq(0, 1, by = 0.2)),
      na.value = "transparent",
      name    = "Utility"
    )
  } +
  {if (!is.null(nsw_boundary))
    geom_sf(data = nsw_boundary, fill = NA, colour = "black",
            linewidth = 0.4)
  } +
  annotation_scale(
    location = "bl", width_hint = 0.25,
    pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
    style = "ticks", text_cex = 0.7
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    height = unit(1.2, "cm"), width = unit(1.0, "cm"),
    pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
    style = north_arrow_fancy_orienteering(text_size = 7)
  ) +
  labs(
    title    = paste0(soil_property, " ", mode_label, " Utility Map"),
    subtitle = paste0("Depth: ", gsub("X", "", depth_layer)),
    x = "Longitude", y = "Latitude"
  ) +
  coord_sf(expand = FALSE) +
  theme_pub +
  theme(
    legend.position   = "right",
    legend.key.height = unit(1.8, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    axis.text         = element_text(size = 8)
  )

map_png <- file.path(output_base,
                     paste0(soil_property, "_", tolower(mode_label),
                            "_utility_map_", depth_layer, ".png"))
ggsave(map_png, p_map, width = 10, height = 8, dpi = 300)
cat("  Map saved      : ", map_png, "\n", sep = "")

# =============================================================================
# 4. Utility Distribution Histogram
# =============================================================================
cat("  Generating histogram...\n")

hist_df <- data.frame(utility = util_vals)
med_val <- median(util_vals)

p_hist <- ggplot(hist_df, aes(x = utility)) +
  geom_histogram(bins = 60, fill = "#2166AC", colour = "white",
                 linewidth = 0.2, alpha = 0.85) +
  geom_vline(xintercept = med_val, linetype = "dashed",
             colour = "#B2182B", linewidth = 0.8) +
  {if (utility_mode == "condition" && !is_utility_diff)
    geom_vline(xintercept = 1, linetype = "dotted",
               colour = "forestgreen", linewidth = 0.7)
  } +
  annotate("text", x = med_val + 0.04, y = Inf, vjust = 2,
           label = paste0("Median = ", round(med_val, 3)),
           colour = "#B2182B", size = 3.5, hjust = 0) +
  {if (utility_mode == "condition" && !is_utility_diff)
    scale_x_continuous()
  else
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))
  } +
  labs(
    x     = if (is_utility_diff) "Condition Utility (0\u20131)" else "Utility (0\u20131)",
    y     = "Pixel count",
    title = paste0(soil_property, " ", mode_label, " Utility Distribution"),
    subtitle = paste0("Depth: ", gsub("X", "", depth_layer))
  ) +
  theme_pub

hist_png <- file.path(output_base,
                      paste0(soil_property, "_", tolower(mode_label),
                             "_utility_hist_", depth_layer, ".png"))
ggsave(hist_png, p_hist, width = 8, height = 5, dpi = 300)
cat("  Histogram saved: ", hist_png, "\n", sep = "")

# =============================================================================
# 5. Combined panel (graph + map side by side)
# =============================================================================
cat("  Generating combined panel...\n")

p_combined <- plot_grid(
  p_curve, p_map,
  ncol = 2, rel_widths = c(0.45, 0.55),
  labels = c("(a)", "(b)"), label_size = 12
)

combined_png <- file.path(output_base,
                          paste0(soil_property, "_", tolower(mode_label),
                                 "_utility_combined_", depth_layer, ".png"))
ggsave(combined_png, p_combined, width = 16, height = 7, dpi = 300)
cat("  Combined saved : ", combined_png, "\n\n", sep = "")

# =============================================================================
# Summary
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("DONE — ", soil_property, " ", mode_label, " Utility outputs in: ", output_base, "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
