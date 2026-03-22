# =============================================================================
# Calibrate k for Condition Utility & Sensitivity Analysis
# =============================================================================
# Companion to 11_3_soil_utility_map.R
#
# Part 1: Data-driven k calibration
#   Computes ΔU from per-class genosoil/phenosoil means, then picks k so
#   that the median ΔU maps to a target condition score (default 0.5).
#   k = -ln(target) / median(ΔU)
#
# Part 2: Sensitivity analysis
#   Generates condition utility maps for k = 5, 10, 15, 20 side by side.
# =============================================================================

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
soil_property  <- "AWC"             # <-- "pH" or "Bulk_Density" or "AWC"
depth_layer    <- "X0.30cm"
target_score   <- 0.5              # desired condition score at median ΔU
k_candidates   <- c(5, 10, 15, 20)  # sensitivity analysis values

# --- Capacity function per property (used to compute U_geno / U_pheno) -------
capacity_config <- list(
  pH = list(
    fn     = function(x, params) exp(-(x - params$mu)^2 / 2),
    params = list(mu = 6.5)
  ),
  Bulk_Density = list(
    fn     = function(x, params) 1 / (1 + exp(params$a * (x - params$b))),
    params = list(a = 10, b = 1.5)
  ),
  AWC = list(
    fn     = function(x, params) 1 / (1 + exp(-params$c * (x - params$d))),
    params = list(c = 1.0, d = 12)
  )
)

if (!soil_property %in% names(capacity_config))
  stop("No capacity config for '", soil_property,
       "'. Available: ", paste(names(capacity_config), collapse = ", "))

cap_fn     <- capacity_config[[soil_property]]$fn
cap_params <- capacity_config[[soil_property]]$params

# --- Paths -------------------------------------------------------------------
input_base   <- file.path("Data", "data_out", "pedogenon_zonal_means")
geno_raster  <- file.path(input_base, "Geno",
                          paste0(soil_property, "_geno_pedogenon_mean_",
                                 depth_layer, ".tif"))
pheno_raster <- file.path(input_base, "Pheno",
                          paste0(soil_property, "_pheno_pedogenon_mean_",
                                 depth_layer, ".tif"))
zonal_csv    <- file.path("Data", "data_out", "condition_maps", soil_property,
                          paste0(soil_property, "_condition_per_class_",
                                 depth_layer, ".csv"))

output_base  <- file.path("Data", "data_out",
                          paste0(soil_property, "_utility"),
                          "k_calibration")
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

# NSW boundary
nsw_boundary <- tryCatch({
  if (requireNamespace("ozmaps", quietly = TRUE)) {
    oz <- ozmaps::ozmap_states
    nsw <- oz[oz$NAME == "New South Wales", ]
    st_transform(nsw, 4326)
  } else NULL
}, error = function(e) NULL)

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
    plot.margin      = margin(5, 5, 5, 5)
  )

# =============================================================================
# PART 1 — Data-driven k Calibration from Per-class Means
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART 1: Data-driven k Calibration\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (!file.exists(zonal_csv)) stop("Per-class CSV not found: ", zonal_csv)
zonal_df <- read.csv(zonal_csv)

# Compute per-class utility difference using capacity function
zonal_df$u_geno  <- cap_fn(zonal_df$geno_mean, cap_params)
zonal_df$u_pheno <- cap_fn(zonal_df$pheno_mean, cap_params)
zonal_df$u_diff  <- zonal_df$u_geno - zonal_df$u_pheno

# Only consider degraded classes (ΔU > 0) for calibration
degraded <- zonal_df$u_diff[!is.na(zonal_df$u_diff) & zonal_df$u_diff > 0]

cat("  Total pedogenon classes : ", nrow(zonal_df), "\n", sep = "")
cat("  Degraded classes (ΔU>0): ", length(degraded), "\n", sep = "")
cat("  Improved classes (ΔU<0): ", sum(zonal_df$u_diff < 0, na.rm = TRUE), "\n", sep = "")
cat("  No-change classes (ΔU=0): ", sum(zonal_df$u_diff == 0, na.rm = TRUE), "\n", sep = "")
cat("  NA classes             : ", sum(is.na(zonal_df$u_diff)), "\n\n", sep = "")

cat("  ΔU distribution (degraded classes only):\n")
cat("    Min    : ", round(min(degraded, na.rm = TRUE), 6), "\n", sep = "")
cat("    Q1     : ", round(quantile(degraded, 0.25, na.rm = TRUE), 6), "\n", sep = "")
cat("    Median : ", round(median(degraded, na.rm = TRUE), 6), "\n", sep = "")
cat("    Mean   : ", round(mean(degraded, na.rm = TRUE), 6), "\n", sep = "")
cat("    Q3     : ", round(quantile(degraded, 0.75, na.rm = TRUE), 6), "\n", sep = "")
cat("    Max    : ", round(max(degraded, na.rm = TRUE), 6), "\n\n", sep = "")

# Calibrate k: k = -ln(target) / median(ΔU)
median_du   <- median(degraded, na.rm = TRUE)
k_calibrated <- -log(target_score) / median_du

cat("  Target condition score at median ΔU: ", target_score, "\n", sep = "")
cat("  Calibrated k = -ln(", target_score, ") / ", round(median_du, 6),
    " = ", round(k_calibrated, 2), "\n\n", sep = "")

# Show what each candidate k produces at the median ΔU
cat("  Condition score at median ΔU for candidate k values:\n")
for (k_val in c(k_candidates, round(k_calibrated, 1))) {
  score <- exp(-k_val * median_du)
  marker <- if (abs(k_val - round(k_calibrated, 1)) < 0.01) " <-- calibrated" else ""
  cat("    k = ", sprintf("%5.1f", k_val), " → U_condition = ",
      round(score, 4), marker, "\n", sep = "")
}

# --- ΔU histogram ------------------------------------------------------------
cat("\n  Generating ΔU distribution plot...\n")

p_du_hist <- ggplot(zonal_df, aes(x = u_diff)) +
  geom_histogram(bins = 40, fill = "#2166AC", colour = "white",
                 linewidth = 0.2, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "forestgreen",
             linewidth = 0.7) +
  geom_vline(xintercept = median_du, linetype = "dashed", colour = "#B2182B",
             linewidth = 0.8) +
  annotate("text", x = median_du + 0.005, y = Inf, vjust = 2,
           label = paste0("Median(degraded) = ", round(median_du, 4)),
           colour = "#B2182B", size = 3.5, hjust = 0) +
  labs(
    x        = expression(Delta*"U  ("*U[geno] - U[pheno]*")"),
    y        = "Pedogenon class count",
    title    = expression(Delta*"U Distribution Across Pedogenon Classes"),
    subtitle = paste0("Depth: ", gsub("X", "", depth_layer))
  ) +
  theme_pub

du_hist_png <- file.path(output_base,
                         paste0(soil_property, "_delta_u_distribution_",
                                depth_layer, ".png"))
ggsave(du_hist_png, p_du_hist, width = 8, height = 5, dpi = 300)
cat("  Saved: ", du_hist_png, "\n", sep = "")

# --- Condition curves for different k ----------------------------------------
cat("  Generating k comparison curves...\n")

k_all <- sort(unique(c(k_candidates, round(k_calibrated, 1))))
x_seq <- seq(-0.2, max(degraded, na.rm = TRUE) * 1.5, length.out = 500)

curve_list <- lapply(k_all, function(k_val) {
  data.frame(
    du        = x_seq,
    condition = ifelse(x_seq < 0, 1, exp(-k_val * x_seq)),
    k         = paste0("k = ", k_val)
  )
})
curve_all <- do.call(rbind, curve_list)
curve_all$k <- factor(curve_all$k, levels = paste0("k = ", k_all))

p_curves <- ggplot(curve_all, aes(x = du, y = condition, colour = k)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = target_score, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = median_du, linetype = "dashed",
             colour = "#B2182B", linewidth = 0.5) +
  annotate("text", x = median_du + 0.003, y = 0.05,
           label = paste0("median ΔU = ", round(median_du, 4)),
           colour = "#B2182B", size = 3, hjust = 0) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.05)) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x        = expression(Delta*"U  ("*U[geno] - U[pheno]*")"),
    y        = "Condition Utility",
    title    = "Condition Utility Curves for Different k Values",
    subtitle = paste0("Depth: ", gsub("X", "", depth_layer)),
    colour   = NULL
  ) +
  theme_pub +
  theme(legend.position = c(0.85, 0.80))

curves_png <- file.path(output_base,
                        paste0(soil_property, "_k_comparison_curves_",
                               depth_layer, ".png"))
ggsave(curves_png, p_curves, width = 8, height = 5, dpi = 300)
cat("  Saved: ", curves_png, "\n", sep = "")

# --- Save calibration summary CSV --------------------------------------------
calib_summary <- data.frame(
  k     = k_all,
  score_at_median_du = exp(-k_all * median_du),
  score_at_q1_du     = exp(-k_all * quantile(degraded, 0.25, na.rm = TRUE)),
  score_at_q3_du     = exp(-k_all * quantile(degraded, 0.75, na.rm = TRUE))
)
calib_csv <- file.path(output_base,
                       paste0(soil_property, "_k_calibration_summary_",
                              depth_layer, ".csv"))
write.csv(calib_summary, calib_csv, row.names = FALSE)
cat("  Saved: ", calib_csv, "\n\n", sep = "")

# =============================================================================
# PART 2 — Sensitivity Analysis: Maps for Each k
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART 2: Sensitivity Analysis — Maps for k = ",
    paste(k_candidates, collapse = ", "), "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (!file.exists(geno_raster))  stop("Genosoil raster not found: ",  geno_raster)
if (!file.exists(pheno_raster)) stop("Phenosoil raster not found: ", pheno_raster)

geno_rast  <- rast(geno_raster)
pheno_rast <- rast(pheno_raster)

cat("  Computing U_geno and U_pheno...\n")
u_geno  <- app(geno_rast,  fun = function(x) cap_fn(x, cap_params))
u_pheno <- app(pheno_rast, fun = function(x) cap_fn(x, cap_params))
u_diff_rast <- u_geno - u_pheno

# Aggregate for faster plotting
target_ncell <- 2e6
agg_factor <- max(1, round(sqrt(ncell(u_diff_rast) / target_ncell)))
if (agg_factor > 1) {
  cat("  Aggregating by factor ", agg_factor, "...\n", sep = "")
  u_diff_plot <- aggregate(u_diff_rast, fact = agg_factor, fun = "mean",
                           na.rm = TRUE)
} else {
  u_diff_plot <- u_diff_rast
}

# Generate a map for each k value
map_list <- list()

for (k_val in k_candidates) {
  cat("  Generating map for k = ", k_val, "...\n", sep = "")

  cond_rast <- app(u_diff_plot, fun = function(d) ifelse(d < 0, 1, exp(-k_val * d)))
  names(cond_rast) <- "Utility"

  vals <- values(cond_rast, mat = FALSE)
  vals <- vals[!is.na(vals)]

  p <- ggplot() +
    geom_spatraster(data = cond_rast) +
    scale_fill_whitebox_c(
      palette  = "muted",
      limits   = c(0, 1),
      breaks   = seq(0, 1, by = 0.2),
      labels   = sprintf("%.1f", seq(0, 1, by = 0.2)),
      na.value = "transparent",
      name     = "Utility"
    ) +
    {if (!is.null(nsw_boundary))
      geom_sf(data = nsw_boundary, fill = NA, colour = "black",
              linewidth = 0.3)
    } +
    labs(
      title    = paste0("k = ", k_val),
      subtitle = paste0("Median = ", round(median(vals), 3)),
      x = NULL, y = NULL
    ) +
    coord_sf(expand = FALSE) +
    theme_pub +
    theme(
      legend.position   = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      axis.text         = element_text(size = 6),
      plot.title        = element_text(size = 12),
      plot.subtitle     = element_text(size = 9)
    )

  map_list[[as.character(k_val)]] <- p

  # Save individual map
  indiv_png <- file.path(output_base,
                         paste0(soil_property, "_condition_utility_k",
                                k_val, "_", depth_layer, ".png"))
  ggsave(indiv_png, p, width = 8, height = 7, dpi = 300)

  # Save individual raster (full resolution)
  full_cond <- app(u_diff_rast, fun = function(d) ifelse(d < 0, 1, exp(-k_val * d)))
  indiv_tif <- file.path(output_base,
                         paste0(soil_property, "_condition_utility_k",
                                k_val, "_", depth_layer, ".tif"))
  writeRaster(full_cond, indiv_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  cat("    Saved: ", indiv_png, "\n", sep = "")
}

# --- Combined sensitivity panel ----------------------------------------------
cat("\n  Generating combined sensitivity panel...\n")

p_panel <- plot_grid(
  plotlist = map_list,
  ncol     = 2,
  labels   = paste0("(", letters[seq_along(k_candidates)], ")"),
  label_size = 11
)

panel_png <- file.path(output_base,
                       paste0(soil_property, "_k_sensitivity_panel_",
                              depth_layer, ".png"))
ggsave(panel_png, p_panel, width = 16, height = 14, dpi = 300)
cat("  Saved: ", panel_png, "\n\n", sep = "")

# =============================================================================
# Summary
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("DONE\n")
cat("  Calibrated k   : ", round(k_calibrated, 2),
    "  (median ΔU = ", round(median_du, 6),
    ", target score = ", target_score, ")\n", sep = "")
cat("  Sensitivity maps: k = ", paste(k_candidates, collapse = ", "), "\n", sep = "")
cat("  Outputs in      : ", output_base, "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n")
