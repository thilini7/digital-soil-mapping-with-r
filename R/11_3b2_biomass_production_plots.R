# =============================================================================
# Biomass Production Utility — Plots Only (Parts 3–6)
# =============================================================================
# Memory-optimised script that loads pre-computed rasters from 11_3b and
# generates maps, histograms, categorical maps, bar charts, and panels.
# Run this after Parts 1–2 have completed and saved TIFs.
#
# Memory strategy (175 M-cell rasters):
#   • Never call values() on full rasters — use terra::global() for stats
#   • Sample ~50 k pixels for histograms
#   • Aggregate to ~500 k cells for plotting
#   • Save each plot to PNG immediately, then rm() + gc()
#   • For multi-panel composites, stitch saved PNGs with magick
#
# Reference:
#   References: see project documentation
# =============================================================================

# --- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(tidyterra)
  library(sf)
  library(ggspatial)
  library(cowplot)
})

set.seed(42)
gc(verbose = FALSE)

# --- Terra options -----------------------------------------------------------
temp_dir <- file.path(getwd(), "temp_terra")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
terraOptions(memfrac = 0.3, tempdir = temp_dir)

# =============================================================================
# CONFIGURATION
# =============================================================================
depth_layer   <- "X0.30cm"
depth_label   <- gsub("X", "", depth_layer)
HIST_SAMPLE_N <- 50000L        # pixels sampled for histograms
AGG_TARGET    <- 500000L       # target cells for map plotting

utility_base <- file.path("Data", "data_out")
output_base  <- file.path(utility_base, "Biomass_production_utility")

cap_tif    <- file.path(output_base,
                        paste0("biomass_overall_capacity_", depth_layer, ".tif"))
con_tif    <- file.path(output_base,
                        paste0("biomass_overall_condition_", depth_layer, ".tif"))
change_tif <- file.path(output_base,
                        paste0("biomass_change_", depth_layer, ".tif"))

# Individual utility raster paths (for Part 6 comparison panel)
raster_paths <- list(
  capacity = list(
    pH  = file.path(utility_base, "pH_utility",
                    paste0("pH_capacity_utility_", depth_layer, ".tif")),
    BD  = file.path(utility_base, "Bulk_Density_utility",
                    paste0("Bulk_Density_capacity_utility_", depth_layer, ".tif")),
    AWC = file.path(utility_base, "AWC_utility",
                    paste0("AWC_capacity_utility_", depth_layer, ".tif"))
  ),
  condition = list(
    pH  = file.path(utility_base, "pH_utility",
                    paste0("pH_condition_utility_", depth_layer, ".tif")),
    BD  = file.path(utility_base, "Bulk_Density_utility",
                    paste0("Bulk_Density_condition_utility_", depth_layer, ".tif")),
    AWC = file.path(utility_base, "AWC_utility",
                    paste0("AWC_condition_utility_", depth_layer, ".tif"))
  )
)

# =============================================================================
# Verify inputs exist
# =============================================================================
cat("Checking pre-computed rasters...\n")
needed  <- c(cap_tif, con_tif, change_tif)
missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  cat("MISSING rasters:\n")
  cat(paste0("  ", missing, "\n"))
  stop("Run 11_3b_biomass_production_utility.R (Parts 1-2) first.")
}
cat("  All 3 summary rasters found.\n\n")

# =============================================================================
# Quick stats via terra::global()  (no values() → zero extra memory)
# =============================================================================
cat("Computing raster statistics (disk-based)...\n")
for (info in list(
  list(tif = cap_tif,    label = "Capacity"),
  list(tif = con_tif,    label = "Condition"),
  list(tif = change_tif, label = "Change")
)) {
  r   <- rast(info$tif)
  s1  <- global(r, fun = "min",   na.rm = TRUE)
  s2  <- global(r, fun = "max",   na.rm = TRUE)
  s3  <- global(r, fun = "mean",  na.rm = TRUE)
  nc  <- global(r, fun = "notNA")
  # Median via small sample (avoids values() on 175 M cells)
  samp <- spatSample(r, size = 50000, method = "random", na.rm = TRUE,
                     as.df = TRUE)
  med  <- median(samp[[1]], na.rm = TRUE)
  cat(sprintf("  %-10s: min=%.4f  med~%.4f  mean=%.4f  max=%.4f  (n=%s)\n",
              info$label, s1$min, med, s3$mean, s2$max,
              format(nc$notNA, big.mark = ",")))
  rm(r, s1, s2, s3, nc, samp, med)
}
cat("\n")
gc(verbose = FALSE)

# =============================================================================
# Helper: sample values from a raster (memory-safe)
# =============================================================================
sample_raster <- function(tif_path, n = HIST_SAMPLE_N) {
  r <- rast(tif_path)
  idx <- spatSample(r, size = n, method = "random", na.rm = TRUE,
                    as.df = TRUE)
  vals <- idx[[1]]
  rm(r, idx); gc(verbose = FALSE)
  return(as.numeric(vals))
}

# =============================================================================
# NSW boundary
# =============================================================================
nsw_boundary <- tryCatch({
  if (requireNamespace("ozmaps", quietly = TRUE)) {
    oz  <- ozmaps::ozmap_states
    nsw <- oz[oz$NAME == "New South Wales", ]
    st_transform(nsw, 4326)
  } else {
    NULL
  }
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
    plot.margin      = margin(10, 10, 10, 10)
  )

# =============================================================================
# Helper: make a map (loads from file, aggregates aggressively)
# =============================================================================
make_map <- function(tif_path, title, subtitle, palette, limits, breaks,
                     diverging = FALSE, midpoint = 0) {
  rast_in <- rast(tif_path)
  agg_factor <- max(1, round(sqrt(ncell(rast_in) / AGG_TARGET)))
  cat(sprintf("    Aggregating %s (factor=%d)...\n", basename(tif_path), agg_factor))
  if (agg_factor > 1) {
    plot_rast <- aggregate(rast_in, fact = agg_factor, fun = "mean",
                           na.rm = TRUE)
  } else {
    plot_rast <- rast_in
  }
  rm(rast_in); gc(verbose = FALSE)

  p <- ggplot() +
    geom_spatraster(data = plot_rast)

  if (diverging) {
    p <- p + scale_fill_gradient2(
      low = "#1A9850", mid = "#FFFFBF", high = "#B2182B",
      midpoint = midpoint,
      limits   = limits,
      breaks   = breaks,
      labels   = sprintf("%.2f", breaks),
      na.value = "transparent",
      name     = "\u0394U"
    )
  } else {
    p <- p + scale_fill_whitebox_c(
      palette  = palette,
      limits   = limits,
      breaks   = breaks,
      labels   = sprintf("%.1f", breaks),
      na.value = "transparent",
      name     = "Utility"
    )
  }

  p <- p +
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
    labs(title = title, subtitle = subtitle,
         x = "Longitude", y = "Latitude") +
    coord_sf(expand = FALSE) +
    theme_pub +
    theme(
      legend.position   = "right",
      legend.key.height = unit(1.8, "cm"),
      legend.key.width  = unit(0.4, "cm"),
      axis.text         = element_text(size = 8)
    )

  return(p)
}

# =============================================================================
# Helper: make a histogram
# =============================================================================
make_hist <- function(vals, title, subtitle, x_lab, fill_col = "#2166AC",
                      xlimits = c(0, 1)) {
  df <- data.frame(val = as.numeric(vals))
  med_val <- median(df$val, na.rm = TRUE)

  ggplot(df, aes(x = val)) +
    geom_histogram(bins = 60, fill = fill_col, colour = "white",
                   linewidth = 0.2, alpha = 0.85) +
    geom_vline(xintercept = med_val, linetype = "dashed",
               colour = "#B2182B", linewidth = 0.8) +
    annotate("text", x = med_val + 0.02, y = Inf, vjust = 2,
             label = paste0("Median = ", round(med_val, 3)),
             colour = "#B2182B", size = 3.5, hjust = 0) +
    scale_x_continuous(limits = xlimits,
                       breaks = seq(xlimits[1], xlimits[2], length.out = 6)) +
    labs(x = x_lab, y = "Pixel count", title = title, subtitle = subtitle) +
    theme_pub
}

# =============================================================================
# Helper: categorical map (from TIF path)
# =============================================================================
make_cat_map <- function(tif_path, title, subtitle, cat_labels, cat_colours) {
  rast_in <- rast(tif_path)
  agg_factor <- max(1, round(sqrt(ncell(rast_in) / AGG_TARGET)))
  if (agg_factor > 1) {
    plot_rast <- aggregate(rast_in, fact = agg_factor, fun = "modal",
                           na.rm = TRUE)
    levels(plot_rast) <- levels(rast_in)
  } else {
    plot_rast <- rast_in
  }
  rm(rast_in); gc(verbose = FALSE)

  p <- ggplot() +
    geom_spatraster(data = plot_rast) +
    scale_fill_manual(
      values   = setNames(cat_colours, cat_labels),
      na.value = "transparent",
      name     = "Category",
      drop     = FALSE
    ) +
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
    labs(title = title, subtitle = subtitle,
         x = "Longitude", y = "Latitude") +
    coord_sf(expand = FALSE) +
    theme_pub +
    theme(
      legend.position   = "right",
      legend.key.height = unit(0.6, "cm"),
      legend.key.width  = unit(0.6, "cm"),
      axis.text         = element_text(size = 8)
    )

  return(p)
}

# =============================================================================
# Helper: bar chart
# =============================================================================
make_bar <- function(cat_table, cat_pct, cat_labels, cat_colours,
                     title, subtitle) {
  df <- data.frame(
    category = factor(cat_labels, levels = cat_labels),
    count    = as.numeric(cat_table),
    pct      = as.numeric(cat_pct)
  )
  ggplot(df, aes(x = category, y = pct, fill = category)) +
    geom_col(colour = "grey30", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = paste0(pct, "%")), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = setNames(cat_colours, cat_labels),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Area (%)", title = title, subtitle = subtitle) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9))
}

# =============================================================================
# PART 3: Generate individual maps (one at a time, save & discard)
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART 3: Generating maps\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# --- 3a. Capacity map --------------------------------------------------------
cat("  Capacity map...\n")
cap_map_png <- file.path(output_base,
                         paste0("biomass_capacity_map_", depth_layer, ".png"))
p <- make_map(
  cap_tif,
  title    = "Capacity - Biomass Production Function",
  subtitle = paste0("Mean(U_pH, U_BD, U_AWC) | Depth: ", depth_label),
  palette  = "muted",
  limits   = c(0, 1),
  breaks   = seq(0, 1, by = 0.2)
)
ggsave(cap_map_png, p, width = 10, height = 8, dpi = 300)
cat("  Saved: ", cap_map_png, "\n")
rm(p); gc(verbose = FALSE)

# --- 3b. Condition map --------------------------------------------------------
cat("  Condition map...\n")
con_map_png <- file.path(output_base,
                         paste0("biomass_condition_map_", depth_layer, ".png"))
p <- make_map(
  con_tif,
  title    = "Condition - Biomass Production Function",
  subtitle = paste0("Mean(U*_pH, U*_BD, U*_AWC) | Depth: ", depth_label),
  palette  = "muted",
  limits   = c(0, 1),
  breaks   = seq(0, 1, by = 0.2)
)
ggsave(con_map_png, p, width = 10, height = 8, dpi = 300)
cat("  Saved: ", con_map_png, "\n")
rm(p); gc(verbose = FALSE)

# --- 3c. Change map (mean ΔU) -------------------------------------------------
cat("  Change map...\n")
chg_r    <- rast(change_tif)
chg_min  <- global(chg_r, fun = "min", na.rm = TRUE)$min
chg_max  <- global(chg_r, fun = "max", na.rm = TRUE)$max
chg_lim  <- max(abs(chg_min), abs(chg_max))
rm(chg_r); gc(verbose = FALSE)
chg_breaks <- pretty(c(-chg_lim, chg_lim), n = 6)

chg_map_png <- file.path(output_base,
                         paste0("biomass_change_map_", depth_layer, ".png"))
p <- make_map(
  change_tif,
  title     = "Biomass Production Change",
  subtitle  = paste0("Mean \u0394U (U_geno \u2212 U_pheno) | Depth: ", depth_label),
  palette   = NULL,
  limits    = range(chg_breaks),
  breaks    = chg_breaks,
  diverging = TRUE,
  midpoint  = 0
)
ggsave(chg_map_png, p, width = 10, height = 8, dpi = 300)
cat("  Saved: ", chg_map_png, "\n\n")
rm(p); gc(verbose = FALSE)

# =============================================================================
# PART 4: Histograms (sampled — 50 k pixels, not 175 M)
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART 4: Generating histograms (sampled)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Capacity
cap_samp <- sample_raster(cap_tif)
p <- make_hist(cap_samp,
               title = "Capacity Distribution",
               subtitle = paste0("Depth: ", depth_label, "  (n=",
                                 format(length(cap_samp), big.mark = ","), " sample)"),
               x_lab = "Capacity Utility (0\u20131)")
cap_hist_png <- file.path(output_base,
                          paste0("biomass_capacity_hist_", depth_layer, ".png"))
ggsave(cap_hist_png, p, width = 8, height = 5, dpi = 300)
rm(p, cap_samp); gc(verbose = FALSE)

# Condition
con_samp <- sample_raster(con_tif)
p <- make_hist(con_samp,
               title = "Condition Distribution",
               subtitle = paste0("Depth: ", depth_label, "  (n=",
                                 format(length(con_samp), big.mark = ","), " sample)"),
               x_lab = "Condition Utility (0\u20131)")
con_hist_png <- file.path(output_base,
                          paste0("biomass_condition_hist_", depth_layer, ".png"))
ggsave(con_hist_png, p, width = 8, height = 5, dpi = 300)
rm(p, con_samp); gc(verbose = FALSE)

# Change
chg_samp <- sample_raster(change_tif)
p <- make_hist(chg_samp,
               title = "Biomass Production Change Distribution",
               subtitle = paste0("Depth: ", depth_label, "  (n=",
                                 format(length(chg_samp), big.mark = ","), " sample)"),
               x_lab = "Mean \u0394U  (<0 = Improved | >0 = Degraded)",
               fill_col = "#7570B3",
               xlimits  = range(chg_breaks))
chg_hist_png <- file.path(output_base,
                          paste0("biomass_change_hist_", depth_layer, ".png"))
ggsave(chg_hist_png, p, width = 8, height = 5, dpi = 300)
rm(p, chg_samp); gc(verbose = FALSE)
cat("  Saved histograms.\n\n")

# =============================================================================
# PART 5: Categorical classification
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART 5: Categorical classification\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Helper: classify raster, save, return freq table (no values() call)
classify_and_save <- function(tif_in, tif_out, brks, labels) {
  r <- rast(tif_in)
  rcl <- matrix(c(brks[-length(brks)], brks[-1], seq_along(labels)), ncol = 3)
  cat_r <- classify(r, rcl = rcl, include.lowest = TRUE)
  levels(cat_r) <- data.frame(id = seq_along(labels), category = labels)
  writeRaster(cat_r, tif_out, overwrite = TRUE)
  rm(r); gc(verbose = FALSE)

  # freq() returns label strings when levels are set
  fr <- freq(cat_r)
  counts <- setNames(rep(0L, length(labels)), labels)
  for (i in seq_len(nrow(fr))) {
    lbl <- as.character(fr$value[i])
    if (lbl %in% labels) counts[lbl] <- fr$count[i]
  }
  rm(cat_r); gc(verbose = FALSE)
  total <- sum(counts)
  pcts  <- round(counts / total * 100, 1)
  return(list(table = counts, pct = pcts, tif = tif_out))
}

# --- 5a. Condition categories -------------------------------------------------
condition_breaks  <- c(0, 0.25, 0.50, 0.75, 0.90, 1.0)
condition_labels  <- c("Highly Degraded", "Degraded",
                       "Slightly Degraded", "No Change", "Improved")
condition_colours <- c("#B2182B", "#EF8A62", "#FDDBC7", "#B8E186", "#1B7837")

cat("  Classifying condition...\n")
con_cat_tif <- file.path(output_base,
                         paste0("biomass_condition_category_", depth_layer, ".tif"))
con_cat <- classify_and_save(con_tif, con_cat_tif, condition_breaks,
                             condition_labels)
cat("  Condition categories:\n")
for (i in seq_along(condition_labels)) {
  cat(sprintf("    %-20s: %8s px  (%5.1f%%)\n",
              condition_labels[i],
              format(con_cat$table[i], big.mark = ","),
              con_cat$pct[i]))
}
cat("\n")

# --- 5b. Capacity categories --------------------------------------------------
capacity_breaks  <- c(0, 0.25, 0.50, 0.75, 1.0)
capacity_labels  <- c("Low", "Moderate", "High", "Very High")
capacity_colours <- c("#D73027", "#FDAE61", "#A6D96A", "#1A9850")

cat("  Classifying capacity...\n")
cap_cat_tif <- file.path(output_base,
                         paste0("biomass_capacity_category_", depth_layer, ".tif"))
cap_cat <- classify_and_save(cap_tif, cap_cat_tif, capacity_breaks,
                             capacity_labels)
cat("  Capacity categories:\n")
for (i in seq_along(capacity_labels)) {
  cat(sprintf("    %-20s: %8s px  (%5.1f%%)\n",
              capacity_labels[i],
              format(cap_cat$table[i], big.mark = ","),
              cap_cat$pct[i]))
}
cat("\n")

# --- 5c. Change categories (Improved / No Change / Degraded) -----------------
change_cat_breaks  <- c(-Inf, -0.05, -0.01, 0.01, 0.05, Inf)
change_cat_labels  <- c("Improved", "Slightly Improved", "No Change",
                        "Slightly Degraded", "Degraded")
change_cat_colours <- c("#1A9850", "#A6D96A", "#FFFFBF", "#FDAE61", "#D73027")

cat("  Classifying change...\n")
chg_cat_tif <- file.path(output_base,
                         paste0("biomass_change_category_", depth_layer, ".tif"))
chg_cat <- classify_and_save(change_tif, chg_cat_tif, change_cat_breaks,
                             change_cat_labels)
cat("  Change categories:\n")
for (i in seq_along(change_cat_labels)) {
  cat(sprintf("    %-22s: %8s px  (%5.1f%%)\n",
              change_cat_labels[i],
              format(chg_cat$table[i], big.mark = ","),
              chg_cat$pct[i]))
}
cat("\n")

# --- 5e. Generate categorical maps (one at a time) ----------------------------
cat("  Generating categorical maps...\n")

con_cat_png <- file.path(output_base,
                         paste0("biomass_condition_category_map_",
                                depth_layer, ".png"))
p <- make_cat_map(con_cat_tif,
                  title = "Biomass Production Condition Classification",
                  subtitle = paste0("Depth: ", depth_label),
                  cat_labels = condition_labels,
                  cat_colours = condition_colours)
ggsave(con_cat_png, p, width = 10, height = 8, dpi = 300)
cat("  Saved: ", con_cat_png, "\n")
rm(p); gc(verbose = FALSE)

cap_cat_png <- file.path(output_base,
                         paste0("biomass_capacity_category_map_",
                                depth_layer, ".png"))
p <- make_cat_map(cap_cat_tif,
                  title = "Biomass Production Capacity Classification",
                  subtitle = paste0("Depth: ", depth_label),
                  cat_labels = capacity_labels,
                  cat_colours = capacity_colours)
ggsave(cap_cat_png, p, width = 10, height = 8, dpi = 300)
cat("  Saved: ", cap_cat_png, "\n")
rm(p); gc(verbose = FALSE)

chg_cat_png <- file.path(output_base,
                         paste0("biomass_change_category_map_",
                                depth_layer, ".png"))
p <- make_cat_map(chg_cat_tif,
                  title = "Biomass Production Change Classification",
                  subtitle = paste0("Mean \u0394U | Depth: ", depth_label),
                  cat_labels = change_cat_labels,
                  cat_colours = change_cat_colours)
ggsave(chg_cat_png, p, width = 10, height = 8, dpi = 300)
cat("  Saved: ", chg_cat_png, "\n")
rm(p); gc(verbose = FALSE)

# --- 5f. Categorical bar charts -----------------------------------------------
cat("  Generating category bar charts...\n")

p <- make_bar(con_cat$table, con_cat$pct, condition_labels, condition_colours,
              title = "Condition Category Distribution",
              subtitle = paste0("Depth: ", depth_label))
con_bar_png <- file.path(output_base,
                         paste0("biomass_condition_category_bar_",
                                depth_layer, ".png"))
ggsave(con_bar_png, p, width = 8, height = 5, dpi = 300)
rm(p); gc(verbose = FALSE)

p <- make_bar(cap_cat$table, cap_cat$pct, capacity_labels, capacity_colours,
              title = "Capacity Category Distribution",
              subtitle = paste0("Depth: ", depth_label))
cap_bar_png <- file.path(output_base,
                         paste0("biomass_capacity_category_bar_",
                                depth_layer, ".png"))
ggsave(cap_bar_png, p, width = 8, height = 5, dpi = 300)
rm(p); gc(verbose = FALSE)

p <- make_bar(chg_cat$table, chg_cat$pct, change_cat_labels, change_cat_colours,
              title = "Change Category Distribution",
              subtitle = paste0("Depth: ", depth_label))
chg_bar_png <- file.path(output_base,
                         paste0("biomass_change_category_bar_",
                                depth_layer, ".png"))
ggsave(chg_bar_png, p, width = 8, height = 5, dpi = 300)
rm(p); gc(verbose = FALSE)
cat("  Saved bar charts.\n\n")

# --- 5g. Categorical combined panel via magick (PNG stitching) ----------------
cat("  Generating categorical combined panel...\n")
if (requireNamespace("magick", quietly = TRUE)) {
  library(magick)

  stitch_row <- function(left_png, right_png) {
    l <- image_read(left_png)  |> image_scale("3000x2400")
    r <- image_read(right_png) |> image_scale("2000x2400")
    image_append(c(l, r))
  }

  row1 <- stitch_row(cap_cat_png, cap_bar_png)
  row2 <- stitch_row(con_cat_png, con_bar_png)
  row3 <- stitch_row(chg_cat_png, chg_bar_png)
  panel <- image_append(c(row1, row2, row3), stack = TRUE)

  cat_panel_png <- file.path(output_base,
                             paste0("biomass_category_panel_",
                                    depth_layer, ".png"))
  image_write(panel, cat_panel_png, format = "png", quality = 95)
  cat("  Saved: ", cat_panel_png, "\n")
  rm(row1, row2, row3, panel); gc(verbose = FALSE)
} else {
  cat("  [magick not installed — skipping combined category panel]\n")
}
cat("\n")

# =============================================================================
# PART 6: Combined panels (PNG stitching via magick)
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART 6: Combined panels\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (requireNamespace("magick", quietly = TRUE)) {
  library(magick)

  # --- 6a. Capacity + Condition (2-panel) ------------------------------------
  l <- image_read(cap_map_png) |> image_scale("3000x2400")
  r <- image_read(con_map_png) |> image_scale("3000x2400")
  panel2 <- image_append(c(l, r))

  cap_con_png <- file.path(output_base,
                           paste0("biomass_capacity_condition_panel_",
                                  depth_layer, ".png"))
  image_write(panel2, cap_con_png, format = "png", quality = 95)
  cat("  Saved: ", cap_con_png, "\n")
  rm(l, r, panel2); gc(verbose = FALSE)

  # --- 6b. Capacity + Condition + Change (3-panel) ---------------------------
  a <- image_read(cap_map_png) |> image_scale("2600x2080")
  b <- image_read(con_map_png) |> image_scale("2600x2080")
  c <- image_read(chg_map_png) |> image_scale("2600x2080")
  panel3 <- image_append(c(a, b, c))

  three_png <- file.path(output_base,
                         paste0("biomass_three_panel_", depth_layer, ".png"))
  image_write(panel3, three_png, format = "png", quality = 95)
  cat("  Saved: ", three_png, "\n")
  rm(a, b, c, panel3); gc(verbose = FALSE)
} else {
  cat("  [magick not installed — skipping combined panels]\n")
}

# --- 6c. Individual property comparison (3x2 grid) ---------------------------
cat("  Individual property comparison panel...\n")
all_indiv <- unlist(raster_paths)
indiv_missing <- all_indiv[!file.exists(all_indiv)]
if (length(indiv_missing) == 0) {
  indiv_pngs <- c()
  for (mode in c("capacity", "condition")) {
    for (prop in c("pH", "BD", "AWC")) {
      prop_label <- switch(prop, pH = "pH", BD = "Bulk Density", AWC = "AWC")
      mode_label <- tools::toTitleCase(mode)
      tag <- paste0(prop, "_", mode)
      tmp_png <- file.path(output_base,
                           paste0("_tmp_", tag, "_", depth_layer, ".png"))
      p <- make_map(
        raster_paths[[mode]][[prop]],
        title    = paste0(prop_label, " ", mode_label),
        subtitle = paste0("Depth: ", depth_label),
        palette  = "muted",
        limits   = c(0, 1),
        breaks   = seq(0, 1, by = 0.2)
      ) + theme(plot.title = element_text(size = 11),
                plot.subtitle = element_text(size = 8),
                axis.text = element_text(size = 6),
                axis.title = element_text(size = 8))
      ggsave(tmp_png, p, width = 8, height = 7, dpi = 200)
      indiv_pngs <- c(indiv_pngs, tmp_png)
      rm(p); gc(verbose = FALSE)
    }
  }

  if (requireNamespace("magick", quietly = TRUE)) {
    imgs <- lapply(indiv_pngs, function(f) image_read(f) |> image_scale("1600x1400"))
    row1 <- image_append(c(imgs[[1]], imgs[[2]], imgs[[3]]))
    row2 <- image_append(c(imgs[[4]], imgs[[5]], imgs[[6]]))
    grid_panel <- image_append(c(row1, row2), stack = TRUE)

    grid_png <- file.path(output_base,
                          paste0("biomass_individual_panel_", depth_layer, ".png"))
    image_write(grid_panel, grid_png, format = "png", quality = 95)
    cat("  Saved: ", grid_png, "\n")
    rm(imgs, row1, row2, grid_panel); gc(verbose = FALSE)
  }

  # Clean up temp files
  file.remove(indiv_pngs)
} else {
  cat("  [Individual utility rasters missing — skipping comparison panel]\n")
}
cat("\n")

# =============================================================================
# Summary
# =============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("DONE — Biomass production plots in: ", output_base, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  Rasters:\n")
cat("    ", cap_tif, "\n")
cat("    ", con_tif, "\n")
cat("    ", change_tif, "\n")
cat("  Category rasters:\n")
cat("    ", cap_cat_tif, "\n")
cat("    ", con_cat_tif, "\n")
cat("    ", chg_cat_tif, "\n")
cat("  Maps:  continuous + categorical (capacity, condition, change)\n")
cat("  Charts: histograms + category bar charts\n")
cat("  Panels: combined continuous, categorical, and individual\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
