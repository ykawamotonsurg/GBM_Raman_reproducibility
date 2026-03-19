#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# =========================================================
# Fig1D public script
# - Reads public processed CSV only
# - Recreates area-averaged spectra (900–1800 cm^-1)
# =========================================================

IN_DIR  <- file.path("processed_data", "fig1")
OUT_DIR <- file.path("output", "fig1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILE_SPEC <- file.path(IN_DIR, "fig1d_area_averaged_spectra.csv")

spec_df <- readr::read_csv(FILE_SPEC, show_col_types = FALSE)

stopifnot(all(c("area4", "raman_shift", "mean_intensity") %in% names(spec_df)))

x_min <- 900
x_max <- 1800

# plotted dashed lines = nearest sampled positions from the original code
ref_plot  <- c(1002, 1154, 1438, 1512, 1656)
# labels shown = conventional peak labels
ref_label <- c(1002, 1156, 1440, 1514, 1654)

area_levels <- c("CWM","ED","IM","CE")
area_colors <- c(
  "CWM" = "#F8766D",
  "ED"  = "#7CAE00",
  "IM"  = "#00BFC4",
  "CE"  = "#C77CFF"
)

spec_df <- spec_df %>%
  dplyr::filter(!is.na(area4), !is.na(raman_shift), !is.na(mean_intensity)) %>%
  dplyr::mutate(
    area4 = factor(area4, levels = area_levels),
    raman_shift = as.numeric(raman_shift),
    mean_intensity = as.numeric(mean_intensity)
  ) %>%
  dplyr::filter(raman_shift >= x_min, raman_shift <= x_max)

y_max <- max(spec_df$mean_intensity, na.rm = TRUE)
y_min <- min(spec_df$mean_intensity, na.rm = TRUE)
y_span <- y_max - y_min
y_top   <- y_max + 0.18 * y_span
label_y <- y_max + 0.06 * y_span

ref_df <- data.frame(
  shift_plot = ref_plot,
  shift_lab  = ref_label,
  label_show = as.character(ref_label)
)

p <- ggplot(spec_df, aes(x = raman_shift, y = mean_intensity, color = area4)) +
  geom_segment(
    data = ref_df,
    aes(x = shift_plot, xend = shift_plot, y = y_min, yend = y_max),
    inherit.aes = FALSE,
    linetype = 2,
    linewidth = 0.6,
    colour = "grey40"
  ) +
  geom_line(linewidth = 1.2) +
  geom_text(
    data = ref_df,
    aes(x = shift_lab, y = label_y, label = label_show),
    inherit.aes = FALSE,
    angle = 90,
    vjust = 0.15,
    size = 5.0,
    fontface = "bold",
    colour = "grey15"
  ) +
  scale_color_manual(values = area_colors, drop = FALSE) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_top)) +
  labs(
    x = expression("Raman shift (cm"^-1*")"),
    y = "Relative intensity (a.u.)",
    color = "Area"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 16, 10, 10)
  )

ggsave(
  filename = file.path(OUT_DIR, "Fig1D_area_averaged_spectra.png"),
  plot = p,
  width = 12.5,
  height = 7.5,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(OUT_DIR, "Fig1D_area_averaged_spectra.pdf"),
  plot = p,
  width = 12.5,
  height = 7.5,
  bg = "white"
)

message("Done: Fig1D was written to ", OUT_DIR)
