#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# =========================================================
# Fig2A public script
# - Reads GitHub public processed CSV only
# - Recreates average spectra: Carotenoid-low vs Carotenoid-high
# =========================================================

IN_DIR  <- file.path("processed_data", "fig2")
OUT_DIR <- file.path("output", "fig2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILE_SPEC <- file.path(IN_DIR, "fig2a_average_spectra_carotenoidHigh_vs_Low.csv")
FILE_THR  <- file.path(IN_DIR, "fig2a_threshold_info.csv")

spec_df <- readr::read_csv(FILE_SPEC, show_col_types = FALSE)

if (file.exists(FILE_THR)) {
  thr_df <- readr::read_csv(FILE_THR, show_col_types = FALSE)
}

stopifnot(all(c("group", "shift", "mean") %in% names(spec_df)))

spec_df <- spec_df %>%
  dplyr::filter(!is.na(group), !is.na(shift), !is.na(mean)) %>%
  dplyr::mutate(
    group = factor(as.character(group), levels = c("Carotenoid-low", "Carotenoid-high")),
    shift = as.numeric(shift),
    mean  = as.numeric(mean)
  ) %>%
  dplyr::filter(shift >= 900, shift <= 1800)

ref_df <- data.frame(
  shift_plot = c(1156, 1514),
  label_show = c("1156", "1514")
)

line_colors <- c(
  "Carotenoid-low"  = "#F8766D",
  "Carotenoid-high" = "#00BFC4"
)

y_max <- max(spec_df$mean, na.rm = TRUE)
y_min <- min(spec_df$mean, na.rm = TRUE)
y_top <- y_max * 1.20
label_y <- y_max * 1.03

p <- ggplot(spec_df, aes(x = shift, y = mean, color = group)) +
  geom_segment(
    data = ref_df,
    aes(x = shift_plot, xend = shift_plot, y = y_min, yend = y_max * 0.95),
    inherit.aes = FALSE,
    linetype = 2,
    linewidth = 0.8,
    colour = "grey40"
  ) +
  geom_line(linewidth = 1.6) +
  geom_text(
    data = ref_df,
    aes(x = shift_plot, y = label_y, label = label_show),
    inherit.aes = FALSE,
    angle = 90,
    vjust = -0.05,
    size = 6.5,
    fontface = "bold",
    colour = "grey15"
  ) +
  scale_color_manual(values = line_colors, drop = FALSE) +
  coord_cartesian(xlim = c(900, 1800), ylim = c(y_min, y_top)) +
  labs(
    x = expression("Raman shift (cm"^-1*")"),
    y = "Intensity (a.u.)",
    color = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(colour = "#D0D0D0", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.ticks   = element_line(linewidth = 0.8),
    legend.position = "right",
    legend.text = element_text(size = 18),
    plot.margin = margin(10, 18, 10, 10)
  )

ggsave(
  filename = file.path(OUT_DIR, "Fig2A_average_spectra_carotenoidHigh_vs_Low.png"),
  plot = p,
  width = 11.5,
  height = 6.7,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(OUT_DIR, "Fig2A_average_spectra_carotenoidHigh_vs_Low.pdf"),
  plot = p,
  width = 11.5,
  height = 6.7,
  bg = "white"
)

message("Done: Fig2A was written to ", OUT_DIR)
