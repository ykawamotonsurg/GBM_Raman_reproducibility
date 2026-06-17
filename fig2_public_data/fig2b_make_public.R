#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# =========================================================
# Fig2B public script (point + error bar style)
# - Reads GitHub public processed CSV only
# - Fixes clipping / spacing / label placement
# =========================================================

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub(file_arg, "", hit[1]))))
  getwd()
}

IN_DIR  <- script_dir()
OUT_DIR <- file.path(IN_DIR, "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILE_FRAC <- file.path(IN_DIR, "fig2b_regionwise_carotenoidHigh_fraction.csv")

df <- readr::read_csv(FILE_FRAC, show_col_types = FALSE)

stopifnot(all(c("region", "frac_pct") %in% names(df)))

region_levels <- c("CWM", "ED", "IM", "CE")

if (!("n_cells" %in% names(df))) df$n_cells <- NA_integer_
if (!("region_label" %in% names(df))) {
  if (all(is.finite(df$n_cells))) {
    df$region_label <- paste0(df$region, "\n(n=", df$n_cells, ")")
  } else {
    df$region_label <- df$region
  }
}
if (!("pct_label" %in% names(df))) {
  df$pct_label <- sprintf("%.1f%%", df$frac_pct)
}
if (!("ci_lo_pct" %in% names(df))) df$ci_lo_pct <- NA_real_
if (!("ci_hi_pct" %in% names(df))) df$ci_hi_pct <- NA_real_

df <- df %>%
  dplyr::mutate(
    region = factor(as.character(region), levels = region_levels),
    region_label = as.character(region_label),
    frac_pct = as.numeric(frac_pct),
    ci_lo_pct = as.numeric(ci_lo_pct),
    ci_hi_pct = as.numeric(ci_hi_pct),
    pct_label = as.character(pct_label)
  ) %>%
  dplyr::arrange(region)

ymax_data <- max(c(df$ci_hi_pct, df$frac_pct), na.rm = TRUE)
if (!is.finite(ymax_data)) ymax_data <- max(df$frac_pct, na.rm = TRUE)
y_upper <- max(30, ymax_data + 5)

label_df <- df %>%
  dplyr::mutate(
    y_label = ifelse(is.finite(ci_hi_pct), ci_hi_pct + 1.2, frac_pct + 1.2)
  )

p <- ggplot(df, aes(x = region, y = frac_pct)) +
  {
    if (all(is.finite(df$ci_lo_pct)) && all(is.finite(df$ci_hi_pct))) {
      list(
        geom_errorbar(
          aes(ymin = ci_lo_pct, ymax = ci_hi_pct),
          width = 0,
          linewidth = 1.2,
          colour = "black"
        ),
        geom_point(size = 5.0, colour = "black")
      )
    } else {
      geom_point(size = 5.0, colour = "black")
    }
  } +
  geom_text(
    data = label_df,
    aes(y = y_label, label = pct_label),
    size = 6.4
  ) +
  scale_x_discrete(labels = df$region_label, drop = FALSE, expand = expansion(add = 0.35)) +
  scale_y_continuous(
    limits = c(0, y_upper),
    breaks = c(0, 10, 20, 30),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Cell fraction (%)"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.ticks   = element_line(linewidth = 0.8),
    plot.margin = margin(10, 26, 12, 26)
  )

ggsave(
  filename = file.path(OUT_DIR, "Fig2B_regionwise_carotenoidHigh_fraction.png"),
  plot = p,
  width = 8.8,
  height = 6.3,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(OUT_DIR, "Fig2B_regionwise_carotenoidHigh_fraction.pdf"),
  plot = p,
  width = 8.8,
  height = 6.3,
  bg = "white"
)

message("Done: Fig2B was written to ", OUT_DIR)
