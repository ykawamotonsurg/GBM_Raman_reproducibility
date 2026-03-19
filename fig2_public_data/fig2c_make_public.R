#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# =========================================================
# Fig2C final public script (matched more closely to original)
# - Reads joined public CSV only
# - Uses probability scale for Spearman stats
# - Uses logit scale for x display / regression line
# =========================================================

IN_DIR  <- file.path("processed_data", "fig2")
OUT_DIR <- file.path("output", "fig2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILE_BRIDGE <- file.path(IN_DIR, "fig2c_bridge_plot_data.csv")

fmt_num <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
fmt_num_p_plain <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}
inv_logit <- function(x) 1 / (1 + exp(-x))

spearman_summary <- function(df, x_col, y_col) {
  sub <- df %>% dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  if (nrow(sub) < 4) return(list(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(sub[[x_col]], sub[[y_col]], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value)
}

df <- readr::read_csv(FILE_BRIDGE, show_col_types = FALSE)

stopifnot(all(c("region", "anchored_score_logit", "carotenoid_high_fraction_pct") %in% names(df)))

if ("bag_score" %in% names(df)) {
  df$bag_score_prob <- as.numeric(df$bag_score)
} else {
  df$bag_score_prob <- inv_logit(as.numeric(df$anchored_score_logit))
}

COLORS_REGION <- c(
  "CWM" = "#BDBDBD",
  "ED"  = "#6BAED6",
  "IM"  = "#FDAE6B",
  "CE"  = "#EF3B2C"
)

df <- df %>%
  dplyr::filter(!is.na(.data$region), !is.na(.data$anchored_score_logit), !is.na(.data$carotenoid_high_fraction_pct)) %>%
  dplyr::mutate(
    area4 = factor(as.character(.data$region), levels = c("CWM", "ED", "IM", "CE")),
    bag_score_logit = as.numeric(.data$anchored_score_logit),
    frac_caro_high = as.numeric(.data$carotenoid_high_fraction_pct) / 100,
    bag_score_prob = as.numeric(.data$bag_score_prob)
  )

st_all <- spearman_summary(df, "bag_score_prob", "frac_caro_high")
st_sub <- spearman_summary(df %>% dplyr::filter(area4 %in% c("ED", "IM", "CE")), "bag_score_prob", "frac_caro_high")

ann_frac <- paste0(
  "All regions: \u03c1 = ", fmt_num(st_all$rho), ", P = ", fmt_num_p_plain(st_all$p),
  "\nED/IM/CE only: \u03c1 = ", fmt_num(st_sub$rho), ", P = ", fmt_num_p_plain(st_sub$p)
)

x_break_probs <- c(0.01, 0.10, 0.50, 0.90, 0.99)
x_breaks <- qlogis(x_break_probs)
x_labels <- percent_format(accuracy = 1)(x_break_probs)

frac_y_top <- max(df$frac_caro_high, na.rm = TRUE) * 1.10
frac_y_breaks <- pretty(c(0, frac_y_top), n = 5)
frac_y_breaks <- frac_y_breaks[frac_y_breaks >= 0 & frac_y_breaks <= frac_y_top]

p <- ggplot(df, aes(x = bag_score_logit, y = frac_caro_high, color = area4)) +
  geom_point(size = 3.0, alpha = 0.92, position = position_jitter(width = 0.05, height = 0)) +
  geom_smooth(
    aes(x = bag_score_logit, y = frac_caro_high),
    inherit.aes = FALSE,
    method = "lm", se = FALSE,
    linewidth = 0.7, color = "black"
  ) +
  annotate(
    "text",
    x = min(df$bag_score_logit, na.rm = TRUE) + 0.02 * diff(range(df$bag_score_logit, na.rm = TRUE)),
    y = frac_y_top * 0.985,
    label = ann_frac,
    hjust = 0, vjust = 1, size = 6.0
  ) +
  scale_color_manual(values = COLORS_REGION, drop = FALSE) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = frac_y_breaks,
    limits = c(NA, frac_y_top),
    expand = expansion(add = c(0.004, 0))
  ) +
  labs(
    x = "Anchored Raman score (logit display)",
    y = "Carotenoid-high fraction",
    color = "Region"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = grid::unit(2.5, "pt"),
    legend.title = element_text(face = "bold", color = "black", size = 18),
    legend.text  = element_text(color = "black", size = 14),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "bottom",
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  filename = file.path(OUT_DIR, "Fig2C_anchoredScore_vs_carotenoidFraction.png"),
  plot = p, width = 12.3, height = 5.5, dpi = 500, bg = "white"
)
ggsave(
  filename = file.path(OUT_DIR, "Fig2C_anchoredScore_vs_carotenoidFraction.pdf"),
  plot = p, width = 12.3, height = 5.5, bg = "white"
)

message("Done: Fig2C was written to ", OUT_DIR)
