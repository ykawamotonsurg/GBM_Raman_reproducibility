#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

infile <- file.path(".", "fig3b_mri_correlation.csv")
outdir <- file.path(".", "outputs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- readr::read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    enhancement_fraction = as.numeric(enhancement_fraction),
    carotenoid_high_cells_in_CE_percent = as.numeric(carotenoid_high_cells_in_CE_percent),
    ce_profiled_cells_n = as.numeric(ce_profiled_cells_n)
  ) %>%
  filter(
    is.finite(enhancement_fraction),
    is.finite(carotenoid_high_cells_in_CE_percent),
    is.finite(ce_profiled_cells_n)
  )

lm_fit <- lm(carotenoid_high_cells_in_CE_percent ~ enhancement_fraction, data = dat)
line_df <- data.frame(
  enhancement_fraction = seq(
    max(0, min(dat$enhancement_fraction) - 0.02),
    min(1, max(dat$enhancement_fraction) + 0.02),
    length.out = 100
  )
)
line_df$carotenoid_high_cells_in_CE_percent <- predict(lm_fit, newdata = line_df)

p <- ggplot(dat, aes(x = enhancement_fraction, y = carotenoid_high_cells_in_CE_percent)) +
  geom_point(size = 3.2, shape = 16, alpha = 0.95) +
  geom_line(
    data = line_df,
    aes(x = enhancement_fraction, y = carotenoid_high_cells_in_CE_percent),
    inherit.aes = FALSE,
    linewidth = 1.05,
    color = "black"
  ) +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10), expand = c(0.01, 0.01)) +
  labs(
    x = "Enhancement fraction",
    y = "Carotenoid-high cells in CE (%)"
  ) +
  theme_classic(base_size = 12.8) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.margin = margin(7, 8, 6, 6, unit = "pt")
  )

ggsave(file.path(outdir, "Fig3B_mri_correlation.pdf"), p, width = 108, height = 88, units = "mm", device = cairo_pdf)
ggsave(file.path(outdir, "Fig3B_mri_correlation.png"), p, width = 108, height = 88, units = "mm", dpi = 600)
