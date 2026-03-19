#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
})

infile <- file.path(".", "fig3c_program_scores_ir_msr1_state.csv")
outdir <- file.path(".", "outputs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

program_order <- c(
  "Hypoxia–glycolysis",
  "HO-1 axis (HMOX1)",
  "Opsonin sensing",
  "Lysosomal degradation",
  "Receptor uptake",
  "Lipid handling",
  "Efflux/export"
)

dat <- readr::read_csv(infile, show_col_types = FALSE) %>%
  filter(niche %in% c("I", "R")) %>%
  group_by(niche_display, program_label, group_display) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    program_label = factor(program_label, levels = program_order),
    program_label = forcats::fct_rev(program_label),
    niche_display = factor(niche_display, levels = c("Interstitial", "Rim"))
  )

wide <- dat %>%
  tidyr::pivot_wider(names_from = group_display, values_from = score) %>%
  rename(score_neg = `MSR1−`, score_pos = `MSR1+`)

pts <- wide %>%
  select(niche_display, program_label, score_neg, score_pos) %>%
  pivot_longer(cols = c(score_neg, score_pos), names_to = "which", values_to = "value") %>%
  mutate(which = recode(which, score_neg = "MSR1−", score_pos = "MSR1+"))

p <- ggplot() +
  geom_segment(
    data = wide,
    aes(x = score_neg, xend = score_pos, y = program_label, yend = program_label),
    linewidth = 1.1,
    color = "grey25"
  ) +
  geom_point(
    data = pts,
    aes(x = value, y = program_label, color = which),
    size = 3.2,
    alpha = 0.95
  ) +
  scale_color_manual(values = c("MSR1−" = "#2C6DB2", "MSR1+" = "#E6862E"), name = NULL) +
  facet_wrap(~niche_display, nrow = 1) +
  labs(x = "Program score", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 18, 10, 10)
  )

ggsave(file.path(outdir, "Fig3C_state_IR_program7.pdf"), p, width = 13.8, height = 5.2, units = "in", device = cairo_pdf)
ggsave(file.path(outdir, "Fig3C_state_IR_program7.png"), p, width = 13.8, height = 5.2, units = "in", dpi = 300)
