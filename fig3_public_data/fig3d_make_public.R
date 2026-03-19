#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

infile <- file.path(".", "fig3d_msr1pos_programs_ir_percase.csv")
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
  mutate(
    program_label = factor(program_label, levels = program_order),
    niche_display = factor(niche_display, levels = c("Interstitial", "Rim"))
  ) %>%
  group_by(case_public_id, program_label, niche_display) %>%
  summarise(score = mean(score, na.rm = TRUE), n_cells = sum(n_cells), .groups = "drop")

p <- ggplot(dat, aes(x = niche_display, y = score)) +
  geom_line(aes(group = case_public_id), linewidth = 0.95, alpha = 0.80, color = "grey35") +
  geom_point(aes(size = n_cells), alpha = 0.90, color = "grey10") +
  facet_wrap(~program_label, nrow = 1, scales = "free_y") +
  scale_size_continuous(name = "Cells (n)", range = c(2.0, 5.5)) +
  labs(x = NULL, y = "Program score") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey45"),
    plot.margin = margin(10, 18, 10, 10),
    panel.spacing.x = unit(1.0, "lines")
  )

ggsave(file.path(outdir, "Fig3D_MSR1pos_programs_IR_program7.pdf"), p, width = 18.0, height = 5.0, units = "in", device = cairo_pdf)
ggsave(file.path(outdir, "Fig3D_MSR1pos_programs_IR_program7.png"), p, width = 18.0, height = 5.0, units = "in", dpi = 300)
