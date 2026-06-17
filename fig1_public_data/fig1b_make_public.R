#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# =========================================================
# Fig1B public script
# - Reads public processed CSV only
# - Recreates the final manuscript-style bag score plot
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

FILE_BAG <- file.path(IN_DIR, "fig1b_bag_scores.csv")
FILE_THR <- file.path(IN_DIR, "fig1c_threshold_metrics.csv")

safe_qlogis <- function(p, eps = 1e-6) {
  qlogis(pmin(pmax(p, eps), 1 - eps))
}

bag_df <- readr::read_csv(FILE_BAG, show_col_types = FALSE)
thr_df <- readr::read_csv(FILE_THR, show_col_types = FALSE)

stopifnot(all(c("case_id", "area4", "bag_score") %in% names(bag_df)))

region_levels <- c("CWM", "ED", "IM", "CE")

bag_df <- bag_df %>%
  dplyr::filter(!is.na(case_id), !is.na(area4), !is.na(bag_score)) %>%
  dplyr::mutate(
    area4 = factor(area4, levels = region_levels),
    anchored_score_logit = safe_qlogis(bag_score)
  )

thr_logit <- NA_real_

if (all(c("threshold_label", "threshold_logit") %in% names(thr_df))) {
  thr_logit <- thr_df %>%
    dplyr::filter(.data$threshold_label == "manuscript_threshold") %>%
    dplyr::slice(1) %>%
    dplyr::pull(.data$threshold_logit)
}

if ((length(thr_logit) == 0 || is.na(thr_logit)) &&
    all(c("comparison", "threshold_logit") %in% names(thr_df))) {
  thr_logit <- thr_df %>%
    dplyr::filter(.data$comparison %in% c("IM_vs_ED", "IM vs ED")) %>%
    dplyr::slice(1) %>%
    dplyr::pull(.data$threshold_logit)
}

if (length(thr_logit) == 0) thr_logit <- NA_real_

# Use a fixed seed for deterministic jitter placement
set.seed(20260319)

p <- ggplot(bag_df, aes(x = area4, y = anchored_score_logit)) +
  geom_boxplot(
    width = 0.72,
    outlier.shape = NA,
    linewidth = 0.75,
    colour = "black",
    fill = "white"
  ) +
  geom_jitter(
    width = 0.13,
    height = 0,
    size = 3.2,
    alpha = 0.6,
    colour = "black"
  ) +
  {
    if (length(thr_logit) == 1 && is.finite(thr_logit)) {
      geom_hline(yintercept = thr_logit, linetype = "dashed", linewidth = 0.75)
    }
  } +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "Region",
    y = "Anchored Raman score (logit)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(colour = "#D0D0D0", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.ticks   = element_line(linewidth = 0.8)
  )

ggsave(
  filename = file.path(OUT_DIR, "Fig1B_bag_scores_logit.png"),
  plot = p,
  width = 8.0,
  height = 5.2,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(OUT_DIR, "Fig1B_bag_scores_logit.pdf"),
  plot = p,
  width = 8.0,
  height = 5.2,
  bg = "white"
)

message("Done: Fig1B was written to ", OUT_DIR)
