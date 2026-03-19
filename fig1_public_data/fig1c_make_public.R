#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# =========================================================
# Fig1C public script
# - Reads public processed CSV only
# - Recreates CE vs CWM and IM vs ED ROC plots as separate panels
# =========================================================

IN_DIR  <- file.path("processed_data", "fig1")
OUT_DIR <- file.path("output", "fig1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILE_ROC_CE  <- file.path(IN_DIR, "fig1c_roc_ce_cwm.csv")
FILE_ROC_IM  <- file.path(IN_DIR, "fig1c_roc_im_ed.csv")
FILE_AUC     <- file.path(IN_DIR, "fig1c_auc_summary.csv")

roc_ce <- readr::read_csv(FILE_ROC_CE, show_col_types = FALSE)
roc_im <- readr::read_csv(FILE_ROC_IM, show_col_types = FALSE)
auc_df <- readr::read_csv(FILE_AUC, show_col_types = FALSE)

prep_roc <- function(df) {
  stopifnot(all(c("fpr", "sensitivity") %in% names(df)))
  df %>%
    dplyr::filter(is.finite(.data$fpr), is.finite(.data$sensitivity)) %>%
    dplyr::mutate(
      fpr = as.numeric(.data$fpr),
      sensitivity = as.numeric(.data$sensitivity)
    ) %>%
    dplyr::arrange(.data$fpr, .data$sensitivity)
}

get_auc <- function(auc_df, comp_patterns = character()) {
  # accept either comparison-style or positive/negative-area style columns
  if (all(c("comparison", "auc") %in% names(auc_df))) {
    out <- auc_df %>%
      dplyr::filter(.data$comparison %in% comp_patterns) %>%
      dplyr::slice(1) %>%
      dplyr::pull(.data$auc)
    if (length(out) > 0) return(out[[1]])
  }
  if (all(c("positive_area", "negative_area", "auc") %in% names(auc_df))) {
    # CE vs CWM
    if (any(comp_patterns %in% c("CE_vs_CWM", "CE vs CWM"))) {
      out <- auc_df %>%
        dplyr::filter(.data$positive_area == "CE", .data$negative_area == "CWM") %>%
        dplyr::slice(1) %>%
        dplyr::pull(.data$auc)
      if (length(out) > 0) return(out[[1]])
    }
    # IM vs ED
    if (any(comp_patterns %in% c("IM_vs_ED", "IM vs ED"))) {
      out <- auc_df %>%
        dplyr::filter(.data$positive_area == "IM", .data$negative_area == "ED") %>%
        dplyr::slice(1) %>%
        dplyr::pull(.data$auc)
      if (length(out) > 0) return(out[[1]])
    }
  }
  NA_real_
}

make_single_roc <- function(df, auc_value, out_png, out_pdf) {
  p <- ggplot(df, aes(x = fpr, y = sensitivity)) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", linewidth = 0.8, colour = "black"
    ) +
    geom_path(
      linewidth = 1.1,
      colour = "black",
      lineend = "butt",
      linejoin = "mitre"
    ) +
    annotate(
      "text",
      x = 0.73, y = 0.10,
      label = sprintf("AUC = %.3f", auc_value),
      size = 6.2,
      fontface = "plain"
    ) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05), expand = FALSE) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
    ) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
    ) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
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

  ggsave(out_png, p, width = 6.3, height = 6.0, dpi = 400, bg = "white")
  ggsave(out_pdf, p, width = 6.3, height = 6.0, bg = "white")
}

roc_ce <- prep_roc(roc_ce)
roc_im <- prep_roc(roc_im)

auc_ce <- get_auc(auc_df, c("CE_vs_CWM", "CE vs CWM"))
auc_im <- get_auc(auc_df, c("IM_vs_ED", "IM vs ED"))

make_single_roc(
  df = roc_ce,
  auc_value = auc_ce,
  out_png = file.path(OUT_DIR, "Fig1C_ROC_CE_vs_CWM.png"),
  out_pdf = file.path(OUT_DIR, "Fig1C_ROC_CE_vs_CWM.pdf")
)

make_single_roc(
  df = roc_im,
  auc_value = auc_im,
  out_png = file.path(OUT_DIR, "Fig1C_ROC_IM_vs_ED.png"),
  out_pdf = file.path(OUT_DIR, "Fig1C_ROC_IM_vs_ED.pdf")
)

message("Done: Fig1C ROC plots were written to ", OUT_DIR)
