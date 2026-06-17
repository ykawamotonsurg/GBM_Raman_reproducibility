#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# =========================================================
# Fig1D public script
# - Reads public processed CSV only
# - Recreates patient-level peak comparison panels
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

FILE_PEAK <- file.path(IN_DIR, "fig1d_peak_patient_means.csv")
FILE_PVAL <- file.path(IN_DIR, "fig1d_peak_stats_ed_vs_im.csv")

peak_df <- readr::read_csv(FILE_PEAK, show_col_types = FALSE)
pval_df <- readr::read_csv(FILE_PVAL, show_col_types = FALSE)

stopifnot(all(c("case_id", "Area", "peak_cm1", "value") %in% names(peak_df)))

area_levels <- c("CWM","ED","IM","CE")
area_colors <- c(
  "CWM" = "#F8766D",
  "ED"  = "#7CAE00",
  "IM"  = "#00BFC4",
  "CE"  = "#C77CFF"
)

peak_df <- peak_df %>%
  dplyr::filter(!is.na(case_id), !is.na(Area), !is.na(peak_cm1), !is.na(value)) %>%
  dplyr::mutate(
    Area = factor(as.character(Area), levels = area_levels),
    peak_cm1_int = as.integer(as.character(peak_cm1)),
    value = as.numeric(value)
  )

pval_df <- pval_df %>%
  dplyr::mutate(peak_cm1_int = as.integer(as.character(peak_cm1)))

make_peak_plot <- function(peaks, out_png, out_pdf) {
  dfp <- peak_df %>%
    dplyr::filter(.data$peak_cm1_int %in% peaks) %>%
    dplyr::mutate(
      peak_cm1_fac = factor(as.character(.data$peak_cm1_int),
                            levels = as.character(peaks))
    )

  ymax <- dfp %>%
    dplyr::group_by(.data$peak_cm1_int) %>%
    dplyr::summarise(y = max(.data$value, na.rm = TRUE), .groups = "drop")

  anno <- NULL
  if (all(c("peak_cm1_int", "contrast_label") %in% names(pval_df))) {
    anno <- pval_df %>%
      dplyr::filter(.data$peak_cm1_int %in% peaks) %>%
      dplyr::left_join(ymax, by = "peak_cm1_int") %>%
      dplyr::mutate(
        peak_cm1_fac = factor(as.character(.data$peak_cm1_int),
                              levels = as.character(peaks)),
        x = 2.5,
        y = .data$y * 1.10
      )
  }

  p <- ggplot(dfp, aes(x = Area, y = value, color = Area)) +
    geom_line(aes(group = case_id), colour = "grey80", linewidth = 0.7, alpha = 0.8) +
    geom_boxplot(
      aes(group = Area),
      colour = "black",
      outlier.shape = NA,
      fill = "white",
      linewidth = 1.0
    ) +
    geom_point(size = 3.1, alpha = 0.9) +
    facet_wrap(~ peak_cm1_fac, scales = "free_y", nrow = 1) +
    scale_color_manual(values = area_colors, drop = FALSE) +
    labs(
      x = "Area",
      y = "Relative intensity",
      color = "Area"
    ) +
    theme_bw(base_size = 16) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = "grey50"),
      strip.text = element_text(face = "bold", size = 16),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 20)
    )

  if (!is.null(anno) && nrow(anno) > 0) {
    p <- p + geom_text(
      data = anno,
      aes(x = x, y = y, label = contrast_label),
      inherit.aes = FALSE,
      size = 5
    )
  }

  ggsave(out_png, p, width = 14, height = 6, dpi = 400, bg = "white")
  ggsave(out_pdf, p, width = 14, height = 6, bg = "white")
}

make_peak_plot(
  peaks   = c(1440, 1654),
  out_png = file.path(OUT_DIR, "Fig1D_LipidProtein_1440_1654.png"),
  out_pdf = file.path(OUT_DIR, "Fig1D_LipidProtein_1440_1654.pdf")
)

make_peak_plot(
  peaks   = c(1156, 1514),
  out_png = file.path(OUT_DIR, "Fig1D_Carotenoids_1156_1514.png"),
  out_pdf = file.path(OUT_DIR, "Fig1D_Carotenoids_1156_1514.pdf")
)

message("Done: Fig1D was written to ", OUT_DIR)
