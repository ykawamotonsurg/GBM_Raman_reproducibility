#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

IN_DIR  <- file.path("processed_data", "fig2")
OUT_DIR <- file.path("output", "fig2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cand_main <- c(
  file.path(IN_DIR, "fig2d_samecell_raman_if_public.csv"),
  file.path(IN_DIR, "fig2e_samecell_raman_if_public.csv")
)
cand_sum <- c(
  file.path(IN_DIR, "fig2d_samecell_summary_p99.csv"),
  file.path(IN_DIR, "fig2e_samecell_summary_p99.csv")
)

FILE_MAIN <- cand_main[file.exists(cand_main)][1]
FILE_SUM  <- cand_sum[file.exists(cand_sum)][1]

if (is.na(FILE_MAIN)) stop("Could not find same-cell public CSV in processed_data/fig2.")
dat_raw <- readr::read_csv(FILE_MAIN, show_col_types = FALSE)
if (!is.na(FILE_SUM)) sum_raw <- readr::read_csv(FILE_SUM, show_col_types = FALSE)

pick_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Missing required column. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x))
  if (is.character(x)) {
    y <- trimws(tolower(x))
    y[y %in% c("true","t","p","pos","positive","1","msr1+","high")] <- "1"
    y[y %in% c("false","f","n","neg","negative","0","msr1-","low")] <- "0"
    y[y %in% c("", "na", "nan")] <- NA_character_
    suppressWarnings(out <- as.integer(y))
    return(out)
  }
  suppressWarnings(as.integer(x))
}

col_case <- pick_col(dat_raw, c("case_id", "case", "case_label", "surgery_date"))
col_cell <- pick_col(dat_raw, c("cell_id", "cell", "cell_local_id"), required = FALSE)
col_caro <- pick_col(dat_raw, c("carotenoid_high", "carotenoid01", "carotenoid"))
col_msr  <- pick_col(dat_raw, c("MSR1_pos_p99", "msr1_pos_p99", "msr01", "msr"))

dat <- dat_raw %>%
  mutate(.tmp_row_id = row_number()) %>%
  transmute(
    case_raw = as.character(.data[[col_case]]),
    cell_id = if (!is.na(col_cell)) suppressWarnings(as.integer(.data[[col_cell]])) else .data$.tmp_row_id,
    carotenoid01 = to01(.data[[col_caro]]),
    msr01 = to01(.data[[col_msr]])
  ) %>%
  filter(!is.na(carotenoid01), !is.na(msr01)) %>%
  filter(carotenoid01 %in% c(0, 1), msr01 %in% c(0, 1))

# Desired display order: current Case2, Case3, Case1 -> display as Case 1, Case 2, Case 3
case_levels_raw <- unique(dat$case_raw)
case_num <- suppressWarnings(as.integer(gsub("^.*?(\\d+)$", "\\1", case_levels_raw)))
if (all(!is.na(case_num))) {
  case_levels_raw <- case_levels_raw[order(case_num)]
}
# rotate order: 2,3,1 when there are exactly 3 cases; otherwise keep existing order
if (length(case_levels_raw) == 3) {
  case_levels_raw <- case_levels_raw[c(2, 3, 1)]
}
case_map <- tibble(
  case_raw = case_levels_raw,
  case_display = paste("Case", seq_along(case_levels_raw))
)

dat <- dat %>%
  left_join(case_map, by = "case_raw") %>%
  mutate(case = factor(case_display, levels = case_map$case_display))

case_2x2 <- dat %>%
  group_by(case) %>%
  summarise(
    A = sum(carotenoid01 == 1 & msr01 == 1),
    B = sum(carotenoid01 == 1 & msr01 == 0),
    C = sum(carotenoid01 == 0 & msr01 == 1),
    D = sum(carotenoid01 == 0 & msr01 == 0),
    .groups = "drop"
  ) %>%
  mutate(
    n_MSR1_pos = A + C,
    n_MSR1_neg = B + D
  )

readr::write_csv(case_2x2, file.path(OUT_DIR, "Fig2D_case_2x2_raw_counts_p99.csv"))
readr::write_csv(case_map, file.path(OUT_DIR, "Fig2D_case_relabel_map.csv"))

N_COL_TILE <- 6
TILE_W <- 0.92
TILE_H <- 0.92

make_tile_positions <- function(df_panel, ncol_tile = 6) {
  if (nrow(df_panel) == 0) {
    return(tibble(order_in_panel = integer(0), x_tile = numeric(0), y_tile = numeric(0)))
  }
  idx <- seq_len(nrow(df_panel))
  tibble(
    order_in_panel = idx,
    x_tile = ((idx - 1) %% ncol_tile) + 1,
    y_tile = -(((idx - 1) %/% ncol_tile) + 1)
  )
}

plot_dat <- dat %>%
  mutate(
    msr_class = factor(ifelse(msr01 == 1, "MSR1+", "MSR1−"), levels = c("MSR1+", "MSR1−")),
    carotenoid_class = factor(
      ifelse(carotenoid01 == 1, "Raman carotenoid-high", "Raman carotenoid-low"),
      levels = c("Raman carotenoid-high", "Raman carotenoid-low")
    )
  ) %>%
  arrange(case, msr_class, desc(carotenoid01), cell_id) %>%
  group_by(case, msr_class) %>%
  mutate(order_in_panel = row_number()) %>%
  group_modify(~ bind_cols(.x, make_tile_positions(.x, ncol_tile = N_COL_TILE))) %>%
  ungroup()

readr::write_csv(plot_dat, file.path(OUT_DIR, "Fig2D_unitplot_source_p99.csv"))

panel_stats <- dat %>%
  mutate(
    msr_class = factor(ifelse(msr01 == 1, "MSR1+", "MSR1−"), levels = c("MSR1+", "MSR1−"))
  ) %>%
  group_by(case, msr_class) %>%
  summarise(
    n_total = n(),
    n_caro_high = sum(carotenoid01 == 1),
    .groups = "drop"
  ) %>%
  complete(
    case = factor(case_map$case_display, levels = case_map$case_display),
    msr_class = factor(c("MSR1+", "MSR1−"), levels = c("MSR1+", "MSR1−")),
    fill = list(n_total = 0, n_caro_high = 0)
  ) %>%
  mutate(label_frac = ifelse(n_total > 0, paste0("carotenoid-high: ", n_caro_high, "/", n_total), NA_character_))

label_df <- panel_stats %>%
  filter(n_total > 0) %>%
  transmute(case, msr_class, x = (N_COL_TILE + 1) / 2, y = 1.35, label = label_frac)

empty_note_df <- panel_stats %>%
  filter(n_total == 0) %>%
  transmute(case, msr_class, x = (N_COL_TILE + 1) / 2, y = -1.8, label = "No evaluable cells")

fill_values <- c(
  "Raman carotenoid-high" = "black",
  "Raman carotenoid-low"  = "white"
)

p <- ggplot() +
  geom_tile(
    data = plot_dat,
    aes(x = x_tile, y = y_tile, fill = carotenoid_class),
    width = TILE_W, height = TILE_H,
    color = "black", linewidth = 0.6
  ) +
  geom_text(data = label_df, aes(x = x, y = y, label = label), size = 4.0) +
  geom_text(data = empty_note_df, aes(x = x, y = y, label = label), size = 4.0, color = "grey25") +
  facet_grid(rows = vars(case), cols = vars(msr_class), switch = "y") +
  scale_fill_manual(values = fill_values, drop = FALSE) +
  scale_x_continuous(limits = c(0.2, N_COL_TILE + 0.8), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3.4, 1.9), breaks = NULL, expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                             override.aes = list(color = "black", linewidth = 0.6))) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 12),
    strip.text.y.left = element_text(face = "bold", angle = 0, size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(1.0, "lines"),
    plot.margin = margin(10, 12, 10, 12)
  )

ggsave(file.path(OUT_DIR, "Fig2D_samecell_unitplot_by_MSR1_p99.png"),
       p, width = 7.2, height = 7.8, units = "in", dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Fig2D_samecell_unitplot_by_MSR1_p99.pdf"),
       p, width = 7.2, height = 7.8, units = "in", bg = "white")

message("Done: Fig2D was written to ", OUT_DIR)
