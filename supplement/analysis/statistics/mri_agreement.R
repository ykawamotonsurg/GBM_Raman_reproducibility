# MRI rater agreement from the same measured volumes used for association.
# Pearson, psych::ICC and Bland-Altman calculations retain the saved method.
options(warn=1,digits=17)
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble);library(psych)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop('Usage: Rscript mri_agreement.R REPOSITORY_ROOT NEW_OUTPUT_DIRECTORY')
B<-normalizePath(args[1],winslash='/',mustWork=TRUE);O<-args[2]
if(dir.exists(O))stop('Choose a new output directory')
dir.create(O,recursive=TRUE)
z<-as.data.frame(read_csv(file.path(B,'supplement/data/MRI_CaseData.csv'),show_col_types=FALSE))
stopifnot(nrow(z)==18,!anyDuplicated(z$patient_id),!anyNA(z$patient_id))
df<-data.frame(case=z$patient_id,TC_1=z$rater1_total_volume_mm3,TC_2=z$rater2_total_volume_mm3,VNEC_1=z$rater1_necrotic_volume_mm3,VNEC_2=z$rater2_necrotic_volume_mm3)
stopifnot(all(is.finite(as.matrix(df[,-1]))),all(df$TC_1>0),all(df$TC_2>0),all(df$VNEC_1>=0 & df$VNEC_1<df$TC_1),all(df$VNEC_2>=0 & df$VNEC_2<df$TC_2))
df<-df %>% mutate(EV_1=TC_1-VNEC_1,EV_2=TC_2-VNEC_2,EF_1=(TC_1-VNEC_1)/TC_1,EF_2=(TC_2-VNEC_2)/TC_2)
get_cor_summary <- function(x, y, var_name) {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  pear <- cor.test(x, y, method = "pearson")

  tibble(
    variable = var_name,
    n = length(x),
    pearson_r = unname(pear$estimate),
    pearson_p = pear$p.value
  )
}

get_icc_summary <- function(x, y, var_name) {
  keep <- complete.cases(x, y)
  dat <- data.frame(rater1 = x[keep], rater2 = y[keep])

  icc_res <- psych::ICC(dat)
  out <- as.data.frame(icc_res$results)

  if (!"type" %in% colnames(out)) {
    stop("Column 'type' not found in psych::ICC results.")
  }

  icc21 <- out[out$type == "ICC2", , drop = FALSE]
  icc2k <- out[out$type == "ICC2k", , drop = FALSE]

  if (nrow(icc21) != 1) stop("ICC2 row not found in psych::ICC results.")
  if (nrow(icc2k) != 1) stop("ICC2k row not found in psych::ICC results.")

  tibble(
    variable = var_name,
    n = sum(keep),
    ICC2_1 = as.numeric(icc21$ICC),
    ICC2_1_lwr = as.numeric(icc21$`lower bound`),
    ICC2_1_upr = as.numeric(icc21$`upper bound`),
    ICC2_1_p = as.numeric(icc21$p),
    ICC2_k = as.numeric(icc2k$ICC),
    ICC2_k_lwr = as.numeric(icc2k$`lower bound`),
    ICC2_k_upr = as.numeric(icc2k$`upper bound`),
    ICC2_k_p = as.numeric(icc2k$p)
  )
}

make_ba_df <- function(x, y, case_id = NULL) {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  if (!is.null(case_id)) case_id <- case_id[keep]

  mean_xy <- (x + y) / 2
  diff_xy <- y - x
  bias <- mean(diff_xy)
  sd_diff <- sd(diff_xy)
  loa_upper <- bias + 1.96 * sd_diff
  loa_lower <- bias - 1.96 * sd_diff

  tibble(
    case = if (is.null(case_id)) seq_along(x) else case_id,
    mean_xy = mean_xy,
    diff_xy = diff_xy,
    bias = bias,
    loa_upper = loa_upper,
    loa_lower = loa_lower
  )
}

cor_tc   <- get_cor_summary(df$TC_1,   df$TC_2,   "Tumor volume")
cor_vnec <- get_cor_summary(df$VNEC_1, df$VNEC_2, "Necrotic volume")
cor_ev   <- get_cor_summary(df$EV_1,   df$EV_2,   "Enhancing volume")
cor_ef   <- get_cor_summary(df$EF_1,   df$EF_2,   "Enhancing fraction")

cor_tbl <- bind_rows(cor_tc, cor_vnec, cor_ev, cor_ef)

icc_tc   <- get_icc_summary(df$TC_1,   df$TC_2,   "Tumor volume")
icc_vnec <- get_icc_summary(df$VNEC_1, df$VNEC_2, "Necrotic volume")
icc_ev   <- get_icc_summary(df$EV_1,   df$EV_2,   "Enhancing volume")
icc_ef   <- get_icc_summary(df$EF_1,   df$EF_2,   "Enhancing fraction")

icc_tbl <- bind_rows(icc_tc, icc_vnec, icc_ev, icc_ef)

agreement_tbl <- cor_tbl %>%
  left_join(icc_tbl, by = c("variable", "n")) %>%
  mutate(
    ICC2_1_CI = paste0(sprintf("%.3f", ICC2_1), " [",
                       sprintf("%.3f", ICC2_1_lwr), ", ",
                       sprintf("%.3f", ICC2_1_upr), "]")
  ) %>%
  select(
    variable, n,
    pearson_r, pearson_p,
    ICC2_1, ICC2_1_lwr, ICC2_1_upr, ICC2_1_p,
    ICC2_k, ICC2_k_lwr, ICC2_k_upr, ICC2_k_p,
    ICC2_1_CI
  )


write_csv(agreement_tbl,file.path(O,'mri_agreement.csv'))
for(metric in c('tumor','necrotic')){
 x<-if(metric=='tumor')df$TC_1 else df$VNEC_1
 y<-if(metric=='tumor')df$TC_2 else df$VNEC_2
 out<-make_ba_df(x,y,df$case);names(out)[names(out)=='case']<-'patient_id'
 write_csv(out,file.path(O,paste0('mri_',metric,'_volume.csv')))
}
writeLines(capture.output(sessionInfo()),file.path(O,'sessionInfo.txt'))
writeLines(c('ICC(2,1): psych::ICC type ICC2','Pearson: cor.test method pearson','Bland-Altman: rater 2 minus rater 1; bias +/- 1.96 * sample SD'),file.path(O,'methods.txt'))
