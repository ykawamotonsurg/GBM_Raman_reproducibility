# Temporal evaluation methods adapted to the de-identified input schema.
# Function definitions only. Sourcing this file does not train or select a threshold.
# The fixed outputs are in input.rds and the two s2b_temporal CSV files.
# Required packages: dplyr, recipes, tibble, purrr, xgboost.
# group_id retains the recorded split/CV grouping; CWM is the display label for NAWM.
# The fit and threshold functions retain the saved method, not a guarantee of identical
# newly trained predictions in a different software environment.

suppressPackageStartupMessages({library(dplyr);library(recipes);library(tibble);library(purrr);library(xgboost)})

make_train_styleA <- function(train_full) {
    train_full %>% dplyr::filter(region %in% c("CE", "CWM"))
}

prep_and_bake_pca <- function(train_df, test_df, pca_k = 5) {
    rec <- recipes::recipe(y ~ ., data = train_df) %>% recipes::update_role(group_id, patient_id, spectrum_id, partition, fold, region, new_role = "id") %>% recipes::step_impute_mean(recipes::all_predictors()) %>% recipes::step_zv(recipes::all_predictors()) %>% 
        recipes::step_normalize(recipes::all_predictors()) %>% recipes::step_pca(recipes::all_predictors(), num_comp = pca_k)
    rec_prep <- recipes::prep(rec, training = train_df, retain = TRUE)
    train_baked <- recipes::bake(rec_prep, new_data = NULL)
    test_baked <- recipes::bake(rec_prep, new_data = test_df)
    list(train = train_baked, test = test_baked, recipe = rec_prep)
}

get_x_y <- function(df_baked) {
    x <- df_baked %>% dplyr::select(dplyr::starts_with("PC")) %>% as.matrix()
    y <- df_baked$y
    list(x = x, y = y)
}

predict_xgb_pca <- function(train_baked, test_baked, trees = 300, depth = 3, eta = 0.1) {
    tr <- get_x_y(train_baked)
    te <- get_x_y(test_baked)
    y_num <- ifelse(tr$y == "tumor", 1, 0)
    dtrain <- xgb.DMatrix(data = tr$x, label = y_num)
    dtest <- xgb.DMatrix(data = te$x)
    params <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = depth, eta = eta, subsample = 1, 
        colsample_bytree = 1)
    fit <- xgb.train(params = params, data = dtrain, nrounds = trees, verbose = 0)
    as.numeric(predict(fit, dtest))
}

bag_score_from_vec <- function(v, agg_name) {
    x <- v[!is.na(v)]
    if (length(x) == 0) 
        return(NA_real_)
    if (agg_name == "mean") 
        return(mean(x))
    if (agg_name == "max") 
        return(max(x))
    if (agg_name == "q90") 
        return(as.numeric(stats::quantile(x, probs = 0.9, na.rm = TRUE, names = FALSE)))
    if (grepl("^top\\d+mean$", agg_name)) {
        k <- as.integer(sub("^top(\\d+)mean$", "\\1", agg_name))
        k <- min(k, length(x))
        return(mean(sort(x, decreasing = TRUE)[1:k]))
    }
    stop("Unknown agg_name: ", agg_name)
}

make_bag_scores <- function(point_df, agg_name) {
    point_df %>% dplyr::mutate(group_id = as.character(group_id), region = factor(as.character(region), levels = c("CWM", 
        "ED", "IM", "CE"))) %>% dplyr::group_by(group_id, region) %>% dplyr::summarise(n_points = sum(!is.na(.pred_tumor)), 
        bag_score = bag_score_from_vec(.pred_tumor, agg_name), .groups = "drop") %>% dplyr::mutate(truth_bin = factor(ifelse(region %in% 
        c("CE", "IM"), "tumor", "non_tumor"), levels = c("tumor", "non_tumor")))
}

metrics_im_ed_at_thr <- function(bag_df, thr) {
    sub <- bag_df %>% dplyr::filter(region %in% c("IM", "ED")) %>% dplyr::filter(!is.na(bag_score))
    pred <- ifelse(sub$bag_score >= thr, "tumor", "non_tumor")
    IM_sens <- mean(pred[sub$region == "IM"] == "tumor", na.rm = TRUE)
    ED_spec <- mean(pred[sub$region == "ED"] == "non_tumor", na.rm = TRUE)
    tibble(IM_sens = IM_sens, ED_spec = ED_spec, bal = (IM_sens + ED_spec)/2, minv = pmin(IM_sens, ED_spec), prod = IM_sens * 
        ED_spec)
}

optimize_thr_im_ed <- function(bag_df) {
    s <- bag_df %>% dplyr::filter(region %in% c("IM", "ED")) %>% dplyr::pull(bag_score)
    s <- s[!is.na(s)]
    if (length(s) < 2) 
        stop("Not enough IM/ED bag scores to optimize threshold.")
    u <- sort(unique(s))
    cand <- c(u[1] - 1e-09, u, u[length(u)] + 1e-09)
    grid <- purrr::map_dfr(cand, ~{
        tibble(thr = .x) %>% dplyr::bind_cols(metrics_im_ed_at_thr(bag_df, .x))
    })
    best <- grid %>% dplyr::arrange(dplyr::desc(bal), dplyr::desc(minv), dplyr::desc(prod)) %>% dplyr::slice(1)
    list(thr_best = best$thr[1], best = best, grid = grid)
}

