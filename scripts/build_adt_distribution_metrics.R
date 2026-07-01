#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(mclust)
  library(mixtools)
  library(ggplot2)
})

options(warn = 1)

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

has_flag <- function(flag) {
  flag %in% commandArgs(trailingOnly = TRUE)
}

print_help <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/build_adt_distribution_metrics.R \\\n",
    "    --seurat-rds tests/data/test_seurat_qc.rds \\\n",
    "    --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \\\n",
    "    --prior-review ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv \\\n",
    "    --outdir ADTriage_output/06_distribution_metrics \\\n",
    "    --group-by Adult_cHSPCs_Ultima_majority_voting\n\n",
    "Options:\n",
    "  --seurat-rds       Input Seurat RDS.\n",
    "  --mapping          Step 1 feature mapping TSV.\n",
    "  --prior-review     Step 5 resolved prior review TSV for per-feature target groups.\n",
    "  --outdir           Output directory.\n",
    "  --assay            ADT assay name. Default: ADT.\n",
    "  --group-by         Metadata column used for target cells. Default: Adult_cHSPCs_Ultima_majority_voting.\n",
    "  --max-k            Maximum mixture components evaluated by BIC. Default: 3.\n",
    "  --min-bic-delta    Minimum BIC improvement to move from k=1 to k=2 or k=2 to k=3. Default: 10.\n",
    "  --fit-sample-size  Max cells used for each GMM fit. Default: 3000.\n",
    "  --min-fit-cells    Minimum positive CLR cells for mixture fitting. Default: 50.\n",
    "  --min-unique       Minimum unique positive CLR values for mixture fitting. Default: 20.\n",
    "  --no-positive-only Fit all finite CLR values instead of CLR > 0. Default: positive-only.\n",
    "  --threshold-method Selected threshold method. Default: noise_3sd.\n",
    "  --plot-positive-only Plot only CLR > 0 values in density PDF. Default: TRUE.\n",
    "  --plot-all-values    Plot all finite CLR values in density PDF.\n",
    "  --margin-den       Density mismatch margin for noise sigma correction. Default: 0.1.\n",
    "  --min-quality-positive-rate  Positive-rate lower bound for pass quality. Default: 5.\n",
    "  --max-quality-positive-rate  Positive-rate upper bound for pass quality. Default: 95.\n",
    "  --max-gmm-overlap-area       Maximum adjacent-component overlap for pass quality. Default: 0.2.\n",
    "  --plot-width       PDF page width in inches. Default: 8.\n",
    "  --plot-height      PDF page height in inches. Default: 6.\n",
    "  --help             Show this help.\n",
    sep = ""
  )
}

if (has_flag("--help")) {
  print_help()
  quit(save = "no", status = 0)
}

seurat_rds <- get_arg("--seurat-rds")
mapping_path <- get_arg("--mapping")
prior_review_path <- get_arg("--prior-review")
outdir <- get_arg("--outdir")
assay_name <- get_arg("--assay", "ADT")
group_by <- get_arg("--group-by", "Adult_cHSPCs_Ultima_majority_voting")
max_k <- as.integer(get_arg("--max-k", "3"))
min_bic_delta <- as.numeric(get_arg("--min-bic-delta", "10"))
fit_sample_size <- as.integer(get_arg("--fit-sample-size", "3000"))
min_fit_cells <- as.integer(get_arg("--min-fit-cells", "50"))
min_unique <- as.integer(get_arg("--min-unique", "20"))
fit_positive_only <- !has_flag("--no-positive-only")
threshold_method <- get_arg("--threshold-method", "noise_3sd")
plot_positive_only <- !has_flag("--plot-all-values")
margin_den <- as.numeric(get_arg("--margin-den", "0.1"))
min_quality_positive_rate <- as.numeric(get_arg("--min-quality-positive-rate", "5"))
max_quality_positive_rate <- as.numeric(get_arg("--max-quality-positive-rate", "95"))
max_gmm_overlap_area <- as.numeric(get_arg("--max-gmm-overlap-area", "0.2"))
plot_width <- as.numeric(get_arg("--plot-width", "8"))
plot_height <- as.numeric(get_arg("--plot-height", "6"))

if (is.null(seurat_rds) || is.null(mapping_path) || is.null(prior_review_path) || is.null(outdir)) {
  print_help()
  stop("--seurat-rds, --mapping, --prior-review, and --outdir are required.", call. = FALSE)
}
if (!max_k %in% c(2L, 3L)) {
  stop("--max-k must be 2 or 3.", call. = FALSE)
}
if (!threshold_method %in% c("noise_3sd", "intersection_1_2")) {
  stop("--threshold-method must be noise_3sd or intersection_1_2.", call. = FALSE)
}

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

get_layer_data <- function(srt, assay, layer = "data") {
  tryCatch(
    GetAssayData(srt, assay = assay, layer = layer),
    error = function(e) GetAssayData(srt, assay = assay, slot = layer)
  )
}

split_groups <- function(value) {
  if (length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(character())
  }
  groups <- trimws(strsplit(value, ";", fixed = TRUE)[[1]])
  groups[nzchar(groups)]
}

join_groups <- function(groups) {
  paste(sort(unique(groups[nzchar(groups)])), collapse = ";")
}

optional_value <- function(df, col, default = "") {
  if (is.null(df) || nrow(df) == 0 || !col %in% colnames(df)) {
    return(default)
  }
  value <- df[[col]][[1]]
  if (length(value) == 0 || is.na(value)) {
    return(default)
  }
  as.character(value)
}

optional_logical <- function(df, col, default = FALSE) {
  value <- optional_value(df, col, default = ifelse(default, "TRUE", "FALSE"))
  toupper(value) %in% c("TRUE", "T", "1", "YES", "Y")
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) NA_real_ else sd(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 7))
}

ratio_safe <- function(numerator, denominator, floor = 0.01) {
  if (!is.finite(numerator)) {
    return(NA_real_)
  }
  numerator / max(denominator, floor, na.rm = TRUE)
}

validate_columns <- function(df, required, name) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(name, " missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

select_fit_values <- function(expr_norm, positive_only = TRUE) {
  x <- expr_norm[is.finite(expr_norm)]
  if (positive_only) {
    x <- x[x > 0]
  }
  x
}

downsample_fit_values <- function(x, n_max = 3000) {
  if (length(x) > n_max) {
    set.seed(42)
    sample(x, n_max)
  } else {
    x
  }
}

compute_bic_values <- function(x, max_k = 3) {
  result <- rep(NA_real_, max_k)
  names(result) <- paste0("k", seq_len(max_k))
  for (k in seq_len(max_k)) {
    fit <- tryCatch(
      mclust::Mclust(x, G = k, modelNames = c("E", "V"), verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$bic)) {
      result[[k]] <- fit$bic
    }
  }
  result
}

select_k_from_bic <- function(bic_values, min_delta = 10) {
  if (!is.finite(bic_values[["k1"]])) {
    return(NA_integer_)
  }
  selected <- 1L
  if (length(bic_values) >= 2 && is.finite(bic_values[["k2"]]) &&
      bic_values[["k2"]] - bic_values[["k1"]] >= min_delta) {
    selected <- 2L
  }
  if (length(bic_values) >= 3 && selected >= 2L && is.finite(bic_values[["k3"]]) &&
      bic_values[["k3"]] - bic_values[["k2"]] >= min_delta) {
    selected <- 3L
  }
  selected
}

fit_normal_mix <- function(x, k) {
  fit <- NULL
  fit_error <- NULL
  invisible(capture.output({
    fit <- tryCatch(
      mixtools::normalmixEM(
        x,
        k = k,
        epsilon = 0.01,
        maxit = 1000,
        maxrestarts = 10,
        verb = FALSE
      ),
      error = function(e) {
        fit_error <<- e
        NULL
      }
    )
  }))
  if (!is.null(fit_error) || is.null(fit)) {
    msg <- if (!is.null(fit_error)) conditionMessage(fit_error) else "normalmixEM returned NULL."
    return(list(status = "fit_failed", message = msg))
  }
  ord <- order(fit$mu)
  list(
    status = "ok",
    method = "normalmixEM",
    message = "",
    k = k,
    mu = as.numeric(fit$mu[ord]),
    sigma = as.numeric(fit$sigma[ord]),
    lambda = as.numeric(fit$lambda[ord]),
    loglik = fit$loglik,
    n_iter = fit$restarts
  )
}

fit_mclust_fallback <- function(x, k) {
  fit <- tryCatch(
    mclust::Mclust(x, G = k, modelNames = c("E", "V"), verbose = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error") || is.null(fit)) {
    msg <- if (inherits(fit, "error")) conditionMessage(fit) else "Mclust returned NULL."
    return(list(status = "fit_failed", message = msg))
  }
  mu <- as.numeric(fit$parameters$mean)
  sigma <- sqrt(as.numeric(fit$parameters$variance$sigmasq))
  if (length(sigma) == 1) {
    sigma <- rep(sigma, k)
  }
  lambda <- as.numeric(fit$parameters$pro)
  if (length(mu) != k || length(sigma) != k || length(lambda) != k ||
      any(!is.finite(mu)) || any(!is.finite(sigma)) || any(sigma <= 0)) {
    return(list(status = "fit_failed", message = "Mclust returned invalid parameters."))
  }
  ord <- order(mu)
  list(
    status = "ok",
    method = paste0("mclust_", fit$modelName),
    message = "",
    k = k,
    mu = mu[ord],
    sigma = sigma[ord],
    lambda = lambda[ord],
    loglik = fit$loglik,
    n_iter = NA_integer_
  )
}

correct_noise_sigma <- function(fit, x, margin_den = 0.1) {
  if (!identical(fit$status, "ok") || length(x) < 2) {
    fit$sigma_corrected <- FALSE
    return(fit)
  }
  kde <- tryCatch(density(x, from = 0), error = function(e) NULL)
  if (is.null(kde)) {
    fit$sigma_corrected <- FALSE
    return(fit)
  }
  mu1 <- fit$mu[[1]]
  sigma1 <- fit$sigma[[1]]
  lambda1 <- fit$lambda[[1]]
  empirical_y <- approx(kde$x, kde$y, xout = mu1, rule = 2)$y
  model_y <- dnorm(mu1, mean = mu1, sd = sigma1) * lambda1
  if (!is.finite(empirical_y) || empirical_y <= 0 || !is.finite(model_y)) {
    fit$sigma_corrected <- FALSE
    return(fit)
  }
  rel_diff <- abs(model_y - empirical_y) / empirical_y
  if (rel_diff <= margin_den) {
    fit$sigma_corrected <- FALSE
    return(fit)
  }
  objective <- function(sd1) {
    pred <- dnorm(mu1, mean = mu1, sd = sd1) * lambda1
    (pred - empirical_y)^2
  }
  opt <- tryCatch(
    optim(par = sigma1, fn = objective, method = "L-BFGS-B", lower = 1e-4, upper = 10),
    error = function(e) NULL
  )
  if (!is.null(opt) && is.finite(opt$par) && opt$par > 0) {
    fit$sigma[[1]] <- as.numeric(opt$par)
    fit$sigma_corrected <- TRUE
  } else {
    fit$sigma_corrected <- FALSE
  }
  fit
}

weighted_intersection <- function(mu1, sd1, w1, mu2, sd2, w2) {
  f <- function(z) {
    w1 * dnorm(z, mean = mu1, sd = sd1) - w2 * dnorm(z, mean = mu2, sd = sd2)
  }
  threshold <- tryCatch({
    if (f(mu1) * f(mu2) <= 0) {
      uniroot(f, interval = c(mu1, mu2))$root
    } else {
      (mu1 * sd2 + mu2 * sd1) / (sd1 + sd2)
    }
  }, error = function(e) {
    (mu1 * sd2 + mu2 * sd1) / (sd1 + sd2)
  })
  as.numeric(threshold)
}

component_overlap <- function(mu1, sd1, w1, mu2, sd2, w2) {
  x_min <- min(mu1 - 8 * sd1, mu2 - 8 * sd2)
  x_max <- max(mu1 + 8 * sd1, mu2 + 8 * sd2)
  grid <- seq(x_min, x_max, length.out = 4000)
  dx <- diff(grid)[[1]]
  sum(pmin(w1 * dnorm(grid, mu1, sd1), w2 * dnorm(grid, mu2, sd2))) * dx
}

summarize_fit <- function(fit, threshold_method = "noise_3sd") {
  if (!identical(fit$status, "ok")) {
    return(fit)
  }
  k <- fit$k
  noise_3sd <- fit$mu[[1]] + 3 * fit$sigma[[1]]
  mid_3sd <- if (k >= 3) fit$mu[[2]] + 3 * fit$sigma[[2]] else NA_real_
  intersection_1_2 <- weighted_intersection(
    fit$mu[[1]], fit$sigma[[1]], fit$lambda[[1]],
    fit$mu[[2]], fit$sigma[[2]], fit$lambda[[2]]
  )
  intersection_2_3 <- if (k >= 3) {
    weighted_intersection(
      fit$mu[[2]], fit$sigma[[2]], fit$lambda[[2]],
      fit$mu[[3]], fit$sigma[[3]], fit$lambda[[3]]
    )
  } else {
    NA_real_
  }
  overlap_1_2 <- component_overlap(
    fit$mu[[1]], fit$sigma[[1]], fit$lambda[[1]],
    fit$mu[[2]], fit$sigma[[2]], fit$lambda[[2]]
  )
  overlap_2_3 <- if (k >= 3) {
    component_overlap(
      fit$mu[[2]], fit$sigma[[2]], fit$lambda[[2]],
      fit$mu[[3]], fit$sigma[[3]], fit$lambda[[3]]
    )
  } else {
    NA_real_
  }
  selected_threshold <- switch(
    threshold_method,
    noise_3sd = noise_3sd,
    intersection_1_2 = intersection_1_2
  )
  fit$threshold_noise_3sd <- noise_3sd
  fit$threshold_mid_3sd <- mid_3sd
  fit$threshold_intersection_1_2 <- intersection_1_2
  fit$threshold_intersection_2_3 <- intersection_2_3
  fit$selected_threshold <- selected_threshold
  fit$selected_threshold_method <- threshold_method
  fit$overlap_1_2 <- overlap_1_2
  fit$overlap_2_3 <- overlap_2_3
  fit$separation_1_2 <- (fit$mu[[2]] - fit$mu[[1]]) / (fit$sigma[[1]] + fit$sigma[[2]])
  fit$separation_2_3 <- if (k >= 3) {
    (fit$mu[[3]] - fit$mu[[2]]) / (fit$sigma[[2]] + fit$sigma[[3]])
  } else {
    NA_real_
  }
  fit
}

fit_marker_distribution <- function(expr_norm) {
  fit_values_raw <- select_fit_values(expr_norm, fit_positive_only)
  if (length(fit_values_raw) < min_fit_cells) {
    return(list(
      status = "fit_failed",
      message = "Too few cells in fit input.",
      fit_input = ifelse(fit_positive_only, "clr_positive_only", "all_finite_clr"),
      fit_n_cells = length(fit_values_raw),
      fit_n_unique = length(unique(fit_values_raw)),
      selected_k = NA_integer_
    ))
  }
  if (length(unique(fit_values_raw)) < min_unique) {
    return(list(
      status = "fit_failed",
      message = "Too few unique values in fit input.",
      fit_input = ifelse(fit_positive_only, "clr_positive_only", "all_finite_clr"),
      fit_n_cells = length(fit_values_raw),
      fit_n_unique = length(unique(fit_values_raw)),
      selected_k = NA_integer_
    ))
  }

  fit_values <- downsample_fit_values(fit_values_raw, fit_sample_size)
  bic_values <- compute_bic_values(fit_values, max_k = max_k)
  selected_k <- select_k_from_bic(bic_values, min_bic_delta)
  base <- list(
    fit_input = ifelse(fit_positive_only, "clr_positive_only", "all_finite_clr"),
    fit_n_cells = length(fit_values),
    fit_n_cells_full = length(fit_values_raw),
    fit_n_unique = length(unique(fit_values_raw)),
    bic_k1 = unname(ifelse("k1" %in% names(bic_values), bic_values[["k1"]], NA_real_)),
    bic_k2 = unname(ifelse("k2" %in% names(bic_values), bic_values[["k2"]], NA_real_)),
    bic_k3 = unname(ifelse("k3" %in% names(bic_values), bic_values[["k3"]], NA_real_)),
    bic_delta_2_vs_1 = unname(ifelse(all(c("k1", "k2") %in% names(bic_values)), bic_values[["k2"]] - bic_values[["k1"]], NA_real_)),
    bic_delta_3_vs_2 = unname(ifelse(all(c("k2", "k3") %in% names(bic_values)), bic_values[["k3"]] - bic_values[["k2"]], NA_real_)),
    selected_k = selected_k
  )

  if (!is.finite(selected_k) || selected_k < 2) {
    return(c(base, list(status = "unimodal_like", message = "BIC did not support k >= 2.")))
  }

  fit <- fit_normal_mix(fit_values, selected_k)
  if (!identical(fit$status, "ok")) {
    fallback <- fit_mclust_fallback(fit_values, selected_k)
    if (!identical(fallback$status, "ok")) {
      return(c(base, list(
        status = "fit_failed",
        message = paste(fit$message, fallback$message, sep = "; ")
      )))
    }
    fit <- fallback
  }
  fit <- correct_noise_sigma(fit, fit_values, margin_den)
  fit <- summarize_fit(fit, threshold_method)
  c(base, fit)
}

compute_target_info <- function(meta, group_by, expected_groups) {
  if (length(expected_groups) == 0) {
    return(list(
      target_idx = integer(),
      status = "no_target_groups",
      message = "No expected_positive_groups available for TNR."
    ))
  }
  values <- as.character(meta[[group_by]])
  values[is.na(values)] <- "NA"
  target_idx <- which(values %in% expected_groups)
  if (length(target_idx) == 0) {
    return(list(
      target_idx = integer(),
      status = "no_target_cells",
      message = "No cells matched expected_positive_groups."
    ))
  }
  list(target_idx = target_idx, status = "ok", message = "")
}

compute_target_raw_metrics <- function(expr_raw, target_idx) {
  result <- list(
    tnr_raw = NA_real_,
    target_dropout_pct = NA_real_,
    target_n_cells = length(target_idx),
    non_target_n_cells = length(expr_raw) - length(target_idx)
  )
  if (length(target_idx) > 0 && length(target_idx) < length(expr_raw)) {
    non_target_idx <- setdiff(seq_along(expr_raw), target_idx)
    result$tnr_raw <- ratio_safe(mean(expr_raw[target_idx]), mean(expr_raw[non_target_idx]))
    result$target_dropout_pct <- mean(expr_raw[target_idx] == 0) * 100
  }
  result
}

compute_threshold_metrics <- function(expr_norm, expr_raw, fit, target_idx, target_raw_metrics) {
  result <- list(
    positive_rate = NA_real_,
    positive_n_cells = NA_integer_,
    negative_n_cells = NA_integer_,
    target_positive_rate = NA_real_,
    non_target_positive_rate = NA_real_,
    snr_raw = NA_real_,
    positive_cv_raw = NA_real_,
    si_raw = NA_real_,
    rd_raw = NA_real_,
    dr_raw = NA_real_
  )
  if (!identical(fit$status, "ok") || !is.finite(fit$selected_threshold)) {
    return(c(result, target_raw_metrics))
  }
  labels <- ifelse(expr_norm > fit$selected_threshold, "Positive", "Negative")
  pos_idx <- which(labels == "Positive")
  neg_idx <- which(labels == "Negative")
  result$positive_rate <- length(pos_idx) / length(expr_norm) * 100
  result$positive_n_cells <- length(pos_idx)
  result$negative_n_cells <- length(neg_idx)
  if (length(pos_idx) > 1 && length(neg_idx) > 1) {
    pos_clr <- expr_norm[pos_idx]
    neg_clr <- expr_norm[neg_idx]
    result$snr_raw <- ratio_safe(mean(expr_raw[pos_idx]), mean(expr_raw[neg_idx]))
    result$positive_cv_raw <- safe_sd(pos_clr) / safe_mean(pos_clr) * 100
    result$si_raw <- (safe_mean(pos_clr) - safe_mean(neg_clr)) / (2 * safe_sd(neg_clr))
    result$rd_raw <- (safe_mean(pos_clr) - safe_mean(neg_clr)) / (safe_sd(pos_clr) + safe_sd(neg_clr))
    result$dr_raw <- safe_mean(pos_clr) - safe_mean(neg_clr)
  }
  if (length(target_idx) > 0 && length(target_idx) < length(expr_norm)) {
    non_target_idx <- setdiff(seq_along(expr_norm), target_idx)
    result$target_positive_rate <- mean(labels[target_idx] == "Positive") * 100
    result$non_target_positive_rate <- mean(labels[non_target_idx] == "Positive") * 100
  }
  c(result, target_raw_metrics)
}

classify_quality <- function(fit, metrics) {
  if (identical(fit$status, "fit_failed")) {
    return(list(status = "fit_failed", message = fit$message))
  }
  if (identical(fit$status, "unimodal_like")) {
    return(list(status = "unimodal_like", message = fit$message))
  }
  if (!identical(fit$status, "ok")) {
    return(list(status = fit$status, message = fit$message))
  }
  if (!is.finite(metrics$positive_rate)) {
    return(list(status = "threshold_uncertain", message = "Selected threshold did not produce a finite positive rate."))
  }
  if (is.finite(fit$overlap_1_2) && fit$overlap_1_2 > max_gmm_overlap_area) {
    return(list(status = "threshold_uncertain", message = "First two components have high overlap."))
  }
  if (metrics$positive_rate < min_quality_positive_rate) {
    return(list(status = "rare_positive", message = "Positive rate is below pass threshold."))
  }
  if (metrics$positive_rate > max_quality_positive_rate) {
    return(list(status = "majority_positive", message = "Most cells are above the selected threshold."))
  }
  if (fit$selected_k == 3) {
    return(list(status = "trimodal_candidate", message = "BIC selected k=3; inspect both threshold candidates."))
  }
  list(status = "bimodal_candidate", message = "")
}

is_finite_number <- function(x) {
  length(x) == 1 && is.finite(x)
}

classify_adt_signal <- function(expected_groups, metrics, raw_mean_umi, dynamic_range_clr, gmm_quality_status, step5_info) {
  positive_rate <- metrics$positive_rate
  positive_n_cells <- metrics$positive_n_cells
  target_positive_rate <- metrics$target_positive_rate
  non_target_positive_rate <- metrics$non_target_positive_rate
  tnr_raw <- metrics$tnr_raw

  has_target <- length(expected_groups) > 0
  has_threshold <- is_finite_number(positive_rate)
  has_positive_n <- is_finite_number(positive_n_cells)
  has_target_positive_rates <- is_finite_number(target_positive_rate) && is_finite_number(non_target_positive_rate)
  has_tnr <- is_finite_number(tnr_raw)
  target_gap <- if (has_target_positive_rates) target_positive_rate - non_target_positive_rate else NA_real_
  relative_target_enrichment <- has_target_positive_rates &&
    target_positive_rate >= 3 * max(non_target_positive_rate, 0.1)
  absolute_target_enrichment <- has_target_positive_rates && target_gap > 5
  target_enriched <- has_tnr && tnr_raw >= 1.5 &&
    (absolute_target_enrichment || relative_target_enrichment)
  step5_has_observed_adt <- length(step5_info$observed_adt_positive_groups) > 0
  step5_no_expected_adt_overlap <- has_target && step5_has_observed_adt &&
    length(step5_info$expected_adt_overlap_groups) == 0
  threshold_off_target <- has_target_positive_rates &&
    target_positive_rate <= non_target_positive_rate
  detected_signal <- has_threshold && is_finite_number(positive_n_cells) && positive_n_cells > 0
  low_positive_count <- has_positive_n && positive_n_cells < 20
  low_global_signal <- is_finite_number(raw_mean_umi) && is_finite_number(dynamic_range_clr) &&
    raw_mean_umi < 0.05 && dynamic_range_clr < 0.2

  if (!has_target) {
    if (!has_threshold) {
      if (low_global_signal) {
        return(list(class = "no_target_no_signal", reason = "No expected target groups and global raw/CLR signal is low."))
      }
      return(list(class = "no_target_threshold_unavailable", reason = "No expected target groups and no usable threshold was selected."))
    }
    if (positive_rate < 5 && low_global_signal) {
      return(list(class = "no_target_no_signal", reason = "No expected target groups, low positive rate, and low global signal."))
    }
    if (positive_rate > 80) {
      return(list(class = "no_target_broad_signal", reason = "No expected target groups but most cells are above threshold."))
    }
    return(list(class = "no_target_detected_signal", reason = "No expected target groups, but thresholded ADT signal is detected."))
  }

  if (!has_threshold) {
    if (has_tnr && tnr_raw >= 1.5) {
      return(list(class = "target_raw_enriched_no_threshold", reason = "Target raw UMI is enriched, but no usable threshold was selected."))
    }
    return(list(class = "no_target_enrichment", reason = "No threshold and no clear target raw enrichment."))
  }

  if (!has_tnr && !has_target_positive_rates) {
    return(list(class = "target_signal_unresolved", reason = "Expected target groups exist, but target enrichment metrics are unavailable."))
  }

  if (target_enriched) {
    if (positive_rate > 80) {
      return(list(class = "target_enriched_broad", reason = "Broad thresholded signal still shows target enrichment."))
    }
    if (low_positive_count) {
      return(list(class = "rare_target_enriched_low_count", reason = "Very few positive cells, but positives are target-enriched; treat as biologically plausible but unstable."))
    }
    if (positive_rate < 5) {
      return(list(class = "rare_target_enriched", reason = "Low global positive rate but target cells are enriched."))
    }
    return(list(class = "target_enriched", reason = "Target cells are enriched by raw TNR and thresholded positive rate."))
  }

  if (detected_signal && (step5_no_expected_adt_overlap || threshold_off_target)) {
    return(list(class = "distribution_detected_but_off_target", reason = "A thresholded ADT-positive population exists, but it is not enriched in the expected target groups."))
  }

  if (has_tnr && has_target_positive_rates &&
      tnr_raw < 1.2 && target_gap <= 5) {
    return(list(class = "no_target_enrichment", reason = "Target and non-target signal are not clearly separated."))
  }

  if (positive_rate < 5) {
    if (low_global_signal) {
      return(list(class = "no_signal", reason = "Low positive rate with low global raw/CLR signal."))
    }
    return(list(class = "rare_or_ambiguous", reason = "Low positive rate without clear target enrichment."))
  }

  if (positive_rate > 80) {
    return(list(class = "broad_positive", reason = "Most cells are above threshold without clear target specificity."))
  }

  if (identical(gmm_quality_status, "trimodal_candidate")) {
    return(list(class = "ambiguous_trimodal", reason = "Trimodal distribution without clear target enrichment under selected threshold."))
  }

  list(class = "ambiguous", reason = "Signal does not meet target-enriched or no-signal rules.")
}

component_value <- function(fit, field, idx) {
  if (!identical(fit$status, "ok") || length(fit[[field]]) < idx) {
    return(NA_real_)
  }
  fit[[field]][[idx]]
}

make_density_plot <- function(expr_norm, fit, metric_row, gmm_row, title) {
  plot_values <- expr_norm[is.finite(expr_norm)]
  if (plot_positive_only) {
    plot_values <- plot_values[plot_values > 0]
  }
  if (length(plot_values) < 2) {
    plot_values <- expr_norm[is.finite(expr_norm)]
  }
  df <- data.frame(CLR = plot_values)
  base <- ggplot(df, aes(x = CLR)) +
    geom_density(fill = "#2C7FB8", color = "#2C7FB8", alpha = 0.42, linewidth = 0.7) +
    labs(title = title, x = "ADT CLR expression", y = "Density") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9, color = "black")
    )

  label_lines <- c(
    sprintf("fit: %s, k=%s, quality=%s", gmm_row$fit_method, gmm_row$selected_k, gmm_row$gmm_quality_status),
    sprintf("plot: %s (n=%s)", ifelse(plot_positive_only, "CLR>0", "all CLR"), length(plot_values)),
    sprintf("threshold(%s)=%.2f", gmm_row$selected_threshold_method, gmm_row$selected_threshold_clr),
    sprintf("pos=%.2f%%, TNR=%s", metric_row$positive_rate, ifelse(is.na(metric_row$tnr_raw), "NA", sprintf("%.2f", metric_row$tnr_raw))),
    sprintf("class=%s", metric_row$adt_signal_class),
    sprintf("target=%s", ifelse(nzchar(metric_row$expected_positive_groups), metric_row$expected_positive_groups, "none"))
  )

  if (identical(fit$status, "ok")) {
    x_min <- min(plot_values, fit$mu - 6 * fit$sigma, na.rm = TRUE)
    x_max <- max(plot_values, fit$mu + 6 * fit$sigma, na.rm = TRUE)
    x_grid <- seq(x_min, x_max, length.out = 1000)
    curve_rows <- lapply(seq_len(fit$selected_k), function(i) {
      data.frame(
        x = x_grid,
        y = fit$lambda[[i]] * dnorm(x_grid, fit$mu[[i]], fit$sigma[[i]]),
        component = paste0("C", i)
      )
    })
    curves <- do.call(rbind, curve_rows)
    base <- base +
      geom_line(data = curves, aes(x = x, y = y, color = component),
                inherit.aes = FALSE, linetype = "dashed", linewidth = 0.7) +
      geom_vline(xintercept = fit$selected_threshold, color = "#B2182B", linetype = "solid", linewidth = 0.65) +
      geom_vline(xintercept = fit$threshold_intersection_1_2, color = "#6A51A3", linetype = "dotdash", linewidth = 0.55) +
      scale_color_manual(values = c("C1" = "#555555", "C2" = "#D95F0E", "C3" = "#238B45"))
    if (is.finite(fit$threshold_intersection_2_3)) {
      base <- base +
        geom_vline(xintercept = fit$threshold_intersection_2_3, color = "#238B45", linetype = "dotdash", linewidth = 0.55)
    }
  } else {
    label_lines <- c(label_lines, paste("status=", gmm_row$gmm_status))
  }

  base +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = paste(label_lines, collapse = "\n"),
      hjust = 1.02,
      vjust = 1.05,
      size = 3.0,
      fill = "white"
    ) +
    theme(legend.position = "none")
}

message("Reading Seurat object: ", seurat_rds)
bm <- readRDS(seurat_rds)

if (!assay_name %in% Assays(bm)) {
  stop("Assay not found: ", assay_name, call. = FALSE)
}
if (!group_by %in% colnames(bm@meta.data)) {
  stop("Metadata column not found: ", group_by, call. = FALSE)
}

message("Normalizing ", assay_name, " assay with CLR margin=2")
bm <- NormalizeData(bm, assay = assay_name, normalization.method = "CLR", margin = 2, verbose = FALSE)

mapping <- read_tsv(mapping_path)
validate_columns(mapping, c("feature_id", "feature_name", "target_class"), "Mapping table")

prior_review <- read_tsv(prior_review_path)
validate_columns(
  prior_review,
  c("feature_id", "feature_name", "group_by", "expected_positive_groups", "prior_method"),
  "Prior review table"
)
prior_review <- prior_review[prior_review$group_by == group_by, , drop = FALSE]

norm_mat <- get_layer_data(bm, assay_name, "data")
raw_mat <- get_layer_data(bm, assay_name, "counts")
assay_features <- rownames(norm_mat)
mapping <- mapping[mapping$feature_name %in% assay_features, , drop = FALSE]
if (nrow(mapping) == 0) {
  stop("No mapping rows matched features in assay: ", assay_name, call. = FALSE)
}

prior_by_id <- split(prior_review, prior_review$feature_id)
meta <- bm@meta.data
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

distribution_rows <- list()
gmm_rows <- list()
plots <- list()

message("Processing ", nrow(mapping), " mapping rows against ", length(assay_features), " ", assay_name, " features")
for (i in seq_len(nrow(mapping))) {
  row <- mapping[i, , drop = FALSE]
  feature_id <- row$feature_id[[1]]
  feature_name <- row$feature_name[[1]]
  target_class <- row$target_class[[1]]

  expr_norm <- as.numeric(norm_mat[feature_name, ])
  expr_raw <- as.numeric(raw_mat[feature_name, ])

  prior_row <- prior_by_id[[feature_id]]
  expected_groups <- character()
  prior_method <- "missing_prior_review"
  prior_confidence <- NA_real_
  observed_adt_positive_groups <- ""
  observed_validation_high_groups <- ""
  adt_validation_overlap_groups <- ""
  expected_adt_overlap_groups <- ""
  expected_validation_overlap_groups <- ""
  step5_needs_manual_review <- FALSE
  step5_review_reason <- ""
  if (!is.null(prior_row) && nrow(prior_row) > 0) {
    expected_groups <- split_groups(prior_row$expected_positive_groups[[1]])
    prior_method <- prior_row$prior_method[[1]]
    prior_confidence <- suppressWarnings(as.numeric(optional_value(prior_row, "prior_confidence", NA_character_)))
    observed_adt_positive_groups <- optional_value(prior_row, "observed_adt_positive_groups")
    observed_validation_high_groups <- optional_value(prior_row, "observed_validation_high_groups")
    adt_validation_overlap_groups <- optional_value(prior_row, "adt_validation_overlap_groups")
    expected_adt_overlap_groups <- optional_value(prior_row, "expected_adt_overlap_groups")
    expected_validation_overlap_groups <- optional_value(prior_row, "expected_validation_overlap_groups")
    step5_needs_manual_review <- optional_logical(prior_row, "needs_manual_review")
    step5_review_reason <- optional_value(prior_row, "review_reason")
  }
  step5_info <- list(
    observed_adt_positive_groups = split_groups(observed_adt_positive_groups),
    observed_validation_high_groups = split_groups(observed_validation_high_groups),
    adt_validation_overlap_groups = split_groups(adt_validation_overlap_groups),
    expected_adt_overlap_groups = split_groups(expected_adt_overlap_groups),
    expected_validation_overlap_groups = split_groups(expected_validation_overlap_groups)
  )

  target_info <- compute_target_info(meta, group_by, expected_groups)
  target_raw_metrics <- compute_target_raw_metrics(expr_raw, target_info$target_idx)
  fit <- fit_marker_distribution(expr_norm)
  threshold_metrics <- compute_threshold_metrics(expr_norm, expr_raw, fit, target_info$target_idx, target_raw_metrics)
  quality <- classify_quality(fit, threshold_metrics)

  q05 <- safe_quantile(expr_norm, 0.05)
  q25 <- safe_quantile(expr_norm, 0.25)
  q75 <- safe_quantile(expr_norm, 0.75)
  q95 <- safe_quantile(expr_norm, 0.95)
  zero_rate <- mean(expr_norm <= 0, na.rm = TRUE) * 100
  nonzero_rate <- mean(expr_norm > 0, na.rm = TRUE) * 100
  raw_mean_umi_value <- safe_mean(expr_raw)
  dynamic_range_value <- q95 - q05
  signal_class <- classify_adt_signal(
    expected_groups,
    threshold_metrics,
    raw_mean_umi_value,
    dynamic_range_value,
    quality$status,
    step5_info
  )

  dist_row <- data.frame(
    feature_id = feature_id,
    feature_name = feature_name,
    target_class = target_class,
    assay = assay_name,
    group_by = group_by,
    expected_positive_groups = join_groups(expected_groups),
    prior_method = prior_method,
    prior_confidence = prior_confidence,
    observed_adt_positive_groups = observed_adt_positive_groups,
    observed_validation_high_groups = observed_validation_high_groups,
    adt_validation_overlap_groups = adt_validation_overlap_groups,
    expected_adt_overlap_groups = expected_adt_overlap_groups,
    expected_validation_overlap_groups = expected_validation_overlap_groups,
    step5_needs_manual_review = step5_needs_manual_review,
    step5_review_reason = step5_review_reason,
    target_status = target_info$status,
    n_cells = length(expr_norm),
    n_nonzero_raw = sum(expr_raw > 0),
    zero_rate = zero_rate,
    nonzero_rate = nonzero_rate,
    mean_clr = safe_mean(expr_norm),
    median_clr = stats::median(expr_norm, na.rm = TRUE),
    sd_clr = safe_sd(expr_norm),
    min_clr = min(expr_norm, na.rm = TRUE),
    q05_clr = q05,
    q25_clr = q25,
    q75_clr = q75,
    q95_clr = q95,
    max_clr = max(expr_norm, na.rm = TRUE),
    dynamic_range_clr = dynamic_range_value,
    raw_mean_umi = raw_mean_umi_value,
    raw_median_umi = stats::median(expr_raw, na.rm = TRUE),
    raw_total_umi = sum(expr_raw, na.rm = TRUE),
    selected_threshold_clr = ifelse(is.null(fit$selected_threshold), NA_real_, fit$selected_threshold),
    selected_threshold_method = ifelse(is.null(fit$selected_threshold_method), "", fit$selected_threshold_method),
    positive_rate = threshold_metrics$positive_rate,
    positive_n_cells = threshold_metrics$positive_n_cells,
    negative_n_cells = threshold_metrics$negative_n_cells,
    target_n_cells = target_raw_metrics$target_n_cells,
    non_target_n_cells = target_raw_metrics$non_target_n_cells,
    target_positive_rate = threshold_metrics$target_positive_rate,
    non_target_positive_rate = threshold_metrics$non_target_positive_rate,
    snr_raw = threshold_metrics$snr_raw,
    tnr_raw = target_raw_metrics$tnr_raw,
    target_dropout_pct = target_raw_metrics$target_dropout_pct,
    positive_cv_raw = threshold_metrics$positive_cv_raw,
    distribution_status = ifelse(identical(fit$status, "ok"), "ok", fit$status),
    gmm_quality_status = quality$status,
    gmm_quality_message = quality$message,
    adt_signal_class = signal_class$class,
    adt_signal_reason = signal_class$reason,
    message = paste(c(fit$message, target_info$message), collapse = "; "),
    stringsAsFactors = FALSE
  )

  gmm_row <- data.frame(
    feature_id = feature_id,
    feature_name = feature_name,
    target_class = target_class,
    assay = assay_name,
    group_by = group_by,
    expected_positive_groups = join_groups(expected_groups),
    prior_method = prior_method,
    prior_confidence = prior_confidence,
    observed_adt_positive_groups = observed_adt_positive_groups,
    observed_validation_high_groups = observed_validation_high_groups,
    adt_validation_overlap_groups = adt_validation_overlap_groups,
    expected_adt_overlap_groups = expected_adt_overlap_groups,
    expected_validation_overlap_groups = expected_validation_overlap_groups,
    step5_needs_manual_review = step5_needs_manual_review,
    step5_review_reason = step5_review_reason,
    fit_input = ifelse(is.null(fit$fit_input), ifelse(fit_positive_only, "clr_positive_only", "all_finite_clr"), fit$fit_input),
    gmm_status = ifelse(identical(fit$status, "ok"), "ok", fit$status),
    gmm_message = fit$message,
    gmm_quality_status = quality$status,
    gmm_quality_message = quality$message,
    adt_signal_class = signal_class$class,
    adt_signal_reason = signal_class$reason,
    selected_k = ifelse(is.null(fit$selected_k), NA_integer_, fit$selected_k),
    fit_method = ifelse(is.null(fit$method), "", fit$method),
    fit_n_cells = ifelse(is.null(fit$fit_n_cells), NA_integer_, fit$fit_n_cells),
    fit_n_cells_full = ifelse(is.null(fit$fit_n_cells_full), NA_integer_, fit$fit_n_cells_full),
    fit_n_unique = ifelse(is.null(fit$fit_n_unique), NA_integer_, fit$fit_n_unique),
    bic_k1 = ifelse(is.null(fit$bic_k1), NA_real_, fit$bic_k1),
    bic_k2 = ifelse(is.null(fit$bic_k2), NA_real_, fit$bic_k2),
    bic_k3 = ifelse(is.null(fit$bic_k3), NA_real_, fit$bic_k3),
    bic_delta_2_vs_1 = ifelse(is.null(fit$bic_delta_2_vs_1), NA_real_, fit$bic_delta_2_vs_1),
    bic_delta_3_vs_2 = ifelse(is.null(fit$bic_delta_3_vs_2), NA_real_, fit$bic_delta_3_vs_2),
    sigma_corrected = ifelse(is.null(fit$sigma_corrected), FALSE, fit$sigma_corrected),
    selected_threshold_clr = ifelse(is.null(fit$selected_threshold), NA_real_, fit$selected_threshold),
    selected_threshold_method = ifelse(is.null(fit$selected_threshold_method), "", fit$selected_threshold_method),
    threshold_noise_3sd = ifelse(is.null(fit$threshold_noise_3sd), NA_real_, fit$threshold_noise_3sd),
    threshold_mid_3sd = ifelse(is.null(fit$threshold_mid_3sd), NA_real_, fit$threshold_mid_3sd),
    threshold_intersection_1_2 = ifelse(is.null(fit$threshold_intersection_1_2), NA_real_, fit$threshold_intersection_1_2),
    threshold_intersection_2_3 = ifelse(is.null(fit$threshold_intersection_2_3), NA_real_, fit$threshold_intersection_2_3),
    component1_mean = component_value(fit, "mu", 1),
    component2_mean = component_value(fit, "mu", 2),
    component3_mean = component_value(fit, "mu", 3),
    component1_sd = component_value(fit, "sigma", 1),
    component2_sd = component_value(fit, "sigma", 2),
    component3_sd = component_value(fit, "sigma", 3),
    component1_weight = component_value(fit, "lambda", 1),
    component2_weight = component_value(fit, "lambda", 2),
    component3_weight = component_value(fit, "lambda", 3),
    separation_1_2 = ifelse(is.null(fit$separation_1_2), NA_real_, fit$separation_1_2),
    separation_2_3 = ifelse(is.null(fit$separation_2_3), NA_real_, fit$separation_2_3),
    overlap_1_2 = ifelse(is.null(fit$overlap_1_2), NA_real_, fit$overlap_1_2),
    overlap_2_3 = ifelse(is.null(fit$overlap_2_3), NA_real_, fit$overlap_2_3),
    positive_rate = threshold_metrics$positive_rate,
    positive_n_cells = threshold_metrics$positive_n_cells,
    negative_n_cells = threshold_metrics$negative_n_cells,
    snr_raw = threshold_metrics$snr_raw,
    tnr_raw = target_raw_metrics$tnr_raw,
    positive_cv_raw = threshold_metrics$positive_cv_raw,
    si_raw = threshold_metrics$si_raw,
    rd_raw = threshold_metrics$rd_raw,
    dr_raw = threshold_metrics$dr_raw,
    target_n_cells = target_raw_metrics$target_n_cells,
    target_positive_rate = threshold_metrics$target_positive_rate,
    non_target_positive_rate = threshold_metrics$non_target_positive_rate,
    target_dropout_pct = target_raw_metrics$target_dropout_pct,
    stringsAsFactors = FALSE
  )

  distribution_rows[[length(distribution_rows) + 1]] <- dist_row
  gmm_rows[[length(gmm_rows) + 1]] <- gmm_row
  title <- sprintf("[%03d] %s / %s", i, feature_id, feature_name)
  plots[[length(plots) + 1]] <- make_density_plot(expr_norm, fit, dist_row, gmm_row, title)
}

distribution_df <- do.call(rbind, distribution_rows)
gmm_df <- do.call(rbind, gmm_rows)

round_numeric_df <- function(df) {
  for (col in colnames(df)) {
    if (is.numeric(df[[col]])) {
      df[[col]] <- round(df[[col]], 6)
    }
  }
  df
}
distribution_df <- round_numeric_df(distribution_df)
gmm_df <- round_numeric_df(gmm_df)

distribution_path <- file.path(outdir, "adt_distribution_metrics.tsv")
gmm_path <- file.path(outdir, "adt_gmm_metrics.tsv")
pdf_path <- file.path(outdir, "adt_density_plots.pdf")

write.table(distribution_df, distribution_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gmm_df, gmm_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("Writing density PDF: ", pdf_path)
pdf(pdf_path, width = plot_width, height = plot_height, onefile = TRUE)
for (p in plots) {
  print(p)
}
dev.off()

message("Wrote: ", distribution_path, " rows=", nrow(distribution_df))
message("Wrote: ", gmm_path, " rows=", nrow(gmm_df))
message("Wrote: ", pdf_path, " pages=", length(plots))
