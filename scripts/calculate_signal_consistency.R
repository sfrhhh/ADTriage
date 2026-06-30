#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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
    "  Rscript scripts/calculate_signal_consistency.R \\\n",
    "    --group-summary-dir ADTriage_output/03_group_summary \\\n",
    "    --outdir ADTriage_output/04_signal_consistency\n\n",
    "Options:\n",
    "  --group-summary-dir  Directory containing Step 3 group summary TSV files.\n",
    "  --outdir             Output directory for Step 4.\n",
    "  --min-groups         Minimum complete groups required for correlation. Default: 3.\n",
    "  --help               Show this help.\n",
    sep = ""
  )
}

if (has_flag("--help")) {
  print_help()
  quit(save = "no", status = 0)
}

summary_dir <- get_arg("--group-summary-dir")
outdir <- get_arg("--outdir")
min_groups <- as.integer(get_arg("--min-groups", "3"))

if (is.null(summary_dir) || is.null(outdir)) {
  print_help()
  stop("--group-summary-dir and --outdir are required.", call. = FALSE)
}
if (is.na(min_groups) || min_groups < 2) {
  stop("--min-groups must be >= 2.", call. = FALSE)
}

required_columns <- function(data, columns, label) {
  missing <- setdiff(columns, colnames(data))
  if (length(missing) > 0) {
    stop(label, " missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

safe_cor <- function(x, y, method) {
  if (length(x) < 2 || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = method, use = "complete.obs"))
}

safe_lm <- function(x, y) {
  if (length(x) < 2 || stats::sd(x) == 0) {
    return(list(slope = NA_real_, r_squared = NA_real_))
  }
  fit <- stats::lm(y ~ x)
  list(
    slope = unname(stats::coef(fit)[["x"]]),
    r_squared = unname(summary(fit)$r.squared)
  )
}

format_stat <- function(value) {
  if (is.na(value)) {
    "NA"
  } else {
    sprintf("%.3f", value)
  }
}

clean_file_part <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

calculate_summary <- function(data, validation_col, validation_label, min_groups) {
  keys <- unique(data[, c("group_by", "feature_id", "feature_name", "target_class", "validation_assay", "validation_feature")])
  records <- vector("list", nrow(keys))

  for (i in seq_len(nrow(keys))) {
    key <- keys[i, , drop = FALSE]
    subset_idx <- data$group_by == key$group_by &
      data$feature_id == key$feature_id &
      data$feature_name == key$feature_name &
      data$validation_feature == key$validation_feature
    sub <- data[subset_idx, , drop = FALSE]
    complete <- sub[sub$status == "ok" &
      is.finite(sub$mean_ADT_CLR) &
      is.finite(sub[[validation_col]]), , drop = FALSE]

    n_groups <- nrow(complete)
    min_cells <- if (n_groups > 0) min(complete$n_cells) else NA_integer_
    pearson <- if (n_groups >= min_groups) safe_cor(complete[[validation_col]], complete$mean_ADT_CLR, "pearson") else NA_real_
    spearman <- if (n_groups >= min_groups) safe_cor(complete[[validation_col]], complete$mean_ADT_CLR, "spearman") else NA_real_
    lm_stats <- if (n_groups >= min_groups) {
      safe_lm(complete[[validation_col]], complete$mean_ADT_CLR)
    } else {
      list(slope = NA_real_, r_squared = NA_real_)
    }

    records[[i]] <- data.frame(
      group_by = key$group_by,
      feature_id = key$feature_id,
      feature_name = key$feature_name,
      target_class = key$target_class,
      validation_assay = key$validation_assay,
      validation_feature = key$validation_feature,
      validation_signal = validation_label,
      pearson_correlation = pearson,
      spearman_correlation = spearman,
      linear_model_slope = lm_stats$slope,
      linear_model_r_squared = lm_stats$r_squared,
      n_groups = n_groups,
      min_cells_per_group = min_cells,
      status = if (n_groups >= min_groups) "ok" else "insufficient_groups",
      message = if (n_groups >= min_groups) "" else paste0("Need at least ", min_groups, " complete groups"),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, records)
}

plot_feature <- function(data, summary_row, validation_col, x_label) {
  sub <- data[data$group_by == summary_row$group_by &
    data$feature_id == summary_row$feature_id &
    data$feature_name == summary_row$feature_name &
    data$validation_feature == summary_row$validation_feature, , drop = FALSE]
  complete <- sub[sub$status == "ok" &
    is.finite(sub$mean_ADT_CLR) &
    is.finite(sub[[validation_col]]), , drop = FALSE]

  title <- sprintf("%s | %s vs %s", summary_row$feature_id, summary_row$feature_name, summary_row$validation_feature)
  subtitle <- sprintf(
    "Spearman r=%s; Pearson r=%s; R2=%s; n=%s",
    format_stat(summary_row$spearman_correlation),
    format_stat(summary_row$pearson_correlation),
    format_stat(summary_row$linear_model_r_squared),
    summary_row$n_groups
  )

  p <- ggplot(complete, aes(x = .data[[validation_col]], y = mean_ADT_CLR)) +
    geom_point(aes(size = n_cells), alpha = 0.85, color = "#2C7FB8") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#D95F0E", linewidth = 0.45) +
    geom_text(aes(label = group_level), vjust = -0.7, size = 2.8, check_overlap = TRUE) +
    scale_size_continuous(name = "n cells", range = c(1.6, 4.2)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = "ADT CLR group mean"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  if (nrow(complete) < 2) {
    p <- p + annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "Insufficient complete groups",
      color = "#D62728",
      size = 4,
      fontface = "bold"
    )
  }
  p
}

write_plot_pdfs <- function(data, summary, validation_col, x_label, prefix, outdir) {
  for (group_by in unique(summary$group_by)) {
    summary_sub <- summary[summary$group_by == group_by, , drop = FALSE]
    pdf_path <- file.path(outdir, sprintf("%s_%s_consistency_scatter_plots.pdf", prefix, clean_file_part(group_by)))
    grDevices::pdf(pdf_path, width = 6.5, height = 5.2, onefile = TRUE)
    for (i in seq_len(nrow(summary_sub))) {
      print(plot_feature(data, summary_sub[i, , drop = FALSE], validation_col, x_label))
    }
    grDevices::dev.off()
    message("Wrote: ", pdf_path, " pages=", nrow(summary_sub))
  }
}

read_tsv <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

rna_path <- file.path(summary_dir, "rna_group_summary.tsv")
cd45_path <- file.path(summary_dir, "cd45isoform_group_summary.tsv")

rna <- read_tsv(rna_path, "RNA group summary")
cd45 <- read_tsv(cd45_path, "CD45 isoform group summary")

required_columns(
  rna,
  c(
    "group_by", "group_level", "n_cells", "feature_id", "feature_name", "target_class",
    "validation_assay", "validation_feature", "mean_ADT_CLR", "status", "mean_RNA_lognormalized"
  ),
  "RNA group summary"
)
required_columns(
  cd45,
  c(
    "group_by", "group_level", "n_cells", "feature_id", "feature_name", "target_class",
    "validation_assay", "validation_feature", "mean_ADT_CLR", "status", "mean_CD45isoform_signal"
  ),
  "CD45 isoform group summary"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rna_summary <- calculate_summary(rna, "mean_RNA_lognormalized", "RNA log-normalized group mean", min_groups)
cd45_summary <- calculate_summary(cd45, "mean_CD45isoform_signal", "CD45isoforms group mean", min_groups)

rna_out <- file.path(outdir, "rna_adt_correlation_summary.tsv")
cd45_out <- file.path(outdir, "cd45isoform_adt_correlation_summary.tsv")
write.table(rna_summary, rna_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cd45_summary, cd45_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", rna_out, " rows=", nrow(rna_summary))
message("Wrote: ", cd45_out, " rows=", nrow(cd45_summary))

write_plot_pdfs(
  rna,
  rna_summary,
  "mean_RNA_lognormalized",
  "RNA log-normalized group mean",
  "rna_adt",
  outdir
)
write_plot_pdfs(
  cd45,
  cd45_summary,
  "mean_CD45isoform_signal",
  "CD45isoforms group mean",
  "cd45isoform_adt",
  outdir
)
