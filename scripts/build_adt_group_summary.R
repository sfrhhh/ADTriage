#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
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
    "  Rscript scripts/build_adt_group_summary.R \\\n",
    "    --seurat-rds tests/data/test_seurat_qc.rds \\\n",
    "    --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \\\n",
    "    --outdir ADTriage_output/03_group_summary \\\n",
    "    --group-by Adult_cHSPCs_Ultima_majority_voting,Standardpcaclusters\n\n",
    "Options:\n",
    "  --seurat-rds  Input Seurat RDS.\n",
    "  --mapping     Step 1 llm_review mapping TSV.\n",
    "  --outdir      Output directory.\n",
    "  --group-by    Comma-separated metadata columns.\n",
    "  --help        Show this help.\n",
    sep = ""
  )
}

if (has_flag("--help")) {
  print_help()
  quit(save = "no", status = 0)
}

seurat_rds <- get_arg("--seurat-rds")
mapping_path <- get_arg("--mapping")
outdir <- get_arg("--outdir")
group_by_arg <- get_arg("--group-by")

if (is.null(seurat_rds) || is.null(mapping_path) || is.null(outdir) || is.null(group_by_arg)) {
  print_help()
  stop("--seurat-rds, --mapping, --outdir, and --group-by are required.", call. = FALSE)
}

group_by_cols <- trimws(strsplit(group_by_arg, ",", fixed = TRUE)[[1]])
group_by_cols <- group_by_cols[nzchar(group_by_cols)]
if (length(group_by_cols) == 0) {
  stop("--group-by did not contain any metadata columns.", call. = FALSE)
}

read_mapping <- function(path) {
  mapping <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c(
    "feature_id", "feature_name", "target_class", "human_gene_symbol",
    "validation_assay", "validation_feature", "mapping_method", "mapping_confidence"
  )
  missing <- setdiff(required, colnames(mapping))
  if (length(missing) > 0) {
    stop("Mapping table missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  mapping
}

get_layer_data <- function(srt, assay, layer = "data") {
  tryCatch(
    GetAssayData(srt, assay = assay, layer = layer),
    error = function(e) {
      GetAssayData(srt, assay = assay, slot = layer)
    }
  )
}

ordered_group_levels <- function(values) {
  if (is.factor(values)) {
    levels(values)
  } else {
    unique(as.character(values))
  }
}

make_group_info <- function(meta, group_by) {
  values <- meta[[group_by]]
  values_chr <- as.character(values)
  values_chr[is.na(values_chr)] <- "NA"
  levels_use <- ordered_group_levels(values)
  levels_use <- as.character(levels_use)
  if (any(is.na(values))) {
    levels_use <- unique(c(levels_use, "NA"))
  }
  levels_use <- levels_use[levels_use %in% values_chr]

  do.call(
    rbind,
    lapply(levels_use, function(level) {
      cells <- rownames(meta)[values_chr == level]
      data.frame(
        group_by = group_by,
        group_level = level,
        n_cells = length(cells),
        stringsAsFactors = FALSE
      )
    })
  )
}

mean_matrix_by_group <- function(mat, meta, group_by, features) {
  features <- unique(features[nzchar(features)])
  present <- intersect(features, rownames(mat))
  group_info <- make_group_info(meta, group_by)
  result <- list(group_info = group_info, means = list(), present = present)

  for (i in seq_len(nrow(group_info))) {
    level <- group_info$group_level[[i]]
    values_chr <- as.character(meta[[group_by]])
    values_chr[is.na(values_chr)] <- "NA"
    cells <- rownames(meta)[values_chr == level]
    if (length(present) == 0 || length(cells) == 0) {
      result$means[[level]] <- numeric(0)
    } else {
      result$means[[level]] <- Matrix::rowMeans(mat[present, cells, drop = FALSE])
    }
  }
  result
}

primary_assay_for_row <- function(target_class) {
  if (identical(target_class, "hashtag")) {
    "HTO"
  } else {
    "ADT"
  }
}

build_adt_summary <- function(mapping, grouped_means) {
  records <- list()
  for (group_by in names(grouped_means)) {
    group_info <- grouped_means[[group_by]]$group_info
    for (i in seq_len(nrow(mapping))) {
      row <- mapping[i, , drop = FALSE]
      assay <- primary_assay_for_row(row$target_class)
      means_obj <- grouped_means[[group_by]][[assay]]
      feature <- row$feature_name
      present <- feature %in% means_obj$present
      status <- if (present) "ok" else "missing_feature"
      message <- if (present) "" else paste0("Feature not found in ", assay, ": ", feature)

      for (j in seq_len(nrow(group_info))) {
        level <- group_info$group_level[[j]]
        value <- if (present) unname(means_obj$means[[level]][[feature]]) else NA_real_
        records[[length(records) + 1]] <- data.frame(
          group_by = group_by,
          group_level = level,
          n_cells = group_info$n_cells[[j]],
          feature_id = row$feature_id,
          feature_name = row$feature_name,
          target_class = row$target_class,
          assay = assay,
          feature = feature,
          mean_CLR = value,
          mean_ADT_CLR = if (identical(assay, "ADT")) value else NA_real_,
          mean_HTO_CLR = if (identical(assay, "HTO")) value else NA_real_,
          status = status,
          message = message,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, records)
}

build_validation_summary <- function(mapping, grouped_means, validation_assay, mean_col) {
  rows <- mapping[mapping$validation_assay == validation_assay & nzchar(mapping$validation_feature), , drop = FALSE]
  records <- list()
  for (group_by in names(grouped_means)) {
    group_info <- grouped_means[[group_by]]$group_info
    validation_obj <- grouped_means[[group_by]][[validation_assay]]
    adt_obj <- grouped_means[[group_by]][["ADT"]]

    for (i in seq_len(nrow(rows))) {
      row <- rows[i, , drop = FALSE]
      validation_feature <- row$validation_feature
      adt_feature <- row$feature_name
      validation_present <- validation_feature %in% validation_obj$present
      adt_present <- adt_feature %in% adt_obj$present
      status <- if (validation_present && adt_present) {
        "ok"
      } else if (!validation_present && !adt_present) {
        "missing_adt_and_validation_feature"
      } else if (!validation_present) {
        "missing_validation_feature"
      } else {
        "missing_adt_feature"
      }
      message <- paste(
        c(
          if (!adt_present) paste0("Feature not found in ADT: ", adt_feature),
          if (!validation_present) paste0("Feature not found in ", validation_assay, ": ", validation_feature)
        ),
        collapse = "; "
      )

      for (j in seq_len(nrow(group_info))) {
        level <- group_info$group_level[[j]]
        adt_value <- if (adt_present) unname(adt_obj$means[[level]][[adt_feature]]) else NA_real_
        validation_value <- if (validation_present) {
          unname(validation_obj$means[[level]][[validation_feature]])
        } else {
          NA_real_
        }
        record <- data.frame(
          group_by = group_by,
          group_level = level,
          n_cells = group_info$n_cells[[j]],
          feature_id = row$feature_id,
          feature_name = row$feature_name,
          target_class = row$target_class,
          validation_assay = row$validation_assay,
          validation_feature = row$validation_feature,
          mean_ADT_CLR = adt_value,
          status = status,
          message = message,
          stringsAsFactors = FALSE
        )
        record[[mean_col]] <- validation_value
        records[[length(records) + 1]] <- record
      }
    }
  }
  if (length(records) == 0) {
    return(data.frame())
  }
  do.call(rbind, records)
}

message("Reading Seurat object: ", seurat_rds)
bm <- readRDS(seurat_rds)

missing_group_by <- setdiff(group_by_cols, colnames(bm@meta.data))
if (length(missing_group_by) > 0) {
  stop("Missing group_by metadata columns: ", paste(missing_group_by, collapse = ", "), call. = FALSE)
}

required_assays <- c("RNA", "ADT", "CD45isoforms")
missing_assays <- setdiff(required_assays, Assays(bm))
if (length(missing_assays) > 0) {
  stop("Missing required assays: ", paste(missing_assays, collapse = ", "), call. = FALSE)
}

message("Normalizing ADT assay with CLR margin=2")
bm <- NormalizeData(bm, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE)

message("Normalizing CD45isoforms assay with CLR")
bm <- NormalizeData(bm, assay = "CD45isoforms", normalization.method = "CLR", verbose = FALSE)

if ("HTO" %in% Assays(bm)) {
  message("Normalizing HTO assay with CLR")
  bm <- NormalizeData(bm, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
}

mapping <- read_mapping(mapping_path)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

adt_features <- unique(mapping$feature_name)
hto_features <- unique(mapping$feature_name[mapping$target_class == "hashtag"])
rna_features <- unique(mapping$validation_feature[mapping$validation_assay == "RNA"])
cd45_features <- unique(mapping$validation_feature[mapping$validation_assay == "CD45isoforms"])

assay_mats <- list(
  ADT = get_layer_data(bm, "ADT", "data"),
  RNA = get_layer_data(bm, "RNA", "data"),
  CD45isoforms = get_layer_data(bm, "CD45isoforms", "data")
)
if ("HTO" %in% Assays(bm)) {
  assay_mats[["HTO"]] <- get_layer_data(bm, "HTO", "data")
}

grouped_means <- list()
for (group_by in group_by_cols) {
  message("Summarizing group_by: ", group_by)
  grouped_means[[group_by]] <- list(
    group_info = make_group_info(bm@meta.data, group_by),
    ADT = mean_matrix_by_group(assay_mats[["ADT"]], bm@meta.data, group_by, adt_features),
    RNA = mean_matrix_by_group(assay_mats[["RNA"]], bm@meta.data, group_by, rna_features),
    CD45isoforms = mean_matrix_by_group(assay_mats[["CD45isoforms"]], bm@meta.data, group_by, cd45_features)
  )
  if ("HTO" %in% names(assay_mats)) {
    grouped_means[[group_by]][["HTO"]] <- mean_matrix_by_group(
      assay_mats[["HTO"]], bm@meta.data, group_by, hto_features
    )
  } else {
    grouped_means[[group_by]][["HTO"]] <- list(
      group_info = grouped_means[[group_by]]$group_info,
      means = list(),
      present = character()
    )
  }
}

adt_summary <- build_adt_summary(mapping, grouped_means)
rna_summary <- build_validation_summary(mapping, grouped_means, "RNA", "mean_RNA_lognormalized")
cd45_summary <- build_validation_summary(
  mapping, grouped_means, "CD45isoforms", "mean_CD45isoform_signal"
)

adt_path <- file.path(outdir, "adt_group_summary.tsv")
rna_path <- file.path(outdir, "rna_group_summary.tsv")
cd45_path <- file.path(outdir, "cd45isoform_group_summary.tsv")

write.table(adt_summary, adt_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rna_summary, rna_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cd45_summary, cd45_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote: ", adt_path, " rows=", nrow(adt_summary))
message("Wrote: ", rna_path, " rows=", nrow(rna_summary))
message("Wrote: ", cd45_path, " rows=", nrow(cd45_summary))
