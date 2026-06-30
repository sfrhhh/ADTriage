#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scop)
  library(Seurat)
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
    "  Rscript scripts/plot_adt_validation_features.R \\\n",
    "    --seurat-rds tests/data/test_seurat_qc.rds \\\n",
    "    --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \\\n",
    "    --outdir ADTriage_output/02_feature_plots \\\n",
    "    --n-features 10\n\n",
    "Options:\n",
    "  --seurat-rds       Input Seurat RDS.\n",
    "  --mapping          Feature mapping TSV from Step 1.\n",
    "  --outdir           Output directory.\n",
    "  --n-features       Number of mapping rows to plot. Default: 10.\n",
    "  --reduction        Reduction name or auto. Default: auto.\n",
    "  --pairs-per-row    Feature pairs per row in final montage. Default: 2.\n",
    "  --panel-px         Square panel PNG size in pixels. Default: 1000.\n",
    "  --panel-dpi        Panel PNG DPI. Default: 200.\n",
    "  --pt-size          Point size passed to FeatureDimPlot. Default: 0.5.\n",
    "  --max-page-height-px  Max montage page height in pixels. Default: 9600.\n",
    "  --magick           ImageMagick executable for montage.\n",
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
outdir <- get_arg("--outdir")
n_features <- as.integer(get_arg("--n-features", "10"))
reduction_arg <- get_arg("--reduction", "auto")
pairs_per_row <- as.integer(get_arg("--pairs-per-row", "2"))
panel_px <- as.integer(get_arg("--panel-px", "1000"))
panel_dpi <- as.integer(get_arg("--panel-dpi", "200"))
point_size <- as.numeric(get_arg("--pt-size", "0.5"))
max_page_height_px <- as.integer(get_arg("--max-page-height-px", "9600"))
magick_path <- get_arg("--magick", "/data/project/yanhuixu/frsui/software/bin/magick")

if (is.null(seurat_rds) || is.null(mapping_path) || is.null(outdir)) {
  print_help()
  stop("--seurat-rds, --mapping, and --outdir are required.", call. = FALSE)
}
if (pairs_per_row < 1) {
  stop("--pairs-per-row must be >= 1.", call. = FALSE)
}
rows_per_page <- max(1, floor(max_page_height_px / panel_px))

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

find_reduction <- function(srt, requested = "auto") {
  reductions <- Reductions(srt)
  if (length(reductions) == 0) {
    stop("No reductions found in Seurat object.", call. = FALSE)
  }
  if (!identical(tolower(requested), "auto")) {
    hit <- reductions[tolower(reductions) == tolower(requested)]
    if (length(hit) == 0) {
      stop(sprintf(
        "Requested reduction '%s' not found. Available reductions: %s",
        requested,
        paste(reductions, collapse = ", ")
      ), call. = FALSE)
    }
    return(hit[[1]])
  }

  priority <- c("standardpcaumap2d", "standardumap2d", "harmony.umap", "cca.umap", "umap")
  for (candidate in priority) {
    hit <- reductions[tolower(reductions) == candidate]
    if (length(hit) > 0) {
      return(hit[[1]])
    }
  }
  reductions[[1]]
}

clean_file_part <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

theme_panel <- function(p) {
  p +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      aspect.ratio = 1,
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.key.height = unit(0.55, "cm"),
      legend.key.width = unit(0.45, "cm"),
      plot.background = element_rect(colour = "black", fill = "white", linewidth = 0.2),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

placeholder_plot <- function(title, message) {
  message <- paste(strwrap(message, width = 24), collapse = "\n")
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message,
      color = "#D62728",
      size = 3.2,
      fontface = "bold",
      lineheight = 0.95
    ) +
    labs(title = title) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      aspect.ratio = 1,
      plot.background = element_rect(colour = "black", fill = "white", linewidth = 0.2),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

feature_exists <- function(srt, assay, feature) {
  assay %in% Assays(srt) && feature %in% rownames(srt[[assay]])
}

make_panel_plot <- function(srt, assay, feature, reduction, title, density = FALSE) {
  tryCatch({
    if (!assay %in% Assays(srt)) {
      stop(sprintf("Assay not found: %s", assay), call. = FALSE)
    }
    if (!feature %in% rownames(srt[[assay]])) {
      stop(sprintf("Feature not found in %s: %s", assay, feature), call. = FALSE)
    }
    p <- FeatureDimPlot(
      srt = srt,
      features = feature,
      reduction = paste0("^", escape_regex(reduction), "$"),
      assay = assay,
      layer = "data",
      raster = FALSE,
      pt.size = point_size,
      add_density = density,
      density_filled = density,
      density_filled_palette = "Greys",
      density_color = "grey35",
      pt.alpha = ifelse(density, 0.35, 0.8),
      aspect.ratio = 1,
      title = title,
      combine = TRUE,
      force = TRUE,
      verbose = FALSE
    )
    list(plot = theme_panel(p), status = "ok", message = "")
  }, error = function(e) {
    msg <- conditionMessage(e)
    list(
      plot = placeholder_plot(title, paste("Plot error", msg, sep = "\n")),
      status = "error",
      message = msg
    )
  })
}

left_spec_for_row <- function(row) {
  validation_assay <- row[["validation_assay"]]
  validation_feature <- row[["validation_feature"]]
  target_class <- row[["target_class"]]

  if (target_class %in% c("hashtag", "spike_control", "control") ||
      validation_assay %in% c("", "none", "manual", "HTO") ||
      validation_feature %in% c("", NA)) {
    return(list(
      assay = "",
      feature = "",
      panel_feature = "",
      title_suffix = "No validation feature",
      status = "placeholder",
      message = "No paired validation feature"
    ))
  }

  list(
    assay = validation_assay,
    feature = validation_feature,
    panel_feature = validation_feature,
    title_suffix = paste(validation_assay, validation_feature, sep = ":"),
    status = "plot",
    message = ""
  )
}

right_spec_for_row <- function(row) {
  target_class <- row[["target_class"]]
  assay <- if (identical(target_class, "hashtag")) "HTO" else "ADT"
  list(
    assay = assay,
    feature = row[["feature_name"]],
    panel_feature = row[["feature_name"]],
    title_suffix = paste(assay, row[["feature_name"]], sep = ":"),
    status = "plot",
    message = ""
  )
}

panel_title <- function(pair_index, side, title_suffix) {
  sprintf("[%02d%s] %s", pair_index, side, title_suffix)
}

save_panel <- function(plot, filename) {
  ggsave(
    filename = filename,
    plot = plot,
    width = panel_px / panel_dpi,
    height = panel_px / panel_dpi,
    dpi = panel_dpi,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )
}

make_panel_record <- function(plot_type, pair_index, side, spec, result, panel_png) {
  pair_row <- ceiling(pair_index / pairs_per_row)
  pair_col <- ((pair_index - 1) %% pairs_per_row) + 1
  page_index <- ceiling(pair_row / rows_per_page)
  page_pair_row <- ((pair_row - 1) %% rows_per_page) + 1
  side_offset <- ifelse(side == "L", 0, 1)
  grid_row <- pair_row
  grid_col <- (pair_col - 1) * 2 + side_offset + 1
  panel_index <- (grid_row - 1) * (pairs_per_row * 2) + grid_col
  page_grid_row <- page_pair_row
  page_grid_col <- grid_col
  page_panel_index <- (page_grid_row - 1) * (pairs_per_row * 2) + page_grid_col

  data.frame(
    plot_type = plot_type,
    panel_index = panel_index,
    page_index = page_index,
    page_panel_index = page_panel_index,
    pair_index = pair_index,
    pair_row = pair_row,
    pair_col = pair_col,
    page_pair_row = page_pair_row,
    panel_side = side,
    grid_row = grid_row,
    grid_col = grid_col,
    page_grid_row = page_grid_row,
    page_grid_col = page_grid_col,
    feature_id = current_row[["feature_id"]],
    feature_name = current_row[["feature_name"]],
    target_class = current_row[["target_class"]],
    assay = spec$assay,
    feature = spec$panel_feature,
    validation_assay = current_row[["validation_assay"]],
    validation_feature = current_row[["validation_feature"]],
    panel_title = panel_title(pair_index, side, spec$title_suffix),
    status = result$status,
    message = result$message,
    panel_png = panel_png,
    montage_png = "",
    stringsAsFactors = FALSE
  )
}

compose_montage <- function(panel_files, output_prefix) {
  if (!file.exists(magick_path)) {
    stop("magick executable not found: ", magick_path, call. = FALSE)
  }
  tile <- sprintf("%dx", pairs_per_row * 2)
  panels_per_page <- rows_per_page * pairs_per_row * 2
  pages <- split(panel_files, ceiling(seq_along(panel_files) / panels_per_page))
  output_files <- character()

  for (page_idx in seq_along(pages)) {
    page_suffix <- if (length(pages) == 1) {
      ".png"
    } else {
      sprintf("_page%03d.png", page_idx)
    }
    output_png <- paste0(output_prefix, page_suffix)
    status <- system2(
      magick_path,
      args = c(
        "montage",
        pages[[page_idx]],
        "-tile", tile,
        "-geometry", "+0+0",
        "-background", "white",
        output_png
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    exit_status <- attr(status, "status")
    if (!is.null(exit_status) && exit_status != 0) {
      stop("magick montage failed: ", paste(status, collapse = "\n"), call. = FALSE)
    }
    output_files <- c(output_files, output_png)
  }
  output_files
}

message("Reading Seurat object: ", seurat_rds)
bm <- readRDS(seurat_rds)

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

reduction <- find_reduction(bm, reduction_arg)
message("Using reduction: ", reduction)

mapping <- read.delim(mapping_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(mapping) == 0) {
  stop("Mapping table has no rows.", call. = FALSE)
}
selected <- head(mapping, n_features)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
panel_root <- file.path(outdir, sprintf("panel_png_top%s", nrow(selected)))
feature_panel_dir <- file.path(panel_root, "plot_feature")
dir.create(feature_panel_dir, recursive = TRUE, showWarnings = FALSE)

plot_types <- c("plot_feature")
layout_records <- list()

for (plot_type in plot_types) {
  density <- identical(plot_type, "plot_feature_density")
  panel_dir <- feature_panel_dir
  panel_files <- character()

  message("Building panel PNGs for ", plot_type)
  for (i in seq_len(nrow(selected))) {
    current_row <<- selected[i, , drop = FALSE]
    specs <- list(L = left_spec_for_row(current_row), R = right_spec_for_row(current_row))

    for (side in c("L", "R")) {
      spec <- specs[[side]]
      title <- panel_title(i, side, spec$title_suffix)
      panel_file <- file.path(
        panel_dir,
        sprintf(
          "%03d_%s_%s_%s.png",
          i,
          side,
          clean_file_part(current_row[["feature_id"]]),
          clean_file_part(spec$title_suffix)
        )
      )

      if (identical(spec$status, "placeholder")) {
        result <- list(
          plot = placeholder_plot(title, spec$message),
          status = "placeholder",
          message = spec$message
        )
      } else {
        result <- make_panel_plot(
          srt = bm,
          assay = spec$assay,
          feature = spec$panel_feature,
          reduction = reduction,
          title = title,
          density = density
        )
      }

      save_panel(result$plot, panel_file)
      panel_files <- c(panel_files, panel_file)
      layout_records[[length(layout_records) + 1]] <- make_panel_record(
        plot_type = plot_type,
        pair_index = i,
        side = side,
        spec = spec,
        result = result,
        panel_png = panel_file
      )
    }
  }

  output_prefix <- file.path(outdir, sprintf("adt_validation_%s_top%s", plot_type, nrow(selected)))
  output_pngs <- compose_montage(panel_files, output_prefix)
  for (output_png in output_pngs) {
    message("Wrote: ", output_png)
  }
  for (page_idx in seq_along(output_pngs)) {
    rows <- vapply(layout_records, function(x) {
      identical(x$plot_type, plot_type) && x$page_index == page_idx
    }, logical(1))
    for (record_idx in which(rows)) {
      layout_records[[record_idx]]$montage_png <- output_pngs[[page_idx]]
    }
  }
}

layout_index <- do.call(rbind, layout_records)
layout_path <- file.path(outdir, sprintf("adt_validation_plot_layout_index_top%s.tsv", nrow(selected)))
write.table(layout_index, layout_path, sep = "\t", quote = FALSE, row.names = FALSE)

status_path <- file.path(outdir, sprintf("adt_validation_plot_status_top%s.tsv", nrow(selected)))
status_cols <- c(
  "plot_type", "pair_index", "feature_id", "feature_name", "target_class",
  "validation_assay", "validation_feature", "panel_side", "assay", "feature",
  "status", "message", "page_index", "grid_row", "grid_col", "page_grid_row",
  "page_grid_col", "montage_png"
)
write.table(layout_index[, status_cols], status_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote: ", layout_path)
message("Wrote: ", status_path)
