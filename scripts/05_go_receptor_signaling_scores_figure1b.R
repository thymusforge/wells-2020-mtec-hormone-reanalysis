# 05_go_receptor_signaling_scores_figure1b.R
# Score two GO receptor signaling pathways on the Figure 1b UMAP:
#   - GOBP_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY / GO:0030520
#   - GOBP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY / GO:0030521
#
# Output:
#   output/05_go_receptor_signaling_scores/
#     05A_featureplot_<term>.png/pdf
#     05B_boxplot_<term>_by_stage.png/pdf
#     05C_panel_<term>.png/pdf
#     05D_featureplot_estrogen_androgen_combined.png/pdf
#     05E_boxplot_estrogen_androgen_by_stage.png/pdf
#     05_gene_set_summary.csv
#     05_cell_scores.csv
#     05_stage_score_summary.csv
#     05_sessionInfo.txt

rm(list = ls())
gc()

# -------------------------------------------------------------------------
# 0. Locate project root
# -------------------------------------------------------------------------

required_relative_file <- file.path("output", "control_seurat_pca_umap.rds")

get_project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in 1:5) {
    if (dir.exists(file.path(current, "data")) &&
        dir.exists(file.path(current, "scripts"))) {
      return(current)
    }
    parent <- dirname(current)
    if (parent == current) break
    current <- parent
  }
  stop("Could not find project root containing both data/ and scripts/")
}

project_root <- get_project_root()

if (!file.exists(file.path(project_root, required_relative_file))) {
  stop(
    "Expected file:\n",
    required_relative_file,
    "\n\nCurrent working directory:\n",
    getwd(),
    "\n\nPlease run 02_run_pca_umap.R first."
  )
}

output_dir   <- file.path(project_root, "output")
input_file   <- file.path(output_dir, "control_seurat_pca_umap.rds")
result_dir   <- file.path(output_dir, "05_go_receptor_signaling_scores")
result_file  <- file.path(output_dir, "control_seurat_go_receptor_signaling_scores.rds")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project root:\n", project_root, "\n\n", sep = "")
cat("Input file:\n", input_file, "\n\n", sep = "")
cat("Output directory:\n", result_dir, "\n\n", sep = "")

# -------------------------------------------------------------------------
# 1. Packages
# -------------------------------------------------------------------------

required_packages <- c("Seurat", "SeuratObject", "ggplot2", "Matrix", "msigdbr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\n\nInstall them first with:\n",
    'install.packages(c("Seurat", "SeuratObject", "ggplot2", "Matrix", "msigdbr"))'
  )
}

cat("Package versions:\n")
for (pkg in required_packages) {
  cat("  ", pkg, ": ", as.character(utils::packageVersion(pkg)), "\n", sep = "")
}
cat("\n")

use_patchwork <- requireNamespace("patchwork", quietly = TRUE)

# -------------------------------------------------------------------------
# 2. Helper functions
# -------------------------------------------------------------------------

get_data_layer <- function(obj, assay = "RNA") {
  mat <- tryCatch(
    SeuratObject::LayerData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )

  if (is.null(mat)) {
    mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = "data"),
      error = function(e) NULL
    )
  }

  if (is.null(mat)) {
    stop("Could not retrieve RNA data layer from the Seurat object.")
  }

  mat
}

scale_vector <- function(x) {
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- 0
  z
}

save_plot <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(result_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )

  ggplot2::ggsave(
    filename = file.path(result_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height
  )
}

load_msig_go_bp <- function(species = "Mus musculus") {
  msig_fun <- msigdbr::msigdbr
  msig_args <- names(formals(msig_fun))

  tryCatch(
    {
      if ("collection" %in% msig_args) {
        msigdbr::msigdbr(
          species = species,
          collection = "C5",
          subcollection = "GO:BP"
        )
      } else {
        msigdbr::msigdbr(
          species = species,
          category = "C5",
          subcategory = "GO:BP"
        )
      }
    },
    error = function(e) {
      stop("Could not load msigdbr GO:BP gene sets: ", conditionMessage(e))
    }
  )
}

get_go_term_genes <- function(msig_go_bp, term_name, go_id) {
  required_cols <- c("gs_name", "gene_symbol")
  if (!all(required_cols %in% colnames(msig_go_bp))) {
    stop("msigdbr GO:BP table does not contain expected columns: gs_name and gene_symbol.")
  }

  exact_source <- if ("gs_exact_source" %in% colnames(msig_go_bp)) {
    msig_go_bp$gs_exact_source
  } else {
    rep(NA_character_, nrow(msig_go_bp))
  }

  keep <- msig_go_bp$gs_name == term_name | exact_source == go_id
  genes <- sort(unique(msig_go_bp$gene_symbol[keep]))

  if (length(genes) == 0) {
    stop("No genes found for ", term_name, " / ", go_id, ".")
  }

  genes
}

calculate_mean_expression_score <- function(obj, data_mat, genes, score_name, min_genes = 3) {
  present_genes <- intersect(unique(genes), rownames(data_mat))

  if (length(present_genes) < min_genes) {
    stop(
      "Gene set '", score_name, "' has only ", length(present_genes),
      " genes present in the dataset."
    )
  }

  raw_score <- Matrix::colMeans(data_mat[present_genes, colnames(obj), drop = FALSE])
  z_score <- scale_vector(raw_score)

  obj[[paste0(score_name, "_raw")]] <- raw_score
  obj[[paste0(score_name, "_z")]] <- z_score

  list(
    obj = obj,
    present_genes = present_genes,
    raw_score = raw_score,
    z_score = z_score
  )
}

make_feature_plot <- function(plot_df, score_col, title, score_limits) {
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = UMAP_1, y = UMAP_2, color = .data[[score_col]])
  ) +
    ggplot2::geom_point(size = 0.42, alpha = 0.9) +
    ggplot2::scale_color_gradient(
      low = "#D9D9D9",
      high = "#D7301F",
      limits = score_limits,
      oob = scales::squish,
      name = "Score\n(z)"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

make_box_plot <- function(score_long_one, title) {
  stage_medians <- stats::aggregate(
    score_z ~ stage_plot,
    data = score_long_one,
    FUN = stats::median
  )
  stage_order <- stage_medians$stage_plot[order(stage_medians$score_z)]

  score_long_one$stage_plot <- factor(score_long_one$stage_plot, levels = stage_order)

  ggplot2::ggplot(
    score_long_one,
    ggplot2::aes(x = stage_plot, y = score_z)
  ) +
    ggplot2::geom_boxplot(
      width = 0.58,
      outlier.size = 0.25,
      linewidth = 0.25,
      fill = "#7FB3D5",
      color = "#333333"
    ) +
    ggplot2::geom_jitter(
      width = 0.18,
      size = 0.25,
      alpha = 0.35,
      color = "#1F1F1F"
    ) +
    ggplot2::labs(title = title, x = NULL, y = "Mean-expression score (z)") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

# -------------------------------------------------------------------------
# 3. Load object and GO gene sets
# -------------------------------------------------------------------------

cat("===== 1. Loading Figure 1b UMAP object =====\n")

control_obj <- readRDS(input_file)
SeuratObject::DefaultAssay(control_obj) <- "RNA"

required_meta <- c("stage_figure1b", "stage_plot", "exp")
missing_meta <- setdiff(required_meta, colnames(control_obj[[]]))

if (length(missing_meta) > 0) {
  stop(
    "Required metadata column(s) missing: ",
    paste(missing_meta, collapse = ", "),
    "\nPlease make sure this object was generated by 02_run_pca_umap.R."
  )
}

reduction_name <- if ("umap_reanalysis" %in% SeuratObject::Reductions(control_obj)) {
  "umap_reanalysis"
} else if ("umap" %in% SeuratObject::Reductions(control_obj)) {
  "umap"
} else {
  stop("No UMAP reduction found. Expected 'umap_reanalysis' from 02_run_pca_umap.R.")
}

figure1b_obj <- subset(control_obj, subset = stage_figure1b != "unknown")
figure1b_obj$stage_plot <- droplevels(figure1b_obj$stage_plot)

cat("\nCells retained for Figure 1b states:\n")
print(sort(table(figure1b_obj$stage_plot), decreasing = TRUE))
cat("\nUMAP reduction used: ", reduction_name, "\n\n", sep = "")

cat("===== 2. Loading GO:BP receptor signaling gene sets =====\n")

msig_go_bp <- load_msig_go_bp()

target_terms <- data.frame(
  score_label = c("Estrogen_receptor_signaling", "Androgen_receptor_signaling"),
  term_name = c(
    "GOBP_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY"
  ),
  go_id = c("GO:0030520", "GO:0030521"),
  plot_title = c(
    "GO:0030520 Estrogen receptor signaling",
    "GO:0030521 Androgen receptor signaling"
  ),
  stringsAsFactors = FALSE
)

go_gene_sets <- setNames(
  lapply(
    seq_len(nrow(target_terms)),
    function(i) {
      get_go_term_genes(
        msig_go_bp = msig_go_bp,
        term_name = target_terms$term_name[i],
        go_id = target_terms$go_id[i]
      )
    }
  ),
  target_terms$score_label
)

data_mat <- get_data_layer(control_obj, assay = "RNA")

gene_set_summary <- data.frame()
score_columns_z <- character(0)

for (score_label in names(go_gene_sets)) {
  score_result <- calculate_mean_expression_score(
    obj = control_obj,
    data_mat = data_mat,
    genes = go_gene_sets[[score_label]],
    score_name = score_label,
    min_genes = 3
  )

  control_obj <- score_result$obj
  score_columns_z <- c(score_columns_z, paste0(score_label, "_z"))

  gene_set_summary <- rbind(
    gene_set_summary,
    data.frame(
      score_label = score_label,
      term_name = target_terms$term_name[target_terms$score_label == score_label],
      go_id = target_terms$go_id[target_terms$score_label == score_label],
      input_genes = length(unique(go_gene_sets[[score_label]])),
      present_genes = length(score_result$present_genes),
      present_gene_symbols = paste(score_result$present_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )

  utils::write.table(
    score_result$present_genes,
    file = file.path(result_dir, paste0("05_genes_used_", score_label, ".txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

utils::write.csv(
  gene_set_summary,
  file = file.path(result_dir, "05_gene_set_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

figure1b_obj <- subset(control_obj, subset = stage_figure1b != "unknown")
figure1b_obj$stage_plot <- droplevels(figure1b_obj$stage_plot)

# -------------------------------------------------------------------------
# 4. Build plotting data
# -------------------------------------------------------------------------

cat("\n===== 3. Building plots and summaries =====\n")

emb <- as.data.frame(SeuratObject::Embeddings(figure1b_obj, reduction = reduction_name))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
emb$cell <- rownames(emb)

meta <- figure1b_obj[[]]
meta$cell <- rownames(meta)

plot_df <- merge(
  emb,
  meta[, c("cell", "stage_plot", "exp", score_columns_z), drop = FALSE],
  by = "cell",
  sort = FALSE
)

utils::write.csv(
  plot_df[, c("cell", "stage_plot", "exp", score_columns_z), drop = FALSE],
  file = file.path(result_dir, "05_cell_scores.csv"),
  row.names = FALSE,
  quote = FALSE
)

score_long <- do.call(
  rbind,
  lapply(
    seq_len(nrow(target_terms)),
    function(i) {
      score_label <- target_terms$score_label[i]
      z_col <- paste0(score_label, "_z")
      raw_col <- paste0(score_label, "_raw")
      data.frame(
        cell = rownames(meta),
        stage_plot = meta$stage_plot,
        exp = meta$exp,
        score_label = score_label,
        term_name = target_terms$term_name[i],
        go_id = target_terms$go_id[i],
        score_z = meta[[z_col]],
        score_raw = meta[[raw_col]],
        stringsAsFactors = FALSE
      )
    }
  )
)

score_long$stage_plot <- factor(score_long$stage_plot, levels = levels(figure1b_obj$stage_plot))
score_long$score_label <- factor(score_long$score_label, levels = target_terms$score_label)

stage_score_summary <- stats::aggregate(
  cbind(score_z, score_raw) ~ stage_plot + score_label + term_name + go_id,
  data = score_long,
  FUN = function(x) mean(x, na.rm = TRUE)
)

stage_score_median <- stats::aggregate(
  score_z ~ stage_plot + score_label,
  data = score_long,
  FUN = function(x) stats::median(x, na.rm = TRUE)
)
colnames(stage_score_median)[colnames(stage_score_median) == "score_z"] <- "median_score_z"

stage_score_summary <- merge(
  stage_score_summary,
  stage_score_median,
  by = c("stage_plot", "score_label"),
  all.x = TRUE,
  sort = FALSE
)

utils::write.csv(
  stage_score_summary,
  file = file.path(result_dir, "05_stage_score_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

sample_score_summary <- stats::aggregate(
  cbind(score_z, score_raw) ~ stage_plot + exp + score_label + term_name + go_id,
  data = score_long,
  FUN = function(x) mean(x, na.rm = TRUE)
)

utils::write.csv(
  sample_score_summary,
  file = file.path(result_dir, "05_sample_score_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

score_limits <- range(plot_df[, score_columns_z, drop = TRUE], finite = TRUE)

feature_plots <- list()
box_plots <- list()

for (i in seq_len(nrow(target_terms))) {
  score_label <- target_terms$score_label[i]
  z_col <- paste0(score_label, "_z")
  file_safe <- tolower(score_label)
  plot_title <- target_terms$plot_title[i]

  p_feature <- make_feature_plot(
    plot_df = plot_df,
    score_col = z_col,
    title = plot_title,
    score_limits = score_limits
  )

  p_box <- make_box_plot(
    score_long_one = score_long[score_long$score_label == score_label, , drop = FALSE],
    title = paste0(plot_title, " by Figure 1b state")
  )

  feature_plots[[score_label]] <- p_feature
  box_plots[[score_label]] <- p_box

  save_plot(
    p_feature,
    paste0("05A_featureplot_", file_safe),
    width = 5.2,
    height = 4.6
  )

  save_plot(
    p_box,
    paste0("05B_boxplot_", file_safe, "_by_stage"),
    width = 7.2,
    height = 4.8
  )

  if (use_patchwork) {
    p_panel <- p_feature / p_box + patchwork::plot_layout(heights = c(1, 1.05))
    save_plot(
      p_panel,
      paste0("05C_panel_", file_safe),
      width = 7.2,
      height = 8.2
    )
  }
}

if (use_patchwork) {
  p_feature_combined <- feature_plots[[1]] + feature_plots[[2]] +
    patchwork::plot_layout(ncol = 2, guides = "collect")
  save_plot(
    p_feature_combined,
    "05D_featureplot_estrogen_androgen_combined",
    width = 9.8,
    height = 4.8
  )
}

p_box_combined <- ggplot2::ggplot(
  score_long,
  ggplot2::aes(x = stage_plot, y = score_z)
) +
  ggplot2::geom_boxplot(
    width = 0.58,
    outlier.size = 0.25,
    linewidth = 0.25,
    fill = "#7FB3D5",
    color = "#333333"
  ) +
  ggplot2::geom_jitter(
    width = 0.18,
    size = 0.22,
    alpha = 0.32,
    color = "#1F1F1F"
  ) +
  ggplot2::facet_wrap(
    ~ score_label,
    ncol = 1,
    scales = "free_y",
    labeller = ggplot2::as_labeller(c(
      Estrogen_receptor_signaling = "GO:0030520 Estrogen receptor signaling",
      Androgen_receptor_signaling = "GO:0030521 Androgen receptor signaling"
    ))
  ) +
  ggplot2::labs(
    title = "GO receptor signaling scores across Figure 1b states",
    x = NULL,
    y = "Mean-expression score (z)"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(hjust = 0.5),
    strip.background = ggplot2::element_blank()
  )

save_plot(
  p_box_combined,
  "05E_boxplot_estrogen_androgen_by_stage",
  width = 7.4,
  height = 7.8
)

# -------------------------------------------------------------------------
# 5. Save object and session information
# -------------------------------------------------------------------------

saveRDS(control_obj, file = result_file)

summary_lines <- c(
  paste0("Input object: ", input_file),
  paste0("Output object: ", result_file),
  paste0("Output directory: ", result_dir),
  "",
  paste0("UMAP reduction used: ", reduction_name),
  paste0("Cells scored after excluding unknown Figure 1b states: ", ncol(figure1b_obj)),
  "",
  "GO terms scored:",
  paste0("  ", target_terms$go_id, " / ", target_terms$term_name),
  "",
  "Score definition:",
  "  Raw score = mean log-normalized RNA expression across present genes.",
  "  Z score = raw score scaled across all cells in the full control object.",
  "",
  "Main output files:",
  "  05A_featureplot_<term>.png/pdf",
  "  05B_boxplot_<term>_by_stage.png/pdf",
  "  05C_panel_<term>.png/pdf if patchwork is installed",
  "  05D_featureplot_estrogen_androgen_combined.png/pdf if patchwork is installed",
  "  05E_boxplot_estrogen_androgen_by_stage.png/pdf",
  "  05_gene_set_summary.csv",
  "  05_cell_scores.csv",
  "  05_stage_score_summary.csv"
)

writeLines(summary_lines, file.path(result_dir, "05_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "05_sessionInfo.txt"))

cat("\nDone.\n")
cat("Saved object:\n", result_file, "\n", sep = "")
cat("Saved plots and tables under:\n", result_dir, "\n", sep = "")
