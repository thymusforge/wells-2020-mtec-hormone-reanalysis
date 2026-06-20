# 05_specific_genes_dotplot_go_receptor_signaling.R
# DotPlot individual genes from two GO receptor signaling terms across Figure 1b states:
#   - GOBP_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY / GO:0030520
#   - GOBP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY / GO:0030521
#
# Output:
#   output/05_go_receptor_signaling_scores/
#     05F_dotplot_<term>_genes_by_stage.png/pdf
#     05G_dotplot_shared_genes_by_stage.png/pdf
#     05H_dotplot_term_unique_genes_by_stage.png/pdf
#     05_specific_gene_expression_summary.csv
#     05_specific_gene_membership.csv
#     05_specific_genes_sessionInfo.txt

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

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
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
  sort(unique(msig_go_bp$gene_symbol[keep]))
}

make_dotplot <- function(obj, genes, title, max_genes_per_plot = 70) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) {
    stop("No requested genes are present in the Seurat object for plot: ", title)
  }

  if (length(genes) > max_genes_per_plot) {
    warning(
      "Plot '", title, "' has ", length(genes), " genes. ",
      "Only the first ", max_genes_per_plot, " genes will be plotted."
    )
    genes <- genes[seq_len(max_genes_per_plot)]
  }

  Seurat::DotPlot(
    object = obj,
    features = genes,
    group.by = "stage_plot",
    assay = "RNA",
    dot.scale = 5.2
  ) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Average\nexpression\n(scaled)"
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL,
      size = "% cells\nexpressing"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    )
}

# -------------------------------------------------------------------------
# 3. Load object and GO genes
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

figure1b_obj <- subset(control_obj, subset = stage_figure1b != "unknown")
figure1b_obj$stage_plot <- droplevels(figure1b_obj$stage_plot)

cat("\nCells retained for Figure 1b states:\n")
print(sort(table(figure1b_obj$stage_plot), decreasing = TRUE))

target_terms <- data.frame(
  term_group = c("Estrogen_receptor_signaling", "Androgen_receptor_signaling"),
  term_name = c(
    "GOBP_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY"
  ),
  go_id = c("GO:0030520", "GO:0030521"),
  plot_title = c(
    "GO:0030520 estrogen receptor signaling genes",
    "GO:0030521 androgen receptor signaling genes"
  ),
  stringsAsFactors = FALSE
)

cat("\n===== 2. Loading GO:BP term genes =====\n")

msig_go_bp <- load_msig_go_bp()

term_genes <- setNames(
  lapply(
    seq_len(nrow(target_terms)),
    function(i) {
      genes <- get_go_term_genes(
        msig_go_bp = msig_go_bp,
        term_name = target_terms$term_name[i],
        go_id = target_terms$go_id[i]
      )
      intersect(genes, rownames(figure1b_obj))
    }
  ),
  target_terms$term_group
)

shared_genes <- Reduce(intersect, term_genes)
unique_genes <- list(
  Estrogen_unique = setdiff(term_genes$Estrogen_receptor_signaling, term_genes$Androgen_receptor_signaling),
  Androgen_unique = setdiff(term_genes$Androgen_receptor_signaling, term_genes$Estrogen_receptor_signaling)
)

membership <- do.call(
  rbind,
  lapply(
    names(term_genes),
    function(term_group) {
      data.frame(
        term_group = term_group,
        term_name = target_terms$term_name[target_terms$term_group == term_group],
        go_id = target_terms$go_id[target_terms$term_group == term_group],
        gene = term_genes[[term_group]],
        also_in_other_term = term_genes[[term_group]] %in% shared_genes,
        stringsAsFactors = FALSE
      )
    }
  )
)

utils::write.csv(
  membership,
  file = file.path(result_dir, "05_specific_gene_membership.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("Present genes per GO term:\n")
print(table(membership$term_group))
cat("\nShared genes:\n")
print(shared_genes)
cat("\n")

# -------------------------------------------------------------------------
# 4. Expression summaries
# -------------------------------------------------------------------------

cat("===== 3. Calculating gene expression summaries =====\n")

data_mat <- get_data_layer(figure1b_obj, assay = "RNA")
all_genes <- sort(unique(membership$gene))
stage_levels <- levels(figure1b_obj$stage_plot)
stage_vector <- as.character(figure1b_obj$stage_plot)

summary_df <- do.call(
  rbind,
  lapply(
    all_genes,
    function(gene) {
      do.call(
        rbind,
        lapply(
          stage_levels,
          function(stage) {
            cells <- colnames(figure1b_obj)[stage_vector == stage]
            expr <- as.numeric(data_mat[gene, cells, drop = TRUE])

            data.frame(
              gene = gene,
              stage_plot = stage,
              average_expression = mean(expr),
              percent_expressing = mean(expr > 0) * 100,
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

summary_df <- merge(
  summary_df,
  unique(membership[, c("gene", "term_group", "go_id", "also_in_other_term")]),
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

utils::write.csv(
  summary_df,
  file = file.path(result_dir, "05_specific_gene_expression_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -------------------------------------------------------------------------
# 5. DotPlots
# -------------------------------------------------------------------------

cat("===== 4. Plotting gene dotplots =====\n")

for (i in seq_len(nrow(target_terms))) {
  term_group <- target_terms$term_group[i]
  file_safe <- tolower(term_group)
  genes <- term_genes[[term_group]]

  p <- make_dotplot(
    obj = figure1b_obj,
    genes = genes,
    title = target_terms$plot_title[i],
    max_genes_per_plot = 70
  )

  height <- max(5.6, 0.18 * length(genes) + 2.2)
  save_plot(
    p,
    paste0("05F_dotplot_", file_safe, "_genes_by_stage"),
    width = 8.8,
    height = height
  )
}

if (length(shared_genes) > 0) {
  p_shared <- make_dotplot(
    obj = figure1b_obj,
    genes = shared_genes,
    title = "Genes shared by estrogen and androgen receptor signaling GO terms",
    max_genes_per_plot = 70
  )

  save_plot(
    p_shared,
    "05G_dotplot_shared_genes_by_stage",
    width = 7.8,
    height = max(4.8, 0.2 * length(shared_genes) + 2.2)
  )
}

for (nm in names(unique_genes)) {
  genes <- unique_genes[[nm]]
  if (length(genes) == 0) {
    next
  }

  p_unique <- make_dotplot(
    obj = figure1b_obj,
    genes = genes,
    title = paste0(gsub("_", " ", nm), " GO term genes"),
    max_genes_per_plot = 70
  )

  save_plot(
    p_unique,
    paste0("05H_dotplot_", tolower(nm), "_genes_by_stage"),
    width = 8.8,
    height = max(5.6, 0.18 * length(genes) + 2.2)
  )
}

# -------------------------------------------------------------------------
# 6. Save session information
# -------------------------------------------------------------------------

summary_lines <- c(
  paste0("Input object: ", input_file),
  paste0("Output directory: ", result_dir),
  "",
  "Purpose:",
  "  DotPlot individual genes from GO:0030520 and GO:0030521 across Figure 1b states.",
  "",
  "GO terms:",
  paste0("  ", target_terms$go_id, " / ", target_terms$term_name),
  "",
  paste0("Estrogen receptor signaling genes present: ", length(term_genes$Estrogen_receptor_signaling)),
  paste0("Androgen receptor signaling genes present: ", length(term_genes$Androgen_receptor_signaling)),
  paste0("Shared genes present: ", length(shared_genes)),
  "",
  "Main output files:",
  "  05F_dotplot_<term>_genes_by_stage.png/pdf",
  "  05G_dotplot_shared_genes_by_stage.png/pdf",
  "  05H_dotplot_<term_unique>_genes_by_stage.png/pdf",
  "  05_specific_gene_membership.csv",
  "  05_specific_gene_expression_summary.csv"
)

writeLines(summary_lines, file.path(result_dir, "05_specific_genes_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "05_specific_genes_sessionInfo.txt"))

cat("\nDone.\n")
cat("Saved plots and tables under:\n", result_dir, "\n", sep = "")
