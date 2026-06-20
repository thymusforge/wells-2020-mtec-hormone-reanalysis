# 02_run_pca_umap.R
# Wells et al. 2020 Figure 1b reanalysis
#
# 输入：
#   output/control_seurat_initial.rds
#
# 主要步骤：
#   1. 使用作者旧流程对应的 mean.var.plot 阈值选择高变基因
#   2. ScaleData（不回归 nUMI、线粒体比例或细胞周期）
#   3. PCA（计算 20 PCs；作者最终使用 PC 1:13）
#   4. 邻居图和聚类（k = 30；resolution = 0.6）
#   5. UMAP（PC 1:13；n.neighbors = 30；min.dist = 0.3）
#   6. 按作者 stage 注释和新聚类分别作图
#
# 重要：
# 输入对象只有已经标准化的 RNA "data" layer。
# 本脚本不会运行 NormalizeData()。

rm(list = ls())
gc()

# -------------------------------------------------------------------------
# 0. 自动寻找项目根目录
# -------------------------------------------------------------------------

required_relative_file <- file.path(
  "output",
  "control_seurat_initial.rds"
)

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
    getwd()
  )
}

output_dir  <- file.path(project_root, "output")
input_file  <- file.path(output_dir, "control_seurat_initial.rds")
result_file <- file.path(output_dir, "control_seurat_pca_umap.rds")

cat("Project root:\n", project_root, "\n\n", sep = "")

# -------------------------------------------------------------------------
# 1. 检查 packages
# -------------------------------------------------------------------------

required_packages <- c("Seurat", "SeuratObject", "ggplot2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\n\nInstall them first with:\n",
    'install.packages(c("Seurat", "SeuratObject", "ggplot2"))'
  )
}

cat("Package versions:\n")
for (pkg in required_packages) {
  cat(
    "  ", pkg, ": ",
    as.character(utils::packageVersion(pkg)),
    "\n",
    sep = ""
  )
}
cat("\n")

# -------------------------------------------------------------------------
# 2. 读取初始 Seurat 对象
# -------------------------------------------------------------------------

cat("===== 1. Loading initial Seurat object =====\n")

control_obj <- readRDS(input_file)
SeuratObject::DefaultAssay(control_obj) <- "RNA"

cat("\nInput object:\n")
print(control_obj)

if (!"data" %in% SeuratObject::Layers(control_obj[["RNA"]])) {
  stop("The RNA assay does not contain a 'data' layer.")
}

if ("counts" %in% SeuratObject::Layers(control_obj[["RNA"]])) {
  warning(
    "A counts layer is present unexpectedly. ",
    "This script will still skip NormalizeData()."
  )
}

if (!identical(dim(control_obj), c(20309L, 2434L))) {
  warning(
    "Unexpected object dimensions: ",
    paste(dim(control_obj), collapse = " x ")
  )
}

cat("\nNormalizeData() will NOT be run.\n")

# -------------------------------------------------------------------------
# 3. 设置作者 Figure 1b 的显示名称和颜色
# -------------------------------------------------------------------------

required_meta <- c("stage_figure1b", "res.0.6", "exp")

missing_meta <- setdiff(required_meta, colnames(control_obj[[]]))

if (length(missing_meta) > 0) {
  stop(
    "Required metadata column(s) missing: ",
    paste(missing_meta, collapse = ", ")
  )
}

stage_label_map <- c(
  "cTEC"          = "cTEC",
  "Ccl21a_high"   = "Ccl21a-high",
  "TAC_TEC"       = "TAC-TEC",
  "Aire_positive" = "Aire-positive",
  "Late_Aire"     = "Late-Aire",
  "Tuft"          = "Tuft",
  "unknown"       = "unknown"
)

stage_levels <- c(
  "cTEC",
  "Ccl21a-high",
  "TAC-TEC",
  "Aire-positive",
  "Late-Aire",
  "Tuft",
  "unknown"
)

stage_colors <- c(
  "cTEC"          = "#CC6600",
  "Ccl21a-high"   = "#009933",
  "TAC-TEC"       = "#0066CC",
  "Aire-positive" = "#660099",
  "Late-Aire"     = "#FF0000",
  "Tuft"          = "#990000",
  "unknown"       = "#D9D9D9"
)

stage_values <- unname(
  stage_label_map[as.character(control_obj$stage_figure1b)]
)

if (anyNA(stage_values)) {
  bad_values <- unique(
    as.character(control_obj$stage_figure1b)[is.na(stage_values)]
  )

  stop(
    "Unrecognized stage_figure1b value(s): ",
    paste(bad_values, collapse = ", ")
  )
}

control_obj$stage_plot <- factor(
  stage_values,
  levels = stage_levels
)

# -------------------------------------------------------------------------
# 4. 选择高变基因
#
# 作者旧版 Seurat：
#   FindVariableGenes(
#     x.low.cutoff = 0.0125,
#     x.high.cutoff = 3,
#     y.cutoff = 0.5
#   )
#
# 当前 Seurat 的近似对应：
#   selection.method = "mean.var.plot"
#   mean.cutoff = c(0.0125, 3)
#   dispersion.cutoff = c(0.5, Inf)
# -------------------------------------------------------------------------

cat("\n===== 2. Finding variable features =====\n")

control_obj <- Seurat::FindVariableFeatures(
  object = control_obj,
  assay = "RNA",
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.0125, 3),
  dispersion.cutoff = c(0.5, Inf),
  verbose = TRUE
)

variable_features <- SeuratObject::VariableFeatures(control_obj)

cat(
  "\nNumber of variable features selected: ",
  length(variable_features),
  "\n",
  sep = ""
)

if (length(variable_features) < 50) {
  stop(
    "Too few variable features were selected (",
    length(variable_features),
    "). Stop before PCA."
  )
}

utils::write.table(
  variable_features,
  file = file.path(output_dir, "02_variable_features.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# 保存高变基因图
variable_plot <- Seurat::VariableFeaturePlot(control_obj) +
  ggplot2::ggtitle(
    paste0(
      "Variable features: mean.var.plot (n = ",
      length(variable_features),
      ")"
    )
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  filename = file.path(output_dir, "02_variable_features.png"),
  plot = variable_plot,
  width = 7,
  height = 5,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir, "02_variable_features.pdf"),
  plot = variable_plot,
  width = 7,
  height = 5
)

# -------------------------------------------------------------------------
# 5. ScaleData
#
# 只缩放高变基因。
# 不回归 nUMI、percent_mito、exp 或 cell cycle。
# -------------------------------------------------------------------------

cat("\n===== 3. Scaling variable features =====\n")

control_obj <- Seurat::ScaleData(
  object = control_obj,
  assay = "RNA",
  features = variable_features,
  do.center = TRUE,
  do.scale = TRUE,
  verbose = TRUE
)

# -------------------------------------------------------------------------
# 6. PCA
#
# 作者 RunPCA 的 seed.use = 42。
# 作者检查 PC 1:20，最终使用 PC 1:13。
# -------------------------------------------------------------------------

cat("\n===== 4. Running PCA =====\n")

control_obj <- Seurat::RunPCA(
  object = control_obj,
  assay = "RNA",
  features = variable_features,
  npcs = 20,
  seed.use = 42,
  approx = TRUE,
  verbose = TRUE
)

if (!"pca" %in% SeuratObject::Reductions(control_obj)) {
  stop("PCA reduction was not created.")
}

pca_embeddings <- SeuratObject::Embeddings(
  control_obj,
  reduction = "pca"
)

cat(
  "\nPCA dimensions: ",
  nrow(pca_embeddings), " cells x ",
  ncol(pca_embeddings), " PCs\n",
  sep = ""
)

# 保存 elbow plot
elbow_plot <- Seurat::ElbowPlot(
  object = control_obj,
  ndims = 20
) +
  ggplot2::geom_vline(
    xintercept = 13,
    linetype = "dashed"
  ) +
  ggplot2::ggtitle("PCA elbow plot; PC 1:13 used downstream") +
  ggplot2::theme_classic()

ggplot2::ggsave(
  filename = file.path(output_dir, "02_pca_elbow.png"),
  plot = elbow_plot,
  width = 7,
  height = 5,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir, "02_pca_elbow.pdf"),
  plot = elbow_plot,
  width = 7,
  height = 5
)

# -------------------------------------------------------------------------
# 7. 构建邻居图
#
# 作者 Seurat 2.3.4 的 FindClusters 默认：
#   k.param = 30
#   prune.SNN = 1/15
#   nn.eps = 0
#
# 当前 Seurat 把该过程拆为 FindNeighbors + FindClusters。
# 使用 RANN 精确邻居搜索，比当前默认 Annoy 更接近旧流程。
# -------------------------------------------------------------------------

cat("\n===== 5. Building nearest-neighbor and SNN graphs =====\n")

control_obj <- Seurat::FindNeighbors(
  object = control_obj,
  reduction = "pca",
  dims = 1:13,
  k.param = 30,
  compute.SNN = TRUE,
  prune.SNN = 1 / 15,
  nn.method = "rann",
  nn.eps = 0,
  verbose = TRUE
)

graph_names <- SeuratObject::Graphs(control_obj)

cat("\nGraphs stored in object:\n")
print(graph_names)

snn_graphs <- grep("_snn$", graph_names, value = TRUE)

if (length(snn_graphs) == 0) {
  stop("No SNN graph was found after FindNeighbors().")
}

snn_graph <- snn_graphs[1]
cat("SNN graph used for clustering: ", snn_graph, "\n", sep = "")

# -------------------------------------------------------------------------
# 8. 聚类
#
# 作者参数：
#   resolution = 0.6
#   algorithm = 1
#   n.start = 100
#   n.iter = 10
#   random.seed = 0
# -------------------------------------------------------------------------

cat("\n===== 6. Clustering cells =====\n")

control_obj <- Seurat::FindClusters(
  object = control_obj,
  graph.name = snn_graph,
  resolution = 0.6,
  algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  verbose = TRUE
)

control_obj$recluster_res0.6 <- as.character(
  SeuratObject::Idents(control_obj)
)

cat("\nNew cluster sizes:\n")
print(
  sort(
    table(control_obj$recluster_res0.6),
    decreasing = TRUE
  )
)

# 比较作者原 cluster 和新 cluster
author_cluster <- as.character(
  control_obj[[]][["res.0.6"]]
)

new_cluster <- as.character(
  control_obj$recluster_res0.6
)

comparison_table <- table(
  author_cluster = author_cluster,
  new_cluster = new_cluster
)

utils::write.csv(
  as.data.frame.matrix(comparison_table),
  file = file.path(
    output_dir,
    "02_author_vs_reclustered_clusters.csv"
  ),
  quote = FALSE
)

stage_cluster_table <- table(
  author_stage = control_obj$stage_plot,
  new_cluster = new_cluster
)

utils::write.csv(
  as.data.frame.matrix(stage_cluster_table),
  file = file.path(
    output_dir,
    "02_author_stage_vs_reclustered_clusters.csv"
  ),
  quote = FALSE
)

# -------------------------------------------------------------------------
# 9. UMAP
#
# 作者 Seurat 2.3.4 默认：
#   n_neighbors = 30
#   min_dist = 0.3
#   metric = "correlation"
#   seed.use 默认 42
#
# 原论文使用 Python umap-learn；这里用当前 Seurat 的 R-native uwot，
# 但显式保留上述数值参数。图形会接近，但不会逐点相同。
# -------------------------------------------------------------------------

cat("\n===== 7. Running UMAP =====\n")

control_obj <- Seurat::RunUMAP(
  object = control_obj,
  reduction = "pca",
  dims = 1:13,
  umap.method = "uwot",
  n.neighbors = 30L,
  min.dist = 0.3,
  metric = "correlation",
  seed.use = 42L,
  reduction.name = "umap_reanalysis",
  reduction.key = "UMAPR_",
  verbose = TRUE
)

if (!"umap_reanalysis" %in% SeuratObject::Reductions(control_obj)) {
  stop("UMAP reduction was not created.")
}

umap_embeddings <- SeuratObject::Embeddings(
  control_obj,
  reduction = "umap_reanalysis"
)

cat(
  "\nUMAP dimensions: ",
  nrow(umap_embeddings), " cells x ",
  ncol(umap_embeddings), " dimensions\n",
  sep = ""
)

# -------------------------------------------------------------------------
# 10. 绘图
# -------------------------------------------------------------------------

cat("\n===== 8. Creating UMAP plots =====\n")

# Figure 1b 最终绘图脚本删除 unknown
figure1b_obj <- subset(
  control_obj,
  subset = stage_figure1b != "unknown"
)

# 按作者细胞类型着色，带图例
plot_stage_legend <- Seurat::DimPlot(
  object = figure1b_obj,
  reduction = "umap_reanalysis",
  group.by = "stage_plot",
  cols = stage_colors[stage_levels[stage_levels != "unknown"]],
  pt.size = 0.75,
  shuffle = TRUE,
  seed = 42,
  raster = FALSE
) +
  ggplot2::labs(
    title = "Figure 1b reanalysis",
    color = NULL
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

# 按作者细胞类型着色，不带图例
plot_stage_no_legend <- plot_stage_legend +
  Seurat::NoLegend()

# 按我们重新计算的 cluster 着色
plot_new_clusters <- Seurat::DimPlot(
  object = control_obj,
  reduction = "umap_reanalysis",
  group.by = "recluster_res0.6",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.75,
  shuffle = TRUE,
  seed = 42,
  raster = FALSE
) +
  ggplot2::labs(
    title = "Reclustered cells: resolution 0.6",
    color = "Cluster"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

# 按两个control样本着色，用于检查batch mixing
plot_experiment <- Seurat::DimPlot(
  object = control_obj,
  reduction = "umap_reanalysis",
  group.by = "exp",
  pt.size = 0.75,
  shuffle = TRUE,
  seed = 42,
  raster = FALSE
) +
  ggplot2::labs(
    title = "Control sample origin",
    color = NULL
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

# 保存 Figure 1b-style 图
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "02_figure1b_reanalysis_with_legend.png"
  ),
  plot = plot_stage_legend,
  width = 8,
  height = 6,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "02_figure1b_reanalysis_with_legend.pdf"
  ),
  plot = plot_stage_legend,
  width = 8,
  height = 6
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "02_figure1b_reanalysis_no_legend.png"
  ),
  plot = plot_stage_no_legend,
  width = 7,
  height = 6,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "02_figure1b_reanalysis_no_legend.pdf"
  ),
  plot = plot_stage_no_legend,
  width = 7,
  height = 6
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "02_umap_reclustered.png"
  ),
  plot = plot_new_clusters,
  width = 8,
  height = 6,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "02_umap_by_experiment.png"
  ),
  plot = plot_experiment,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 11. 保存分析后的对象和运行摘要
# -------------------------------------------------------------------------

cat("\n===== 9. Saving results =====\n")

saveRDS(
  control_obj,
  file = result_file
)

summary_file <- file.path(
  output_dir,
  "02_pca_umap_summary.txt"
)

summary_lines <- c(
  paste0("Input object: ", input_file),
  paste0("Output object: ", result_file),
  "",
  paste0("Genes: ", nrow(control_obj)),
  paste0("Cells: ", ncol(control_obj)),
  paste0(
    "Variable features: ",
    length(SeuratObject::VariableFeatures(control_obj))
  ),
  "NormalizeData: NOT RUN",
  "ScaleData features: variable features only",
  "PCA seed: 42",
  "PCA computed: 20",
  "PCs used downstream: 1:13",
  "FindNeighbors k.param: 30",
  "FindNeighbors nn.method: rann",
  "FindNeighbors prune.SNN: 1/15",
  "FindClusters resolution: 0.6",
  "FindClusters algorithm: 1",
  "FindClusters n.start: 100",
  "FindClusters n.iter: 10",
  "FindClusters random.seed: 0",
  "UMAP implementation: uwot",
  "UMAP n.neighbors: 30",
  "UMAP min.dist: 0.3",
  "UMAP metric: correlation",
  "UMAP seed: 42",
  "",
  "New cluster sizes:",
  paste(
    capture.output(
      print(
        sort(
          table(control_obj$recluster_res0.6),
          decreasing = TRUE
        )
      )
    ),
    collapse = "\n"
  ),
  "",
  "Author stage counts (unknown retained in object):",
  paste(
    capture.output(
      print(
        sort(
          table(control_obj$stage_plot),
          decreasing = TRUE
        )
      )
    ),
    collapse = "\n"
  )
)

writeLines(summary_lines, summary_file)

session_file <- file.path(
  output_dir,
  "02_sessionInfo.txt"
)

writeLines(
  capture.output(sessionInfo()),
  session_file
)

cat("\nDone.\n")
cat("Saved processed Seurat object:\n", result_file, "\n", sep = "")
cat("Saved Figure 1b-style plot:\n",
    file.path(output_dir, "02_figure1b_reanalysis_with_legend.png"),
    "\n",
    sep = "")
cat("Saved PCA elbow plot:\n",
    file.path(output_dir, "02_pca_elbow.png"),
    "\n",
    sep = "")
cat("Saved run summary:\n", summary_file, "\n", sep = "")
