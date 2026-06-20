# 01_create_seurat_object.R
# Wells et al. 2020 Figure 1b reanalysis
#
# 目的：
# 1. 读取作者从 control Seurat 对象导出的 log-normalized expression matrix
# 2. 读取并匹配 controls_meta.csv
# 3. 创建一个只有 RNA "data" layer、没有 raw counts layer 的 Seurat v5 对象
# 4. 保存为 output/control_seurat_initial.rds
#
# 重要：
# GSE137699_combinedControl.csv 已经是标准化表达矩阵。
# 本脚本不会运行 NormalizeData()，避免二次标准化。
rm(list = ls())
gc()

# -------------------------------------------------------------------------
# 0. 自动寻找项目根目录
# -------------------------------------------------------------------------

required_relative_files <- c(
  file.path("data", "GSE137699_combinedControl.csv"),
  file.path("data", "controls_meta.csv")
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
missing_required <- required_relative_files[
  !file.exists(file.path(project_root, required_relative_files))
]

if (length(missing_required) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_required, collapse = "\n"),
    "\n\nSee README.md Data source and data/README.md for data preparation instructions."
  )
}

data_dir   <- file.path(project_root, "data")
output_dir <- file.path(project_root, "output")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

control_file <- file.path(data_dir, "GSE137699_combinedControl.csv")
meta_file    <- file.path(data_dir, "controls_meta.csv")
output_file  <- file.path(output_dir, "control_seurat_initial.rds")

cat("Project root:\n", project_root, "\n\n", sep = "")

# -------------------------------------------------------------------------
# 1. 检查所需 R packages
# -------------------------------------------------------------------------

required_packages <- c("data.table", "Matrix", "SeuratObject")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\n\nInstall them first with:\n",
    'install.packages(c("data.table", "Matrix", "SeuratObject"))'
  )
}

cat("Package versions:\n")
for (pkg in required_packages) {
  cat("  ", pkg, ": ", as.character(utils::packageVersion(pkg)), "\n", sep = "")
}
cat("\n")

# -------------------------------------------------------------------------
# 2. 读取 metadata
# -------------------------------------------------------------------------

cat("===== 1. Reading controls metadata =====\n")

meta <- utils::read.csv(
  meta_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (anyDuplicated(rownames(meta))) {
  stop("Duplicate cell barcodes were found in controls_meta.csv.")
}

cat("Metadata dimensions: ",
    nrow(meta), " cells x ", ncol(meta), " columns\n", sep = "")

cat("Metadata columns:\n")
print(colnames(meta))

# 保存原始 stage，并增加论文 Figure 1b 使用的名称
if (!"stage" %in% colnames(meta)) {
  stop("The metadata does not contain a 'stage' column.")
}

meta$stage_original <- meta$stage
meta$stage_figure1b <- meta$stage

meta$stage_figure1b[
  meta$stage_figure1b == "Early_Aire"
] <- "TAC_TEC"

meta$stage_figure1b[
  meta$stage_figure1b == "Cortico_medullary"
] <- "cTEC"

cat("\nFigure 1b stage counts (unknown retained for now):\n")
print(sort(table(meta$stage_figure1b), decreasing = TRUE))

# -------------------------------------------------------------------------
# 3. 读取第一行，取得 2,434 个 cell barcodes
# -------------------------------------------------------------------------

cat("\n===== 2. Reading expression-matrix header =====\n")

header_line <- readLines(
  control_file,
  n = 1,
  warn = FALSE,
  encoding = "UTF-8"
)

if (length(header_line) != 1) {
  stop("Could not read the header of the expression matrix.")
}

cell_names <- strsplit(header_line, ",", fixed = TRUE)[[1]]
cell_names <- sub('^"', "", cell_names)
cell_names <- sub('"$', "", cell_names)

if (anyDuplicated(cell_names)) {
  stop("Duplicate cell barcodes were found in the expression-matrix header.")
}

cat("Expression matrix contains ", length(cell_names),
    " cell columns.\n", sep = "")

# -------------------------------------------------------------------------
# 4. 读取表达矩阵
#
# 该 CSV 的 header 只有 cell names；
# 每个数据行的第一列是 gene name，后面是 expression values。
# 因此这里跳过 header，以 header = FALSE 读取所有数据行。
# -------------------------------------------------------------------------

cat("\n===== 3. Reading normalized expression matrix =====\n")
cat("This file is large and may take several minutes to load.\n\n")

expr_dt <- data.table::fread(
  file = control_file,
  skip = 1,
  header = FALSE,
  data.table = TRUE,
  showProgress = TRUE,
  check.names = FALSE
)

expected_columns <- length(cell_names) + 1L

if (ncol(expr_dt) != expected_columns) {
  stop(
    "Unexpected number of columns in the expression matrix.\n",
    "Expected: 1 gene-name column + ", length(cell_names),
    " cell columns = ", expected_columns, "\n",
    "Observed: ", ncol(expr_dt)
  )
}

gene_names <- as.character(expr_dt[[1]])
expr_dt[[1]] <- NULL

if (anyNA(gene_names) || any(gene_names == "")) {
  stop("Missing or blank gene names were found.")
}

if (anyDuplicated(gene_names)) {
  duplicate_n <- sum(duplicated(gene_names))
  warning(
    duplicate_n,
    " duplicated gene name(s) were found. ",
    "make.unique() will be used to preserve every row."
  )
  gene_names <- make.unique(gene_names)
}

non_numeric_columns <- which(
  !vapply(expr_dt, is.numeric, logical(1))
)

if (length(non_numeric_columns) > 0) {
  stop(
    "Non-numeric expression columns were detected: ",
    paste(head(non_numeric_columns, 20), collapse = ", ")
  )
}

cat("Converting expression table to a numeric matrix...\n")
expr_matrix <- as.matrix(expr_dt)

rm(expr_dt)
gc()

rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- cell_names

cat(
  "Dense expression matrix dimensions: ",
  nrow(expr_matrix), " genes x ",
  ncol(expr_matrix), " cells\n",
  sep = ""
)

if (any(!is.finite(expr_matrix))) {
  stop("The expression matrix contains NA, NaN, or infinite values.")
}

if (any(expr_matrix < 0)) {
  stop("Negative expression values were detected unexpectedly.")
}

# 该数据应包含非整数的小数值，因为它是标准化表达矩阵
non_integer_fraction <- mean(
  abs(expr_matrix - round(expr_matrix)) > .Machine$double.eps^0.5
)

cat(
  "Fraction of non-integer values: ",
  signif(non_integer_fraction, 4), "\n",
  sep = ""
)

if (non_integer_fraction == 0) {
  warning(
    "All values appear to be integers. ",
    "Please verify whether this file is truly the normalized GEO export."
  )
}

# -------------------------------------------------------------------------
# 5. 检查 metadata 与表达矩阵 barcode
# -------------------------------------------------------------------------

cat("\n===== 4. Matching expression data to metadata =====\n")

if (!setequal(colnames(expr_matrix), rownames(meta))) {
  missing_in_meta <- setdiff(colnames(expr_matrix), rownames(meta))
  missing_in_expr <- setdiff(rownames(meta), colnames(expr_matrix))

  stop(
    "Expression matrix and metadata do not contain the same cells.\n",
    "Cells missing from metadata: ", length(missing_in_meta), "\n",
    "Cells missing from expression matrix: ", length(missing_in_expr)
  )
}

# 即使原始顺序不同，也强制按照表达矩阵列顺序重排 metadata
meta <- meta[colnames(expr_matrix), , drop = FALSE]

if (!identical(colnames(expr_matrix), rownames(meta))) {
  stop("Cell order could not be aligned between expression matrix and metadata.")
}

cat("All ", ncol(expr_matrix),
    " cell barcodes match exactly and are in the same order.\n", sep = "")

# -------------------------------------------------------------------------
# 6. 转成稀疏矩阵
# -------------------------------------------------------------------------

cat("\n===== 5. Converting to a sparse matrix =====\n")

normalized_sparse <- Matrix::Matrix(
  expr_matrix,
  sparse = TRUE
)

rm(expr_matrix)
gc()

cat("Sparse matrix class: ", class(normalized_sparse)[1], "\n", sep = "")
cat(
  "Sparse matrix dimensions: ",
  nrow(normalized_sparse), " genes x ",
  ncol(normalized_sparse), " cells\n",
  sep = ""
)

cat(
  "Non-zero entries: ",
  length(normalized_sparse@x), "\n",
  sep = ""
)

# -------------------------------------------------------------------------
# 7. 创建 data-only Seurat v5 object
#
# 注意：
# normalized_sparse 放入 RNA assay 的 "data" layer。
# 不创建假的 counts layer，也不运行 NormalizeData()。
# -------------------------------------------------------------------------

cat("\n===== 6. Creating a data-only Seurat v5 object =====\n")

rna_assay <- SeuratObject::CreateAssay5Object(
  data = normalized_sparse
)

control_obj <- SeuratObject::CreateSeuratObject(
  counts = rna_assay,
  assay = "RNA",
  project = "Wells2020_Figure1b",
  meta.data = meta
)

SeuratObject::DefaultAssay(control_obj) <- "RNA"

cat("\nSeurat object:\n")
print(control_obj)

cat("\nRNA assay layers:\n")
print(SeuratObject::Layers(control_obj[["RNA"]]))

if (!"data" %in% SeuratObject::Layers(control_obj[["RNA"]])) {
  stop("The RNA assay does not contain the expected 'data' layer.")
}

if ("counts" %in% SeuratObject::Layers(control_obj[["RNA"]])) {
  warning(
    "A counts layer was unexpectedly created. ",
    "Inspect the object before continuing."
  )
}

if (!identical(colnames(control_obj), rownames(control_obj[[]]))) {
  stop("Seurat object cell order and metadata row order do not match.")
}

expected_dimensions <- c(20309L, 2434L)

if (!identical(as.integer(dim(control_obj)), expected_dimensions)) {
  warning(
    "Object dimensions differ from the expected 20,309 genes x 2,434 cells.\n",
    "Observed: ",
    paste(dim(control_obj), collapse = " x ")
  )
}

# -------------------------------------------------------------------------
# 8. 保存对象和检查报告
# -------------------------------------------------------------------------

cat("\n===== 7. Saving the initial Seurat object =====\n")

saveRDS(
  control_obj,
  file = output_file
)

summary_file <- file.path(
  output_dir,
  "01_control_seurat_initial_summary.txt"
)

summary_lines <- c(
  paste0("Project root: ", project_root),
  paste0("Input expression file: ", control_file),
  paste0("Input metadata file: ", meta_file),
  paste0("Output Seurat object: ", output_file),
  "",
  paste0("Genes: ", nrow(control_obj)),
  paste0("Cells: ", ncol(control_obj)),
  paste0(
    "RNA layers: ",
    paste(SeuratObject::Layers(control_obj[["RNA"]]), collapse = ", ")
  ),
  "",
  "Cells per experiment:",
  paste(capture.output(print(table(control_obj$exp))), collapse = "\n"),
  "",
  "Cells per Figure 1b stage:",
  paste(
    capture.output(
      print(sort(table(control_obj$stage_figure1b), decreasing = TRUE))
    ),
    collapse = "\n"
  )
)

writeLines(summary_lines, summary_file)

session_file <- file.path(output_dir, "01_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)

cat("\nDone.\n")
cat("Saved Seurat object:\n", output_file, "\n", sep = "")
cat("Saved summary:\n", summary_file, "\n", sep = "")
cat("Saved session information:\n", session_file, "\n", sep = "")

cat("\nFinal validation:\n")
cat("  Dimensions: ",
    nrow(control_obj), " genes x ",
    ncol(control_obj), " cells\n", sep = "")
cat("  Cell/metadata order identical: ",
    identical(colnames(control_obj), rownames(control_obj[[]])),
    "\n", sep = "")
cat("  NormalizeData() was NOT run.\n")
