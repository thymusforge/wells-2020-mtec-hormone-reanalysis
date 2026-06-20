# 00_inspect_geo_files.R
# 目的：
# 1) 检查 combinedControl、allSamples 和 controls_meta 的基本结构
# 2) 判断 Figure 1b 应该使用哪一个表达矩阵
# 3) 不修改任何原始数据，也不创建 Seurat 对象

rm(list = ls())

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
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

control_file <- file.path(data_dir, "GSE137699_combinedControl.csv")
all_file     <- file.path(data_dir, "GSE137699_allSamples.csv")
meta_file    <- file.path(data_dir, "controls_meta.csv")

required <- c(control_file, meta_file)
missing_required <- required[!file.exists(required)]

if (length(missing_required) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_required, collapse = "\n"),
    "\n\nPlease run this script from the project root or one of its subdirectories.",
    "\nIf these data files are missing, see README.md Data source and data/README.md."
  )
}

# 只读取CSV第一行，取得细胞条形码；不会把整个大矩阵读入内存
read_cell_header <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  first_line <- readLines(con, n = 1, warn = FALSE)
  if (length(first_line) != 1) {
    stop("Could not read the header from: ", path)
  }

  cells <- strsplit(first_line, ",", fixed = TRUE)[[1]]
  cells <- sub('^"', "", cells)
  cells <- sub('"$', "", cells)
  cells
}

cat("===== 1. Read controls metadata =====\n")
meta <- read.csv(
  meta_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Metadata dimensions:", nrow(meta), "cells x", ncol(meta), "columns\n")
cat("Metadata columns:\n")
print(colnames(meta))

if ("exp" %in% colnames(meta)) {
  cat("\nCells per control sample:\n")
  print(table(meta$exp))
}

if ("stage" %in% colnames(meta)) {
  cat("\nCells per annotated stage:\n")
  print(sort(table(meta$stage), decreasing = TRUE))
}

cat("\n===== 2. Inspect combinedControl header =====\n")
control_cells <- read_cell_header(control_file)

cat("Number of expression-matrix columns:", length(control_cells), "\n")
cat("First 5 cell names:\n")
print(head(control_cells, 5))

cat("\nDoes combinedControl exactly match metadata row names?\n")
exact_match <- identical(control_cells, rownames(meta))
print(exact_match)

cat("\nBarcode overlap:\n")
cat(
  "Metadata cells found in combinedControl:",
  sum(rownames(meta) %in% control_cells),
  "/", nrow(meta), "\n"
)
cat(
  "combinedControl cells found in metadata:",
  sum(control_cells %in% rownames(meta)),
  "/", length(control_cells), "\n"
)

cat("\nCell-name prefixes in combinedControl:\n")
prefix <- sub("_(.*)$", "", control_cells)
print(table(prefix))

cat("\n===== 3. Inspect allSamples header, if present =====\n")
if (file.exists(all_file)) {
  all_cells <- read_cell_header(all_file)

  cat("Number of allSamples expression-matrix columns:", length(all_cells), "\n")
  cat("First 5 cell names:\n")
  print(head(all_cells, 5))

  cat("\nHow many Figure 1b control cells are also present in allSamples?\n")
  cat(sum(control_cells %in% all_cells), "/", length(control_cells), "\n")

  cat("\nCell-name prefixes in allSamples:\n")
  all_prefix <- sub("_(.*)$", "", all_cells)
  print(sort(table(all_prefix), decreasing = TRUE))
} else {
  cat(
    "allSamples file was not found at:\n", all_file,
    "\nThis does not prevent Figure 1b reanalysis.\n"
  )
}

cat("\n===== 4. Interpretation =====\n")
if (exact_match && length(control_cells) == 2434 && nrow(meta) == 2434) {
  cat(
    "combinedControl and controls_meta describe the same 2,434 control cells.\n",
    "This is the appropriate processed dataset for Figure 1b reanalysis.\n",
    "allSamples is a broader combined dataset and is not needed for Figure 1b.\n",
    sep = ""
  )
} else {
  cat(
    "The files are not an exact one-to-one match.\n",
    "Do not create a Seurat object until the barcode mismatch is resolved.\n",
    sep = ""
  )
}

cat(
  "\nImportant: these CSV files are expression/metadata exports, not the original\n",
  "seurat_controls_merged.rda object. They do not contain the original UMAP coordinates.\n",
  sep = ""
)
