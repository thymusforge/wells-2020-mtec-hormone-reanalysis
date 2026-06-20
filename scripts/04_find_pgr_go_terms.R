# 04_find_pgr_go_terms.R
# Find GO terms/pathways that contain the mouse Pgr gene in MSigDB.
#
# Output:
#   output/04_kegg_go_sex_hormone_terms/
#     04B_pgr_go_term_summary.csv
#     04B_pgr_go_terms.csv
#     04B_pgr_go_geneset_members.csv
#     04B_sessionInfo.txt

rm(list = ls())
gc()

# -------------------------------------------------------------------------
# 0. Locate project root
# -------------------------------------------------------------------------

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
output_dir   <- file.path(project_root, "output")
result_dir   <- file.path(output_dir, "04_kegg_go_sex_hormone_terms")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project root:\n", project_root, "\n\n", sep = "")
cat("Output directory:\n", result_dir, "\n\n", sep = "")

# -------------------------------------------------------------------------
# 1. Packages
# -------------------------------------------------------------------------

required_packages <- c("msigdbr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\n\nInstall first with:\n",
    'install.packages("msigdbr")'
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

load_msig_collection <- function(collection, subcollection = NULL, species = "Mus musculus") {
  msig_fun <- msigdbr::msigdbr
  msig_args <- names(formals(msig_fun))

  tryCatch(
    {
      if ("collection" %in% msig_args) {
        args <- list(species = species, collection = collection)
        if (!is.null(subcollection)) {
          args$subcollection <- subcollection
        }
        do.call(msigdbr::msigdbr, args)
      } else {
        args <- list(species = species, category = collection)
        if (!is.null(subcollection)) {
          args$subcategory <- subcollection
        }
        do.call(msigdbr::msigdbr, args)
      }
    },
    error = function(e) {
      warning(
        "Could not load MSigDB collection ",
        collection,
        ifelse(is.null(subcollection), "", paste0("/", subcollection)),
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

discover_go_subcollections <- function() {
  fallback_subcollections <- c("GO:BP", "GO:CC", "GO:MF")

  if (!exists("msigdbr_collections", where = asNamespace("msigdbr"), mode = "function")) {
    return(fallback_subcollections)
  }

  collections <- tryCatch(
    msigdbr::msigdbr_collections(),
    error = function(e) NULL
  )

  if (is.null(collections) || nrow(collections) == 0) {
    return(fallback_subcollections)
  }

  collection_col <- intersect(c("gs_collection", "collection", "category"), colnames(collections))[1]
  subcollection_col <- intersect(
    c("gs_subcollection", "subcollection", "subcategory"),
    colnames(collections)
  )[1]

  if (is.na(collection_col) || is.na(subcollection_col)) {
    return(fallback_subcollections)
  }

  keep <- collections[[collection_col]] == "C5" &
    grepl("^GO:", collections[[subcollection_col]], ignore.case = TRUE)

  found <- sort(unique(collections[[subcollection_col]][keep]))
  if (length(found) == 0) {
    fallback_subcollections
  } else {
    found
  }
}

standardize_go <- function(msig, subcollection) {
  if (is.null(msig) || nrow(msig) == 0) {
    return(data.frame())
  }

  required_cols <- c("gs_name", "gene_symbol")
  if (!all(required_cols %in% colnames(msig))) {
    stop(
      "msigdbr output for ",
      subcollection,
      " does not contain expected columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  optional_cols <- c("gs_id", "gs_exact_source", "gs_description", "gs_collection")
  for (col in optional_cols) {
    if (!col %in% colnames(msig)) {
      msig[[col]] <- NA_character_
    }
  }

  msig$source <- "GO"
  msig$gs_subcollection <- subcollection
  msig[, c(
    "source",
    "gs_collection",
    "gs_subcollection",
    "gs_name",
    "gs_id",
    "gs_exact_source",
    "gs_description",
    "gene_symbol"
  )]
}

# -------------------------------------------------------------------------
# 3. Load GO terms and find Pgr-containing gene sets
# -------------------------------------------------------------------------

target_gene <- "Pgr"

cat("===== Loading MSigDB GO terms =====\n")

go_subcollections <- discover_go_subcollections()
cat("GO subcollections queried:\n")
print(go_subcollections)
cat("\n")

go_tables <- lapply(
  go_subcollections,
  function(subcollection) {
    msig <- load_msig_collection(collection = "C5", subcollection = subcollection)
    standardize_go(msig, subcollection)
  }
)

go_tables <- go_tables[vapply(go_tables, nrow, integer(1)) > 0]
if (length(go_tables) == 0) {
  stop("No GO terms were loaded from msigdbr.")
}

go_all <- do.call(rbind, go_tables)
rownames(go_all) <- NULL

cat("GO terms loaded:\n")
print(table(go_all$gs_subcollection))
cat("\n")

pgr_rows <- go_all[toupper(go_all$gene_symbol) == toupper(target_gene), , drop = FALSE]

pgr_terms <- unique(pgr_rows[, c(
  "source",
  "gs_collection",
  "gs_subcollection",
  "gs_name",
  "gs_id",
  "gs_exact_source",
  "gs_description"
)])

pgr_terms$total_genes_in_set <- as.integer(vapply(
  pgr_terms$gs_name,
  function(term) length(unique(go_all$gene_symbol[go_all$gs_name == term])),
  integer(1)
))

pgr_terms$gene_found <- target_gene
pgr_terms <- pgr_terms[order(pgr_terms$gs_subcollection, pgr_terms$gs_name), ]
rownames(pgr_terms) <- NULL

pgr_members <- go_all[go_all$gs_name %in% pgr_terms$gs_name, , drop = FALSE]
pgr_members <- pgr_members[order(pgr_members$gs_subcollection, pgr_members$gs_name, pgr_members$gene_symbol), ]
rownames(pgr_members) <- NULL

summary_counts <- aggregate(
  gs_name ~ gs_subcollection,
  data = pgr_terms,
  FUN = function(x) length(unique(x))
)
colnames(summary_counts)[colnames(summary_counts) == "gs_name"] <- "pgr_go_term_count"

summary_grid <- data.frame(
  gs_subcollection = go_subcollections,
  stringsAsFactors = FALSE
)

summary_counts <- merge(summary_grid, summary_counts, by = "gs_subcollection", all.x = TRUE, sort = FALSE)
summary_counts$pgr_go_term_count[is.na(summary_counts$pgr_go_term_count)] <- 0L

total_go_counts <- aggregate(
  gs_name ~ gs_subcollection,
  data = unique(go_all[, c("gs_subcollection", "gs_name")]),
  FUN = length
)
colnames(total_go_counts)[colnames(total_go_counts) == "gs_name"] <- "total_go_terms"

summary_counts <- merge(summary_counts, total_go_counts, by = "gs_subcollection", all.x = TRUE, sort = FALSE)
summary_counts$target_gene <- target_gene
summary_counts <- summary_counts[, c("target_gene", "gs_subcollection", "pgr_go_term_count", "total_go_terms")]

cat("Pgr-containing GO term counts:\n")
print(summary_counts)
cat("\n")

# -------------------------------------------------------------------------
# 4. Save outputs
# -------------------------------------------------------------------------

utils::write.csv(
  summary_counts,
  file = file.path(result_dir, "04B_pgr_go_term_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

utils::write.csv(
  pgr_terms,
  file = file.path(result_dir, "04B_pgr_go_terms.csv"),
  row.names = FALSE,
  quote = TRUE
)

utils::write.csv(
  pgr_members,
  file = file.path(result_dir, "04B_pgr_go_geneset_members.csv"),
  row.names = FALSE,
  quote = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(result_dir, "04B_sessionInfo.txt")
)

cat("Done.\n")
cat("Summary file:\n", file.path(result_dir, "04B_pgr_go_term_summary.csv"), "\n", sep = "")
cat("Pgr GO terms:\n", file.path(result_dir, "04B_pgr_go_terms.csv"), "\n", sep = "")
