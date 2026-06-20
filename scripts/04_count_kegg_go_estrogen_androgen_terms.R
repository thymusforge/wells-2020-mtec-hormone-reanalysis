# 04_count_kegg_go_estrogen_androgen_terms.R
# Count KEGG and GO terms related to estrogen/estradiol or androgen.
#
# Output:
#   output/04_kegg_go_sex_hormone_terms/
#     04_term_count_summary.csv
#     04_matching_terms_all.csv
#     04_matching_terms_GO.csv
#     04_matching_terms_KEGG.csv
#     04_sessionInfo.txt

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

discover_subcollections <- function(collection, pattern, fallback_subcollections) {
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

  keep <- collections[[collection_col]] == collection &
    grepl(pattern, collections[[subcollection_col]], ignore.case = TRUE)

  found <- sort(unique(collections[[subcollection_col]][keep]))
  if (length(found) == 0) {
    fallback_subcollections
  } else {
    found
  }
}

load_msig_subcollections <- function(collection, subcollections, source_name_prefix) {
  tables <- lapply(
    subcollections,
    function(subcollection) {
      msig <- load_msig_collection(collection = collection, subcollection = subcollection)
      standardized <- standardize_terms(msig, paste0(source_name_prefix, "_", subcollection))
      if (nrow(standardized) > 0) {
        standardized$source <- source_name_prefix
        standardized$gs_subcollection <- subcollection
      }
      standardized
    }
  )

  tables <- tables[vapply(tables, nrow, integer(1)) > 0]
  if (length(tables) == 0) {
    return(data.frame())
  }

  unique(do.call(rbind, tables))
}

standardize_terms <- function(msig, source_name) {
  if (is.null(msig) || nrow(msig) == 0) {
    return(data.frame())
  }

  msig <- as.data.frame(msig, stringsAsFactors = FALSE)

  required_cols <- c("gs_name", "gene_symbol")
  if (!all(required_cols %in% colnames(msig))) {
    stop(
      "msigdbr output for ",
      source_name,
      " does not contain expected columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  optional_cols <- c("gs_id", "gs_exact_source", "gs_description", "gs_collection", "gs_subcollection")
  for (col in optional_cols) {
    if (!col %in% colnames(msig)) {
      msig[[col]] <- NA_character_
    }
  }

  term_info <- unique(msig[, c(optional_cols, "gs_name")])
  term_info$source <- source_name

  gene_counts <- stats::aggregate(
    gene_symbol ~ gs_name,
    data = msig,
    FUN = function(x) length(unique(x))
  )
  colnames(gene_counts)[colnames(gene_counts) == "gene_symbol"] <- "total_genes_in_set"

  if ("db_gene_symbol" %in% colnames(msig)) {
    db_gene_counts <- stats::aggregate(
      db_gene_symbol ~ gs_name,
      data = msig,
      FUN = function(x) length(unique(x))
    )
    colnames(db_gene_counts)[colnames(db_gene_counts) == "db_gene_symbol"] <- "total_db_genes_in_set"
    gene_counts <- merge(gene_counts, db_gene_counts, by = "gs_name", all.x = TRUE, sort = FALSE)
  } else {
    gene_counts$total_db_genes_in_set <- NA_integer_
  }

  term_info <- merge(term_info, gene_counts, by = "gs_name", all.x = TRUE, sort = FALSE)

  term_info[, c(
    "source",
    "gs_collection",
    "gs_subcollection",
    "gs_name",
    "gs_id",
    "gs_exact_source",
    "gs_description",
    "total_genes_in_set",
    "total_db_genes_in_set"
  )]
}

match_terms <- function(term_df) {
  if (nrow(term_df) == 0) {
    term_df$keyword_group <- character(0)
    term_df$matched_keyword <- character(0)
    return(term_df)
  }

  searchable_text <- paste(
    term_df$gs_name,
    term_df$gs_description,
    term_df$gs_exact_source,
    sep = " "
  )

  estrogen_hit <- grepl("ESTROGEN|ESTRADIOL", searchable_text, ignore.case = TRUE)
  androgen_hit <- grepl("ANDROGEN", searchable_text, ignore.case = TRUE)

  out <- rbind(
    transform(
      term_df[estrogen_hit, , drop = FALSE],
      keyword_group = "Estrogen_or_estradiol",
      matched_keyword = "ESTROGEN|ESTRADIOL"
    ),
    transform(
      term_df[androgen_hit, , drop = FALSE],
      keyword_group = "Androgen",
      matched_keyword = "ANDROGEN"
    )
  )

  out <- unique(out)
  rownames(out) <- NULL
  out
}

# -------------------------------------------------------------------------
# 3. Load KEGG and GO terms from MSigDB
# -------------------------------------------------------------------------

cat("===== Loading MSigDB KEGG and GO terms =====\n")

kegg_subcollections <- discover_subcollections(
  collection = "C2",
  pattern = "KEGG",
  fallback_subcollections = c("CP:KEGG", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS")
)

go_subcollections <- discover_subcollections(
  collection = "C5",
  pattern = "^GO:",
  fallback_subcollections = c("GO:BP", "GO:CC", "GO:MF")
)

cat("KEGG subcollections queried:\n")
print(kegg_subcollections)
cat("\nGO subcollections queried:\n")
print(go_subcollections)
cat("\n")

kegg_terms <- load_msig_subcollections(
  collection = "C2",
  subcollections = kegg_subcollections,
  source_name_prefix = "KEGG"
)

go_terms <- load_msig_subcollections(
  collection = "C5",
  subcollections = go_subcollections,
  source_name_prefix = "GO"
)

term_tables <- list(
  KEGG = kegg_terms,
  GO = go_terms
)

term_tables <- term_tables[vapply(term_tables, nrow, integer(1)) > 0]
if (length(term_tables) == 0) {
  stop("No KEGG or GO terms were loaded from msigdbr.")
}

all_terms <- unique(do.call(rbind, term_tables))
rownames(all_terms) <- NULL

cat("Terms loaded:\n")
print(table(all_terms$source))
cat("\n")

# -------------------------------------------------------------------------
# 4. Count estrogen/androgen-related terms
# -------------------------------------------------------------------------

matching_terms <- match_terms(all_terms)

all_sources <- unique(all_terms$source)
all_keyword_groups <- c("Estrogen_or_estradiol", "Androgen")
summary_grid <- expand.grid(
  source = all_sources,
  keyword_group = all_keyword_groups,
  stringsAsFactors = FALSE
)

if (nrow(matching_terms) > 0) {
  summary_counts <- aggregate(
    gs_name ~ source + keyword_group,
    data = matching_terms,
    FUN = function(x) length(unique(x))
  )
  colnames(summary_counts)[colnames(summary_counts) == "gs_name"] <- "matching_term_count"

  summary_counts <- merge(
    summary_grid,
    summary_counts,
    by = c("source", "keyword_group"),
    all.x = TRUE,
    sort = FALSE
  )
} else {
  summary_counts <- summary_grid
  summary_counts$matching_term_count <- 0L
}
summary_counts$matching_term_count[is.na(summary_counts$matching_term_count)] <- 0L

total_by_source <- aggregate(
  gs_name ~ source,
  data = all_terms,
  FUN = function(x) length(unique(x))
)
colnames(total_by_source)[colnames(total_by_source) == "gs_name"] <- "total_terms_in_source"

summary_counts <- merge(summary_counts, total_by_source, by = "source", all.x = TRUE, sort = FALSE)
summary_counts <- summary_counts[order(summary_counts$source, summary_counts$keyword_group), ]
rownames(summary_counts) <- NULL

cat("Matching term counts:\n")
print(summary_counts)
cat("\n")

# -------------------------------------------------------------------------
# 5. Save outputs
# -------------------------------------------------------------------------

utils::write.csv(
  summary_counts,
  file = file.path(result_dir, "04_term_count_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

utils::write.csv(
  matching_terms,
  file = file.path(result_dir, "04_matching_terms_all.csv"),
  row.names = FALSE,
  quote = TRUE
)

utils::write.csv(
  matching_terms[grepl("^GO", matching_terms$source), , drop = FALSE],
  file = file.path(result_dir, "04_matching_terms_GO.csv"),
  row.names = FALSE,
  quote = TRUE
)

utils::write.csv(
  matching_terms[matching_terms$source == "KEGG", , drop = FALSE],
  file = file.path(result_dir, "04_matching_terms_KEGG.csv"),
  row.names = FALSE,
  quote = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(result_dir, "04_sessionInfo.txt")
)

cat("Done.\n")
cat("Summary file:\n", file.path(result_dir, "04_term_count_summary.csv"), "\n", sep = "")
cat("All matching terms:\n", file.path(result_dir, "04_matching_terms_all.csv"), "\n", sep = "")
