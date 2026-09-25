## DuckDB-backed store access. One connection is opened at app start; the two small index
## tables are held in memory so that a single-gene query returns ~22k rows and every
## downstream filter is plain dplyr.

suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(arrow); library(dplyr)
})

#' DuckDB driver, tolerant of the version that happens to be installed.
#'
#' `shared_home` was added in duckdb 1.5.5 and silences the extension/secret-home notice.
#' It does not exist in earlier builds — rocker/shiny-verse:4.5.3 ships duckdb 1.5.2, where
#' passing it fails with "... must be empty". Probe the formals rather than the version
#' string so this keeps working in both directions.
duckdb_driver <- function() {
  if ("shared_home" %in% names(formals(duckdb::duckdb))) {
    duckdb::duckdb(shared_home = FALSE)
  } else {
    duckdb::duckdb()
  }
}

#' Open the store. Returns a handle used by query_gene().
store_open <- function(store_dir, memory_limit = "4GB") {
  gi_f <- file.path(store_dir, "gene_index.parquet")
  sm_f <- file.path(store_dir, "samples.parquet")
  ex_d <- file.path(store_dir, "expression")

  present <- file.exists(c(gi_f, sm_f, ex_d))
  if (!all(present)) {
    missing <- c(gi_f, sm_f, ex_d)[!present]
    stop(sprintf(paste0("store incomplete at '%s' (missing: %s).\n",
                        "Build it with:\n",
                        "  Rscript bin/01_build_store.R\n",
                        "  Rscript bin/02_repartition_store.R\n",
                        "  Rscript bin/03_build_indexes.R"),
                 store_dir, paste(basename(missing), collapse = ", ")), call. = FALSE)
  }

  con <- dbConnect(duckdb_driver())
  dbExecute(con, sprintf("SET memory_limit='%s'", memory_limit))

  list(con        = con,
       store_dir  = store_dir,
       expr_glob  = file.path(store_dir, "expression", "**", "*.parquet"),
       samples    = as.data.frame(arrow::read_parquet(sm_f)),
       gene_index = as.data.frame(arrow::read_parquet(gi_f)))
}

store_close <- function(st) {
  try(dbDisconnect(st$con, shutdown = TRUE), silent = TRUE)
  invisible(NULL)
}

#' Bucket a gene symbol the same way the ETL did.
gene_bucket <- function(gene) {
  b <- toupper(substr(gene, 1, 1))
  if (grepl("^[A-Z]$", b)) b else "OTHER"
}

#' All measurements for one gene, with sample metadata joined.
#'
#' The gene_bucket predicate is supplied explicitly so partition pruning fires; relying on
#' DuckDB to derive it from gene_symbol would scan every bucket.
query_gene <- function(st, gene, datasets = NULL) {
  gene <- toupper(trimws(as.character(gene)))
  if (!length(gene) || is.na(gene) || !nzchar(gene)) return(.empty_gene_frame())

  sql <- sprintf(
    "SELECT gene_symbol, dataset, sample_id, value
     FROM read_parquet('%s', hive_partitioning = 1)
     WHERE gene_bucket = ? AND gene_symbol = ?",
    st$expr_glob)

  d <- dbGetQuery(st$con, sql, params = list(gene_bucket(gene), gene))
  if (!nrow(d)) return(.empty_gene_frame())

  out <- dplyr::left_join(d, st$samples, by = c("sample_id", "dataset"))
  if (!is.null(datasets)) out <- out[out$dataset %in% datasets, , drop = FALSE]
  out
}

#' Append explicit not-detected rows for samples that have no value for this gene.
#'
#' Mass spectrometry reports a value only where a protein was detected. Wachtel's raw file
#' carries a row for all 15 samples per gene, but `iTop3` is NA where the protein was below
#' the detection limit, and the ETL drops those. Dropping them from the PLOT is misleading:
#' a protein detected in 2 of 6 FN lines at high abundance reads as "high in FN" when it is
#' actually absent in 4 of 6. In label-free MS non-detection is evidence of low protein, not
#' missing-at-random, so the denominator has to stay visible.
#'
#' Returns `d` with extra rows carrying `value = NA` for every missing sample.
add_undetected <- function(st, d, datasets) {
  roster <- st$samples[st$samples$dataset %in% datasets, , drop = FALSE]
  if (!nrow(roster)) return(d)

  have <- paste(d$dataset, d$sample_id, sep = "\r")
  miss <- roster[!paste(roster$dataset, roster$sample_id, sep = "\r") %in% have, ,
                 drop = FALSE]
  if (!nrow(miss)) return(d)

  gene <- if (nrow(d)) d$gene_symbol[1] else NA_character_
  pad <- data.frame(gene_symbol = gene, dataset = miss$dataset,
                    sample_id = miss$sample_id, value = NA_real_,
                    stringsAsFactors = FALSE)
  pad <- cbind(pad, miss[, setdiff(names(miss), c("dataset", "sample_id")), drop = FALSE])
  rbind(d[, names(pad), drop = FALSE], pad)
}

.empty_gene_frame <- function() {
  data.frame(gene_symbol = character(0), dataset = character(0),
             sample_id = character(0), value = numeric(0),
             assay = character(0), unit = character(0), group = character(0),
             tissue = character(0), subgroup = character(0),
             is_tumour = logical(0), stringsAsFactors = FALSE)
}

#' n / median / IQR per dataset x group. The app always shows this next to the plot:
#' a violin over Wachtel's n=15 or GSE108022's five controls is misleading without it.
summarise_gene <- function(d) {
  if (!nrow(d)) {
    return(data.frame(dataset = character(0), group = character(0), unit = character(0),
                      n = integer(0), median = numeric(0),
                      q25 = numeric(0), q75 = numeric(0)))
  }
  d |>
    dplyr::group_by(.data$dataset, .data$group, .data$unit) |>
    dplyr::summarise(n = dplyr::n(),
                     median = stats::median(.data$value, na.rm = TRUE),
                     q25 = unname(stats::quantile(.data$value, 0.25, na.rm = TRUE)),
                     q75 = unname(stats::quantile(.data$value, 0.75, na.rm = TRUE)),
                     .groups = "drop") |>
    dplyr::arrange(.data$dataset, .data$group) |>
    as.data.frame()
}
