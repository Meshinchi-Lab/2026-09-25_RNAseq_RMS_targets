#!/usr/bin/env Rscript
## Stage 2. DuckDB repartitions the sample-blocked staging files into gene-bucketed
## Parquet, out-of-core with an explicit memory ceiling. This is the step that makes
## single-gene queries fast: a query for one gene opens one gene_bucket directory rather
## than scanning ~411M rows.

suppressPackageStartupMessages(library(duckdb))

root <- rprojroot::find_root(rprojroot::has_file("config.yml"))
setwd(root)
cfg <- config::get(file = "config.yml")

if (!dir.exists(cfg$staging_dir)) {
  stop("no staging directory — run bin/01_build_store.R first", call. = FALSE)
}

out <- file.path(cfg$store_dir, "expression")
dir.create(cfg$store_dir, showWarnings = FALSE, recursive = TRUE)
unlink(out, recursive = TRUE)

source(file.path(root, "shiny", "rms_gene_explorer", "R", "store.R"))
con <- dbConnect(duckdb_driver())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

invisible(dbExecute(con, sprintf("SET memory_limit='%s'", cfg$duckdb_memory_limit)))
tmp <- file.path(cfg$store_dir, "duckdb_tmp")
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
invisible(dbExecute(con, sprintf("SET temp_directory='%s'", tmp)))
invisible(dbExecute(con, "SET preserve_insertion_order=false"))

t0 <- Sys.time()
invisible(dbExecute(con, sprintf("
  COPY (
    SELECT gene_symbol,
           dataset,
           sample_id,
           CAST(value AS FLOAT) AS value,
           CASE WHEN regexp_matches(gene_symbol, '^[A-Za-z]')
                THEN upper(substr(gene_symbol, 1, 1))
                ELSE 'OTHER' END AS gene_bucket
    FROM read_parquet('%s/**/*.parquet', hive_partitioning = 1)
  )
  TO '%s'
  (FORMAT PARQUET, PARTITION_BY (gene_bucket), COMPRESSION zstd, OVERWRITE_OR_IGNORE)",
  cfg$staging_dir, out)))

n <- dbGetQuery(con, sprintf(
  "SELECT count(*) AS n FROM read_parquet('%s/**/*.parquet', hive_partitioning = 1)",
  out))$n

unlink(tmp, recursive = TRUE)

message(sprintf("Repartitioned %s rows in %.1f min",
                format(n, big.mark = ","),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
message("Next: Rscript bin/03_build_indexes.R")
