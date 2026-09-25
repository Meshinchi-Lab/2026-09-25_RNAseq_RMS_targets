#!/usr/bin/env Rscript
## Stage 3. Writes the two small tables the app holds in memory:
##   samples.parquet    one row per sample, all datasets
##   gene_index.parquet one row per gene, backing the server-side selectize

suppressPackageStartupMessages({
  library(duckdb); library(arrow); library(dplyr)
})

root <- rprojroot::find_root(rprojroot::has_file("config.yml"))
setwd(root)
source("R/harmonise.R")
source("R/loaders.R")
cfg <- config::get(file = "config.yml")

meta_files <- list.files("meta_staging", pattern = "[.]parquet$", full.names = TRUE)
if (!length(meta_files)) stop("no meta_staging/*.parquet — run bin/01_build_store.R first")

samples <- do.call(rbind, lapply(meta_files,
                                 function(f) as.data.frame(arrow::read_parquet(f))))

## The key is (dataset, sample_id), NOT sample_id alone. The same RMS cell lines (RD, JR,
## TTC442) are measured in both DepMap and Wachtel, which is exactly the cross-modality
## comparison this app exists to show. Uniqueness is enforced per dataset by each loader.
stopifnot(!anyDuplicated(samples[c("dataset", "sample_id")]))
assert_vocab(samples$group, VALID_GROUPS, "group")
assert_vocab(samples$assay, VALID_ASSAYS, "assay")
arrow::write_parquet(samples, file.path(cfg$store_dir, "samples.parquet"))
message(sprintf("samples.parquet: %s rows, %d datasets",
                format(nrow(samples), big.mark = ","),
                dplyr::n_distinct(samples$dataset)))

con <- dbConnect(duckdb::duckdb(shared_home = FALSE))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
invisible(dbExecute(con, sprintf("SET memory_limit='%s'", cfg$duckdb_memory_limit)))

expr_glob <- file.path(cfg$store_dir, "expression", "**", "*.parquet")
gi <- dbGetQuery(con, sprintf("
  SELECT gene_symbol,
         string_agg(DISTINCT dataset, ',' ORDER BY dataset) AS datasets,
         count(*)::INTEGER AS n_samples_total
  FROM read_parquet('%s', hive_partitioning = 1)
  GROUP BY gene_symbol
  ORDER BY gene_symbol", expr_glob))

arrow::write_parquet(gi, file.path(cfg$store_dir, "gene_index.parquet"))
message(sprintf("gene_index.parquet: %s genes", format(nrow(gi), big.mark = ",")))
