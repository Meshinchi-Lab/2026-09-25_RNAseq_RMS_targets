#!/usr/bin/env Rscript
## Stage 4. Reconciles the built store against the sources and fails loudly on mismatch.
## Run after every rebuild — a silently truncated dataset is otherwise invisible.

suppressPackageStartupMessages({ library(duckdb); library(arrow) })

root <- rprojroot::find_root(rprojroot::has_file("config.yml"))
setwd(root)
source("R/harmonise.R")
cfg <- config::get(file = "config.yml")

con <- dbConnect(duckdb::duckdb(shared_home = FALSE))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
invisible(dbExecute(con, sprintf("SET memory_limit='%s'", cfg$duckdb_memory_limit)))

expr_glob <- file.path(cfg$store_dir, "expression", "**", "*.parquet")
samples_f <- file.path(cfg$store_dir, "samples.parquet")
samples <- as.data.frame(arrow::read_parquet(samples_f))

per_ds <- dbGetQuery(con, sprintf("
  SELECT dataset,
         count(*)                     AS n_rows,
         count(DISTINCT sample_id)    AS n_samples,
         count(DISTINCT gene_symbol)  AS n_genes
  FROM read_parquet('%s', hive_partitioning = 1)
  GROUP BY dataset ORDER BY dataset", expr_glob))

expected <- as.data.frame(table(samples$dataset), stringsAsFactors = FALSE)
names(expected) <- c("dataset", "n_samples_meta")
chk <- merge(per_ds, expected, by = "dataset", all = TRUE)
chk$samples_match <- chk$n_samples == chk$n_samples_meta
print(chk)

fails <- character(0)
if (any(!chk$samples_match, na.rm = TRUE) || anyNA(chk$samples_match)) {
  fails <- c(fails, "sample counts in expression store disagree with samples.parquet")
}

## Join on the full key. sample_id alone is not unique across datasets: RD, JR and TTC442
## are measured in both DepMap and Wachtel.
orphans <- dbGetQuery(con, sprintf("
  SELECT count(*) AS n FROM (
    SELECT DISTINCT dataset, sample_id FROM read_parquet('%s', hive_partitioning = 1)
  ) e LEFT JOIN read_parquet('%s') s USING (dataset, sample_id)
  WHERE s.sample_id IS NULL", expr_glob, samples_f))$n
if (orphans > 0) {
  fails <- c(fails, sprintf("%d sample_id(s) in expression have no metadata row", orphans))
}

bad_bucket <- dbGetQuery(con, sprintf("
  SELECT count(*) AS n FROM read_parquet('%s', hive_partitioning = 1)
  WHERE gene_bucket <> CASE WHEN regexp_matches(gene_symbol, '^[A-Za-z]')
                            THEN upper(substr(gene_symbol,1,1)) ELSE 'OTHER' END",
  expr_glob))$n
if (bad_bucket > 0) {
  fails <- c(fails, sprintf("%d row(s) in the wrong gene_bucket", bad_bucket))
}

if (length(fails)) {
  stop("STORE VALIDATION FAILED:\n  - ", paste(fails, collapse = "\n  - "), call. = FALSE)
}
message("\nStore validation PASSED")
