#!/usr/bin/env Rscript
## Stage 1 of the ETL. Melts each registered dataset into staging Parquet, one column
## block at a time, freeing each source object before moving to the next. Blocks are laid
## out by SAMPLE here, which is the wrong layout for querying — stage 2 repartitions by
## gene. Splitting the work this way is what keeps peak RAM bounded: melt and repartition
## have opposite memory profiles.
##
## Usage:
##   Rscript bin/01_build_store.R            # every registered dataset
##   Rscript bin/01_build_store.R GSE108022 Wachtel_protein

suppressPackageStartupMessages(library(arrow))

root <- rprojroot::find_root(rprojroot::has_file("config.yml"))
setwd(root)
source("R/harmonise.R")
source("R/loaders.R")

cfg <- config::get(file = "config.yml")
cfg$source_dir <- path.expand(cfg$source_dir)

args <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args)) args else names(LOADERS)
unknown <- setdiff(targets, names(LOADERS))
if (length(unknown)) stop("unknown dataset(s): ", paste(unknown, collapse = ", "))

dir.create(cfg$staging_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("meta_staging", showWarnings = FALSE, recursive = TRUE)

for (ds in targets) {
  message("=== ", ds, " ===")
  t0 <- Sys.time()

  ld <- LOADERS[[ds]](cfg)

  assert_vocab(ld$samples$group, VALID_GROUPS, "group")
  assert_vocab(ld$samples$assay, VALID_ASSAYS, "assay")
  stopifnot(setequal(names(ld$samples), SAMPLE_COLS),
            !anyDuplicated(ld$samples$sample_id))

  arrow::write_parquet(ld$samples, file.path("meta_staging", paste0(ds, ".parquet")))

  out_dir <- file.path(cfg$staging_dir, paste0("dataset=", ds))
  unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)

  i <- 1L; n_rows <- 0
  repeat {
    b <- ld$blocks(i)
    if (is.null(b)) break
    b$value <- as.numeric(b$value)
    arrow::write_parquet(
      b, file.path(out_dir, sprintf("block-%04d.parquet", i)), compression = "zstd")
    n_rows <- n_rows + nrow(b)
    message(sprintf("  block %3d  %12s rows  (running total %s)",
                    i, format(nrow(b), big.mark = ","), format(n_rows, big.mark = ",")))
    rm(b); gc(verbose = FALSE)
    i <- i + 1L
  }

  rm(ld); gc(verbose = FALSE)
  message(sprintf("  done: %s rows in %.1f min\n", format(n_rows, big.mark = ","),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

message("Stage 1 complete. Next: Rscript bin/02_repartition_store.R")
