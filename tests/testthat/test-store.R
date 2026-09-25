root <- rprojroot::find_root(rprojroot::has_file("config.yml"))
source(file.path(root, "shiny", "rms_gene_explorer", "R", "store.R"))
cfg <- config::get(file = file.path(root, "config.yml"))
store <- file.path(root, cfg$store_dir)

skip_if_no_store <- function() {
  skip_if_not(file.exists(file.path(store, "gene_index.parquet")),
              "store not built; run bin/01-03 first")
}

test_that("store_open fails fast with a helpful message when the store is absent", {
  expect_error(store_open("/nonexistent/store"), "bin/01_build_store.R")
})

test_that("gene_bucket matches the ETL's bucketing rule", {
  expect_equal(gene_bucket("GAPDH"), "G")
  expect_equal(gene_bucket("myod1"), "M")
  expect_equal(gene_bucket("7SK"), "OTHER")
})

test_that("store_open returns con, samples and gene_index", {
  skip_if_no_store()
  st <- store_open(store); on.exit(store_close(st))
  expect_true(all(c("con", "samples", "gene_index") %in% names(st)))
  expect_true(nrow(st$samples) > 0)
  expect_true(nrow(st$gene_index) > 0)
  expect_true(all(c("sample_id", "dataset", "assay", "unit", "group",
                    "tissue", "subgroup", "is_tumour") %in% names(st$samples)))
})

test_that("query_gene returns joined sample metadata", {
  skip_if_no_store()
  st <- store_open(store); on.exit(store_close(st))
  d <- query_gene(st, "GAPDH")
  expect_true(nrow(d) > 0)
  expect_true(all(c("gene_symbol","dataset","sample_id","value",
                    "group","tissue","subgroup","assay","unit") %in% names(d)))
  expect_false(anyNA(d$group))
  expect_equal(unique(d$gene_symbol), "GAPDH")
})

test_that("query_gene restricted to one dataset returns only that dataset", {
  skip_if_no_store()
  st <- store_open(store); on.exit(store_close(st))
  d <- query_gene(st, "GAPDH", datasets = "GSE108022")
  expect_equal(unique(d$dataset), "GSE108022")
})

test_that("an absent gene returns zero rows rather than erroring", {
  skip_if_no_store()
  st <- store_open(store); on.exit(store_close(st))
  expect_equal(nrow(query_gene(st, "NOTAREALGENE123")), 0L)
  expect_equal(nrow(query_gene(st, "")), 0L)
})

test_that("bucket pruning gives the same answer as an unpruned scan", {
  skip_if_no_store()
  st <- store_open(store); on.exit(store_close(st))
  pruned <- query_gene(st, "GAPDH")
  full <- DBI::dbGetQuery(st$con, sprintf(
    "SELECT count(*) AS n FROM read_parquet('%s', hive_partitioning=1)
     WHERE gene_symbol = 'GAPDH'", st$expr_glob))$n
  expect_equal(nrow(pruned), full)
})

test_that("summarise_gene reports one row per dataset x group with n", {
  skip_if_no_store()
  st <- store_open(store); on.exit(store_close(st))
  d <- query_gene(st, "GAPDH")
  s <- summarise_gene(d)
  expect_true(all(c("dataset", "group", "unit", "n", "median", "q25", "q75") %in% names(s)))
  expect_equal(sum(s$n), nrow(d))
  expect_true(all(s$q25 <= s$median & s$median <= s$q75))
})

test_that("summarise_gene handles an empty frame", {
  expect_equal(nrow(summarise_gene(.empty_gene_frame())), 0L)
})

test_that("the sample key is (dataset, sample_id), not sample_id alone", {
  skip_if_no_store()
  st <- store_open(store); on.exit(store_close(st))
  ## RD, JR and TTC442 are RMS cell lines measured in BOTH DepMap and Wachtel. That
  ## overlap is the cross-modality comparison the app exists to show, so sample_id is
  ## unique only within a dataset. A join on sample_id alone would fan out these rows.
  expect_false(anyDuplicated(st$samples[c("dataset", "sample_id")]) > 0)
  shared <- intersect(st$samples$sample_id[st$samples$dataset == "DepMap_26Q1"],
                      st$samples$sample_id[st$samples$dataset == "Wachtel_protein"])
  skip_if(length(shared) == 0, "no cross-modality cell lines in this store")
  d <- query_gene(st, "CD276")
  expect_equal(nrow(d), nrow(unique(d[c("dataset", "sample_id")])))
})

test_that("query_gene returns every dataset that measures a ubiquitous gene", {
  skip_if_no_store()
  st <- store_open(store); on.exit(store_close(st))
  d <- query_gene(st, "GAPDH")
  expect_setequal(unique(d$dataset), unique(st$samples$dataset))
})
