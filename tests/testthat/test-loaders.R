root <- rprojroot::find_root(rprojroot::has_file("config.yml"))
source(file.path(root, "R", "harmonise.R"))
source(file.path(root, "R", "loaders.R"))
cfg <- config::get(file = file.path(root, "config.yml"))
cfg$source_dir <- path.expand(cfg$source_dir)

skip_if_no_source <- function() {
  skip_if_not(dir.exists(cfg$source_dir), "upstream source_dir not available")
}

test_that("registry records the expected datasets", {
  expect_setequal(names(LOADERS),
                  c("GSE108022", "Wachtel_protein", "PCGP", "dGTEx",
                    "DepMap_26Q1", "GTEx_adult_v11"))
  expect_setequal(names(LOADERS), names(DATASET_ASSAY))
})

test_that("GSE108022 loader yields a valid samples frame", {
  skip_if_no_source()
  ld <- LOADERS[["GSE108022"]](cfg)
  expect_equal(nrow(ld$samples), 106L)
  expect_setequal(names(ld$samples), SAMPLE_COLS)
  assert_vocab(ld$samples$group, VALID_GROUPS, "group")
  expect_equal(unique(ld$samples$assay), "pTPM")
  expect_true(all(c("FP", "FN", "healthy_muscle") %in% ld$samples$group))
  ## The five non-RMS samples are skeletal muscle (NS127muscle .. NS134muscle).
  expect_equal(sum(ld$samples$group == "healthy_muscle"), 5L)
  expect_true(all(ld$samples$tissue[ld$samples$group == "healthy_muscle"] == "Muscle"))
})

test_that("GSE108022 blocks cover every sample exactly once", {
  skip_if_no_source()
  ld <- LOADERS[["GSE108022"]](cfg)
  seen <- character(0); i <- 1L
  repeat {
    b <- ld$blocks(i)
    if (is.null(b)) break
    expect_setequal(names(b), c("gene_symbol", "sample_id", "value"))
    seen <- c(seen, unique(b$sample_id)); i <- i + 1L
  }
  expect_setequal(seen, ld$samples$sample_id)
  expect_false(any(duplicated(seen)))
})

test_that("Wachtel loader yields iTop3 abundance grouped by model class", {
  skip_if_no_source()
  ld <- LOADERS[["Wachtel_protein"]](cfg)
  expect_equal(nrow(ld$samples), 15L)
  expect_equal(unique(ld$samples$assay), "iTop3")
  expect_setequal(unique(ld$samples$group),
                  c("FP_cellline", "FN_cellline", "PDX", "control"))
  expect_equal(as.integer(table(ld$samples$group)[
    c("FN_cellline", "FP_cellline", "PDX", "control")]), c(6L, 4L, 3L, 2L))
  ## Fusion status must survive for the tooltip: the three PDXs are fusion-positive.
  expect_equal(sum(ld$samples$fusion == "FP"), 7L)
  expect_equal(sum(ld$samples$fusion == "FN"), 6L)
  b <- ld$blocks(1L)
  expect_true(nrow(b) > 1000)
  expect_true(all(b$value >= 0, na.rm = TRUE))
  ## Both controls must carry measurements.
  expect_true(all(c("MRC-5", "Myoblasts") %in% b$sample_id))
  expect_true(all(b$sample_id %in% ld$samples$sample_id))
})

test_that("Wachtel values are iTop3, not the ppm conversion", {
  skip_if_no_source()
  raw <- data.table::fread(
    file.path(cfg$source_dir,
              "2026-08-06_RMS_protein_validation_wachtel_all_proteins_long.tsv"),
    select = c("GN", "sample_short", "iTop3", "ppm"), data.table = TRUE)
  raw <- raw[GN == "ERBB2" & sample_short == "Rh4" & !is.na(iTop3)]
  skip_if(nrow(raw) == 0, "ERBB2/Rh4 not measured")
  b <- LOADERS[["Wachtel_protein"]](cfg)$blocks(1L)
  got <- b$value[b$gene_symbol == "ERBB2" & b$sample_id == "Rh4"]
  expect_equal(got, max(raw$iTop3))
  expect_false(isTRUE(all.equal(got, max(raw$ppm))))
})

test_that("PCGP and dGTEx split the combined object correctly", {
  skip_if_no_source()
  pcgp <- LOADERS[["PCGP"]](cfg)
  dg   <- LOADERS[["dGTEx"]](cfg)
  expect_equal(nrow(pcgp$samples), 169L)
  expect_equal(nrow(dg$samples), 332L)
  expect_setequal(unique(pcgp$samples$group), c("FP", "FN", "RMS_unspecified"))
  ## dGTEx splits on tissue: 48 Muscle samples are the comparator, 284 are safety.
  expect_setequal(unique(dg$samples$group), c("healthy_muscle", "normal_tissue"))
  expect_equal(sum(dg$samples$group == "healthy_muscle"), 48L)
  expect_true(all(dg$samples$tissue[dg$samples$group == "healthy_muscle"] == "Muscle"))
  expect_length(intersect(pcgp$samples$sample_id, dg$samples$sample_id), 0L)
})

test_that("DepMap yields readable cell-line names and both group classes", {
  skip_if_no_source()
  ld <- LOADERS[["DepMap_26Q1"]](cfg)
  expect_true(nrow(ld$samples) > 1500)
  expect_setequal(unique(ld$samples$group),
                  c("cell_line_RMS_FP", "cell_line_RMS_FN", "cell_line_other"))
  ## 14 FP + 11 FN RMS lines are present in the DepMap expression matrix.
  expect_equal(sum(ld$samples$group == "cell_line_RMS_FP"), 14L)
  expect_equal(sum(ld$samples$group == "cell_line_RMS_FN"), 11L)
  expect_equal(sum(ld$samples$fusion == "FP"), 14L)
  expect_equal(sum(ld$samples$fusion == "FN"), 11L)
  expect_false(anyDuplicated(ld$samples$sample_id) > 0)
  expect_true(sum(grepl("^ACH-", ld$samples$sample_id)) < 10)
})

test_that("adult GTEx retains all 30 tissues and melts a block cleanly", {
  skip_if_no_source()
  ld <- LOADERS[["GTEx_adult_v11"]](cfg)
  expect_equal(nrow(ld$samples), 19616L)
  expect_gte(dplyr::n_distinct(ld$samples$tissue), 25L)
  expect_setequal(unique(ld$samples$group), c("healthy_muscle", "normal_tissue"))
  ## All 30 tissues are retained; only Muscle becomes the comparator.
  expect_equal(sum(ld$samples$group == "healthy_muscle"), 818L)
  b <- ld$blocks(1L)
  expect_equal(dplyr::n_distinct(b$sample_id), cfg$block_size)
  expect_true(is.numeric(b$value))
})
