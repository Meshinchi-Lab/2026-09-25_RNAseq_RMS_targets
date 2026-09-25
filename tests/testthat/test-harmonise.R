source(file.path(rprojroot::find_root(rprojroot::has_file("config.yml")), "R", "harmonise.R"))
suppressPackageStartupMessages(library(dplyr))

test_that("vocabularies are closed sets with no duplicates", {
  expect_false(any(duplicated(VALID_GROUPS)))
  expect_setequal(names(ASSAY_UNITS), VALID_ASSAYS)
  expect_setequal(names(ASSAY_MODALITY), VALID_ASSAYS)
  expect_setequal(names(ASSAY_IS_LOG), VALID_ASSAYS)
  expect_setequal(unname(ASSAY_MODALITY), c("RNA", "protein"))
  expect_setequal(names(GROUP_LABELS), VALID_GROUPS)
  expect_true(all(COMPARISON_GROUPS %in% VALID_GROUPS))
})

test_that("normal_tissue is excluded from the RMS comparison panel", {
  ## Pooling 30 GTEx tissues into one "normal" violin would hide that the comparator
  ## for RMS is skeletal muscle.
  expect_false("normal_tissue" %in% COMPARISON_GROUPS)
  expect_true("healthy_muscle" %in% COMPARISON_GROUPS)
})

test_that("PCGP fusion vocabulary maps onto VALID_GROUPS", {
  expect_equal(canonicalise_group(c("FP", "FN", "UNK"), "PCGP"),
               c("FP", "FN", "RMS_unspecified"))
})

test_that("normal cohorts split on tissue, muscle versus everything else", {
  tis <- c("Muscle", "Brain", "Muscle", "Pancreas")
  expect_equal(canonicalise_group(rep(NA_character_, 4), "dGTEx", tissue = tis),
               c("healthy_muscle", "normal_tissue", "healthy_muscle", "normal_tissue"))
  expect_equal(canonicalise_group(rep(NA_character_, 2), "GTEx_adult_v11",
                                  tissue = c("Muscle", "Thyroid")),
               c("healthy_muscle", "normal_tissue"))
})

test_that("a normal cohort without tissue is an error, not a silent pooling", {
  expect_error(canonicalise_group(rep(NA_character_, 2), "dGTEx"), "requires `tissue`")
})

test_that("GSE108022 non-RMS samples are healthy muscle, not generic controls", {
  expect_equal(canonicalise_group(c("FP", "FN", "muscle"), "GSE108022"),
               c("FP", "FN", "healthy_muscle"))
})

test_that("DepMap splits RMS lines on fusion status, not just RMS vs other", {
  expect_equal(canonicalise_group(c("RMS-FP", "RMS-FN", "adult_nonsarc", "osteosarcoma"),
                                  "DepMap_26Q1"),
               c("cell_line_RMS_FP", "cell_line_RMS_FN",
                 "cell_line_other", "cell_line_other"))
})

test_that("DepMap RMS cell lines carry their fusion status into the tooltip", {
  expect_equal(fusion_from_group(c("cell_line_RMS_FP", "cell_line_RMS_FN",
                                   "cell_line_other")),
               c("FP", "FN", "unknown"))
})

test_that("both DepMap RMS classes reach the comparison panel", {
  expect_true(all(c("cell_line_RMS_FP", "cell_line_RMS_FN") %in% COMPARISON_GROUPS))
  expect_false("cell_line_RMS" %in% VALID_GROUPS)
})

test_that("critical tissues are spelled to match the store vocabulary", {
  expect_length(CRITICAL_TISSUES, 9L)
  ## Spaces, not the notebook's underscores — "Blood_Vessel" would silently match nothing.
  expect_true("Blood Vessel" %in% CRITICAL_TISSUES)
  expect_false(any(grepl("_", CRITICAL_TISSUES)))
  ## Muscle is included: it is the RMS lineage comparator.
  expect_true(MUSCLE_TISSUE %in% CRITICAL_TISSUES)
  expect_true("Blood" %in% CRITICAL_TISSUES)
})

test_that("Wachtel maps to model class, not fusion status", {
  expect_equal(canonicalise_group(c("PDX", "FP-cellline", "FN-cellline", "Control"),
                                  "Wachtel_protein"),
               c("PDX", "FP_cellline", "FN_cellline", "control"))
})

test_that("fusion tag recovers FP/FN from Wachtel model classes for the tooltip", {
  expect_equal(fusion_from_group(c("FP_cellline", "FN_cellline", "PDX", "control")),
               c("FP", "FN", "FP", "control"))
  expect_equal(fusion_from_group(c("FP", "FN", "healthy_muscle", "normal_tissue")),
               c("FP", "FN", "normal", "normal"))
  expect_equal(fusion_from_group("cell_line_other"), "unknown")
})

test_that("unmapped values are a hard error, never a silent NA", {
  expect_error(canonicalise_group("WAT", "PCGP"), "unmapped")
  expect_error(canonicalise_group("FP", "NoSuchDataset"), "no group mapping")
})

test_that("assert_vocab rejects out-of-vocabulary values", {
  expect_silent(assert_vocab(c("FP", "healthy_muscle"), VALID_GROUPS, "group"))
  expect_error(assert_vocab("nonsense", VALID_GROUPS, "group"), "group")
})

## ------------------------------------------------------------------ scale handling

test_that("already-log assays are never logged a second time", {
  ## The bug this guards: DepMap ships log2(TPM+1); applying log2(x+1) again produced
  ## values on a log2(log2(TPM+1)+1) axis.
  v <- c(10, 10)
  out <- rescale_values(v, c("pTPM", "log2_TPM1p"), "log2")
  expect_equal(out[1], log2(11))
  expect_equal(out[2], 10)
})

test_that("linear scale back-transforms already-log assays", {
  out <- rescale_values(c(10, 10), c("pTPM", "log2_TPM1p"), "linear")
  expect_equal(out[1], 10)
  expect_equal(out[2], 2^10 - 1)
})

test_that("the two scales agree: log2 output is log2(linear output + 1)", {
  ## rescale_values always consumes the STORED value, so the invariant is between the
  ## two scale settings, not a composition of the function with itself.
  v <- c(5.5, 5.5, 0, 123.4)
  a <- c("pTPM", "log2_TPM1p", "pTPM", "log2_TPM1p")
  expect_equal(rescale_values(v, a, "log2"),
               log2(rescale_values(v, a, "linear") + 1))
})

test_that("axis label never nests log2 inside log2", {
  expect_equal(scale_axis_label(c("pTPM", "log2_TPM1p"), "log2"), "log2(TPM + 1)")
  expect_equal(scale_axis_label(c("pTPM", "log2_TPM1p"), "linear"), "TPM")
  expect_equal(scale_axis_label("log2_TPM1p", "log2"), "log2(TPM + 1)")
  expect_false(grepl("log2\\(log2", scale_axis_label("log2_TPM1p", "log2")))
})

test_that("a protein-only selection keeps its own units", {
  expect_equal(scale_axis_label("iTop3", "linear"), "iTop3 abundance")
  expect_equal(scale_axis_label("iTop3", "log2"), "log2(iTop3 abundance + 1)")
})

test_that("Wachtel is stored as iTop3, not ppm", {
  expect_equal(unname(DATASET_ASSAY[["Wachtel_protein"]]), "iTop3")
  expect_false(any(grepl("ppm", VALID_ASSAYS)))
  expect_equal(unname(ASSAY_UNITS[["iTop3"]]), "iTop3 abundance")
})
