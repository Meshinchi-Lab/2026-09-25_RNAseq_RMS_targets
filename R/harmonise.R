## Controlled vocabularies and cohort-to-canonical mappings for the RMS gene store.
##
## Unlike `canonicalise_fusion()` in 2026-08-25_RMS_multi_cohort_expression_figures.qmd,
## which falls back to "unknown", every function here ERRORS on unmapped input. A silent
## fallback during a one-time ETL produces a store that is quietly wrong, and nothing
## downstream would ever surface it.

VALID_GROUPS <- c(
  ## RMS tumours
  "FP", "FN", "RMS_unspecified",
  ## Healthy reference. dGTEx / adult GTEx samples split on tissue: skeletal muscle is
  ## the comparator for an RMS panel, every other tissue belongs in the safety view.
  "healthy_muscle", "normal_tissue",
  ## DepMap cancer cell lines. RMS lines are split on fusion status so they line up with
  ## the tumour cohorts rather than collapsing 25 lines into one violin.
  "cell_line_RMS_FP", "cell_line_RMS_FN", "cell_line_other",
  ## Wachtel model classes — the panel is model systems, so the colour carries the
  ## model type and the tooltip carries fusion status.
  "FP_cellline", "FN_cellline", "PDX", "control"
)

## Groups belonging in the RMS-vs-healthy-muscle comparison panel. `normal_tissue` is
## deliberately absent: pooling 30 GTEx tissues into one "normal" violin would hide the
## fact that the comparator for RMS is skeletal muscle.
COMPARISON_GROUPS <- c("FP", "FN", "RMS_unspecified", "healthy_muscle",
                       "cell_line_RMS_FP", "cell_line_RMS_FN", "cell_line_other")

## Display labels. Kept separate from the storage vocabulary so the store never has to be
## rebuilt to reword an axis.
GROUP_LABELS <- c(
  FP                = "FP",
  FN                = "FN",
  RMS_unspecified   = "RMS (fusion unknown)",
  healthy_muscle    = "healthy muscle",
  normal_tissue     = "other normal tissue",
  cell_line_RMS_FP  = "RMS FP cell line",
  cell_line_RMS_FN  = "RMS FN cell line",
  cell_line_other   = "other cancer cell line",
  FP_cellline       = "FP cell line",
  FN_cellline       = "FN cell line",
  PDX               = "PDX",
  control           = "control"
)

VALID_ASSAYS <- c("pTPM", "log2_TPM1p", "iTop3")

ASSAY_UNITS <- c(pTPM       = "pTPM",
                 log2_TPM1p = "log2(TPM+1)",
                 iTop3      = "iTop3 abundance")

ASSAY_MODALITY <- c(pTPM = "RNA", log2_TPM1p = "RNA", iTop3 = "protein")

## Whether the stored values are ALREADY log-transformed. DepMap ships
## OmicsExpressionTPMLogp1 — applying log2(x+1) to it again would double-log the values
## and label the axis log2(log2(TPM+1)+1).
ASSAY_IS_LOG <- c(pTPM = FALSE, log2_TPM1p = TRUE, iTop3 = FALSE)

VALID_FUSION <- c("FP", "FN", "normal", "control", "unknown")

DATASET_ASSAY <- c(PCGP            = "pTPM",
                   dGTEx           = "pTPM",
                   GTEx_adult_v11  = "pTPM",
                   GSE108022       = "pTPM",
                   DepMap_26Q1     = "log2_TPM1p",
                   Wachtel_protein = "iTop3")

## The tissue label that counts as skeletal muscle in dGTEx and adult GTEx. Both use
## `Muscle` as the broad label (detail levels are Muscle - Quadriceps, - Gastrocnemius,
## - Abdominal, - Infrahyoid).
MUSCLE_TISSUE <- "Muscle"

## CAR-T safety-critical tissues, per the project convention in notes.md step 2 and the
## ped_critical_tissues / adult_critical_tissues vectors in
## 2026-08-25_RMS_multi_cohort_expression_figures.qmd. Skeletal muscle is included because
## it is the RMS lineage comparator, not because it is a safety concern.
##
## Spelled with spaces to match the store's tissue vocabulary ("Blood Vessel"), not the
## underscored form the notebook used. `Blood` resolves in adult GTEx (9 tissues) and is
## simply absent from dGTEx (8 tissues) — hence "blood and blood vessel where available".
CRITICAL_TISSUES <- c("Blood", "Blood Vessel", "Brain", "Heart", "Kidney",
                      "Liver", "Lung", "Muscle", "Nerve")

#' Abort unless every value is in `valid`.
assert_vocab <- function(x, valid, what) {
  bad <- setdiff(unique(x[!is.na(x)]), valid)
  if (length(bad)) {
    stop(sprintf("invalid %s value(s): %s\n  allowed: %s",
                 what, paste(bad, collapse = ", "), paste(valid, collapse = ", ")),
         call. = FALSE)
  }
  invisible(x)
}

.group_maps <- list(
  PCGP = c(FP = "FP", FN = "FN", UNK = "RMS_unspecified"),
  ## GSE108022's five non-RMS samples are skeletal muscle (NS127muscle .. NS134muscle),
  ## not a generic control.
  GSE108022 = c(FP = "FP", FN = "FN", muscle = "healthy_muscle"),
  Wachtel_protein = c(`FP-cellline` = "FP_cellline",
                      `FN-cellline` = "FN_cellline",
                      PDX           = "PDX",
                      Control       = "control")
)

#' Map a cohort's native group vocabulary onto VALID_GROUPS.
#'
#' @param x character vector of native values.
#' @param dataset one of names(DATASET_ASSAY).
#' @param tissue optional tissue vector, required for dGTEx / GTEx_adult_v11.
#' @return character vector of the same length, all in VALID_GROUPS.
canonicalise_group <- function(x, dataset, tissue = NULL) {
  stopifnot(length(dataset) == 1L)

  ## Normal-tissue cohorts split on tissue, not on any native group column.
  if (dataset %in% c("dGTEx", "GTEx_adult_v11")) {
    if (is.null(tissue)) {
      stop(sprintf("dataset '%s' requires `tissue` to split muscle from other tissues",
                   dataset), call. = FALSE)
    }
    return(ifelse(!is.na(tissue) & tissue == MUSCLE_TISSUE,
                  "healthy_muscle", "normal_tissue"))
  }

  ## DepMap is rule-based rather than a lookup: the disease_group vocabulary is open, so
  ## anything that is not an RMS line is "some other cancer line". RMS lines carry their
  ## fusion class in the same field (RMS-FP / RMS-FN), so the split comes for free.
  if (identical(dataset, "DepMap_26Q1")) {
    xu <- toupper(as.character(x))
    return(dplyr::case_when(
      grepl("^RMS", xu) & grepl("FP$", xu) ~ "cell_line_RMS_FP",
      grepl("^RMS", xu) & grepl("FN$", xu) ~ "cell_line_RMS_FN",
      TRUE                                 ~ "cell_line_other"))
  }

  map <- .group_maps[[dataset]]
  if (is.null(map)) {
    stop(sprintf("no group mapping defined for dataset '%s'", dataset), call. = FALSE)
  }

  out <- unname(map[as.character(x)])
  if (anyNA(out)) {
    stop(sprintf("unmapped %s group value(s): %s",
                 dataset, paste(unique(x[is.na(out)]), collapse = ", ")),
         call. = FALSE)
  }
  assert_vocab(out, VALID_GROUPS, "group")
  out
}

#' Fusion tag for tooltips, derived from the canonical group.
#'
#' The Wachtel panel colours by model class (FP cell line / FN cell line / PDX), so the
#' underlying fusion status has to travel separately to stay visible on hover. Wachtel's
#' three PDXs are all fusion-positive.
fusion_from_group <- function(group) {
  out <- dplyr::case_when(
    group %in% c("FP", "FP_cellline", "PDX", "cell_line_RMS_FP") ~ "FP",
    group %in% c("FN", "FN_cellline", "cell_line_RMS_FN")        ~ "FN",
    group %in% c("healthy_muscle", "normal_tissue") ~ "normal",
    group == "control"                            ~ "control",
    TRUE                                          ~ "unknown"
  )
  assert_vocab(out, VALID_FUSION, "fusion")
  out
}

#' Put stored values onto a common scale for plotting.
#'
#' `scale = "log2"` returns log2(TPM + 1) for every RNA assay: linear assays are
#' transformed, already-log assays are passed through untouched. `scale = "linear"`
#' back-transforms already-log assays with 2^x - 1. Without this, DepMap's
#' log2(TPM+1) values would be logged a second time and plotted against PCGP's raw pTPM
#' on a shared axis.
rescale_values <- function(value, assay, scale = c("log2", "linear")) {
  scale <- match.arg(scale)
  is_log <- unname(ASSAY_IS_LOG[assay])
  is_log[is.na(is_log)] <- FALSE
  if (identical(scale, "log2")) {
    ifelse(is_log, value, log2(value + 1))
  } else {
    ifelse(is_log, 2^value - 1, value)
  }
}

#' Axis label for a set of assays under a given scale.
scale_axis_label <- function(assays, scale = c("log2", "linear")) {
  scale  <- match.arg(scale)
  assays <- unique(assays[!is.na(assays)])
  if (!length(assays)) return("expression")

  ## A single non-RNA assay (protein) keeps its own units; rescale_values() leaves it
  ## alone on the linear path and log2s it on the log path.
  if (length(assays) == 1L && !identical(unname(ASSAY_MODALITY[assays]), "RNA")) {
    u <- unname(ASSAY_UNITS[assays])
    return(if (identical(scale, "log2")) paste0("log2(", u, " + 1)") else u)
  }

  if (identical(scale, "log2")) "log2(TPM + 1)" else "TPM"
}
