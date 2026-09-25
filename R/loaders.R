## Loader registry. Each loader is `function(cfg)` returning:
##   samples : data.frame with exactly SAMPLE_COLS
##   blocks  : function(i) -> data.frame(gene_symbol, sample_id, value), or NULL past the end
##
## Blocking exists because the adult GTEx source is 18,683 x 19,616; melting it in one pass
## would need tens of GB. Every loader is blocked for uniformity even when it need not be.
##
## Requires R/harmonise.R to be sourced first.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

SAMPLE_COLS <- c("sample_id", "dataset", "assay", "unit",
                 "group", "fusion", "tissue", "subgroup", "is_tumour")

LOADERS <- list()

register_loader <- function(name, fn) {
  LOADERS[[name]] <<- fn
  invisible(NULL)
}

#' Column extractor appropriate to the source object.
#'
#' Sources are a mix of true matrices (GSE108022, DepMap) and data.frames (adult GTEx,
#' PCGP/dGTEx). `[[` works only on the data.frame kind and `[, j]` is wasteful on it, so
#' the right accessor is chosen once per source rather than per column.
.column_getter <- function(mat) {
  if (is.data.frame(mat)) function(j) mat[[j]] else function(j) mat[, j]
}

#' Melt a gene x sample source into long form for one column block.
#'
#' `as.numeric(unlist(...))` covers the adult-GTEx case, where the source is a data.frame
#' whose columns each carry a redundant 18,683-element `names` attribute. It strips that
#' without the whole-object copy `as.matrix()` would cost.
.melt_block <- function(getter, cols, genes) {
  vals <- lapply(cols, function(j) as.numeric(unlist(getter(j), use.names = FALSE)))
  dt <- data.table::data.table(
    gene_symbol = rep(genes, times = length(cols)),
    sample_id   = rep(cols, each = length(genes)),
    value       = unlist(vals, use.names = FALSE)
  )
  dt[!is.na(value)]
}

#' Build a blocks() closure over a wide source.
.blocker <- function(mat, sample_ids, block_size) {
  genes <- rownames(mat)
  if (is.null(genes)) stop("source matrix has no rownames (gene symbols)", call. = FALSE)
  getter <- .column_getter(mat)
  starts <- seq(1L, length(sample_ids), by = block_size)
  function(i) {
    if (i > length(starts)) return(NULL)
    from <- starts[i]
    to   <- min(from + block_size - 1L, length(sample_ids))
    as.data.frame(.melt_block(getter, sample_ids[from:to], genes))
  }
}

.blank_samples <- function(sample_id, dataset, group, tissue, subgroup, is_tumour) {
  assay <- unname(DATASET_ASSAY[[dataset]])
  group <- as.character(group)
  data.frame(
    sample_id = as.character(sample_id),
    dataset   = dataset,
    assay     = assay,
    unit      = unname(ASSAY_UNITS[[assay]]),
    group     = group,
    fusion    = fusion_from_group(group),
    tissue    = as.character(tissue),
    subgroup  = as.character(subgroup),
    is_tumour = as.logical(is_tumour),
    stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------- GSE108022 (Shern 2014)
register_loader("GSE108022", function(cfg) {
  f <- file.path(cfg$source_dir, "2026-08-17_RNAseq_RMS_GSE108022_EDA_pTPM.RDS")
  o <- readRDS(f)
  meta <- as.data.frame(o$meta)

  ## fusion_group is FP/FN for tumours. The five non-RMS samples are skeletal muscle —
  ## their sample_ids are literally NS127muscle .. NS134muscle — so they map to
  ## healthy_muscle, not a generic control.
  native <- ifelse(!meta$is_rms, "muscle", as.character(meta$fusion_group))
  if (anyNA(native)) {
    stop("GSE108022: NA fusion_group on an is_rms sample", call. = FALSE)
  }

  samples <- .blank_samples(
    sample_id = meta$sample_id,
    dataset   = "GSE108022",
    group     = canonicalise_group(native, "GSE108022"),
    tissue    = ifelse(meta$is_rms, NA_character_, MUSCLE_TISSUE),
    subgroup  = meta$histology,
    is_tumour = meta$is_rms)

  list(samples = samples,
       blocks  = .blocker(o$mat, colnames(o$mat), cfg$block_size))
})

## ------------------------------------------------------- Wachtel 2023 RMS surfaceome (MS)
register_loader("Wachtel_protein", function(cfg) {
  expr_f <- file.path(cfg$source_dir,
                      "2026-08-06_RMS_protein_validation_wachtel_all_proteins_long.tsv")
  meta_f <- file.path(cfg$source_dir,
                      "2026-08-06_RMS_protein_validation_wachtel_sample_metadata.tsv")

  meta <- data.table::fread(meta_f, data.table = FALSE)

  ## This panel is model systems, so `group` carries the model class (FP cell line /
  ## FN cell line / PDX / control) and drives the colour. Fusion status travels in
  ## `fusion` for the tooltip — Wachtel's three PDXs are all fusion-positive.
  samples <- .blank_samples(
    sample_id = meta$sample_short,
    dataset   = "Wachtel_protein",
    group     = canonicalise_group(meta$sample_class, "Wachtel_protein"),
    tissue    = NA_character_,
    subgroup  = meta$sample_type,
    is_tumour = meta$group == "RMS")

  ## Already long. Values are iTop3 abundance as published, NOT the ppm conversion —
  ## nothing here is compared against another protein source, so the ppm normalisation
  ## would only add a transform with no consumer. Collapse protein isoforms to the gene
  ## level by max iTop3 and serve as a single block.
  long <- data.table::fread(expr_f, select = c("GN", "sample_short", "iTop3"),
                            data.table = TRUE)
  long <- long[!is.na(iTop3) & !is.na(GN) & GN != ""]
  long <- long[, .(value = max(iTop3, na.rm = TRUE)),
               by = .(gene_symbol = GN, sample_id = sample_short)]
  long <- as.data.frame(long[, .(gene_symbol, sample_id, value)])

  list(samples = samples,
       blocks  = function(i) if (i == 1L) long else NULL)
})

## ------------------------------------------------ PCGP + dGTEx (one object, two cohorts)
.read_combined <- local({
  cache <- NULL
  function(cfg) {
    if (is.null(cache)) {
      f <- file.path(cfg$source_dir,
                     "2026-07-10_RNAseq_RMS_target_discovery_combined_pTPM.RDS")
      cache <<- readRDS(f)
    }
    cache
  }
})

.combined_loader <- function(which_cohort) {
  force(which_cohort)
  function(cfg) {
    o <- .read_combined(cfg)
    meta <- as.data.frame(o$meta)
    meta <- meta[meta$cohort == which_cohort, , drop = FALSE]

    native <- if (identical(which_cohort, "PCGP")) as.character(meta$fusion_positive)
              else rep(NA_character_, nrow(meta))

    ## dGTEx splits on tissue: its 48 Muscle samples are the healthy-muscle comparator,
    ## the other 284 across 23 tissues are the safety reference.
    samples <- .blank_samples(
      sample_id = meta$sample_id,
      dataset   = which_cohort,
      group     = canonicalise_group(native, which_cohort, tissue = meta$tissue_broad),
      tissue    = meta$tissue_broad,
      ## PCGP's informative secondary axis is Diagnosis/Relapse/Xenograft; dGTEx's is the
      ## detailed tissue (Muscle - Quadriceps vs - Gastrocnemius).
      subgroup  = if (identical(which_cohort, "PCGP")) meta$sample_type
                  else meta$tissue_detail,
      is_tumour = identical(which_cohort, "PCGP"))

    list(samples = samples,
         blocks  = .blocker(o$mat, meta$sample_id, cfg$block_size))
  }
}

register_loader("PCGP",  .combined_loader("PCGP"))
register_loader("dGTEx", .combined_loader("dGTEx"))

## ------------------------------------------------------------------ DepMap Public 26Q1
register_loader("DepMap_26Q1", function(cfg) {
  cache_f <- file.path(cfg$source_dir,
                       "2026-09-14_depmap_expression_and_dependency_depmap_26Q1_cache.rds")
  coh_f   <- file.path(cfg$source_dir,
                       "2026-09-14_depmap_expression_and_dependency_cohort_assignment.tsv")

  o <- readRDS(cache_f)
  expr <- o$expr_log2tpm
  coh  <- data.table::fread(coh_f, data.table = FALSE)

  ids <- intersect(colnames(expr), coh$ModelID)
  coh <- coh[match(ids, coh$ModelID), , drop = FALSE]

  samples <- .blank_samples(
    sample_id = ids,
    dataset   = "DepMap_26Q1",
    group     = canonicalise_group(coh$disease_group, "DepMap_26Q1"),
    tissue    = NA_character_,
    subgroup  = coh$OncotreeLineage,
    is_tumour = TRUE)

  ## Label cell lines by their stripped name where available — ACH- IDs are unreadable in
  ## a tooltip. Fall back to the ACH- id if the name is missing or would collide.
  nice <- coh$StrippedCellLineName
  nice <- ifelse(is.na(nice) | !nzchar(nice), ids, nice)
  nice[duplicated(nice)] <- ids[duplicated(nice)]
  samples$sample_id <- nice
  stopifnot(!anyDuplicated(samples$sample_id))

  rename_map <- stats::setNames(samples$sample_id, ids)
  inner <- .blocker(expr, ids, cfg$block_size)

  list(samples = samples,
       blocks = function(i) {
         b <- inner(i)
         if (is.null(b)) return(NULL)
         b$sample_id <- unname(rename_map[b$sample_id])
         b
       })
})

## ------------------------------------------------------- Adult GTEx v11 (all 30 tissues)
## `o$mat` is an 18,683 x 19,616 DATA.FRAME, not a matrix: is.matrix() is FALSE and every
## column carries a redundant 18,683-element names attribute (hence 3.4 GB on disk, ~6 GB
## resident). .melt_block()'s as.numeric(unlist(...)) strips that per block. Do NOT call
## as.matrix() on the whole object — it doubles peak RAM for no benefit.
register_loader("GTEx_adult_v11", function(cfg) {
  f <- file.path(cfg$source_dir,
                 "2026-07-10_RNAseq_RMS_target_discovery_adult_gtex_v11_pTPM.RDS")
  o <- readRDS(f)
  meta <- as.data.frame(o$meta)

  ## Splits on tissue: the 818 Muscle samples are the healthy-muscle comparator, the other
  ## 18,798 across 29 tissues are the on-target/off-tumour safety reference.
  samples <- .blank_samples(
    sample_id = meta$sample_id,
    dataset   = "GTEx_adult_v11",
    group     = canonicalise_group(rep(NA_character_, nrow(meta)), "GTEx_adult_v11",
                                   tissue = meta$tissue_broad_adult),
    tissue    = meta$tissue_broad_adult,
    subgroup  = meta$tissue_detail_adult,
    is_tumour = FALSE)

  list(samples = samples,
       blocks  = .blocker(o$mat, meta$sample_id, cfg$block_size))
})
