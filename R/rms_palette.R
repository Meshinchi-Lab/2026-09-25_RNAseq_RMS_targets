## VENDORED — do not edit here.
## Source: ~/github_repos/2025-12-04_RNAseq_RMS_infant/R/rms_palette.R
## Copied: 2026-09-25
## This is the canonical RMS colour system shared by all eight analysis notebooks.
## Keeping a byte-identical copy is what makes app figures match manuscript figures.
## To update: re-copy from upstream, do not hand-patch.

# =============================================================================
# rms_palette.R — single canonical colour system for the RMS CAR-T project
# -----------------------------------------------------------------------------
# Sourced by every analysis notebook:
#   2025-12-04_RNAseq_RMS_infant_rnaseq_count_nf_EDA.qmd
#   2026-06-03_RNAseq_RMS_PCGP_GTEX_EDA.qmd
#   2026-07-10_RNAseq_RMS_target_discovery.qmd
#   2026-08-06_RMS_RNAseq_microarray_validation.qmd
#   2026-08-06_RMS_protein_validation.qmd
#   2026-08-17_RNAseq_RMS_GSE108022_EDA.qmd
#   2026-08-25_RMS_multi_cohort_expression_figures.qmd
#   2026-09-14_depmap_expression_and_dependency.qmd
#
# DESIGN CONSTRAINT (2026-09-22): blues + greens, yellow only where a third
# nominal category genuinely requires a third hue, and GREY reserved for
# controls / healthy reference / unknown / missing. No reds or oranges —
# the previous Okabe-Ito usage across the repo put the tumour classes in
# warm hues and the healthy reference in green, which inverted the intended
# "grey = reference" reading.
#
# Why this hue set works:
#   * Blue → green → yellow is the viridis axis, so the categorical palette
#     and the sequential ramps live in the same perceptual family.
#   * Blue ↔ yellow is the canonical colour-vision-deficiency-safe diverging
#     pair (deuteranopia / protanopia both preserve it), so the z-score ramp
#     stays readable for ~8 % of male readers who lose red-green contrast.
#   * Reserving grey for reference/missing means a reader can find the
#     control arm in any panel without consulting the legend.
#
# ALIAS HANDLING. The repo carries three incompatible fusion-status
# vocabularies (`FP` / `fusion_positive` / `PAX3-FOXO1`) and several cohort
# spellings. Every accessor below normalises its input to a canonical class
# and returns a vector named with the CALLER's original level strings, so a
# notebook keeps its local vocabulary while the colour is shared.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Core hue ladder
# -----------------------------------------------------------------------------
# Ordered dark-blue -> blue -> teal -> green -> light-green -> amber, plus a
# grey ladder for reference / unknown / missing. Every categorical palette
# below is drawn from these constants so hues never drift between notebooks.
RMS_HUES <- c(
  blue_dk    = "#08306B",   # darkest blue  — primary/discovery class
  blue       = "#08519C",   # strong blue
  blue_md    = "#2171B5",   # mid blue
  blue_lt    = "#6BAED6",   # light blue
  teal       = "#41B6C4",   # cyan-teal     — secondary tumour class
  green_dk   = "#1B4332",   # dark green
  green      = "#238B45",   # green
  green_lt   = "#78C679",   # light green
  green_pale = "#A1DAB4",   # pale green
  amber      = "#FEC44F",   # yellow/amber  — third-category hue only
  amber_dk   = "#8C6D00",   # dark amber    — diverging high pole
  amber_lt   = "#FEE391",   # pale amber
  ink        = "#111111",   # near-black    — dense categorical (cell lines)
  grey_dk    = "#4D4D4D",   # charcoal      — established / known reference
  grey       = "#969696",   # mid grey      — CONTROL / healthy reference
  grey_lt    = "#BDBDBD",   # light grey    — secondary reference
  grey_pale  = "#D9D9D9",   # pale grey     — unknown
  grey_na    = "#E5E5E5"    # near-white    — missing / not measured
)

# Canonical NA / missing fill used by every scale and heatmap in the project.
RMS_NA_COLOUR <- unname(RMS_HUES["grey_na"])

.rms_hue <- function(key) unname(RMS_HUES[key])


# -----------------------------------------------------------------------------
# 1. Fusion status  (FP / FN / normal / unknown)
# -----------------------------------------------------------------------------
# The two tumour classes carry saturated hues; healthy reference and unknown
# are grey so the control arm is identifiable without the legend.
RMS_FUSION_COLOURS <- c(
  FP      = .rms_hue("blue_dk"),
  FN      = .rms_hue("teal"),
  normal  = .rms_hue("grey"),
  unknown = .rms_hue("grey_pale")
)

#' Normalise any fusion-status spelling used anywhere in the repo to one of
#' {FP, FN, normal, unknown}.
rms_fusion_class <- function(x) {
  xu <- toupper(trimws(as.character(x)))
  xu[is.na(xu) | !nzchar(xu)] <- "UNKNOWN"
  out <- rep("unknown", length(xu))
  out[grepl("^(FP|FUSION[_ -]?POSITIVE|POSITIVE)$", xu)]                 <- "FP"
  out[grepl("PAX[37]\\s*(::|--|-|_)\\s*FOXO1", xu)]                      <- "FP"
  out[grepl("PAX3\\s*(::|--|-|_)\\s*(NCOA1|INO80D)", xu)]                <- "FP"
  out[grepl("^(FN|FUSION[_ -]?NEGATIVE|NEGATIVE|FUSIONNEG|NONE)$", xu)]  <- "FN"
  out[grepl("^(NORMAL|HEALTHY([_ ]MUSCLE)?|CONTROL|GTEX|DGTEX|MUSCLE|SKELETAL[_ ]MUSCLE|NON[_ -]?CANCEROUS)$", xu)] <- "normal"
  out[grepl("^(UNKNOWN|NOT[_ ]AVAILABLE|UNKNOWN[_ ]SUBTYPE|NA|OTHER)$", xu)] <- "unknown"
  out
}

#' Named colour vector for a caller-supplied set of fusion levels.
#' Returns a vector named with the ORIGINAL level strings.
rms_fusion_pal <- function(levels) {
  levels <- as.character(levels)
  stats::setNames(unname(RMS_FUSION_COLOURS[rms_fusion_class(levels)]), levels)
}

#' Ordered factor over the canonical 4 classes.
rms_fusion_factor <- function(x) {
  factor(rms_fusion_class(x), levels = names(RMS_FUSION_COLOURS))
}


# -----------------------------------------------------------------------------
# 2. Subcellular localisation  (plasma_membrane / secreted / intracellular)
# -----------------------------------------------------------------------------
# PRECEDENCE (fixed 2026-09-22, per project convention):
#     plasma_membrane  >  secreted  >  intracellular
#
# Rationale. A CAR-T construct requires a surface-accessible epitope, so a
# gene annotated as plasma-membrane is classified there even when it is ALSO
# annotated as secreted — 4 of the 37 candidates carry both flags (FGFR4,
# PTK7, EPHA8, NTN5-class receptors with shed ectodomains). Those genes are
# still membrane targets; the shed form is a separate clinical consideration
# (soluble-antigen decoy sink) captured by `rms_is_shed()` below, NOT by
# demoting them out of the plasma-membrane class.
#
# `secreted` outranks `intracellular` because a secreted antigen is at least
# antibody-accessible, whereas a purely intracellular antigen requires a
# TCR-T / peptide-MHC modality entirely outside the CAR-T scope of this study.
RMS_LOCALIZATION_COLOURS <- c(
  plasma_membrane = .rms_hue("blue"),    # best CAR-T tractability
  secreted        = .rms_hue("green"),   # antibody-accessible, shed risk
  intracellular   = .rms_hue("amber")    # needs TCR-T / peptide-MHC
)

RMS_LOCALIZATION_LABELS <- c(
  plasma_membrane = "Plasma membrane",
  secreted        = "Secreted",
  intracellular   = "Intracellular"
)

#' Derive the 3-level localisation factor from the HPA/COMPARTMENTS flags.
#'
#' @param is_plasma_membrane logical (or coercible) vector
#' @param is_secreted        logical (or coercible) vector
#' @return ordered factor, levels plasma_membrane > secreted > intracellular
rms_localization <- function(is_plasma_membrane, is_secreted) {
  pm  <- as.logical(is_plasma_membrane) %in% TRUE
  sec <- as.logical(is_secreted)        %in% TRUE
  factor(
    ifelse(pm, "plasma_membrane", ifelse(sec, "secreted", "intracellular")),
    levels = names(RMS_LOCALIZATION_COLOURS)
  )
}

#' TRUE for genes that are BOTH plasma-membrane and secreted — i.e. membrane
#' targets with a shed/soluble form that can act as a CAR-T decoy sink.
#' Reported alongside localisation; does not change the localisation class.
rms_is_shed <- function(is_plasma_membrane, is_secreted) {
  (as.logical(is_plasma_membrane) %in% TRUE) & (as.logical(is_secreted) %in% TRUE)
}

rms_localization_pal <- function(labelled = FALSE) {
  p <- RMS_LOCALIZATION_COLOURS
  if (labelled) names(p) <- RMS_LOCALIZATION_LABELS[names(p)]
  p
}

# Shapes carry localisation on scatter/volcano panels where fill is already
# spent on another variable. Stroke-capable pch codes only (21-25).
RMS_LOCALIZATION_SHAPES <- c(
  plasma_membrane = 21L,   # filled circle
  secreted        = 23L,   # filled diamond
  intracellular   = 22L    # filled square
)


# -----------------------------------------------------------------------------
# 3. Target class  (novel vs known)
# -----------------------------------------------------------------------------
# Known targets are the pipeline's positive controls, so they take the
# charcoal "reference" colour; novel candidates — the actual finding — take
# the saturated teal.
RMS_TARGET_CLASS_COLOURS <- c(
  novel = .rms_hue("teal"),
  known = .rms_hue("grey_dk")
)

RMS_TARGET_CLASS_LABELS <- c(
  novel = "Novel candidate",
  known = "Known target"
)

# Shape encoding used on every volcano / effect-size scatter in the project.
RMS_TARGET_CLASS_SHAPES <- c(novel = 21L, known = 24L)   # circle vs triangle-up

#' Derive the 2-level target-class factor. Missing -> "novel" (conservative:
#' unproven is treated as novel, never as an established target).
rms_target_class <- function(is_known_target) {
  factor(ifelse(as.logical(is_known_target) %in% TRUE, "known", "novel"),
         levels = names(RMS_TARGET_CLASS_COLOURS))
}

rms_target_class_pal <- function(labelled = FALSE) {
  p <- RMS_TARGET_CLASS_COLOURS
  if (labelled) names(p) <- RMS_TARGET_CLASS_LABELS[names(p)]
  p
}


# -----------------------------------------------------------------------------
# 4. Cohort / dataset
# -----------------------------------------------------------------------------
# Tumour cohorts run dark-blue -> green -> amber in pipeline order; the two
# healthy-tissue references (dGTEx, adult GTEx) are grey, matching the
# "grey = reference" rule used for fusion status.
RMS_COHORT_COLOURS <- c(
  PCGP           = .rms_hue("blue_dk"),
  OPBG           = .rms_hue("blue_md"),
  GSE108022      = .rms_hue("green"),
  R2_microarray  = .rms_hue("green_lt"),
  DepMap         = .rms_hue("amber"),
  dGTEx          = .rms_hue("grey"),
  adult_GTEx     = .rms_hue("grey_lt"),
  unknown        = .rms_hue("grey_pale")
)

#' Normalise cohort spellings (`R2 microarray`, `R2_microarray`, `AdultGTEx`,
#' `PCGP+dGTEx`, ...) to the canonical keys above.
rms_cohort_key <- function(x) {
  xs <- trimws(as.character(x))
  xu <- toupper(gsub("[^A-Za-z0-9]+", "_", xs))
  out <- rep("unknown", length(xu))
  out[grepl("^PCGP",            xu)] <- "PCGP"
  out[grepl("^OPBG|PAZIENTI",   xu)] <- "OPBG"
  out[grepl("GSE108022|SHERN",  xu)] <- "GSE108022"
  out[grepl("^R2|MICROARRAY",   xu)] <- "R2_microarray"
  out[grepl("DEPMAP|CCLE",      xu)] <- "DepMap"
  out[grepl("ADULT.*GTEX|GTEX.*ADULT", xu)] <- "adult_GTEx"
  out[grepl("^D_?GTEX|DEVELOPMENTAL_GTEX", xu)] <- "dGTEx"
  # `PCGP+dGTEx` names the combined discovery matrix -> read as PCGP.
  out[grepl("^PCGP", xu)] <- "PCGP"
  out
}

rms_cohort_pal <- function(levels) {
  levels <- as.character(levels)
  stats::setNames(unname(RMS_COHORT_COLOURS[rms_cohort_key(levels)]), levels)
}


# -----------------------------------------------------------------------------
# 5. DepMap disease group  /  proteomics cell-line class
# -----------------------------------------------------------------------------
RMS_DISEASE_GROUP_COLOURS <- c(
  `RMS-FP`            = .rms_hue("blue_dk"),
  `RMS-FN`            = .rms_hue("teal"),
  `RMS-unknown`       = .rms_hue("green_pale"),
  `osteosarcoma`      = .rms_hue("green"),
  `other_sarcoma`     = .rms_hue("green_lt"),
  `pediatric_nonsarc` = .rms_hue("grey"),      # reference cohort
  `adult_nonsarc`     = .rms_hue("grey_lt")    # de-emphasised
)

# Wachtel iTop3 sample classes.
#
# Revised 2026-09-22 for separability. The previous assignment (dark blue
# #08306B / near-black #111111 / dark green #1B4332) put all three tumour-model
# classes at the dark end of the ramp: worst-pair CIE-Lab dE was 28.3 and the
# worst-pair greyscale luminance gap was 0.053, i.e. FP- and FN-cellline were
# essentially the same colour once printed or viewed at 2 px point size.
#
# The set below was chosen by maximising the WORST pairwise dE subject to a
# greyscale-separation floor, searched over the palette ladder:
#     worst-pair dE      28.3 -> 54.9
#     worst-pair |dLum| 0.053 -> 0.150
#
# This is the one place the project's "yellow only if absolutely needed" rule
# is invoked: three mutually-distinguishable tumour-model classes plus a grey
# control cannot be built from blues and greens alone without spending two
# slots on near-identical greens (best yellow-free option scored dE 42.0 /
# |dLum| 0.101 using two greens).
#
# Assignment is semantic, not arbitrary:
#   FP-cellline -> the canonical FP navy, so fusion-positive reads the same
#                  here as in every other figure in the project
#   FN-cellline -> green
#   PDX         -> amber, the brightest: PDX is the odd-one-out *model type*
#                  (in-vivo patient-derived, not a cultured line)
#   Control     -> grey, per the project-wide "grey = reference" rule
RMS_SAMPLE_CLASS_COLOURS <- c(
  `PDX`         = .rms_hue("amber"),
  `FP-cellline` = .rms_hue("blue_dk"),
  `FN-cellline` = .rms_hue("green"),
  `Control`     = .rms_hue("grey")
)

# Frejno/Kuster + CCLE RPPA cell-line classes.
RMS_CELL_LINE_COLOURS <- c(
  RMS           = .rms_hue("blue_dk"),
  other_sarcoma = .rms_hue("green_lt"),
  non_cancerous = .rms_hue("grey"),
  other         = .rms_hue("grey_pale")
)


# -----------------------------------------------------------------------------
# 6. Binary flags / provenance
# -----------------------------------------------------------------------------
# Heatmap row-annotation flags (surface, secreted, TF, known, in-trial).
RMS_FLAG_COLOURS <- c(`TRUE` = "#252525", `FALSE` = RMS_NA_COLOUR)

# Validation vote outcome.
RMS_VALIDATED_COLOURS <- c(`TRUE` = .rms_hue("blue"), `FALSE` = .rms_hue("grey_lt"))

# PASS / FAIL for the filter-cascade alluvial.
RMS_PASSFAIL_COLOURS <- c(PASS = .rms_hue("blue"), FAIL = .rms_hue("grey_lt"))

# Why a gene entered the shortlist.
RMS_CANDIDATE_SOURCE_COLOURS <- c(
  transcriptome                    = .rms_hue("blue"),
  clinical_trial                   = .rms_hue("green"),
  `transcriptome + clinical_trial` = .rms_hue("blue_dk"),
  unknown                          = .rms_hue("grey_pale")
)


# -----------------------------------------------------------------------------
# 7. Continuous ramps
# -----------------------------------------------------------------------------
# Sequential = ColorBrewer YlGnBu, i.e. pale-yellow -> green -> teal -> deep
# blue. Chosen over viridis because viridis begins at dark PURPLE (#440154),
# which falls outside the blue/green/yellow constraint this palette is built
# on. YlGnBu is perceptually near-uniform, colour-vision-deficiency-safe, and
# prints legibly in greyscale (monotonic luminance).
RMS_SEQ_ANCHORS <- c("#FFFFD9", "#C7E9B4", "#7FCDBB", "#41B6C4",
                     "#2171B5", "#08306B")

rms_seq_colours <- function(n = 6) {
  grDevices::colorRampPalette(RMS_SEQ_ANCHORS)(n)
}

#' Diverging ramp for row-z / log2FC panels.
#' Blue (low) <-> near-white (0) <-> amber (high). Blue-yellow is the
#' colour-vision-deficiency-safe diverging pair; the previous blue-red ramp
#' was not.
rms_diverging_ramp <- function(limit = 2) {
  circlize::colorRamp2(
    breaks = c(-limit, -limit/2, 0, limit/2, limit),
    colors = c(.rms_hue("blue_dk"), .rms_hue("blue_lt"), "#F7F7F7",
               .rms_hue("amber"),   .rms_hue("amber_dk"))
  )
}

#' Sequential ramp for expression / safety panels (0 .. ceiling).
rms_sequential_ramp <- function(ceiling = 5, n = 6) {
  circlize::colorRamp2(seq(0, ceiling, length.out = n), rms_seq_colours(n))
}

#' CRISPR gene-effect ramp. Essential (negative) reads dark blue; neutral
#' reads near-white; positive (growth-suppressing) reads pale amber.
rms_essentiality_ramp <- function() {
  circlize::colorRamp2(
    breaks = c(-2.0, -1.0, -0.5, 0.0, 0.5),
    colors = c(.rms_hue("blue_dk"), .rms_hue("blue_md"), .rms_hue("blue_lt"),
               "#F0F0F0", .rms_hue("amber_lt"))
  )
}


# -----------------------------------------------------------------------------
# 8. ggplot2 convenience scales
# -----------------------------------------------------------------------------
scale_fill_rms_localization <- function(name = "Localisation", labelled = TRUE, ...) {
  p <- rms_localization_pal(labelled = FALSE)
  ggplot2::scale_fill_manual(
    name = name, values = p, drop = FALSE,
    labels = if (labelled) RMS_LOCALIZATION_LABELS[names(p)] else ggplot2::waiver(),
    na.value = RMS_NA_COLOUR, ...)
}

scale_colour_rms_localization <- function(name = "Localisation", labelled = TRUE, ...) {
  p <- rms_localization_pal(labelled = FALSE)
  ggplot2::scale_colour_manual(
    name = name, values = p, drop = FALSE,
    labels = if (labelled) RMS_LOCALIZATION_LABELS[names(p)] else ggplot2::waiver(),
    na.value = RMS_NA_COLOUR, ...)
}

scale_shape_rms_target_class <- function(name = "Target class", labelled = TRUE, ...) {
  s <- RMS_TARGET_CLASS_SHAPES
  ggplot2::scale_shape_manual(
    name = name, values = s, drop = FALSE,
    labels = if (labelled) RMS_TARGET_CLASS_LABELS[names(s)] else ggplot2::waiver(), ...)
}

scale_fill_rms_target_class <- function(name = "Target class", labelled = TRUE, ...) {
  p <- rms_target_class_pal(labelled = FALSE)
  ggplot2::scale_fill_manual(
    name = name, values = p, drop = FALSE,
    labels = if (labelled) RMS_TARGET_CLASS_LABELS[names(p)] else ggplot2::waiver(),
    na.value = RMS_NA_COLOUR, ...)
}

scale_fill_rms_fusion <- function(name = "Fusion status", ...) {
  ggplot2::scale_fill_manual(name = name, values = RMS_FUSION_COLOURS,
                             drop = FALSE, na.value = RMS_NA_COLOUR, ...)
}

scale_colour_rms_fusion <- function(name = "Fusion status", ...) {
  ggplot2::scale_colour_manual(name = name, values = RMS_FUSION_COLOURS,
                               drop = FALSE, na.value = RMS_NA_COLOUR, ...)
}


# -----------------------------------------------------------------------------
# 9. Shared ggplot theme
# -----------------------------------------------------------------------------
theme_rms_base <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      axis.text        = ggplot2::element_text(colour = "black"),
      legend.position  = "bottom",
      plot.title       = ggplot2::element_text(face = "bold")
    )
}


# -----------------------------------------------------------------------------
# 10. Standard row annotation for candidate heatmaps
# -----------------------------------------------------------------------------
#' Build the project-standard ComplexHeatmap row annotation carrying the two
#' axes every primary figure must show: target class (novel vs known) and
#' subcellular localisation (plasma membrane / secreted / intracellular).
#'
#' @param gene_ids ordered character vector of gene symbols (heatmap row order)
#' @param cand     candidate table containing gene_id, is_known_target,
#'                 is_plasma_membrane, is_secreted
#' @param which    "row" (default) or "column"
#' @param extra    optional named list of additional annotation vectors
#' @param extra_col optional named list of colour maps for `extra`
rms_candidate_annotation <- function(gene_ids, cand, which = "row",
                                      extra = NULL, extra_col = NULL,
                                      show_legend = TRUE) {
  stopifnot(requireNamespace("ComplexHeatmap", quietly = TRUE))
  idx <- match(as.character(gene_ids), as.character(cand$gene_id))
  loc <- rms_localization(cand$is_plasma_membrane[idx], cand$is_secreted[idx])
  tcl <- rms_target_class(cand$is_known_target[idx])

  args <- list(
    `Target class`  = as.character(tcl),
    `Localisation`  = as.character(loc)
  )
  cols <- list(
    `Target class` = RMS_TARGET_CLASS_COLOURS,
    `Localisation` = RMS_LOCALIZATION_COLOURS
  )
  if (!is.null(extra)) {
    args <- c(args, extra)
    if (!is.null(extra_col)) cols <- c(cols, extra_col)
  }
  do.call(ComplexHeatmap::HeatmapAnnotation, c(
    args,
    list(which              = which,
         col                = cols,
         na_col             = RMS_NA_COLOUR,
         annotation_name_gp = grid::gpar(fontsize = 7),
         simple_anno_size   = grid::unit(3, "mm"),
         show_legend        = show_legend)
  ))
}


# -----------------------------------------------------------------------------
# 11. Self-check
# -----------------------------------------------------------------------------
# Every colour must be a valid hex; every categorical palette must have
# unique colours within itself (so no two levels are indistinguishable).
local({
  all_pals <- list(
    RMS_HUES, RMS_FUSION_COLOURS, RMS_LOCALIZATION_COLOURS,
    RMS_TARGET_CLASS_COLOURS, RMS_COHORT_COLOURS,
    RMS_DISEASE_GROUP_COLOURS, RMS_SAMPLE_CLASS_COLOURS,
    RMS_CELL_LINE_COLOURS, RMS_CANDIDATE_SOURCE_COLOURS
  )
  for (p in all_pals) {
    stopifnot(all(grepl("^#[0-9A-Fa-f]{6}$", p)))
    if (anyDuplicated(p)) {
      warning("rms_palette: duplicated colours within a palette: ",
              paste(names(p)[duplicated(p)], collapse = ", "))
    }
  }
})

message(sprintf(
  "rms_palette.R loaded — %d core hues; fusion/localisation/target-class/cohort scales available.",
  length(RMS_HUES)))
