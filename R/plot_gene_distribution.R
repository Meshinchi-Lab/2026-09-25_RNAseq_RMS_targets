## The single plotting contract shared by the Shiny app and
## 2026-08-25_RMS_multi_cohort_expression_figures.qmd. Both callers use this function so
## that app figures and manuscript figures cannot drift apart.
##
## Requires R/rms_palette.R and R/harmonise.R to be sourced first.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

## Minimum samples in a group before a violin or box is drawn for it. Below this the
## points and the median crossbar are shown alone: a density over one or two observations
## is not a distribution, and rendering it as a hairline makes a present group look absent.
MIN_N_DISTRIBUTION <- 3L

#' Colours for the canonical group vocabulary, falling back to the project sequential
#' ramp for free-text values such as GTEx tissue names.
#'
#' `amber` is the project's designated third-nominal-category hue, which is exactly what
#' PDX is alongside FP and FN model systems. Grey is reserved for healthy reference.
rms_group_palette <- function(levels) {
  base <- c(FP              = .rms_hue("blue_dk"),
            FP_cellline     = .rms_hue("blue_dk"),
            FN              = .rms_hue("teal"),
            FN_cellline     = .rms_hue("teal"),
            PDX             = .rms_hue("amber"),
            RMS_unspecified  = .rms_hue("blue_lt"),
            ## DepMap RMS lines share the fusion hues so an FP line reads as FP across
            ## every panel; the mid-blue / green pairing keeps them distinguishable from
            ## the primary tumour cohorts at a glance.
            cell_line_RMS_FP = .rms_hue("blue_md"),
            cell_line_RMS_FN = .rms_hue("green_lt"),
            cell_line_other  = .rms_hue("grey_lt"),
            healthy_muscle  = .rms_hue("grey"),
            normal_tissue   = .rms_hue("grey_lt"),
            control         = .rms_hue("grey_dk"))

  ## Labels reach the scale already humanised, so key on the label too.
  names(base) <- unname(ifelse(names(base) %in% names(GROUP_LABELS),
                               GROUP_LABELS[names(base)], names(base)))

  known <- intersect(levels, names(base))
  other <- setdiff(levels, names(base))
  out <- base[known]
  if (length(other)) {
    out <- c(out, stats::setNames(rms_seq_colours(length(other)), other))
  }
  out
}

#' Humanise a group vector using GROUP_LABELS, leaving free text untouched.
label_groups <- function(x) {
  x <- as.character(x)
  ifelse(x %in% names(GROUP_LABELS), unname(GROUP_LABELS[x]), x)
}

#' Distribution of one gene across cohorts.
#'
#' @param df    data.frame with gene_symbol, dataset, sample_id, value, assay, unit,
#'              group, fusion, tissue, subgroup.
#' @param gene  character(1), used in the plot title.
#' @param facet_by   column to facet columns by.
#' @param group_by   column mapped to x and fill.
#' @param scale      "log2" or "linear"; applied per-assay by rescale_values().
#' @param style      "violin", "box", or "ridge" (ridge renders as a flipped boxplot;
#'                   ggridges is deliberately not a dependency).
#' @param interactive TRUE swaps the jitter layer for ggiraph::geom_point_interactive().
#' @param order_levels optional character vector fixing the x-axis order (raw group
#'                   vocabulary, not labels).
#' @return a ggplot object.
plot_gene_distribution <- function(df,
                                   gene,
                                   facet_by    = c("dataset", "tissue"),
                                   group_by    = c("group", "tissue", "subgroup", "fusion"),
                                   scale       = c("log2", "linear"),
                                   style       = c("violin", "box", "ridge"),
                                   interactive = FALSE,
                                   order_levels = NULL) {
  facet_by <- match.arg(facet_by)
  group_by <- match.arg(group_by)
  scale    <- match.arg(scale)
  style    <- match.arg(style)

  ylab <- scale_axis_label(if (nrow(df)) df$assay else character(0), scale)
  if (!nrow(df)) return(.placeholder_plot(gene, ylab))

  d <- df |>
    mutate(
      .facet     = .data[[facet_by]],
      .raw       = as.character(.data[[group_by]]),
      .tip       = .build_tooltip(df, group_by, ylab),
      .detected  = is.finite(.data$value)
    ) |>
    filter(!is.na(.data$.facet))

  if (!nrow(d) || !any(d$.detected)) return(.placeholder_plot(gene, ylab))

  ## Non-detections are kept as explicit observations rather than dropped, following the
  ## convention in 2026-08-06_RMS_protein_validation.qmd (chunk `wachtel-boxplot`): an NA
  ## iTop3 means "below the MS detection limit", which is evidence of low protein, not a
  ## missing-at-random measurement. Dropping them makes a protein detected in 2 of 6 FN
  ## lines read as "high in FN" when it is absent in 4 of 6.
  ##
  ## Sentinel = 30% of the minimum DETECTED value, computed on the raw scale so the
  ## per-assay transform below places it correctly on either axis. That puts the
  ## non-detections in a visually distinct band instead of colliding with real low values.
  min_detected <- min(d$value[d$.detected], na.rm = TRUE)
  na_sentinel  <- min_detected * 0.30
  separator    <- min_detected * 0.60
  d$value[!d$.detected] <- na_sentinel

  d <- d |>
    mutate(
      ## Per-assay transform: linear assays get log2(x+1), already-log assays pass
      ## through. Without this DepMap's log2(TPM+1) would be logged a second time.
      .value = rescale_values(.data$value, .data$assay, scale)
    ) |>
    filter(is.finite(.data$.value))

  d$.tip[!d$.detected] <- paste0(sub("\n.*$", "", d$.tip[!d$.detected]),
                                 "\nnot detected (below MS floor)")

  ## Per-group n counts DETECTED samples only, and drives whether a distribution mark is
  ## drawn. A violin over Wachtel's 2 controls renders as a hairline and reads as a
  ## missing group; the points are the honest mark at that n.
  d <- d |>
    group_by(.data$.facet, .data$.raw) |>
    mutate(.n_group = sum(.data$.detected)) |>
    ungroup()

  if (!nrow(d)) return(.placeholder_plot(gene, ylab))
  has_undetected <- any(!d$.detected)
  separator_y <- rescale_values(separator, d$assay[1], scale)

  lev <- if (!is.null(order_levels)) {
    label_groups(intersect(order_levels, unique(d$.raw)))
  } else {
    sort(unique(label_groups(d$.raw)))
  }
  d$.group <- factor(label_groups(d$.raw), levels = lev)
  d <- d[!is.na(d$.group), , drop = FALSE]
  if (!nrow(d)) return(.placeholder_plot(gene, ylab))

  pal <- rms_group_palette(levels(d$.group))

  p <- ggplot(d, aes(x = .data$.group, y = .data$.value, fill = .data$.group))

  ## Distribution marks and the median crossbar summarise DETECTED values only, and are
  ## drawn only for groups with enough of them to have a shape.
  big  <- function(x) x[x$.n_group >= MIN_N_DISTRIBUTION & x$.detected, , drop = FALSE]
  seen <- function(x) x[x$.detected, , drop = FALSE]

  p <- p + switch(
    style,
    violin = geom_violin(data = big, scale = "width", trim = TRUE,
                         colour = NA, alpha = 0.55),
    box    = geom_boxplot(data = big, outlier.shape = NA, alpha = 0.55, width = 0.6),
    ridge  = geom_boxplot(data = big, outlier.shape = NA, alpha = 0.55, width = 0.6)
  )

  ## Small panels (the 15-sample Wachtel set) need visible points, since the points ARE
  ## the plot there rather than an overlay on a density.
  pt_size  <- if (nrow(d) <= 60) 2.1 else 0.8
  pt_alpha <- if (nrow(d) <= 60) 0.95 else 0.7
  jit <- position_jitter(width = 0.18, height = 0, seed = 1L)

  ## Where some samples were not detected, shape carries detection status and the group
  ## colour moves to `fill` — pch 21 is a filled circle with an outline, pch 2 an open
  ## triangle that ignores fill, matching the project rule "grey = missing value".
  ## Panels with nothing missing keep the plain dark points, so the RNA tabs are unchanged.
  p <- p + if (has_undetected) {
    if (interactive) {
      ggiraph::geom_point_interactive(
        aes(tooltip = .data$.tip, data_id = .data$sample_id, shape = .data$.detected),
        position = jit, size = pt_size, alpha = pt_alpha,
        colour = "grey15", stroke = 0.45, show.legend = FALSE)
    } else {
      geom_point(aes(shape = .data$.detected), position = jit,
                 size = pt_size, alpha = pt_alpha,
                 colour = "grey15", stroke = 0.45, show.legend = FALSE)
    }
  } else {
    if (interactive) {
      ggiraph::geom_point_interactive(
        aes(tooltip = .data$.tip, data_id = .data$sample_id),
        position = jit, size = pt_size, alpha = pt_alpha,
        colour = "grey20", show.legend = FALSE)
    } else {
      geom_point(position = jit, size = pt_size, alpha = pt_alpha,
                 colour = "grey20", show.legend = FALSE)
    }
  }

  if (has_undetected) {
    p <- p +
      scale_shape_manual(values = c(`TRUE` = 21L, `FALSE` = 2L), guide = "none") +
      ## Separates the detected zone from the not-detected band.
      geom_hline(yintercept = separator_y, linetype = "dotted",
                 colour = "grey40", linewidth = 0.35) +
      annotate("text", x = Inf, y = rescale_values(na_sentinel, d$assay[1], scale),
               label = "not detected\n(below MS floor)", hjust = 1.05, vjust = 0.5,
               size = 2.6, colour = "grey40", fontface = "italic")
  }

  p <- p +
    stat_summary(data = seen, fun = median, geom = "crossbar", width = 0.65,
                 linewidth = 0.35, colour = "black", show.legend = FALSE) +
    facet_grid(cols = vars(.data$.facet), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal, na.value = RMS_NA_COLOUR, name = NULL) +
    labs(title = gene, x = NULL, y = ylab) +
    theme_rms_base(base_size = 10) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1),
          strip.text      = element_text(face = "bold"),
          legend.position = "bottom")

  if (identical(style, "ridge")) p <- p + coord_flip()

  p
}

#' Tooltip text. Always carries fusion status when it is known and is not already the
#' grouping variable — the Wachtel panel colours by model class, so FP/FN would otherwise
#' be invisible.
.build_tooltip <- function(df, group_by, ylab) {
  grp <- label_groups(df[[group_by]])
  base <- sprintf("%s\n%s: %.3g\n%s", df$sample_id, df$unit, df$value, grp)

  show_fusion <- !identical(group_by, "fusion") &&
    "fusion" %in% names(df) &&
    any(df$fusion %in% c("FP", "FN"), na.rm = TRUE)
  if (show_fusion) {
    base <- paste0(base, "\nfusion: ", ifelse(is.na(df$fusion), "unknown", df$fusion))
  }

  has_tissue <- !identical(group_by, "tissue") && "tissue" %in% names(df) &&
    any(!is.na(df$tissue))
  if (has_tissue) {
    base <- ifelse(is.na(df$tissue), base, paste0(base, "\ntissue: ", df$tissue))
  }
  base
}

.placeholder_plot <- function(gene, ylab) {
  ggplot() +
    annotate("text", x = 0, y = 0,
             label = paste0(gene, " not found in the selected datasets"),
             size = 4, colour = .rms_hue("grey")) +
    labs(title = gene, y = ylab, x = NULL) +
    theme_void()
}
