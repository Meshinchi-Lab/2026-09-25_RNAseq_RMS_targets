## RMS-vs-healthy-muscle comparison panel.
##
## Only COMPARISON_GROUPS are shown: dGTEx and adult GTEx contribute their skeletal-muscle
## samples ONLY (48 and 818 respectively). Pooling all 24 / 30 tissues into one "normal"
## violin would hide the fact that the developmental comparator for RMS is skeletal muscle,
## and would drag the normal median toward whatever tissue happens to dominate the cohort.
## The other tissues are the safety question, and live in the normal-tissue panel.

## Left-to-right x order within each facet: tumour classes, then cell lines, then the
## healthy-muscle reference last so the comparator is always the rightmost column.
RNA_GROUP_ORDER <- c("FP", "FN", "RMS_unspecified",
                     "cell_line_RMS_FP", "cell_line_RMS_FN", "cell_line_other",
                     "healthy_muscle")

mod_rna_panel_ui <- function(id) {
  ns <- NS(id)
  ggiraph::girafeOutput(ns("plot"), height = "440px")
}

#' @param data reactive returning the query_gene() frame, already dataset-filtered.
#' @param gene reactive returning the selected symbol.
#' @param opts reactive list(scale=, style=, group_by=).
#' @return reactive returning the STATIC twin of the on-screen plot, for downloads.
mod_rna_panel_server <- function(id, data, gene, opts) {
  moduleServer(id, function(input, output, session) {

    rna <- reactive({
      d <- data()
      d <- d[!is.na(d$assay) & ASSAY_MODALITY[d$assay] == "RNA", , drop = FALSE]
      d[d$group %in% COMPARISON_GROUPS, , drop = FALSE]
    })

    build <- function(interactive) {
      o <- opts()
      plot_gene_distribution(rna(), gene(), group_by = o$group_by,
                             scale = o$scale, style = o$style,
                             interactive = interactive,
                             order_levels = RNA_GROUP_ORDER)
    }


    output$plot <- ggiraph::renderGirafe({
      ggiraph::girafe(
        ggobj = build(TRUE), width_svg = 11, height_svg = 4.4,
        options = list(
          ggiraph::opts_hover(css = "fill:#000000;stroke:none;r:2.5pt;"),
          ggiraph::opts_selection(type = "none"),
          ggiraph::opts_toolbar(saveaspng = FALSE)))
    })

    reactive(build(FALSE))
  })
}
