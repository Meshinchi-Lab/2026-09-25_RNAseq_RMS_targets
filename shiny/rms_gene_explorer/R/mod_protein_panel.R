## Protein modality panel (Wachtel 2023 RMS surfaceome).
##
## A SEPARATE panel with its own y-axis: iTop3 abundance must never be read off an RNA
## pTPM axis. Values are iTop3 as published, not the ppm conversion — nothing here is
## compared against another protein source.
##
## Coloured by model class (FP cell line / FN cell line / PDX / control) because that is
## the axis of variation in a model-system panel; fusion status rides in the tooltip.

PROTEIN_GROUP_ORDER <- c("FP_cellline", "FN_cellline", "PDX", "control")

mod_protein_panel_ui <- function(id) {
  ns <- NS(id)
  ggiraph::girafeOutput(ns("plot"), height = "360px")
}

mod_protein_panel_server <- function(id, store, data, gene, opts) {
  moduleServer(id, function(input, output, session) {

    prot <- reactive({
      d <- data()
      d <- d[!is.na(d$assay) & ASSAY_MODALITY[d$assay] == "protein", , drop = FALSE]
      ## Re-attach the samples the MS did not detect, so all 15 model systems appear and
      ## the detection rate is visible rather than silently absent from the denominator.
      add_undetected(store, d, "Wachtel_protein")
    })

    build <- function(interactive) {
      o <- opts()
      ## group_by is forced to the model class: the sidebar's fusion/tissue/subgroup
      ## choices are meaningless for a 15-sample model panel with no tissue.
      plot_gene_distribution(prot(), gene(), group_by = "group",
                             scale = o$scale, style = o$style,
                             interactive = interactive,
                             order_levels = PROTEIN_GROUP_ORDER)
    }


    output$plot <- ggiraph::renderGirafe({
      ggiraph::girafe(
        ggobj = build(TRUE), width_svg = 9, height_svg = 3.6,
        options = list(
          ggiraph::opts_hover(css = "fill:#000000;stroke:none;r:2.5pt;"),
          ggiraph::opts_selection(type = "none"),
          ggiraph::opts_toolbar(saveaspng = FALSE)))
    })

    reactive(build(FALSE))
  })
}
