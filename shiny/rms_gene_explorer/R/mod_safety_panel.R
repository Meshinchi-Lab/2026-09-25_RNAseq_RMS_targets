## Normal Tissues panel — the on-target/off-tumour view.
##
## This is why all 30 adult GTEx tissues are kept in the store rather than the 9-tissue
## critical subset: expression in pancreas, thyroid, ovary or pituitary is a real CAR-T
## liability, and subsetting the STORE would hide it permanently. The scope toggle here is
## a display filter, so the full data is always one click away.
##
## Faceted by cohort, one x position per tissue, ordered by descending median so the
## highest-risk tissue reads first.

mod_safety_panel_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "d-flex align-items-center gap-3 mb-2",
        radioButtons(ns("scope"), NULL,
                     c("Critical tissues + muscle" = "critical",
                       "All tissues"               = "all"),
                     selected = "critical", inline = TRUE)),
    ggiraph::girafeOutput(ns("plot"), height = "520px")
  )
}

mod_safety_panel_server <- function(id, data, gene, opts) {
  moduleServer(id, function(input, output, session) {

    normals <- reactive({
      d <- data()
      d <- d[d$group %in% c("healthy_muscle", "normal_tissue") & !is.na(d$tissue), ,
             drop = FALSE]
      if (identical(input$scope, "critical")) {
        d <- d[d$tissue %in% CRITICAL_TISSUES, , drop = FALSE]
      }
      d
    })

    build <- function(interactive) {
      d <- normals()
      o <- opts()
      ## Tissue is both the facet-internal x and the colour here, so the ordering has to
      ## be computed from the data rather than a fixed vocabulary.
      ord <- if (nrow(d)) {
        names(sort(tapply(d$value, d$tissue, median, na.rm = TRUE), decreasing = TRUE))
      } else NULL
      plot_gene_distribution(d, gene(), group_by = "tissue",
                             scale = o$scale, style = o$style,
                             interactive = interactive, order_levels = ord)
    }


    output$plot <- ggiraph::renderGirafe({
      ## Widen the canvas when all ~30 tissues are on the x axis, or the labels collide.
      w <- if (identical(input$scope, "critical")) 10 else 15
      ggiraph::girafe(
        ggobj = build(TRUE), width_svg = w, height_svg = 5.2,
        options = list(
          ggiraph::opts_hover(css = "fill:#000000;stroke:none;r:2.5pt;"),
          ggiraph::opts_selection(type = "none"),
          ggiraph::opts_toolbar(saveaspng = FALSE)))
    })

    reactive(build(FALSE))
  })
}
