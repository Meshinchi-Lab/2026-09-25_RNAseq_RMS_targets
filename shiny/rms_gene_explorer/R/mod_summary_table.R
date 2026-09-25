## n / median / IQR per dataset x group. Always rendered — never behind a toggle.
## A violin over Wachtel's n=15 or GSE108022's five healthy controls is actively
## misleading without a visible n.

mod_summary_table_ui <- function(id) {
  ns <- NS(id)
  DT::DTOutput(ns("tbl"))
}

mod_summary_table_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    output$tbl <- DT::renderDT({
      s <- summarise_gene(data())
      if (!nrow(s)) {
        return(DT::datatable(s, rownames = FALSE, selection = "none",
                             options = list(dom = "t",
                                            language = list(emptyTable =
                                              "No measurements for this gene."))))
      }
      s$median <- signif(s$median, 4)
      s$q25    <- signif(s$q25, 4)
      s$q75    <- signif(s$q75, 4)
      ## Humanised labels so the table reads the same as the plot axes — "healthy muscle",
      ## not "healthy_muscle".
      s$group <- label_groups(s$group)
      DT::datatable(
        s, rownames = FALSE, selection = "none",
        colnames = c("Dataset", "Class", "Unit", "n", "Median", "Q25", "Q75"),
        options = list(pageLength = 20, dom = "tip", order = list(list(0, "asc"))))
    })
  })
}
