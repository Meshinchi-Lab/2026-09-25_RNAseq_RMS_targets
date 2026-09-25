## rms_gene_explorer — expression distribution of one gene across RMS-project cohorts.
##
## GOVERNANCE: PCGP is controlled-access under dbGaP phs000720 and the OPBG cohort is
## pre-publication. Run only on the analysis host or behind institutional SSO. Never
## deploy to shinyapps.io and never package with shinylive — shinylive ships the whole
## store into the viewer's browser.

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(ggplot2)
  library(ggiraph); library(dplyr); library(DT)
})

ROOT <- rprojroot::find_root(rprojroot::has_file("config.yml"))

source(file.path(ROOT, "R", "harmonise.R"))
source(file.path(ROOT, "R", "rms_palette.R"))
source(file.path(ROOT, "R", "plot_gene_distribution.R"))
for (f in list.files(file.path(ROOT, "shiny", "rms_gene_explorer", "R"),
                     pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}

CFG   <- config::get(file = file.path(ROOT, "config.yml"))
STORE <- store_open(file.path(ROOT, CFG$store_dir))
onStop(function() store_close(STORE))

DATASETS <- sort(unique(STORE$samples$dataset))
## `else` may not begin a line at top level in R, hence the braces.
DEFAULT_GENE <- {
  if ("CD276" %in% STORE$gene_index$gene_symbol) "CD276" else STORE$gene_index$gene_symbol[1]
}

ui <- page_sidebar(
  title = "RMS gene explorer",
  theme = bs_theme(version = 5, preset = "flatly"),

  sidebar = sidebar(
    width = 320,
    selectizeInput("gene", "Gene symbol", choices = NULL, multiple = FALSE,
                   options = list(placeholder = "e.g. CD276, MYOD1, GAPDH",
                                  maxOptions = 50)),
    checkboxGroupInput("datasets", "Datasets", choices = DATASETS, selected = DATASETS),
    radioButtons("scale", "Scale", c("log2(TPM+1)" = "log2", "Linear TPM" = "linear"),
                 "log2"),
    radioButtons("group_by", "Group by (RMS comparison tab)",
                 c("Class" = "group", "Fusion status" = "fusion",
                   "Subgroup" = "subgroup"), "group"),
    radioButtons("style", "Style",
                 c("Violin" = "violin", "Box" = "box", "Horizontal" = "ridge"), "violin"),
    hr(),
    downloadButton("dl_png", "Plot (PNG)", class = "btn-sm btn-outline-primary w-100"),
    div(class = "mt-2"),
    downloadButton("dl_tsv", "Data (TSV)", class = "btn-sm btn-outline-secondary w-100")
  ),

  navset_card_tab(
    id = "tab",
    nav_panel("RMS vs healthy muscle",  mod_rna_panel_ui("rna")),
    nav_panel("Normal Tissues",         mod_safety_panel_ui("safety")),
    nav_panel("Protein (Wachtel)",      mod_protein_panel_ui("protein"))
  ),

  card(card_header("Sample counts and summary statistics"),
       mod_summary_table_ui("summary"))
)

server <- function(input, output, session) {

  updateSelectizeInput(session, "gene", server = TRUE,
                       choices = STORE$gene_index$gene_symbol,
                       selected = DEFAULT_GENE)

  gene <- reactive(req(input$gene))

  gene_data <- reactive({
    query_gene(STORE, gene(), datasets = input$datasets)
  }) |> bindCache(gene(), input$datasets)

  opts <- reactive(list(scale = input$scale, style = input$style,
                        group_by = input$group_by))

  rna_static    <- mod_rna_panel_server("rna", gene_data, gene, opts)
  safety_static <- mod_safety_panel_server("safety", gene_data, gene, opts)
  prot_static   <- mod_protein_panel_server("protein", STORE, gene_data, gene, opts)
  mod_summary_table_server("summary", gene_data)

  ## The PNG follows whichever tab is in view, so the download always matches the figure
  ## the user is looking at.
  active_static <- reactive({
    switch(input$tab %||% "RMS vs healthy muscle",
           "RMS vs healthy muscle" = rna_static(),
           "Normal Tissues"        = safety_static(),
           "Protein (Wachtel)"     = prot_static(),
           rna_static())
  })

  output$dl_png <- downloadHandler(
    filename = function() {
      sprintf("%s_%s_%s.png", gene(),
              gsub("[^A-Za-z0-9]+", "_", tolower(input$tab %||% "rms")), Sys.Date())
    },
    content = function(file) {
      ggsave(file, plot = active_static(), width = 12, height = 5, dpi = 300, bg = "white")
    })

  output$dl_tsv <- downloadHandler(
    filename = function() sprintf("%s_expression_%s.tsv", gene(), Sys.Date()),
    content = function(file) {
      utils::write.table(gene_data(), file, sep = "\t", quote = FALSE, row.names = FALSE)
    })
}

shinyApp(ui, server)
