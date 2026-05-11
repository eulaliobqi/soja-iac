# ============================================================
# RNASeq Insight Dashboard – Glycine max
# Execute com: Rscript -e "shiny::runApp('dashboard/app.R', port=3838)"
# ============================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(ggplot2)
library(dplyr)
library(readr)
library(pheatmap)
library(scales)

# ── Configuração: diretório de resultados ─────────────────────
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "results")

# ── Funções de carregamento ───────────────────────────────────
load_tsv_safe <- function(path, ...) {
  full <- file.path(RESULTS_DIR, path)
  if (!file.exists(full)) return(NULL)
  tryCatch(read_tsv(full, show_col_types = FALSE, ...), error = function(e) NULL)
}

# ── UI ────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "green",

  dashboardHeader(title = "RNASeq Insight – Soja"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Visão Geral",          tabName = "overview",    icon = icon("chart-bar")),
      menuItem("Expressão Diferencial",tabName = "de",          icon = icon("dna")),
      menuItem("Enriquecimento",        tabName = "enrichment",  icon = icon("project-diagram")),
      menuItem("Splicing",              tabName = "splicing",    icon = icon("code-branch")),
      menuItem("WGCNA",                 tabName = "wgcna",       icon = icon("network-wired")),
      menuItem("Integração",            tabName = "integration", icon = icon("layer-group")),
      menuItem("Dados Brutos",          tabName = "raw",         icon = icon("table"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f4f4; }
      .box { border-radius: 6px; }
      .value-box { border-radius: 6px; }
    "))),

    tabItems(

      # ── Visão Geral ────────────────────────────────────────
      tabItem(tabName = "overview",
        fluidRow(
          valueBoxOutput("vbox_total_de",  width = 3),
          valueBoxOutput("vbox_up",        width = 3),
          valueBoxOutput("vbox_down",      width = 3),
          valueBoxOutput("vbox_splicing",  width = 3)
        ),
        fluidRow(
          box(title = "PCA – Amostras", status = "primary", solidHeader = TRUE,
              width = 6, plotlyOutput("pca_plot", height = "380px")),
          box(title = "Genes DE por Regulação", status = "success", solidHeader = TRUE,
              width = 6, plotlyOutput("de_barplot", height = "380px"))
        )
      ),

      # ── Expressão Diferencial ─────────────────────────────
      tabItem(tabName = "de",
        fluidRow(
          box(width = 3, solidHeader = TRUE, status = "warning", title = "Filtros",
            sliderInput("padj_filter", "FDR máximo", min = 0.001, max = 0.2,
                        value = 0.05, step = 0.001),
            sliderInput("lfc_filter",  "|log2FC| mínimo", min = 0, max = 5,
                        value = 1, step = 0.1),
            radioButtons("reg_filter", "Regulação",
                         choices = c("Todos" = "all", "Up" = "up", "Down" = "down"),
                         selected = "all")
          ),
          box(width = 9, solidHeader = TRUE, status = "primary", title = "Volcano Plot",
            plotlyOutput("volcano_plot", height = "500px"))
        ),
        fluidRow(
          box(width = 12, solidHeader = TRUE, status = "info", title = "MA Plot",
            plotlyOutput("ma_plot", height = "400px"))
        ),
        fluidRow(
          box(width = 12, title = "Tabela de Resultados DESeq2",
            DTOutput("de_table"))
        )
      ),

      # ── Enriquecimento ────────────────────────────────────
      tabItem(tabName = "enrichment",
        tabBox(width = 12,
          tabPanel("GO – Biological Process",
            fluidRow(
              column(6, plotlyOutput("go_bp_plot",   height = "500px")),
              column(6, DTOutput("go_bp_table"))
            )
          ),
          tabPanel("KEGG Pathways",
            fluidRow(
              column(6, plotlyOutput("kegg_plot",    height = "500px")),
              column(6, DTOutput("kegg_table"))
            )
          ),
          tabPanel("GSEA – GO",
            DTOutput("gsea_go_table")
          ),
          tabPanel("GSEA – KEGG",
            DTOutput("gsea_kegg_table")
          )
        )
      ),

      # ── Splicing ──────────────────────────────────────────
      tabItem(tabName = "splicing",
        fluidRow(
          box(width = 12, solidHeader = TRUE, status = "danger",
              title = "Eventos de Splicing Alternativo (rMATS)",
            fluidRow(
              column(4, plotlyOutput("splicing_pie",  height = "350px")),
              column(8, DTOutput("splicing_table"))
            )
          )
        )
      ),

      # ── WGCNA ─────────────────────────────────────────────
      tabItem(tabName = "wgcna",
        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "info", title = "Módulos",
            plotlyOutput("module_size_plot", height = "400px")),
          box(width = 8, solidHeader = TRUE, status = "primary", title = "Hub Genes",
            DTOutput("hub_genes_table"))
        )
      ),

      # ── Integração ────────────────────────────────────────
      tabItem(tabName = "integration",
        fluidRow(
          box(width = 12, solidHeader = TRUE, status = "success",
              title = "Genes Candidatos – Score de Integração",
            plotlyOutput("integration_bubble", height = "550px")
          )
        ),
        fluidRow(
          box(width = 12, title = "Tabela de Candidatos",
            DTOutput("candidates_table"))
        )
      ),

      # ── Dados Brutos ──────────────────────────────────────
      tabItem(tabName = "raw",
        fluidRow(
          box(width = 12, title = "Contagens Normalizadas (VST)",
            DTOutput("norm_counts_table"))
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Carregamento reativo dos dados
  de_all    <- reactive(load_tsv_safe("deseq2/deseq2_results_all.tsv"))
  de_sig    <- reactive(load_tsv_safe("deseq2/deseq2_results_sig.tsv"))
  norm_cnt  <- reactive(load_tsv_safe("deseq2/normalized_counts.tsv"))
  meta_df   <- reactive(load_tsv_safe("counts/sample_metadata.tsv"))
  go_bp     <- reactive(load_tsv_safe("enrichment/go_bp_results.tsv"))
  kegg      <- reactive(load_tsv_safe("enrichment/kegg_results.tsv"))
  gsea_go   <- reactive(load_tsv_safe("enrichment/gsea_go_results.tsv"))
  gsea_kegg <- reactive(load_tsv_safe("enrichment/gsea_kegg_results.tsv"))
  splicing  <- reactive(load_tsv_safe("splicing/splicing_significant.tsv"))
  wgcna_mod <- reactive(load_tsv_safe("wgcna/wgcna_modules.tsv"))
  hub_gene  <- reactive(load_tsv_safe("wgcna/wgcna_hub_genes.tsv"))
  candidates <- reactive(load_tsv_safe("integration/key_candidates.tsv"))

  # DE filtrado
  de_filtered <- reactive({
    req(de_all())
    df <- de_all() %>%
      filter(!is.na(padj), padj < input$padj_filter,
             abs(log2FoldChange) > input$lfc_filter)
    if (input$reg_filter != "all")
      df <- df %>% filter(regulation == input$reg_filter)
    df
  })

  # ── ValueBoxes ─────────────────────────────────────────────
  output$vbox_total_de <- renderValueBox({
    n <- if (!is.null(de_sig())) nrow(de_sig()) else 0
    valueBox(n, "Genes DE Totais", icon = icon("dna"), color = "blue")
  })
  output$vbox_up <- renderValueBox({
    n <- if (!is.null(de_sig())) sum(de_sig()$regulation == "up", na.rm = TRUE) else 0
    valueBox(n, "Up-regulated", icon = icon("arrow-up"), color = "red")
  })
  output$vbox_down <- renderValueBox({
    n <- if (!is.null(de_sig())) sum(de_sig()$regulation == "down", na.rm = TRUE) else 0
    valueBox(n, "Down-regulated", icon = icon("arrow-down"), color = "navy")
  })
  output$vbox_splicing <- renderValueBox({
    n <- if (!is.null(splicing())) nrow(splicing()) else 0
    valueBox(n, "Eventos Splicing", icon = icon("code-branch"), color = "orange")
  })

  # ── PCA ────────────────────────────────────────────────────
  output$pca_plot <- renderPlotly({
    req(norm_cnt(), meta_df())
    mat <- norm_cnt() %>%
      tibble::column_to_rownames("gene_id") %>%
      as.matrix()
    pca <- prcomp(t(mat), scale. = TRUE)
    pct  <- round(summary(pca)$importance[2, 1:2] * 100, 1)
    df_pca <- data.frame(
      PC1 = pca$x[, 1], PC2 = pca$x[, 2],
      sample    = rownames(pca$x),
      condition = meta_df()$condition[match(rownames(pca$x), meta_df()$sample)]
    )
    plot_ly(df_pca, x = ~PC1, y = ~PC2, color = ~condition, text = ~sample,
            type = "scatter", mode = "markers+text",
            textposition = "top center",
            marker = list(size = 12)) %>%
      layout(title = "PCA",
             xaxis = list(title = paste0("PC1 (", pct[1], "%)")),
             yaxis = list(title = paste0("PC2 (", pct[2], "%)")))
  })

  # ── Barplot DE ──────────────────────────────────────────────
  output$de_barplot <- renderPlotly({
    req(de_sig())
    df <- de_sig() %>%
      count(regulation) %>%
      mutate(regulation = factor(regulation, levels = c("up", "down")))
    plot_ly(df, x = ~regulation, y = ~n, type = "bar",
            color = ~regulation,
            colors = c(up = "#D6604D", down = "#2166AC")) %>%
      layout(title = "Genes DE", yaxis = list(title = "Nº de genes"),
             showlegend = FALSE)
  })

  # ── Volcano ─────────────────────────────────────────────────
  output$volcano_plot <- renderPlotly({
    req(de_all())
    df <- de_all() %>%
      filter(!is.na(padj), !is.na(log2FoldChange)) %>%
      mutate(
        color = case_when(
          padj < input$padj_filter & log2FoldChange >  input$lfc_filter ~ "Up",
          padj < input$padj_filter & log2FoldChange < -input$lfc_filter ~ "Down",
          TRUE ~ "NS"
        ),
        log_padj = -log10(padj + 1e-300)
      )
    plot_ly(df, x = ~log2FoldChange, y = ~log_padj,
            color = ~color,
            colors = c(Up = "#D6604D", Down = "#2166AC", NS = "#BBBBBB"),
            text = ~paste("Gene:", gene_id, "<br>LFC:", round(log2FoldChange, 2),
                          "<br>FDR:", formatC(padj, format = "e", digits = 2)),
            type = "scatter", mode = "markers",
            marker = list(size = 4, opacity = 0.7)) %>%
      layout(title = "Volcano Plot",
             xaxis = list(title = "log2 Fold Change"),
             yaxis = list(title = "-log10(FDR)"),
             shapes = list(
               list(type = "line", x0 = input$lfc_filter,  x1 = input$lfc_filter,
                    y0 = 0, y1 = max(df$log_padj), line = list(dash = "dash", color = "grey")),
               list(type = "line", x0 = -input$lfc_filter, x1 = -input$lfc_filter,
                    y0 = 0, y1 = max(df$log_padj), line = list(dash = "dash", color = "grey")),
               list(type = "line", x0 = min(df$log2FoldChange), x1 = max(df$log2FoldChange),
                    y0 = -log10(input$padj_filter), y1 = -log10(input$padj_filter),
                    line = list(dash = "dash", color = "grey"))
             ))
  })

  # ── MA plot ──────────────────────────────────────────────────
  output$ma_plot <- renderPlotly({
    req(de_all())
    df <- de_all() %>%
      filter(!is.na(padj), !is.na(log2FoldChange)) %>%
      mutate(color = ifelse(padj < input$padj_filter & abs(log2FoldChange) > input$lfc_filter,
                            regulation, "ns"))
    plot_ly(df, x = ~log10(baseMean + 1), y = ~log2FoldChange,
            color = ~color,
            colors = c(up = "#D6604D", down = "#2166AC", ns = "#CCCCCC"),
            text = ~paste("Gene:", gene_id),
            type = "scatter", mode = "markers",
            marker = list(size = 3, opacity = 0.6)) %>%
      layout(title = "MA Plot",
             xaxis = list(title = "log10(baseMean + 1)"),
             yaxis = list(title = "log2 Fold Change"))
  })

  # ── Tabelas ──────────────────────────────────────────────────
  output$de_table <- renderDT({
    req(de_filtered())
    de_filtered() %>%
      mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
      datatable(filter = "top", extensions = "Buttons",
                options = list(dom = "Bfrtip",
                               buttons = c("csv", "excel"),
                               pageLength = 15))
  })

  output$go_bp_plot <- renderPlotly({
    req(go_bp())
    df <- go_bp() %>% filter(p.adjust < 0.05) %>% head(20)
    if (nrow(df) == 0) return(NULL)
    plot_ly(df, x = ~Count, y = ~reorder(Description, Count),
            type = "bar", orientation = "h",
            color = ~-log10(p.adjust), colors = "Blues") %>%
      layout(title = "GO – Biological Process",
             xaxis = list(title = "Gene Count"),
             yaxis = list(title = ""))
  })
  output$go_bp_table <- renderDT({
    req(go_bp())
    go_bp() %>% filter(p.adjust < 0.05) %>%
      select(Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
      mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
      datatable(filter = "top", options = list(pageLength = 10))
  })

  output$kegg_plot <- renderPlotly({
    req(kegg())
    df <- kegg() %>% filter(p.adjust < 0.05) %>% head(20)
    if (nrow(df) == 0) return(NULL)
    plot_ly(df, x = ~Count, y = ~reorder(Description, Count),
            type = "bar", orientation = "h",
            color = ~-log10(p.adjust), colors = "Reds") %>%
      layout(title = "KEGG Pathways",
             xaxis = list(title = "Gene Count"),
             yaxis = list(title = ""))
  })
  output$kegg_table <- renderDT({
    req(kegg())
    kegg() %>% filter(p.adjust < 0.05) %>%
      select(Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
      mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
      datatable(filter = "top", options = list(pageLength = 10))
  })

  output$gsea_go_table   <- renderDT({ req(gsea_go());   datatable(gsea_go(),   filter = "top", options = list(pageLength = 15)) })
  output$gsea_kegg_table <- renderDT({ req(gsea_kegg()); datatable(gsea_kegg(), filter = "top", options = list(pageLength = 15)) })

  output$splicing_pie <- renderPlotly({
    req(splicing())
    if (nrow(splicing()) == 0) return(NULL)
    df <- splicing() %>% count(event_type)
    plot_ly(df, labels = ~event_type, values = ~n, type = "pie",
            textinfo = "label+percent+value") %>%
      layout(title = "Eventos por tipo")
  })
  output$splicing_table <- renderDT({
    req(splicing())
    datatable(splicing(), filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  output$module_size_plot <- renderPlotly({
    req(wgcna_mod())
    df <- wgcna_mod() %>% count(module) %>% filter(module != "grey") %>%
      arrange(desc(n))
    plot_ly(df, x = ~reorder(module, n), y = ~n, type = "bar",
            marker = list(color = df$module)) %>%
      layout(title = "Tamanho dos Módulos WGCNA",
             xaxis = list(title = "Módulo"), yaxis = list(title = "Nº de genes"))
  })
  output$hub_genes_table <- renderDT({
    req(hub_gene())
    datatable(hub_gene(), filter = "top",
              options = list(pageLength = 15))
  })

  output$integration_bubble <- renderPlotly({
    req(candidates())
    if (nrow(candidates()) == 0) return(NULL)
    df <- head(candidates(), 50)
    plot_ly(df, x = ~log2FoldChange, y = ~integration_score,
            size = ~-log10(padj + 1e-300),
            color = ~regulation,
            colors = c(up = "#D6604D", down = "#2166AC"),
            text = ~paste("Gene:", gene_id,
                          "<br>LFC:", round(log2FoldChange, 2),
                          "<br>Score:", integration_score),
            type = "scatter", mode = "markers",
            marker = list(sizemode = "diameter", opacity = 0.8)) %>%
      layout(title = "Top 50 – Score de Integração",
             xaxis = list(title = "log2 Fold Change"),
             yaxis = list(title = "Integration Score"))
  })
  output$candidates_table <- renderDT({
    req(candidates())
    datatable(candidates(), filter = "top", extensions = "Buttons",
              options = list(dom = "Bfrtip", buttons = c("csv", "excel"),
                             pageLength = 15))
  })

  output$norm_counts_table <- renderDT({
    req(norm_cnt())
    datatable(norm_cnt(), filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })
}

# ── Launch ───────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
