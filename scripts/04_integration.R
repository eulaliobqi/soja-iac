#!/usr/bin/env Rscript
# ============================================================
# 04_integration.R – Integração multi-ômica e ranking de genes
# DE + Splicing + Enriquecimento + WGCNA
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(tibble)
  library(stringr)
})

# ── Argumentos ───────────────────────────────────────────────
opt_list <- list(
  make_option("--deseq2",      type = "character", help = "DESeq2 significativos (TSV)"),
  make_option("--splicing",    type = "character", help = "Splicing significativos (TSV)"),
  make_option("--go",          type = "character", help = "GO-BP results (TSV)"),
  make_option("--kegg",        type = "character", help = "KEGG results (TSV)"),
  make_option("--wgcna",       type = "character", help = "WGCNA modules (TSV)"),
  make_option("--hub_genes",   type = "character", help = "WGCNA hub genes (TSV)"),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("═══════════════════════════════════════\n")
cat("  Integração Multi-ômica\n")
cat("═══════════════════════════════════════\n")

# Remove sufixo de versão Phytozome para comparações entre datasets
clean_glyma <- function(ids) {
  sub("^(Glyma\\.[0-9A-Za-z]+G[0-9]+).*", "\\1", as.character(ids))
}

theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

save_plot <- function(p, name, w = 10, h = 8) {
  ggsave(file.path(opt$figures_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(opt$figures_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# ── 1. Carrega dados ──────────────────────────────────────────
de_sig <- tryCatch(read_tsv(opt$deseq2, show_col_types = FALSE),
                   error = function(e) data.frame())

splicing_sig <- tryCatch(read_tsv(opt$splicing, show_col_types = FALSE),
                          error = function(e) data.frame())

go_res  <- tryCatch(read_tsv(opt$go,   show_col_types = FALSE),
                    error = function(e) data.frame())
kegg_res <- tryCatch(read_tsv(opt$kegg, show_col_types = FALSE),
                     error = function(e) data.frame())

wgcna_mod <- tryCatch(read_tsv(opt$wgcna,     show_col_types = FALSE),
                      error = function(e) data.frame())
hub_genes  <- tryCatch(read_tsv(opt$hub_genes, show_col_types = FALSE),
                       error = function(e) data.frame())

cat(sprintf("DE significativos: %d\n", nrow(de_sig)))
cat(sprintf("Eventos splicing: %d\n", nrow(splicing_sig)))
cat(sprintf("Termos GO-BP: %d\n", nrow(go_res)))
cat(sprintf("KEGG pathways: %d\n", nrow(kegg_res)))

# ── 2. Genes em vias enriquecidas ─────────────────────────────
get_genes_from_enrichment <- function(enrich_df) {
  if (nrow(enrich_df) == 0 || !"geneID" %in% colnames(enrich_df)) return(character(0))
  genes <- enrich_df$geneID %>%
    str_split("/") %>%
    unlist() %>%
    unique() %>%
    na.omit()
  return(genes)
}

go_genes   <- clean_glyma(get_genes_from_enrichment(go_res))
kegg_genes <- clean_glyma(get_genes_from_enrichment(kegg_res))
enriched_genes <- union(go_genes, kegg_genes)

# Genes com splicing alternativo
splicing_gene_col <- intersect(c("GeneID", "geneSymbol", "gene_id"), colnames(splicing_sig))[1]
splice_genes <- if (!is.na(splicing_gene_col) && nrow(splicing_sig) > 0)
  splicing_sig[[splicing_gene_col]] else character(0)

# Hub genes WGCNA
hub_gene_ids <- if (nrow(hub_genes) > 0) hub_genes$gene_id else character(0)

# ── 3. Integração e scoring ───────────────────────────────────
if (nrow(de_sig) > 0) {
  integrated <- de_sig %>%
    select(gene_id, log2FoldChange, padj, baseMean, regulation) %>%
    mutate(
      # Score normalizado de LFC e baseMean
      lfc_score    = abs(log2FoldChange) / max(abs(log2FoldChange), na.rm = TRUE),
      mean_score   = log10(baseMean + 1) / max(log10(baseMean + 1), na.rm = TRUE),
      sig_score    = -log10(padj + 1e-300) / max(-log10(padj + 1e-300), na.rm = TRUE),

      # Flags de camadas de evidência (usa IDs limpos para comparação cross-dataset)
      gene_id_clean = clean_glyma(gene_id),
      has_splicing  = gene_id_clean %in% clean_glyma(splice_genes),
      in_pathway    = gene_id_clean %in% enriched_genes,
      is_hub        = gene_id_clean %in% clean_glyma(hub_gene_ids),

      # Score composto (0-10)
      integration_score = round(
        (lfc_score * 3 + mean_score * 2 + sig_score * 2 +
         as.numeric(has_splicing) * 1.5 +
         as.numeric(in_pathway)   * 1.0 +
         as.numeric(is_hub)       * 0.5),
        2
      )
    )

  # Adiciona módulo WGCNA
  if (nrow(wgcna_mod) > 0) {
    integrated <- integrated %>%
      left_join(wgcna_mod %>% select(gene_id, module), by = "gene_id")
  }

  # Ranking final
  ranking <- integrated %>%
    arrange(desc(integration_score))

  # Genes candidatos (evidência em 2+ camadas)
  candidates <- integrated %>%
    filter((has_splicing + in_pathway + is_hub) >= 2) %>%
    arrange(desc(integration_score))

} else {
  integrated <- data.frame()
  ranking    <- data.frame()
  candidates <- data.frame()
  cat("AVISO: sem genes DE significativos para integração.\n")
}

# ── 4. Exporta ───────────────────────────────────────────────
write_tsv(integrated, file.path(opt$outdir, "integrated_genes.tsv"))
write_tsv(ranking,    file.path(opt$outdir, "gene_ranking.tsv"))
write_tsv(candidates, file.path(opt$outdir, "key_candidates.tsv"))

cat(sprintf("\nGenes candidatos (2+ camadas de evidência): %d\n", nrow(candidates)))
if (nrow(candidates) > 0) {
  cat("Top 10 candidatos:\n")
  print(head(candidates %>% select(gene_id, log2FoldChange, padj,
                                    has_splicing, in_pathway, is_hub,
                                    integration_score), 10))
}

# ── 5. Figuras ────────────────────────────────────────────────

# 5.1 Diagrama de Venn / UpSet – sobreposição das camadas
if (nrow(integrated) > 0) {

  # Barplot de evidências
  evid_df <- data.frame(
    Camada  = c("Apenas DE", "DE + Splicing", "DE + Via", "DE + Hub", "Multi-camada"),
    N_genes = c(
      sum(!integrated$has_splicing & !integrated$in_pathway & !integrated$is_hub),
      sum( integrated$has_splicing & !integrated$in_pathway & !integrated$is_hub),
      sum(!integrated$has_splicing &  integrated$in_pathway & !integrated$is_hub),
      sum(!integrated$has_splicing & !integrated$in_pathway &  integrated$is_hub),
      sum((integrated$has_splicing + integrated$in_pathway + integrated$is_hub) >= 2)
    )
  )

  p_evid <- ggplot(evid_df, aes(reorder(Camada, N_genes), N_genes)) +
    geom_bar(stat = "identity", fill = "#4393C3", width = 0.6) +
    geom_text(aes(label = N_genes), hjust = -0.2) +
    coord_flip() +
    labs(title = "Genes por camada de evidência",
         x = NULL, y = "Número de genes") +
    theme_pub
  save_plot(p_evid, "evidence_layers", w = 8, h = 5)

  # Bubble plot: integration score
  if (nrow(ranking) > 0) {
    top_ranked <- head(ranking, 30)
    p_bubble <- ggplot(top_ranked,
                       aes(x = log2FoldChange, y = integration_score,
                           size = -log10(padj + 1e-300),
                           color = regulation,
                           label = gene_id)) +
      geom_point(alpha = 0.8) +
      geom_text_repel(size = 3, max.overlaps = 15) +
      scale_color_manual(values = c(up = "#D6604D", down = "#2166AC")) +
      scale_size_continuous(range = c(3, 10), name = "-log10(padj)") +
      labs(title = "Top 30 genes – Score de integração",
           x = "log2 Fold Change",
           y = "Integration Score",
           color = "Regulação") +
      theme_pub
    save_plot(p_bubble, "top_genes_integration", w = 10, h = 8)
  }
}

# 5.2 Tabela sumária de candidatos (se houver)
if (nrow(candidates) > 0) {
  cat_df <- candidates %>%
    mutate(
      Evidence = paste0(
        ifelse(has_splicing, "Splicing|", ""),
        ifelse(in_pathway,   "Pathway|",  ""),
        ifelse(is_hub,       "Hub",       "")
      ) %>% str_remove("\\|$")
    ) %>%
    select(gene_id, log2FoldChange, padj, integration_score, Evidence) %>%
    head(50)

  write_tsv(cat_df, file.path(opt$outdir, "candidates_table.tsv"))
}

cat("\nIntegração concluída. Resultados em:", opt$outdir, "\n")
