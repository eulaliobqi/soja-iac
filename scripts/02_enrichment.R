#!/usr/bin/env Rscript
# ============================================================
# 02_enrichment.R – GO, KEGG e GSEA para Glycine max
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Gmax.eg.db)
  library(AnnotationDbi)
  library(fgsea)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# ── Argumentos ───────────────────────────────────────────────
opt_list <- list(
  make_option("--deseq2",      type = "character", help = "DESeq2 results (TSV, todos os genes)"),
  make_option("--norm_counts", type = "character", help = "Contagens normalizadas (TSV)"),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1.0),
  make_option("--organism",    type = "character", default = "gma"),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("═══════════════════════════════════════\n")
cat("  Enriquecimento Funcional – GO + KEGG + GSEA\n")
cat("═══════════════════════════════════════\n")

# ── 1. Leitura ────────────────────────────────────────────────
res_all <- read_tsv(opt$deseq2, show_col_types = FALSE)

# Genes DE significativos
de_sig <- res_all %>%
  filter(!is.na(padj), padj < opt$padj, abs(log2FoldChange) > opt$lfc)

cat(sprintf("Genes DE: %d | Todos testados: %d\n", nrow(de_sig), nrow(res_all)))

# ── 2. Mapeamento de IDs para ENTREZ ─────────────────────────
# org.Gmax.eg.db usa SYMBOL e ENTREZID
# Os IDs do featureCounts (Glyma.XXX) são mapeados via SYMBOL ou TAIR

map_to_entrez <- function(gene_ids, db = org.Gmax.eg.db) {
  # Tenta múltiplas keytypes na ordem de prioridade
  for (ktype in c("SYMBOL", "TAIR", "ALIAS", "GENENAME")) {
    available <- try(keytypes(db), silent = TRUE)
    if (inherits(available, "try-error")) next
    if (!ktype %in% available) next

    mapped <- tryCatch(
      mapIds(db, keys = gene_ids, column = "ENTREZID",
             keytype = ktype, multiVals = "first"),
      error = function(e) NULL
    )
    if (!is.null(mapped) && sum(!is.na(mapped)) > 0) {
      cat(sprintf("  IDs mapeados via keytype '%s': %d/%d\n",
                  ktype, sum(!is.na(mapped)), length(gene_ids)))
      return(mapped)
    }
  }
  # Fallback: usa o próprio ID numérico se já for ENTREZ
  if (all(grepl("^[0-9]+$", gene_ids[!is.na(gene_ids)]))) {
    cat("  IDs já estão em formato ENTREZID numérico\n")
    return(setNames(gene_ids, gene_ids))
  }
  warning("Mapeamento de IDs falhou. Verifique os gene IDs do featureCounts.")
  return(setNames(rep(NA, length(gene_ids)), gene_ids))
}

# Mapeamento
cat("Mapeando gene IDs para ENTREZID...\n")
all_genes_entrez  <- map_to_entrez(res_all$gene_id)
de_genes_entrez   <- map_to_entrez(de_sig$gene_id)

# Remove NAs
background  <- na.omit(as.character(all_genes_entrez))
de_entrez   <- na.omit(as.character(de_genes_entrez))

cat(sprintf("Background: %d | DE com ENTREZ: %d\n", length(background), length(de_entrez)))

# ── 3. ORA – GO (BP, MF, CC) ─────────────────────────────────
run_go_ora <- function(genes, bg, ont, label) {
  tryCatch({
    res <- enrichGO(
      gene          = genes,
      universe      = bg,
      OrgDb         = org.Gmax.eg.db,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE,
      minGSSize     = 5,
      maxGSSize     = 500
    )
    cat(sprintf("  GO-%s: %d termos enriquecidos\n", ont, nrow(res@result[res@result$p.adjust < 0.05, ])))
    return(res)
  }, error = function(e) {
    message(sprintf("  GO-%s falhou: %s", ont, e$message))
    return(NULL)
  })
}

cat("\nAnálise GO – ORA:\n")
go_bp <- run_go_ora(de_entrez, background, "BP", "BP")
go_mf <- run_go_ora(de_entrez, background, "MF", "MF")
go_cc <- run_go_ora(de_entrez, background, "CC", "CC")

# Exporta
export_enrich <- function(obj, fname) {
  if (is.null(obj)) {
    write_tsv(data.frame(), file.path(opt$outdir, fname))
  } else {
    write_tsv(as.data.frame(obj), file.path(opt$outdir, fname))
  }
}
export_enrich(go_bp, "go_bp_results.tsv")
export_enrich(go_mf, "go_mf_results.tsv")
export_enrich(go_cc, "go_cc_results.tsv")

# ── 4. ORA – KEGG ────────────────────────────────────────────
cat("\nAnálise KEGG – ORA:\n")
kegg_res <- tryCatch({
  res <- enrichKEGG(
    gene          = de_entrez,
    universe      = background,
    organism      = opt$organism,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 5,
    maxGSSize     = 500
  )
  cat(sprintf("  KEGG: %d pathways enriquecidos\n",
              nrow(res@result[res@result$p.adjust < 0.05, ])))
  res
}, error = function(e) {
  message(sprintf("  KEGG falhou: %s", e$message))
  NULL
})
export_enrich(kegg_res, "kegg_results.tsv")

# ── 5. GSEA – GO ─────────────────────────────────────────────
cat("\nGSEA:\n")
# Ranking por sinal (LFC * -log10(padj))
ranked_df <- res_all %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(entrez = as.character(all_genes_entrez[gene_id])) %>%
  filter(!is.na(entrez)) %>%
  mutate(rank_score = log2FoldChange * -log10(padj + 1e-300)) %>%
  arrange(desc(rank_score))

ranked_vec <- setNames(ranked_df$rank_score, ranked_df$entrez)
ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]

gsea_go <- tryCatch({
  res <- gseGO(
    geneList     = ranked_vec,
    OrgDb        = org.Gmax.eg.db,
    ont          = "BP",
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose      = FALSE
  )
  cat(sprintf("  GSEA-GO: %d termos\n", nrow(res@result[res@result$p.adjust < 0.05, ])))
  res
}, error = function(e) {
  message(sprintf("  GSEA-GO falhou: %s", e$message))
  NULL
})

gsea_kegg <- tryCatch({
  res <- gseKEGG(
    geneList     = ranked_vec,
    organism     = opt$organism,
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose      = FALSE
  )
  cat(sprintf("  GSEA-KEGG: %d pathways\n", nrow(res@result[res@result$p.adjust < 0.05, ])))
  res
}, error = function(e) {
  message(sprintf("  GSEA-KEGG falhou: %s", e$message))
  NULL
})

export_enrich(gsea_go,   "gsea_go_results.tsv")
export_enrich(gsea_kegg, "gsea_kegg_results.tsv")

# ── 6. Figuras ────────────────────────────────────────────────
theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

save_plot <- function(p, name, w = 10, h = 8) {
  ggsave(file.path(opt$figures_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(opt$figures_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# GO BP dotplot
if (!is.null(go_bp) && nrow(go_bp) > 0) {
  p_go <- dotplot(go_bp, showCategory = 20, title = "GO – Biological Process") +
    theme_pub
  save_plot(p_go, "go_bp_dotplot", w = 10, h = 9)

  # GO BP enrichment map
  tryCatch({
    go_bp2 <- pairwise_termsim(go_bp)
    p_emap <- emapplot(go_bp2, showCategory = 25) +
      ggtitle("GO BP – Enrichment Map")
    save_plot(p_emap, "go_bp_emap", w = 12, h = 10)
  }, error = function(e) message("emapplot falhou: ", e$message))
}

# KEGG dotplot
if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
  p_kegg <- dotplot(kegg_res, showCategory = 20, title = "KEGG Pathways") +
    theme_pub
  save_plot(p_kegg, "kegg_dotplot", w = 10, h = 9)

  # KEGG barplot
  p_kegg_bar <- barplot(kegg_res, showCategory = 20, title = "KEGG Pathways") +
    theme_pub
  save_plot(p_kegg_bar, "kegg_barplot", w = 10, h = 9)
}

# GSEA GO plot
if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
  p_gsea <- dotplot(gsea_go, showCategory = 20, split = ".sign",
                    title = "GSEA – GO Biological Process") +
    facet_grid(. ~ .sign) +
    theme_pub
  save_plot(p_gsea, "gsea_go_dotplot", w = 12, h = 9)
}

# GSEA KEGG
if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  p_gsea_kegg <- dotplot(gsea_kegg, showCategory = 20, split = ".sign",
                          title = "GSEA – KEGG Pathways") +
    facet_grid(. ~ .sign) +
    theme_pub
  save_plot(p_gsea_kegg, "gsea_kegg_dotplot", w = 12, h = 9)
}

cat("\nEnriquecimento funcional concluído. Resultados em:", opt$outdir, "\n")
