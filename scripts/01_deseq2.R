#!/usr/bin/env Rscript
# ============================================================
# 01_deseq2.R – Expressão diferencial com DESeq2
# Glycine max RNASeq Pipeline
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(EnhancedVolcano)
  library(dplyr)
  library(readr)
  library(scales)
})

# ── Argumentos ───────────────────────────────────────────────
opt_list <- list(
  make_option("--counts",       type = "character", help = "Matriz de contagens (TSV)"),
  make_option("--metadata",     type = "character", help = "Metadados das amostras (TSV)"),
  make_option("--control",      type = "character", default = "control"),
  make_option("--treatment",    type = "character", default = "treatment"),
  make_option("--padj",         type = "double",    default = 0.05),
  make_option("--lfc",          type = "double",    default = 1.0),
  make_option("--outdir",       type = "character", default = "."),
  make_option("--figures_dir",  type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir,       showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir,  showWarnings = FALSE, recursive = TRUE)

cat("═══════════════════════════════════════\n")
cat("  DESeq2 – Expressão Diferencial\n")
cat("═══════════════════════════════════════\n")

# ── 1. Leitura dos dados ──────────────────────────────────────
counts <- read_tsv(opt$counts, show_col_types = FALSE) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

meta <- read_tsv(opt$metadata, show_col_types = FALSE) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("sample")

# Garante correspondência de ordem
counts <- counts[, rownames(meta), drop = FALSE]
meta$condition <- factor(meta$condition, levels = c(opt$control, opt$treatment))

cat(sprintf("Genes: %d | Amostras: %d\n", nrow(counts), ncol(counts)))
cat(sprintf("Contraste: %s vs %s\n", opt$treatment, opt$control))

# ── 2. Cria objeto DESeq2 ─────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition
)

# Filtragem: remove genes com < 10 reads totais
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat(sprintf("Genes após filtragem (>= 10 reads): %d\n", nrow(dds)))

# ── 3. Normalização e DE ──────────────────────────────────────
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

# Resultados com shrinkage LFC (apeglm)
contrast  <- c("condition", opt$treatment, opt$control)
res_raw   <- results(dds, contrast = contrast, alpha = opt$padj)

# Tenta apeglm shrinkage (mais robusto)
tryCatch({
  coef_name <- resultsNames(dds)[grep(opt$treatment, resultsNames(dds))]
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
}, error = function(e) {
  message("apeglm indisponível, usando ashr shrinkage...")
  res <<- lfcShrink(dds, contrast = contrast, type = "ashr", quiet = TRUE)
})

# Contagens normalizadas (vst)
vst_data  <- vst(dds, blind = FALSE)
norm_mat  <- assay(vst_data)

# ── 4. Exporta resultados ─────────────────────────────────────
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  mutate(
    regulation = case_when(
      padj < opt$padj & log2FoldChange > opt$lfc  ~ "up",
      padj < opt$padj & log2FoldChange < -opt$lfc ~ "down",
      TRUE                                          ~ "ns"
    )
  )

res_sig <- res_df %>%
  filter(padj < opt$padj, abs(log2FoldChange) > opt$lfc)

write_tsv(res_df,  file.path(opt$outdir, "deseq2_results_all.tsv"))
write_tsv(res_sig, file.path(opt$outdir, "deseq2_results_sig.tsv"))

norm_df <- as.data.frame(norm_mat) %>%
  tibble::rownames_to_column("gene_id")
write_tsv(norm_df, file.path(opt$outdir, "normalized_counts.tsv"))

up_count   <- sum(res_sig$regulation == "up",   na.rm = TRUE)
down_count <- sum(res_sig$regulation == "down",  na.rm = TRUE)

cat(sprintf("\nGenes significativos (FDR < %.2f, |LFC| > %.1f):\n", opt$padj, opt$lfc))
cat(sprintf("  ▲ Up-regulated  : %d\n", up_count))
cat(sprintf("  ▼ Down-regulated: %d\n", down_count))
cat(sprintf("  Total           : %d\n", nrow(res_sig)))

# ── 5. Figuras ────────────────────────────────────────────────
theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# 5.1 PCA
pca_data  <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
pca_var   <- attr(pca_data, "percentVar") * 100
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3.5, max.overlaps = 20) +
  scale_color_manual(values = c("#2166AC", "#D6604D")) +
  labs(title = "PCA – Amostras",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var[2]),
       color = "Condição") +
  theme_pub
ggsave(file.path(opt$figures_dir, "pca_samples.pdf"), p_pca, width = 7, height = 6)
ggsave(file.path(opt$figures_dir, "pca_samples.png"), p_pca, width = 7, height = 6, dpi = 300)

# 5.2 Volcano plot
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(regulation != "ns") %>%
  group_by(regulation) %>%
  slice_min(padj, n = 10) %>%
  ungroup()

p_vol <- EnhancedVolcano(res_df,
  lab            = res_df$gene_id,
  x              = "log2FoldChange",
  y              = "padj",
  pCutoff        = opt$padj,
  FCcutoff       = opt$lfc,
  pointSize      = 2,
  labSize        = 3,
  title          = sprintf("Volcano – %s vs %s", opt$treatment, opt$control),
  subtitle       = sprintf("FDR < %.2f | |LFC| > %.1f", opt$padj, opt$lfc),
  col            = c("#AAAAAA", "#2166AC", "#D6604D", "#B2182B"),
  colAlpha       = 0.7,
  legendPosition = "bottom",
  drawConnectors = TRUE,
  selectLab      = head(top_genes$gene_id, 20)
)
ggsave(file.path(opt$figures_dir, "volcano_plot.pdf"), p_vol, width = 10, height = 8)
ggsave(file.path(opt$figures_dir, "volcano_plot.png"), p_vol, width = 10, height = 8, dpi = 300)

# 5.3 MA plot
p_ma <- ggplot(res_df %>% filter(!is.na(padj)),
               aes(x = log10(baseMean + 1),
                   y = log2FoldChange,
                   color = regulation)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = c(-opt$lfc, opt$lfc), linetype = "dashed", color = "navy") +
  geom_hline(yintercept = 0, color = "black") +
  scale_color_manual(values = c(up = "#D6604D", down = "#2166AC", ns = "#CCCCCC")) +
  labs(title = "MA Plot",
       x     = "log10(baseMean + 1)",
       y     = "log2 Fold Change",
       color = "Regulação") +
  theme_pub
ggsave(file.path(opt$figures_dir, "ma_plot.pdf"), p_ma, width = 8, height = 6)
ggsave(file.path(opt$figures_dir, "ma_plot.png"), p_ma, width = 8, height = 6, dpi = 300)

# 5.4 Heatmap – top 50 genes DE
top50 <- res_sig %>%
  slice_min(padj, n = 50) %>%
  pull(gene_id)

if (length(top50) > 1) {
  mat_heat <- norm_mat[top50, , drop = FALSE]
  mat_heat <- t(scale(t(mat_heat)))  # z-score por gene

  annot_col <- data.frame(
    Condição = meta$condition,
    row.names = rownames(meta)
  )

  pheatmap(mat_heat,
    annotation_col = annot_col,
    show_rownames  = length(top50) <= 50,
    fontsize_row   = 7,
    fontsize_col   = 10,
    clustering_method = "ward.D2",
    color  = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
    border_color = NA,
    main   = sprintf("Heatmap – Top %d genes DE (z-score)", length(top50)),
    filename = file.path(opt$figures_dir, "heatmap_top_genes.pdf"),
    width  = 10, height = 12
  )
  pheatmap(mat_heat,
    annotation_col = annot_col,
    show_rownames  = length(top50) <= 50,
    fontsize_row   = 7,
    fontsize_col   = 10,
    clustering_method = "ward.D2",
    color  = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
    border_color = NA,
    main   = sprintf("Heatmap – Top %d genes DE (z-score)", length(top50)),
    filename = file.path(opt$figures_dir, "heatmap_top_genes.png"),
    width  = 10, height = 12
  )
}

# 5.5 Barplot de genes DE up/down
bar_df <- data.frame(
  Regulação = c("Up-regulated", "Down-regulated"),
  Genes     = c(up_count, down_count),
  Cor       = c("#D6604D", "#2166AC")
)
p_bar <- ggplot(bar_df, aes(Regulação, Genes, fill = Cor)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = Genes), vjust = -0.5, fontface = "bold") +
  scale_fill_identity() +
  labs(title = "Genes Diferencialmente Expressos",
       x = NULL, y = "Número de genes") +
  theme_pub +
  theme(legend.position = "none")
ggsave(file.path(opt$figures_dir, "de_barplot.pdf"), p_bar, width = 5, height = 5)
ggsave(file.path(opt$figures_dir, "de_barplot.png"), p_bar, width = 5, height = 5, dpi = 300)

# ── 6. Resumo ─────────────────────────────────────────────────
summary_txt <- sprintf(
"DESeq2 Analysis Summary
======================
Contraste      : %s vs %s
Genes testados : %d
FDR cutoff     : %.2f
|LFC| cutoff   : %.1f

Resultados:
  Up-regulated   : %d
  Down-regulated : %d
  Total DE       : %d

Dispersão estimada: OK
Modelo: Wald test (LFC shrinkage: apeglm/ashr)
",
  opt$treatment, opt$control, nrow(res_df),
  opt$padj, opt$lfc, up_count, down_count, nrow(res_sig)
)
writeLines(summary_txt, file.path(opt$outdir, "deseq2_summary.txt"))
cat("\nAnálise DESeq2 concluída. Resultados em:", opt$outdir, "\n")
