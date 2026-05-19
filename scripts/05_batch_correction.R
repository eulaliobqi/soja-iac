#!/usr/bin/env Rscript
# ============================================================
# 05_batch_correction.R – ComBat-Seq + diagnóstico PCA
# Critério: aplicar se PC1 > 40% var AND cor(PC1, batch) > 0.7
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(sva)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

opt_list <- list(
  make_option("--counts",   type = "character", help = "Counts matrix (TSV, gene_id col)"),
  make_option("--metadata", type = "character", help = "Metadados com coluna batch (TSV)"),
  make_option("--outdir",   type = "character", default = "."),
  make_option("--pc1_var",  type = "double",    default = 0.40),
  make_option("--cor_thr",  type = "double",    default = 0.70)
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "counts_corrected.tsv"))
writeLines("", file.path(opt$outdir, "batch_report.txt"))

cat("═══════════════════════════════════════\n")
cat("  ComBat-Seq – Correção de Batch\n")
cat("═══════════════════════════════════════\n")

counts <- read_tsv(opt$counts, show_col_types = FALSE) |>
  tibble::column_to_rownames("gene_id") |>
  as.matrix()
counts <- round(counts)
storage.mode(counts) <- "integer"

meta <- read_tsv(opt$metadata, show_col_types = FALSE) |>
  as.data.frame() |>
  tibble::column_to_rownames("sample")
counts <- counts[, rownames(meta), drop = FALSE]

# ── PCA diagnóstico ───────────────────────────────────────────
pca_diag <- function(mat, meta_df, title) {
  dds_tmp <- DESeqDataSetFromMatrix(mat, colData = meta_df, design = ~1)
  vst_mat  <- assay(vst(dds_tmp, blind = TRUE))
  vst_mat  <- vst_mat[apply(vst_mat, 1, var) > 0, , drop = FALSE]  # Lição 5
  pca      <- prcomp(t(vst_mat), scale. = FALSE)
  pct      <- summary(pca)$importance[2, ] * 100

  df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    condition = meta_df$condition,
    sample    = rownames(meta_df)
  )
  if ("batch" %in% names(meta_df)) df$batch <- meta_df$batch

  p <- ggplot(df, aes(PC1, PC2, color = condition, label = sample)) +
    geom_point(size = 3) + ggrepel::geom_text_repel(size = 3) +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", pct[1]),
         y = sprintf("PC2 (%.1f%%)", pct[2])) +
    theme_bw()

  list(pca = pca, pct = pct, plot = p, vst = vst_mat)
}

pre <- pca_diag(counts, meta, "PCA – Antes da correção")
ggsave(file.path(opt$outdir, "pca_before_batch.pdf"), pre$plot, width = 7, height = 5)

# ── Decisão de aplicar ComBat-Seq ─────────────────────────────
apply_combat <- FALSE
report_lines <- c()

if ("batch" %in% names(meta)) {
  batch_vec   <- meta$batch
  pc1_var_obs <- pre$pct[1] / 100
  cor_val     <- abs(cor(pre$pca$x[, 1], as.numeric(factor(batch_vec))))
  report_lines <- c(
    sprintf("PC1 variância: %.1f%%", pc1_var_obs * 100),
    sprintf("Correlação PC1 × batch: %.3f", cor_val),
    sprintf("Limiar PC1: %.0f%%  |  Limiar cor: %.2f", opt$pc1_var * 100, opt$cor_thr)
  )
  apply_combat <- (pc1_var_obs > opt$pc1_var) && (cor_val > opt$cor_thr)
  report_lines <- c(report_lines,
                    sprintf("Decisão: %s", ifelse(apply_combat, "APLICAR ComBat-Seq", "NÃO aplicar")))
} else {
  report_lines <- "Coluna 'batch' ausente nos metadados — ComBat-Seq não aplicado."
}

writeLines(report_lines, file.path(opt$outdir, "batch_report.txt"))
cat(paste(report_lines, collapse = "\n"), "\n")

# ── ComBat-Seq ────────────────────────────────────────────────
if (apply_combat) {
  corrected <- tryCatch(
    ComBat_seq(counts, batch = batch_vec, group = meta$condition),
    error = function(e) { message("ComBat-Seq falhou: ", e$message); counts }
  )
} else {
  corrected <- counts
}

corrected_df <- as.data.frame(corrected) |> tibble::rownames_to_column("gene_id")
write_tsv(corrected_df, file.path(opt$outdir, "counts_corrected.tsv"))

post <- pca_diag(corrected, meta, "PCA – Após correção")
ggsave(file.path(opt$outdir, "pca_after_batch.pdf"), post$plot, width = 7, height = 5)

cat(sprintf("Contagens corrigidas: %d genes × %d amostras\n",
            nrow(corrected), ncol(corrected)))
cat("Concluído.\n")
