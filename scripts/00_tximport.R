#!/usr/bin/env Rscript
# ============================================================
# 00_tximport.R – Importa quantificações Salmon (tx→gene)
# Glycine max RNASeq Pipeline
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(readr)
  library(dplyr)
})

opt_list <- list(
  make_option("--salmon_dir",  type = "character", help = "Dir com subpastas quant.sf por amostra"),
  make_option("--tx2gene",     type = "character", default = NULL, help = "TSV transcrito→gene"),
  make_option("--gtf",         type = "character", default = NULL, help = "GTF alternativo ao tx2gene"),
  make_option("--samplesheet", type = "character", help = "CSV com colunas sample,condition,path"),
  make_option("--outdir",      type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ─────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "salmon_counts.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "salmon_tpm.tsv"))

cat("═══════════════════════════════════════\n")
cat("  tximport – Salmon → contagens por gene\n")
cat("═══════════════════════════════════════\n")

# ── 1. Amostras ───────────────────────────────────────────────
ss <- read_csv(opt$samplesheet, show_col_types = FALSE)
quant_files <- file.path(opt$salmon_dir, ss$sample, "quant.sf")
names(quant_files) <- ss$sample

missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) stop("quant.sf faltando: ", paste(names(missing), collapse = ", "))

# ── 2. tx2gene ────────────────────────────────────────────────
if (!is.null(opt$tx2gene) && file.exists(opt$tx2gene)) {
  tx2gene <- read_tsv(opt$tx2gene, col_names = c("tx", "gene"), show_col_types = FALSE)
} else {
  # Inferir tx→gene diretamente dos IDs reais do quant.sf
  # Evita mismatch entre GTF txdb e FASTA headers gerados pelo gffread
  cat("Inferindo tx2gene dos IDs do quant.sf...\n")
  quant1  <- read_tsv(quant_files[1], show_col_types = FALSE)
  txids   <- quant1$Name
  # Glycine max Wm82: Glyma.01G000600.1.Wm82.a4.v1 → Glyma.01G000600
  geneids <- gsub("^(Glyma\\.[0-9]+G[0-9]+).*", "\\1", txids)
  # Fallback genérico para IDs que não casaram com o padrão Glyma
  unchanged <- geneids == txids
  geneids[unchanged] <- gsub("\\.[0-9]+$", "", txids[unchanged])
  tx2gene <- data.frame(tx = txids, gene = geneids)
}

cat(sprintf("Transcritos mapeados: %d → genes únicos: %d\n",
            nrow(tx2gene), length(unique(tx2gene$gene))))

# ── 3. tximport ───────────────────────────────────────────────
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene,
                ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)

counts_df <- as.data.frame(txi$counts) |>
  tibble::rownames_to_column("gene_id")
tpm_df    <- as.data.frame(txi$abundance) |>
  tibble::rownames_to_column("gene_id")

# Strip versão nos gene_ids resultantes
counts_df$gene_id <- gsub("\\.[0-9]+$", "", counts_df$gene_id)
tpm_df$gene_id    <- gsub("\\.[0-9]+$", "", tpm_df$gene_id)

write_tsv(counts_df, file.path(opt$outdir, "salmon_counts.tsv"))
write_tsv(tpm_df,    file.path(opt$outdir, "salmon_tpm.tsv"))

cat(sprintf("Genes: %d | Amostras: %d\n", nrow(counts_df), ncol(counts_df) - 1))
cat("Concluído.\n")
