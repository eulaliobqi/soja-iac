#!/usr/bin/env Rscript
# ============================================================
# 10_metanalysis.R – Validação cruzada de DEGs via GEO/SRA
# Glycine max: datasets relevantes GSE99698, GSE107900, GSE143156
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GEOquery)
  library(limma)
  library(readr)
  library(dplyr)
})

opt_list <- list(
  make_option("--deseq2_sig",     type = "character", help = "DEGs sig (TSV)"),
  make_option("--geo_accessions", type = "character",
              default = "GSE99698,GSE107900",
              help = "Accessions GEO separados por vírgula"),
  make_option("--outdir",         type = "character", default = "."),
  make_option("--figures_dir",    type = "character", default = "figures"),
  make_option("--padj_meta",      type = "double",    default = 0.05)
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "metanalysis_validated_genes.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "metanalysis_overlap.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "metanalysis_summary.tsv"))
writeLines("", file.path(opt$outdir, "metanalysis_report.txt"))
writeLines("", file.path(opt$figures_dir, ".meta_placeholder"))

cat("═══════════════════════════════════════\n")
cat("  Metanálise – Validação cruzada GEO\n")
cat("═══════════════════════════════════════\n")

our_degs <- read_tsv(opt$deseq2_sig, show_col_types = FALSE)
our_degs$gene_id <- gsub("\\.[0-9]+$", "", our_degs$gene_id)
our_deg_set <- unique(our_degs$gene_id)
cat(sprintf("Nossos DEGs: %d\n", length(our_deg_set)))

accessions <- strsplit(opt$geo_accessions, ",")[[1]] |> trimws()

validated_by <- list()
summary_rows <- list()

for (acc in accessions) {
  cat(sprintf("\nProcessando %s...\n", acc))

  geo <- tryCatch(
    getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) { message("GEO falhou para ", acc, ": ", e$message); NULL }
  )
  if (is.null(geo)) next

  eset <- geo[[1]]
  expr <- exprs(eset)

  # Normaliza por quartil se necessário
  if (max(expr, na.rm = TRUE) > 100) {
    expr <- normalizeBetweenArrays(expr, method = "quantile")
  }

  pd <- pData(eset)
  if (!"condition" %in% names(pd) && "characteristics_ch1" %in% names(pd)) {
    pd$condition <- pd$characteristics_ch1
  }
  if (!"condition" %in% names(pd)) {
    cat(sprintf("  %s: coluna condition não detectada — pulando\n", acc))
    next
  }

  grps <- as.factor(pd$condition)
  if (nlevels(grps) < 2) {
    cat(sprintf("  %s: menos de 2 grupos — pulando\n", acc))
    next
  }

  design_mat <- model.matrix(~0 + grps)
  fit  <- lmFit(expr, design_mat)
  cont <- makeContrasts(contrasts = paste0(levels(grps)[2], "-", levels(grps)[1]),
                        levels = design_mat)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)

  tt <- topTable(fit2, number = Inf, adjust.method = "BH")
  tt$gene_id <- gsub("\\.[0-9]+$", "", rownames(tt))

  meta_degs  <- unique(tt$gene_id[!is.na(tt$adj.P.Val) & tt$adj.P.Val < opt$padj_meta])
  overlap    <- intersect(our_deg_set, meta_degs)

  cat(sprintf("  DEGs %s: %d | Overlap: %d\n", acc, length(meta_degs), length(overlap)))

  for (g in overlap) validated_by[[g]] <- c(validated_by[[g]], acc)

  summary_rows[[acc]] <- data.frame(
    accession   = acc,
    n_degs_geo  = length(meta_degs),
    n_overlap   = length(overlap),
    pct_overlap = round(length(overlap) / max(length(our_deg_set), 1) * 100, 1)
  )
}

# ── Resultados ────────────────────────────────────────────────
if (length(validated_by) > 0) {
  validated_df <- data.frame(
    gene_id     = names(validated_by),
    n_validated = sapply(validated_by, length),
    datasets    = sapply(validated_by, paste, collapse = ",")
  ) |> arrange(desc(n_validated))
  write_tsv(validated_df, file.path(opt$outdir, "metanalysis_validated_genes.tsv"))

  overlap_df <- data.frame(
    gene_id  = names(validated_by),
    datasets = sapply(validated_by, paste, collapse = ",")
  )
  write_tsv(overlap_df, file.path(opt$outdir, "metanalysis_overlap.tsv"))
  cat(sprintf("\nGenes validados em ≥1 dataset GEO: %d\n", nrow(validated_df)))
}

if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  write_tsv(summary_df, file.path(opt$outdir, "metanalysis_summary.tsv"))
}

writeLines(c(
  sprintf("Datasets processados: %d", length(summary_rows)),
  sprintf("Genes validados: %d", length(validated_by))
), file.path(opt$outdir, "metanalysis_report.txt"))

cat("Concluído.\n")
