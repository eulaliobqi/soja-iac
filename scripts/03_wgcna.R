#!/usr/bin/env Rscript
# ============================================================
# 03_wgcna.R – Co-expressão de genes com WGCNA
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
})

# Permite paralelismo interno do WGCNA
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ── Argumentos ───────────────────────────────────────────────
opt_list <- list(
  make_option("--norm_counts", type = "character", help = "Contagens normalizadas VST (TSV)"),
  make_option("--metadata",    type = "character", help = "Metadados das amostras (TSV)"),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures"),
  make_option("--min_genes",   type = "integer",   default = 5000,
              help = "Mínimo de genes com variância para WGCNA"),
  make_option("--soft_power",  type = "integer",   default = 0,
              help = "Soft thresholding power (0 = auto-detect)")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("═══════════════════════════════════════\n")
cat("  WGCNA – Co-expressão Gênica\n")
cat("═══════════════════════════════════════\n")

# ── 1. Leitura dos dados ──────────────────────────────────────
norm_df <- read_tsv(opt$norm_counts, show_col_types = FALSE) %>%
  column_to_rownames("gene_id")

meta <- read_tsv(opt$metadata, show_col_types = FALSE) %>%
  column_to_rownames("sample")

# Transpõe: amostras = linhas, genes = colunas (formato WGCNA)
expr_mat <- t(as.matrix(norm_df))
cat(sprintf("Expressão: %d amostras x %d genes\n", nrow(expr_mat), ncol(expr_mat)))

# ── 2. Filtragem por variância ────────────────────────────────
gene_vars <- apply(expr_mat, 2, var)
n_keep <- min(opt$min_genes, ncol(expr_mat))
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:n_keep])
expr_mat <- expr_mat[, top_var_genes]
cat(sprintf("Genes selecionados por variância: %d\n", ncol(expr_mat)))

# ── 3. Detecção de outliers de amostras ──────────────────────
sample_tree <- hclust(dist(expr_mat), method = "average")
pdf(file.path(opt$figures_dir, "sample_clustering_tree.pdf"), width = 10, height = 6)
plot(sample_tree, main = "Dendrograma de Amostras", xlab = "", sub = "",
     cex = 0.8, cex.main = 1.2)
abline(h = quantile(dist(expr_mat), 0.95), col = "red", lty = 2)
dev.off()
cat("Dendrograma de amostras gerado.\n")

# Verificação de dados faltantes
gsg <- goodSamplesGenes(expr_mat, verbose = 0)
if (!gsg$allOK) {
  expr_mat <- expr_mat[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("Após limpeza: %d amostras x %d genes\n", nrow(expr_mat), ncol(expr_mat)))
}

# ── 4. Escolha do soft thresholding power ────────────────────
if (opt$soft_power == 0) {
  cat("Detectando soft thresholding power ótimo...\n")
  powers <- c(1:10, seq(12, 30, 2))

  sft <- pickSoftThreshold(
    expr_mat,
    powerVector  = powers,
    networkType  = "signed hybrid",
    verbose      = 0
  )

  # Plota a análise de soft thresholding
  pdf(file.path(opt$figures_dir, "soft_threshold.pdf"), width = 12, height = 5)
  par(mfrow = c(1, 2))

  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit (signed R²)",
       main = "Scale independence",
       type = "n")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = 0.9, col = "red")
  abline(h = 0.8, col = "red", lty = 2)

  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       main = "Mean connectivity", type = "n")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels = powers, cex = 0.9, col = "red")
  dev.off()

  # Seleciona power onde R² > 0.8
  soft_power <- sft$powerEstimate
  if (is.na(soft_power)) {
    soft_power <- 12  # default conservador
    message("Power automático não encontrado, usando power = 12")
  }
} else {
  soft_power <- opt$soft_power
}
cat(sprintf("Soft thresholding power: %d\n", soft_power))

# ── 5. Construção da rede e detecção de módulos ───────────────
cat("Construindo rede de co-expressão (pode demorar)...\n")
net <- blockwiseModules(
  expr_mat,
  power               = soft_power,
  networkType         = "signed hybrid",
  TOMType             = "signed",
  minModuleSize       = 30,
  mergeCutHeight      = 0.25,
  deepSplit           = 2,
  numericLabels       = FALSE,
  pamRespectsDendro   = FALSE,
  maxBlockSize        = ncol(expr_mat),
  saveTOMs            = FALSE,
  verbose             = 0,
  nThreads            = parallel::detectCores() - 1
)

module_colors <- net$colors
gene_modules  <- data.frame(
  gene_id = names(module_colors),
  module  = module_colors,
  stringsAsFactors = FALSE
)

n_modules <- length(unique(module_colors[module_colors != "grey"]))
cat(sprintf("Módulos detectados: %d\n", n_modules))
cat(sprintf("Genes em módulos: %d | Não atribuídos (grey): %d\n",
            sum(module_colors != "grey"), sum(module_colors == "grey")))

# ── 6. Eigengenes e correlação com fenótipo ───────────────────
ME_list <- moduleEigengenes(expr_mat, module_colors)
MEs     <- ME_list$eigengenes
MEs     <- orderMEs(MEs)

# Trait matrix: condition como binário
trait_mat <- data.frame(
  condition_binary = ifelse(meta[rownames(expr_mat), "condition"] ==
                              unique(meta$condition)[2], 1, 0),
  row.names = rownames(expr_mat)
)

# Correlação ME-trait
ME_trait_cor  <- cor(MEs, trait_mat, use = "pairwise.complete.obs")
ME_trait_pval <- corPvalueStudent(ME_trait_cor, nrow(expr_mat))

# Plota heatmap ME-trait
pdf(file.path(opt$figures_dir, "module_trait_heatmap.pdf"), width = 8, height = max(6, n_modules * 0.4))
pheatmap(ME_trait_cor,
         display_numbers = matrix(
           ifelse(ME_trait_pval < 0.05, "*", ""),
           nrow = nrow(ME_trait_cor)
         ),
         color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         border_color = NA,
         cluster_rows = nrow(ME_trait_cor) > 1,
         cluster_cols = ncol(ME_trait_cor) > 1,
         main = "Correlação: Módulos × Fenótipo\n(* p < 0.05)",
         fontsize = 10)
dev.off()

# ── 7. Hub genes ─────────────────────────────────────────────
cat("Identificando hub genes...\n")
gene_module_mem <- cor(expr_mat, MEs, use = "pairwise.complete.obs")

hub_genes_list <- lapply(unique(module_colors[module_colors != "grey"]), function(mod) {
  me_col  <- paste0("ME", mod)
  if (!me_col %in% colnames(gene_module_mem)) return(NULL)
  mod_genes <- names(module_colors)[module_colors == mod]
  mem       <- gene_module_mem[mod_genes, me_col]
  df <- data.frame(
    gene_id    = mod_genes,
    module     = mod,
    membership = mem,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(abs(membership)))
  head(df, 10)
})

hub_genes <- bind_rows(hub_genes_list)

# ── 8. Exporta resultados ─────────────────────────────────────
write_tsv(gene_modules, file.path(opt$outdir, "wgcna_modules.tsv"))
write_tsv(hub_genes,    file.path(opt$outdir, "wgcna_hub_genes.tsv"))
write_tsv(
  as.data.frame(MEs) %>% rownames_to_column("sample"),
  file.path(opt$outdir, "wgcna_eigengenes.tsv")
)

# Sumário de módulos
module_summary <- gene_modules %>%
  group_by(module) %>%
  summarise(n_genes = n()) %>%
  arrange(desc(n_genes))
write_tsv(module_summary, file.path(opt$outdir, "wgcna_module_summary.tsv"))

# Dendrograma de genes e módulos
pdf(file.path(opt$figures_dir, "gene_dendrogram_modules.pdf"), width = 12, height = 5)
plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
                    "Módulo", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Dendrograma de Genes e Módulos")
dev.off()

cat("\nWGCNA concluído. Resultados em:", opt$outdir, "\n")
