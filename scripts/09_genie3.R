#!/usr/bin/env Rscript
# ============================================================
# 09_genie3.R – Inferência de rede regulatória TF→gene (GENIE3)
# Sempre usa run_genie3() com fallback serial (Lição 12)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GENIE3)
  library(readr)
  library(dplyr)
})

opt_list <- list(
  make_option("--norm_counts",   type = "character", help = "Contagens normalizadas (TSV)"),
  make_option("--tf_classified", type = "character", help = "tf_deg_classified.tsv do PlantTFDB"),
  make_option("--n_trees",       type = "integer",   default = 500L),
  make_option("--n_links",       type = "integer",   default = 5000L),
  make_option("--ncores",        type = "integer",   default = 4L),
  make_option("--outdir",        type = "character", default = "."),
  make_option("--figures_dir",   type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs (Lição 13) ─────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "genie3_network.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "genie3_hub_tfs.tsv"))
writeLines("", file.path(opt$outdir, "genie3_summary.txt"))
writeLines("", file.path(opt$figures_dir, ".genie3_placeholder"))

cat("═══════════════════════════════════════\n")
cat("  GENIE3 – Rede Regulatória TF→gene\n")
cat("═══════════════════════════════════════\n")

norm <- read_tsv(opt$norm_counts, show_col_types = FALSE) |>
  tibble::column_to_rownames("gene_id")
rownames(norm) <- gsub("\\.[0-9]+$", "", rownames(norm))

# Reguladores: TFs DEGs presentes na matriz
tf_data    <- read_tsv(opt$tf_classified, show_col_types = FALSE)
tf_data$gene_id <- gsub("\\.[0-9]+$", "", tf_data$gene_id)
regulators <- intersect(unique(tf_data$gene_id), rownames(norm))
targets    <- rownames(norm)

cat(sprintf("TFs na matriz: %d | Genes alvo: %d\n", length(regulators), length(targets)))

if (length(regulators) < 2) {
  writeLines("Reguladores insuficientes para GENIE3.", file.path(opt$outdir, "genie3_summary.txt"))
  cat("Finalizado sem rede.\n")
  quit(status = 0)
}

expr_mat <- as.matrix(norm)

# ── run_genie3 com fallback serial (Lição 12) ─────────────────
run_genie3 <- function(expr_mat, regulators, targets, n_trees, ncores) {
  tryCatch(
    GENIE3(exprMatrix = expr_mat, regulators = regulators, targets = targets,
           treeMethod = "RF", K = "sqrt", nTrees = n_trees, nCores = ncores),
    error = function(e) {
      message(sprintf("GENIE3 (nCores=%d) falhou: %s", ncores, e$message))
      if (ncores > 1L) {
        message("Repetindo com nCores=1...")
        tryCatch(
          GENIE3(exprMatrix = expr_mat, regulators = regulators, targets = targets,
                 treeMethod = "RF", K = "sqrt", nTrees = n_trees, nCores = 1L),
          error = function(e2) { message("GENIE3 serial falhou: ", e2$message); NULL }
        )
      } else NULL
    }
  )
}

weight_mat <- run_genie3(expr_mat, regulators, targets, opt$n_trees, opt$ncores)

if (is.null(weight_mat)) {
  writeLines("GENIE3 falhou em ambos os modos.", file.path(opt$outdir, "genie3_summary.txt"))
  quit(status = 0)
}

# ── Top links ────────────────────────────────────────────────
link_list <- getLinkList(weight_mat, threshold = 0, reportMax = opt$n_links)
names(link_list) <- c("regulator", "target", "weight")
write_tsv(link_list, file.path(opt$outdir, "genie3_network.tsv"))

# Hub TFs por soma de pesos de saída
hub_tfs <- link_list |>
  group_by(regulator) |>
  summarise(total_weight = sum(weight), n_targets = n()) |>
  arrange(desc(total_weight))
write_tsv(hub_tfs, file.path(opt$outdir, "genie3_hub_tfs.tsv"))

writeLines(c(
  sprintf("Links GENIE3: %d", nrow(link_list)),
  sprintf("Reguladores: %d", length(unique(link_list$regulator))),
  sprintf("Top TF: %s (peso=%.3f)", hub_tfs$regulator[1], hub_tfs$total_weight[1])
), file.path(opt$outdir, "genie3_summary.txt"))

cat(sprintf("Rede: %d links | Top TF hub: %s\n",
            nrow(link_list), hub_tfs$regulator[1]))
cat("Concluído.\n")
