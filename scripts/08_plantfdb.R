#!/usr/bin/env Rscript
# ============================================================
# 08_plantfdb.R – Classificação de TFs DEGs via PlantTFDB
# Glycine max: prefixo Gma_ (não Ath_)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

opt_list <- list(
  make_option("--deseq2_sig",    type = "character", help = "DEGs significativos (TSV)"),
  make_option("--deseq2_all",    type = "character", help = "Todos os genes testados (TSV)"),
  make_option("--plantfdb_file", type = "character", default = NULL,
              help = "Gma_TF_list.txt.gz local (opcional; baixa automaticamente se omitido)"),
  make_option("--outdir",        type = "character", default = "."),
  make_option("--figures_dir",   type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "tf_deg_classified.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "tf_family_summary.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "tf_family_enrichment.tsv"))
writeLines("", file.path(opt$outdir, "plantfdb_summary.txt"))
writeLines("", file.path(opt$figures_dir, ".plantfdb_placeholder"))

cat("═══════════════════════════════════════\n")
cat("  PlantTFDB – TFs em Glycine max\n")
cat("═══════════════════════════════════════\n")

# ── 1. Carregar PlantTFDB ─────────────────────────────────────
load_plantfdb <- function(local_file = NULL) {
  PLANTFDB_URL <- "http://planttfdb.gao-lab.org/download/TF_list/Gma_TF_list.txt.gz"
  tmp_gz <- file.path(tempdir(), "Gma_TF_list.txt.gz")

  if (!is.null(local_file) && file.exists(local_file)) {
    cat(sprintf("Carregando PlantTFDB local: %s\n", local_file))
    return(read_tsv(local_file, col_names = c("tf_id", "family", "species"),
                    show_col_types = FALSE, comment = "#"))
  }
  cat("Baixando PlantTFDB Gma...\n")
  result <- tryCatch({
    download.file(PLANTFDB_URL, tmp_gz, quiet = TRUE, method = "auto")
    read_tsv(tmp_gz, col_names = c("tf_id", "family", "species"),
             show_col_types = FALSE, comment = "#")
  }, error = function(e) {
    message("Download PlantTFDB falhou: ", e$message)
    NULL
  })
  result
}

tfdb <- load_plantfdb(opt$plantfdb_file)

if (is.null(tfdb) || nrow(tfdb) == 0) {
  writeLines("PlantTFDB não disponível.", file.path(opt$outdir, "plantfdb_summary.txt"))
  cat("Finalizado sem PlantTFDB.\n")
  quit(status = 0)
}

tfdb$tf_id <- gsub("\\.[0-9]+$", "", tfdb$tf_id)
cat(sprintf("TFs em PlantTFDB Gma: %d (famílias: %d)\n",
            nrow(tfdb), length(unique(tfdb$family))))

# ── 2. DEGs × TFs ────────────────────────────────────────────
degs_sig <- read_tsv(opt$deseq2_sig, show_col_types = FALSE)
degs_all <- read_tsv(opt$deseq2_all, show_col_types = FALSE)
degs_sig$gene_id <- gsub("\\.[0-9]+$", "", degs_sig$gene_id)
degs_all$gene_id <- gsub("\\.[0-9]+$", "", degs_all$gene_id)

tf_degs <- degs_sig |>
  inner_join(tfdb |> select(tf_id, family), by = c("gene_id" = "tf_id"))
cat(sprintf("TFs DEGs: %d\n", nrow(tf_degs)))

write_tsv(tf_degs, file.path(opt$outdir, "tf_deg_classified.tsv"))

if (nrow(tf_degs) == 0) {
  writeLines("Nenhum TF encontrado nos DEGs.", file.path(opt$outdir, "plantfdb_summary.txt"))
  quit(status = 0)
}

# ── 3. Resumo por família ─────────────────────────────────────
fam_sum <- tf_degs |>
  group_by(family) |>
  summarise(
    n_tfs  = n(),
    n_up   = sum(log2FoldChange > 0, na.rm = TRUE),
    n_down = sum(log2FoldChange < 0, na.rm = TRUE),
    mean_lfc = mean(log2FoldChange, na.rm = TRUE)
  ) |>
  arrange(desc(n_tfs))
write_tsv(fam_sum, file.path(opt$outdir, "tf_family_summary.tsv"))

# ── 4. Enriquecimento de famílias (Fisher) ────────────────────
all_tfs     <- unique(tfdb$tf_id)
bg_tfs      <- intersect(degs_all$gene_id, all_tfs)
deg_tfs_ids <- unique(tf_degs$gene_id)

fam_enrichment <- lapply(unique(tfdb$family), function(fam) {
  fam_ids    <- tfdb$tf_id[tfdb$family == fam]
  a  <- sum(deg_tfs_ids %in% fam_ids)
  b  <- length(deg_tfs_ids) - a
  c_ <- sum(bg_tfs %in% fam_ids) - a
  d  <- length(bg_tfs) - a - b - c_
  if ((a + b) == 0 || (a + c_) == 0) return(NULL)
  p  <- fisher.test(matrix(c(a, b, c_, d), 2, 2), alternative = "greater")$p.value
  data.frame(family = fam, n_deg = a, n_bg = a + c_, pvalue = p)
}) |> bind_rows()

if (nrow(fam_enrichment) > 0) {
  fam_enrichment$padj <- p.adjust(fam_enrichment$pvalue, method = "BH")
  write_tsv(fam_enrichment |> arrange(padj),
            file.path(opt$outdir, "tf_family_enrichment.tsv"))
}

# ── 5. Figuras ────────────────────────────────────────────────
top_fam <- fam_sum |> slice_head(n = 15)

p_bar <- ggplot(top_fam, aes(reorder(family, n_tfs), n_tfs)) +
  geom_col(fill = "#1a9641") + coord_flip() +
  labs(title = "Famílias TF – DEGs Glycine max",
       x = NULL, y = "N° TFs DEGs") + theme_bw()
ggsave(file.path(opt$figures_dir, "tf_families_barplot.pdf"), p_bar, width = 8, height = 6)
ggsave(file.path(opt$figures_dir, "tf_families_barplot.png"), p_bar, width = 8, height = 6,
       dpi = 300)

p_updown <- tf_degs |>
  mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down")) |>
  count(family, direction) |>
  filter(family %in% top_fam$family) |>
  ggplot(aes(reorder(family, n), n, fill = direction)) +
  geom_col(position = "dodge") + coord_flip() +
  scale_fill_manual(values = c(Up = "#d7191c", Down = "#2c7bb6")) +
  labs(title = "TFs por direção e família", x = NULL, y = "N°") + theme_bw()
ggsave(file.path(opt$figures_dir, "tf_families_updown.pdf"), p_updown, width = 8, height = 6)
ggsave(file.path(opt$figures_dir, "tf_families_updown.png"), p_updown, width = 8, height = 6,
       dpi = 300)

writeLines(c(
  sprintf("TFs totais Gma PlantTFDB: %d", nrow(tfdb)),
  sprintf("TFs DEGs: %d", nrow(tf_degs)),
  sprintf("Famílias representadas: %d", length(unique(tf_degs$family))),
  sprintf("Top família: %s (%d TFs)", fam_sum$family[1], fam_sum$n_tfs[1])
), file.path(opt$outdir, "plantfdb_summary.txt"))

cat("Concluído.\n")
