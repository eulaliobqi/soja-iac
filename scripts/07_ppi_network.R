#!/usr/bin/env Rscript
# ============================================================
# 07_ppi_network.R – PPI via STRINGdb (espécie 3847 = Glycine max)
# Outputs inicializados antes de tentar conexão (sem internet → ok)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(STRINGdb)
  library(igraph)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

opt_list <- list(
  make_option("--deseq2_sig",  type = "character", help = "DEGs significativos (TSV)"),
  make_option("--score_thr",   type = "integer",   default = 400L),
  make_option("--hub_n",       type = "integer",   default = 20L),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures"),
  make_option("--cache_dir",   type = "character", default = "/tmp/stringdb_gmax")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$cache_dir,   showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs (Lição 13 + 14) ───────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "ppi_edges.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "ppi_nodes.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "hub_genes.tsv"))
writeLines("", file.path(opt$outdir, "network_summary.txt"))
writeLines("", file.path(opt$figures_dir, ".ppi_placeholder"))

cat("═══════════════════════════════════════\n")
cat("  PPI Network – STRINGdb (Glycine max 3847)\n")
cat("═══════════════════════════════════════\n")

degs <- read_tsv(opt$deseq2_sig, show_col_types = FALSE)
degs$gene_id <- gsub("\\.[0-9]+$", "", degs$gene_id)

# ── STRINGdb (pode falhar sem internet) ──────────────────────
ppi_result <- tryCatch({
  string_db <- STRINGdb$new(
    version         = "11.5",
    species         = 3847L,     # Glycine max
    score_threshold = opt$score_thr,
    input_directory = opt$cache_dir
  )

  deg_df <- data.frame(gene = degs$gene_id, lfc = degs$log2FoldChange)
  mapped <- string_db$map(deg_df, "gene", removeUnmappedRows = FALSE)
  mapped <- mapped[!is.na(mapped$STRING_id), ]
  cat(sprintf("DEGs mapeados no STRING: %d / %d\n", nrow(mapped), nrow(degs)))

  if (nrow(mapped) < 2) stop("Genes insuficientes mapeados")

  interactions <- string_db$get_interactions(mapped$STRING_id)
  list(mapped = mapped, interactions = interactions)
}, error = function(e) {
  message("STRINGdb falhou: ", e$message)
  NULL
})

if (is.null(ppi_result)) {
  writeLines("STRINGdb não disponível (sem internet ou mapeamento insuficiente).",
             file.path(opt$outdir, "network_summary.txt"))
  cat("Finalizado sem rede PPI.\n")
  quit(status = 0)
}

mapped       <- ppi_result$mapped
interactions <- ppi_result$interactions

edges_df <- interactions |>
  rename(from = from, to = to, score = combined_score)
write_tsv(edges_df, file.path(opt$outdir, "ppi_edges.tsv"))

# ── Grafo igraph ──────────────────────────────────────────────
g <- graph_from_data_frame(edges_df, directed = FALSE,
                            vertices = mapped[, c("STRING_id", "gene")])
V(g)$degree <- degree(g)

nodes_df <- data.frame(
  string_id = V(g)$name,
  gene      = V(g)$gene,
  degree    = V(g)$degree
) |> arrange(desc(degree))
write_tsv(nodes_df, file.path(opt$outdir, "ppi_nodes.tsv"))

hub_df <- nodes_df |> slice_head(n = opt$hub_n)
write_tsv(hub_df, file.path(opt$outdir, "hub_genes.tsv"))

writeLines(c(
  sprintf("Nodes: %d | Edges: %d", vcount(g), ecount(g)),
  sprintf("Hub top-%d: %s", opt$hub_n, paste(hub_df$gene[1:min(5, nrow(hub_df))], collapse = ", "))
), file.path(opt$outdir, "network_summary.txt"))

# ── Figuras ───────────────────────────────────────────────────
p_deg <- ggplot(nodes_df, aes(degree)) +
  geom_histogram(bins = 30, fill = "#2c7bb6") +
  labs(title = "Distribuição de grau – PPI Glycine max",
       x = "Grau", y = "Frequência") + theme_bw()

ggsave(file.path(opt$figures_dir, "ppi_degree_dist.pdf"), p_deg, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "ppi_degree_dist.png"), p_deg, width = 7, height = 5,
       dpi = 300)

# Sub-rede dos hubs
hub_ids  <- hub_df$string_id
subg     <- induced_subgraph(g, v = hub_ids)
deg_sub  <- degree(subg)

if (vcount(subg) > 1) {
  hub_nodes <- data.frame(
    gene   = V(subg)$gene,
    degree = deg_sub
  ) |> arrange(desc(degree)) |> slice_head(n = 20)

  p_hub <- ggplot(hub_nodes, aes(reorder(gene, degree), degree)) +
    geom_col(fill = "#d7191c") + coord_flip() +
    labs(title = sprintf("Top %d Hub Genes – PPI", min(20, nrow(hub_nodes))),
         x = NULL, y = "Grau") + theme_bw()
  ggsave(file.path(opt$figures_dir, "ppi_hub_genes.pdf"), p_hub, width = 7, height = 6)
  ggsave(file.path(opt$figures_dir, "ppi_hub_genes.png"), p_hub, width = 7, height = 6,
         dpi = 300)
}

cat(sprintf("Rede: %d nós | %d arestas | Hub genes: %d\n",
            vcount(g), ecount(g), nrow(hub_df)))
cat("Concluído.\n")
