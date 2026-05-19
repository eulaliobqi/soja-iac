#!/usr/bin/env Rscript
# ============================================================
# 02b_goseq.R – GO com correção de viés de tamanho (GOseq)
# Glycine max: NUNCA usar genome/id internos (não existem)
# Sempre gene2cat via biomaRt ou arquivo de anotação
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(goseq)
  library(readr)
  library(dplyr)
})

opt_list <- list(
  make_option("--deseq2",   type = "character", help = "DESeq2 results all (TSV)"),
  make_option("--gtf",      type = "character", default = NULL),
  make_option("--go_annot", type = "character", default = NULL,
              help = "TSV gene_id/go_id/namespace (opcional; sem isso usa biomaRt)"),
  make_option("--padj",     type = "double",  default = 0.05),
  make_option("--lfc",      type = "double",  default = 1.0),
  make_option("--outdir",   type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ────────────────────────────────────────
for (f in c("goseq_bp_results.tsv", "goseq_mf_results.tsv", "goseq_cc_results.tsv")) {
  write_tsv(data.frame(), file.path(opt$outdir, f))
}

cat("═══════════════════════════════════════\n")
cat("  GOseq – GO com correção viés de tamanho\n")
cat("═══════════════════════════════════════\n")

# ── 1. Leitura DESeq2 ─────────────────────────────────────────
res <- read_tsv(opt$deseq2, show_col_types = FALSE)
res$gene_id <- gsub("\\.[0-9]+$", "", res$gene_id)

is_deg <- as.integer(!is.na(res$padj) &
                       res$padj < opt$padj &
                       abs(res$log2FoldChange) > opt$lfc)
names(is_deg) <- res$gene_id
cat(sprintf("DEGs: %d / %d\n", sum(is_deg), length(is_deg)))

# ── 2. Comprimentos de gene ───────────────────────────────────
make_txdb <- function(gtf) {
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdbmaker::makeTxDbFromGFF(gtf, format = "GTF")
  } else {
    GenomicFeatures::makeTxDbFromGFF(gtf, format = "GTF")
  }
}

gene_lengths <- setNames(rep(1000L, length(is_deg)), names(is_deg))

if (!is.null(opt$gtf) && file.exists(opt$gtf)) {
  suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(BiocGenerics)
  })
  txdb <- tryCatch(make_txdb(opt$gtf), error = function(e) {
    message("txdb falhou: ", e$message); NULL
  })
  if (!is.null(txdb)) {
    exons_by_gene <- exonsBy(txdb, by = "gene")
    lens <- sum(width(reduce(exons_by_gene)))
    lens_clean <- gsub("\\.[0-9]+$", "", names(lens))
    names(lens) <- lens_clean
    shared <- intersect(names(is_deg), names(lens))
    gene_lengths[shared] <- lens[shared]
    cat(sprintf("Comprimentos reais: %d genes\n", length(shared)))
  }
}

# ── 3. Anotações GO ───────────────────────────────────────────
get_gene2go <- function(all_genes, ontology, go_annot_file = NULL) {
  ont_map <- c(BP = "biological_process", MF = "molecular_function",
               CC = "cellular_component")

  if (!is.null(go_annot_file) && file.exists(go_annot_file)) {
    ann <- read_tsv(go_annot_file, show_col_types = FALSE)
    ann$gene_id <- gsub("\\.[0-9]+$", "", ann$gene_id)
    sub <- ann[ann$namespace == ont_map[ontology] & !is.na(ann$go_id), ]
  } else {
    cat(sprintf("Buscando GO %s via biomaRt...\n", ontology))
    mart <- tryCatch(
      biomaRt::useMart("plants_mart", dataset = "gmax_eg_gene",
                        host = "https://plants.ensembl.org"),
      error = function(e) { message("biomaRt: ", e$message); NULL }
    )
    if (is.null(mart)) return(list())
    go_df <- tryCatch(
      biomaRt::getBM(
        attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
        filters    = "ensembl_gene_id",
        values     = all_genes, mart = mart
      ),
      error = function(e) { message("getBM: ", e$message); data.frame() }
    )
    if (nrow(go_df) == 0) return(list())
    go_df$ensembl_gene_id <- gsub("\\.[0-9]+$", "", go_df$ensembl_gene_id)
    sub <- go_df[go_df$namespace_1003 == ont_map[ontology] & !is.na(go_df$go_id), ]
    names(sub)[1:2] <- c("gene_id", "go_id")
  }
  if (nrow(sub) == 0) return(list())
  split(sub$go_id, sub$gene_id)
}

# ── 4. Executa GOseq por ontologia ───────────────────────────
run_goseq_ont <- function(ont) {
  gene2cat <- get_gene2go(names(is_deg), ont, opt$go_annot)
  if (length(gene2cat) == 0) {
    cat(sprintf("GOseq %s: sem anotações disponíveis\n", ont))
    return(data.frame())
  }

  n_unique <- length(unique(gene_lengths))
  if (n_unique >= 6) {
    pwf <- tryCatch(
      nullp(is_deg, bias.data = gene_lengths, plot.fit = FALSE),
      error = function(e) { message("nullp falhou: ", e$message); NULL }
    )
    method <- "Wallenius"
  } else {
    cat(sprintf("Aviso: %d comprimentos únicos → Hypergeometric\n", n_unique))
    pwf <- data.frame(
      DEgenes   = is_deg,
      bias.data = gene_lengths,
      pwf       = sum(is_deg) / length(is_deg)
    )
    rownames(pwf) <- names(is_deg)
    method <- "Hypergeometric"
  }
  if (is.null(pwf)) return(data.frame())

  res_go <- tryCatch(
    goseq(pwf, gene2cat = gene2cat, method = method),
    error = function(e) { message("goseq ", ont, ": ", e$message); data.frame() }
  )
  if (nrow(res_go) == 0) return(data.frame())
  res_go$padj <- p.adjust(res_go$over_represented_pvalue, method = "BH")
  res_go |>
    filter(padj < 0.05) |>
    arrange(padj) |>
    mutate(ontology = ont)
}

suppressPackageStartupMessages(library(biomaRt))

bp <- run_goseq_ont("BP")
mf <- run_goseq_ont("MF")
cc <- run_goseq_ont("CC")

write_tsv(bp, file.path(opt$outdir, "goseq_bp_results.tsv"))
write_tsv(mf, file.path(opt$outdir, "goseq_mf_results.tsv"))
write_tsv(cc, file.path(opt$outdir, "goseq_cc_results.tsv"))

cat(sprintf("GOseq: BP=%d  MF=%d  CC=%d termos enriquecidos (padj<0.05)\n",
            nrow(bp), nrow(mf), nrow(cc)))
cat("Concluído.\n")
