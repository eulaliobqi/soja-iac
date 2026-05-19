#!/usr/bin/env Rscript
# ============================================================
# 11_lncrna.R – Predição de lncRNAs novos via Biostrings
# Critério: comprimento ≥200 nt + ORF máximo <100 aa
# Input: novel_transcripts.fa gerado pelo LNCRNA_PRED Nextflow
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

opt_list <- list(
  make_option("--fasta",       type = "character", help = "FASTA de transcritos novos"),
  make_option("--min_len_nt",  type = "integer",   default = 200L),
  make_option("--max_orf_aa",  type = "integer",   default = 100L),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "lncrna_all.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "lncrna_candidates.tsv"))
writeLines("", file.path(opt$outdir, "lncrna_summary.txt"))
writeLines("", file.path(opt$figures_dir, ".lncrna_placeholder"))

cat("═══════════════════════════════════════\n")
cat("  lncRNA Prediction – Biostrings\n")
cat("═══════════════════════════════════════\n")

if (!file.exists(opt$fasta)) {
  writeLines("FASTA não encontrado.", file.path(opt$outdir, "lncrna_summary.txt"))
  cat("FASTA ausente — finalizado.\n")
  quit(status = 0)
}

seqs <- readDNAStringSet(opt$fasta)
cat(sprintf("Transcritos novos: %d\n", length(seqs)))

# ── Máximo ORF por sequência ──────────────────────────────────
max_orf_aa <- function(seq) {
  max_len <- 0L
  for (frame in 0:2) {
    s <- subseq(seq, start = frame + 1)
    prot <- suppressWarnings(translate(s, no.init.codon = TRUE))
    parts <- strsplit(as.character(prot), "\\*")[[1]]
    lens  <- nchar(parts)
    if (length(lens) > 0) max_len <- max(max_len, max(lens))
  }
  max_len
}

all_results <- lapply(seq_along(seqs), function(i) {
  id  <- names(seqs)[i]
  seq <- seqs[[i]]
  len <- length(seq)
  orf <- tryCatch(max_orf_aa(seq), error = function(e) NA_integer_)
  data.frame(transcript_id = id, length_nt = len, max_orf_aa = orf)
}) |> bind_rows()

write_tsv(all_results, file.path(opt$outdir, "lncrna_all.tsv"))

candidates <- all_results |>
  filter(length_nt >= opt$min_len_nt, max_orf_aa < opt$max_orf_aa)
write_tsv(candidates, file.path(opt$outdir, "lncrna_candidates.tsv"))

cat(sprintf("lncRNA candidatos: %d / %d\n", nrow(candidates), nrow(all_results)))

# ── Figuras ───────────────────────────────────────────────────
p_len <- ggplot(all_results, aes(length_nt)) +
  geom_histogram(bins = 40, fill = "#7b2d8b") +
  geom_vline(xintercept = opt$min_len_nt, linetype = "dashed", color = "red") +
  labs(title = "Distribuição de comprimento – transcritos novos",
       x = "Comprimento (nt)", y = "Frequência") + theme_bw()
ggsave(file.path(opt$figures_dir, "lncrna_length_dist.pdf"), p_len, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "lncrna_length_dist.png"), p_len, width = 7, height = 5,
       dpi = 300)

all_results$class <- ifelse(
  all_results$length_nt >= opt$min_len_nt & all_results$max_orf_aa < opt$max_orf_aa,
  "lncRNA candidato", "Descartado"
)
p_scatter <- ggplot(all_results, aes(length_nt, max_orf_aa, color = class)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("lncRNA candidato" = "#7b2d8b", Descartado = "grey60")) +
  labs(title = "lncRNA: comprimento vs ORF máximo",
       x = "Comprimento (nt)", y = "ORF máximo (aa)") + theme_bw()
ggsave(file.path(opt$figures_dir, "lncrna_scatter.pdf"), p_scatter, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "lncrna_scatter.png"), p_scatter, width = 7, height = 5,
       dpi = 300)

writeLines(c(
  sprintf("Transcritos novos analisados: %d", nrow(all_results)),
  sprintf("lncRNA candidatos (≥%d nt, <%d aa): %d",
          opt$min_len_nt, opt$max_orf_aa, nrow(candidates))
), file.path(opt$outdir, "lncrna_summary.txt"))

cat("Concluído.\n")
