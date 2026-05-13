#!/usr/bin/env Rscript
# ============================================================
# 02_enrichment.R – GO, KEGG e GSEA para Glycine max
# Estratégia de anotação (em ordem de prioridade):
#   1. Arquivo annotation_info do Phytozome (--go_annot)
#   2. Cache local (gmax_go_cache.rds)
#   3. Ensembl Plants via biomaRt (baixa e armazena cache)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# ── Argumentos ───────────────────────────────────────────────
opt_list <- list(
  make_option("--deseq2",      type = "character", help = "DESeq2 results (TSV, todos os genes)"),
  make_option("--norm_counts", type = "character", help = "Contagens normalizadas (TSV)"),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1.0),
  make_option("--organism",    type = "character", default = "gma"),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures"),
  make_option("--go_annot",    type = "character", default = NULL,
              help = "Arquivo annotation_info do Phytozome ou TSV gene_id/go_id/namespace")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("═══════════════════════════════════════\n")
cat("  Enriquecimento Funcional – GO + KEGG + GSEA\n")
cat("═══════════════════════════════════════\n")

# ── 1. Leitura ────────────────────────────────────────────────
res_all <- read_tsv(opt$deseq2, show_col_types = FALSE)

de_sig <- res_all %>%
  filter(!is.na(padj), padj < opt$padj, abs(log2FoldChange) > opt$lfc)

cat(sprintf("Genes DE: %d | Todos testados: %d\n", nrow(de_sig), nrow(res_all)))

# ── 2. Limpeza de IDs ─────────────────────────────────────────
# Remove sufixo de versão: Glyma.08G200100.Wm82.a2.v1 → Glyma.08G200100
clean_glyma <- function(ids) {
  sub("^(Glyma\\.[0-9A-Za-z]+G[0-9]+).*", "\\1", ids)
}

res_all <- res_all %>% mutate(gene_id_clean = clean_glyma(gene_id))
de_sig  <- de_sig  %>% mutate(gene_id_clean = clean_glyma(gene_id))

bg_ids <- unique(res_all$gene_id_clean)
de_ids <- unique(de_sig$gene_id_clean)

cat(sprintf("Background (únicos): %d | DE (únicos): %d\n", length(bg_ids), length(de_ids)))

# ── 3. Carregamento de anotações GO ───────────────────────────
# Retorna lista com:
#   TERM2GENE_BP/MF/CC: data.frame(go_id, gene_id_clean)
#   TERM2NAME:           data.frame(go_id, term)
#   entrez_map:          data.frame(gene_id_clean, entrezid) – pode ser NULL
CACHE_FILE <- file.path(opt$outdir, "gmax_go_cache.rds")

load_go_annotations <- function(go_annot_file = NULL) {

  # ── 3a. Arquivo fornecido pelo usuário ──────────────────────
  if (!is.null(go_annot_file) && file.exists(go_annot_file)) {
    cat(sprintf("Carregando anotações GO de: %s\n", go_annot_file))

    raw <- read_tsv(go_annot_file, show_col_types = FALSE, comment = "#")

    # Phytozome annotation_info.txt tem colunas fixas; detectamos automaticamente
    if (all(c("locusName", "GO") %in% names(raw))) {
      # Formato Phytozome
      tmp <- raw %>%
        select(gene_id_clean = locusName, go_raw = GO) %>%
        filter(!is.na(go_raw), go_raw != "")
      annot <- do.call(rbind, lapply(seq_len(nrow(tmp)), function(i) {
        gos <- str_trim(strsplit(tmp$go_raw[i], "\\|")[[1]])
        gos <- gos[grepl("^GO:", gos)]
        if (length(gos)) data.frame(gene_id_clean = tmp$gene_id_clean[i],
                                     go_id = gos, stringsAsFactors = FALSE)
      }))
    } else if (all(c("gene_id", "go_id") %in% names(raw))) {
      # Formato simples: gene_id / go_id [/ namespace]
      annot <- raw %>%
        rename(gene_id_clean = gene_id) %>%
        filter(grepl("^GO:", go_id))
    } else {
      stop("Formato de --go_annot não reconhecido. Esperado: locusName/GO (Phytozome) ",
           "ou gene_id/go_id.")
    }

    annot$gene_id_clean <- clean_glyma(annot$gene_id_clean)
    cat(sprintf("  %d anotações GO para %d genes\n",
                nrow(annot), length(unique(annot$gene_id_clean))))
    return(build_term_lists(annot, entrez_map = NULL))
  }

  # ── 3b. Cache local ─────────────────────────────────────────
  if (file.exists(CACHE_FILE)) {
    cat("Carregando cache de anotações GO...\n")
    return(readRDS(CACHE_FILE))
  }

  # ── 3c. Ensembl Plants via biomaRt ──────────────────────────
  cat("Buscando anotações GO via Ensembl Plants (biomaRt)...\n")

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "biomaRt não está instalado. Instale com:\n",
      "  mamba run -n r-analysis R -e 'BiocManager::install(\"biomaRt\")'\n",
      "Ou forneça --go_annot com o arquivo annotation_info do Phytozome."
    )
  }

  mart <- NULL
  for (host in c("https://plants.ensembl.org", "https://asia.ensembl.org",
                  "https://useast.ensembl.org")) {
    mart <- tryCatch(
      biomaRt::useMart("plants_mart", host = host),
      error = function(e) NULL
    )
    if (!is.null(mart)) { cat(sprintf("  Conectado: %s\n", host)); break }
  }
  if (is.null(mart)) stop("Não foi possível conectar ao Ensembl Plants.")

  datasets  <- biomaRt::listDatasets(mart)
  gmax_ds   <- datasets$dataset[grep("gmax|glycine|soy", datasets$dataset,
                                      ignore.case = TRUE)][1]
  if (is.na(gmax_ds)) stop("Dataset Glycine max não encontrado no Ensembl Plants.")
  cat(sprintf("  Dataset: %s\n", gmax_ds))

  gmax_mart <- biomaRt::useDataset(gmax_ds, mart = mart)

  cat("  Baixando anotações GO (pode demorar 2-5 minutos)...\n")
  go_raw <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "go_id", "namespace_1003", "entrezgene_id"),
    mart = gmax_mart
  )
  go_raw <- go_raw[go_raw$go_id != "", ]
  cat(sprintf("  %d anotações obtidas para %d genes Ensembl\n",
              nrow(go_raw), length(unique(go_raw$ensembl_gene_id))))

  # Normaliza IDs para correspondência (remove pontos/underscores, uppercase)
  norm <- function(x) toupper(gsub("[._-]", "", as.character(x)))

  ensembl_df <- data.frame(ensembl_id = unique(go_raw$ensembl_gene_id),
                            norm = norm(unique(go_raw$ensembl_gene_id)),
                            stringsAsFactors = FALSE)
  our_df <- data.frame(gene_id_clean = bg_ids,
                        norm = norm(bg_ids),
                        stringsAsFactors = FALSE)

  matched <- merge(our_df, ensembl_df, by = "norm")[, c("gene_id_clean", "ensembl_id")]
  n_matched <- length(unique(matched$gene_id_clean))
  cat(sprintf("  Genes mapeados: %d/%d (%.1f%%)\n",
              n_matched, length(bg_ids), 100 * n_matched / length(bg_ids)))

  if (n_matched == 0) stop("Nenhum gene Glyma mapeado. Verifique a versão do genoma.")

  # Junta: gene_id_clean → go_id, namespace, entrezid
  annot <- merge(matched, go_raw, by.x = "ensembl_id", by.y = "ensembl_gene_id") %>%
    select(gene_id_clean, go_id, namespace = namespace_1003, entrezid = entrezgene_id)

  entrez_map <- annot %>%
    filter(!is.na(entrezid)) %>%
    select(gene_id_clean, entrezid) %>%
    distinct()

  result <- build_term_lists(annot, entrez_map)
  saveRDS(result, CACHE_FILE)
  cat(sprintf("  Cache salvo em: %s\n", CACHE_FILE))
  return(result)
}

# Constrói as listas TERM2GENE por ontologia e TERM2NAME
build_term_lists <- function(annot, entrez_map) {
  # Adiciona namespace se ausente (via GO.db)
  if (!"namespace" %in% names(annot) || all(is.na(annot$namespace))) {
    if (requireNamespace("GO.db", quietly = TRUE)) {
      ns_map <- suppressMessages(
        AnnotationDbi::select(GO.db::GO.db,
                              keys = unique(annot$go_id),
                              columns = c("GOID", "ONTOLOGY"),
                              keytype = "GOID")
      )
      annot <- merge(annot,
                     ns_map[, c("GOID", "ONTOLOGY")],
                     by.x = "go_id", by.y = "GOID", all.x = TRUE)
      names(annot)[names(annot) == "ONTOLOGY"] <- "namespace"
    }
  }

  # Normaliza namespace para BP/MF/CC
  if ("namespace" %in% names(annot)) {
    annot$namespace <- dplyr::recode(annot$namespace,
      "biological_process" = "BP", "molecular_function" = "MF",
      "cellular_component" = "CC")
  } else {
    annot$namespace <- NA_character_
  }

  # TERM2NAME via GO.db
  TERM2NAME <- data.frame(go_id = character(0), term = character(0))
  if (requireNamespace("GO.db", quietly = TRUE)) {
    go_ids <- unique(annot$go_id)
    go_info <- suppressMessages(
      AnnotationDbi::select(GO.db::GO.db,
                            keys = go_ids, columns = c("GOID", "TERM"),
                            keytype = "GOID")
    )
    TERM2NAME <- go_info[, c("GOID", "TERM")]
    names(TERM2NAME) <- c("go_id", "term")
  }

  make_t2g <- function(ns) {
    sub <- annot[!is.na(annot$namespace) & annot$namespace == ns,
                  c("go_id", "gene_id_clean")]
    unique(sub)
  }

  list(
    TERM2GENE_BP = make_t2g("BP"),
    TERM2GENE_MF = make_t2g("MF"),
    TERM2GENE_CC = make_t2g("CC"),
    TERM2GENE_ALL = unique(annot[, c("go_id", "gene_id_clean")]),
    TERM2NAME     = TERM2NAME,
    entrez_map    = entrez_map
  )
}

go_annots <- tryCatch(
  load_go_annotations(opt$go_annot),
  error = function(e) { message("Anotações GO: ", e$message); NULL }
)

# ── 4. ORA – GO ───────────────────────────────────────────────
run_go_ora <- function(de, bg, TERM2GENE, TERM2NAME, ont_label) {
  tryCatch({
    res <- enricher(
      gene          = de,
      universe      = bg,
      TERM2GENE     = TERM2GENE,
      TERM2NAME     = TERM2NAME,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 5,
      maxGSSize     = 500
    )
    n <- if (!is.null(res)) sum(res@result$p.adjust < 0.05, na.rm = TRUE) else 0
    cat(sprintf("  GO-%s: %d termos enriquecidos\n", ont_label, n))
    res
  }, error = function(e) {
    message(sprintf("  GO-%s falhou: %s", ont_label, e$message))
    NULL
  })
}

export_enrich <- function(obj, fname) {
  df <- if (is.null(obj)) data.frame() else as.data.frame(obj)
  write_tsv(df, file.path(opt$outdir, fname))
}

go_bp <- go_mf <- go_cc <- NULL
cat("\nAnálise GO – ORA:\n")
if (!is.null(go_annots)) {
  go_bp <- run_go_ora(de_ids, bg_ids,
                       go_annots$TERM2GENE_BP, go_annots$TERM2NAME, "BP")
  go_mf <- run_go_ora(de_ids, bg_ids,
                       go_annots$TERM2GENE_MF, go_annots$TERM2NAME, "MF")
  go_cc <- run_go_ora(de_ids, bg_ids,
                       go_annots$TERM2GENE_CC, go_annots$TERM2NAME, "CC")
} else {
  cat("  Sem anotações GO disponíveis – arquivos exportados vazios.\n")
}
export_enrich(go_bp, "go_bp_results.tsv")
export_enrich(go_mf, "go_mf_results.tsv")
export_enrich(go_cc, "go_cc_results.tsv")

# ── 5. ORA – KEGG ────────────────────────────────────────────
cat("\nAnálise KEGG – ORA:\n")
kegg_res <- NULL

# Entrez IDs: primeiro tenta biomaRt, fallback AnnotationHub
get_entrez <- function(gene_ids_clean, entrez_map_df) {
  if (!is.null(entrez_map_df) && nrow(entrez_map_df) > 0) {
    m <- entrez_map_df[entrez_map_df$gene_id_clean %in% gene_ids_clean, ]
    return(na.omit(as.character(unique(m$entrezid))))
  }
  # Fallback AnnotationHub OrgDb
  tryCatch({
    library(AnnotationHub)
    ah <- AnnotationHub(ask = FALSE)
    hits <- query(ah, c("OrgDb", "Glycine"))
    db <- hits[[length(hits)]]
    mapped <- mapIds(db, keys = gene_ids_clean,
                     column = "ENTREZID", keytype = "GID", multiVals = "first")
    na.omit(as.character(mapped))
  }, error = function(e) character(0))
}

entrez_map_df <- if (!is.null(go_annots)) go_annots$entrez_map else NULL
de_entrez  <- get_entrez(de_ids,  entrez_map_df)
bg_entrez  <- get_entrez(bg_ids,  entrez_map_df)

cat(sprintf("  Entrez IDs – background: %d | DE: %d\n",
            length(bg_entrez), length(de_entrez)))

if (length(de_entrez) > 0) {
  kegg_res <- tryCatch({
    res <- enrichKEGG(
      gene          = de_entrez,
      universe      = bg_entrez,
      organism      = opt$organism,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 5,
      maxGSSize     = 500
    )
    cat(sprintf("  KEGG: %d pathways enriquecidos\n",
                sum(res@result$p.adjust < 0.05, na.rm = TRUE)))
    res
  }, error = function(e) {
    message("  KEGG falhou: ", e$message)
    NULL
  })
} else {
  cat("  Sem Entrez IDs disponíveis – KEGG ignorado.\n")
}
export_enrich(kegg_res, "kegg_results.tsv")

# ── 6. GSEA ──────────────────────────────────────────────────
cat("\nGSEA:\n")

ranked_df <- res_all %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(rank_score = log2FoldChange * -log10(padj + 1e-300)) %>%
  arrange(desc(rank_score))

# GSEA-GO: usa Glyma IDs diretamente com TERM2GENE_ALL
gsea_go <- NULL
if (!is.null(go_annots) && nrow(go_annots$TERM2GENE_BP) > 0) {
  ranked_go <- setNames(ranked_df$rank_score, ranked_df$gene_id_clean)
  ranked_go <- ranked_go[!duplicated(names(ranked_go))]

  gsea_go <- tryCatch({
    res <- GSEA(
      geneList     = ranked_go,
      TERM2GENE    = go_annots$TERM2GENE_BP,
      TERM2NAME    = go_annots$TERM2NAME,
      minGSSize    = 15,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose      = FALSE
    )
    cat(sprintf("  GSEA-GO: %d termos\n",
                sum(res@result$p.adjust < 0.05, na.rm = TRUE)))
    res
  }, error = function(e) {
    message("  GSEA-GO falhou: ", e$message)
    NULL
  })
}

# GSEA-KEGG: usa Entrez IDs (só quando temos mapa completo Glyma→ENTREZ)
gsea_kegg <- NULL
if (!is.null(entrez_map_df) && nrow(entrez_map_df) > 0 && length(de_entrez) > 0) {
  ranked_kegg_df <- ranked_df %>%
    left_join(entrez_map_df, by = "gene_id_clean") %>%
    filter(!is.na(entrezid)) %>%
    arrange(desc(rank_score))

  ranked_kegg <- setNames(ranked_kegg_df$rank_score,
                           as.character(ranked_kegg_df$entrezid))
  ranked_kegg <- ranked_kegg[!duplicated(names(ranked_kegg))]

  gsea_kegg <- tryCatch({
    res <- gseKEGG(
      geneList     = ranked_kegg,
      organism     = opt$organism,
      minGSSize    = 15,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose      = FALSE
    )
    cat(sprintf("  GSEA-KEGG: %d pathways\n",
                sum(res@result$p.adjust < 0.05, na.rm = TRUE)))
    res
  }, error = function(e) {
    message("  GSEA-KEGG falhou: ", e$message)
    NULL
  })
}

export_enrich(gsea_go,   "gsea_go_results.tsv")
export_enrich(gsea_kegg, "gsea_kegg_results.tsv")

# ── 7. Figuras ────────────────────────────────────────────────
theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

save_plot <- function(p, name, w = 10, h = 8) {
  ggsave(file.path(opt$figures_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(opt$figures_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300)
}

if (!is.null(go_bp) && nrow(go_bp) > 0) {
  p_go <- dotplot(go_bp, showCategory = 20, title = "GO – Biological Process") +
    theme_pub
  save_plot(p_go, "go_bp_dotplot", w = 10, h = 9)

  tryCatch({
    go_bp2 <- pairwise_termsim(go_bp)
    p_emap <- emapplot(go_bp2, showCategory = 25) +
      ggtitle("GO BP – Enrichment Map")
    save_plot(p_emap, "go_bp_emap", w = 12, h = 10)
  }, error = function(e) message("emapplot falhou: ", e$message))
}

if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
  p_kegg <- dotplot(kegg_res, showCategory = 20, title = "KEGG Pathways") + theme_pub
  save_plot(p_kegg, "kegg_dotplot", w = 10, h = 9)
  p_kegg_bar <- barplot(kegg_res, showCategory = 20, title = "KEGG Pathways") + theme_pub
  save_plot(p_kegg_bar, "kegg_barplot", w = 10, h = 9)
}

if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
  p_gsea <- dotplot(gsea_go, showCategory = 20, split = ".sign",
                    title = "GSEA – GO Biological Process") +
    facet_grid(. ~ .sign) + theme_pub
  save_plot(p_gsea, "gsea_go_dotplot", w = 12, h = 9)
}

if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  p_gsea_kegg <- dotplot(gsea_kegg, showCategory = 20, split = ".sign",
                          title = "GSEA – KEGG Pathways") +
    facet_grid(. ~ .sign) + theme_pub
  save_plot(p_gsea_kegg, "gsea_kegg_dotplot", w = 12, h = 9)
}

cat("\nEnriquecimento funcional concluído. Resultados em:", opt$outdir, "\n")
