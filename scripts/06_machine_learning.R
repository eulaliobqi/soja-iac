#!/usr/bin/env Rscript
# ============================================================
# 06_machine_learning.R – Biomarcadores: RF + SVM + ElasticNet
# Nota: com n<10 amostras, AUC=1 é trivial (LOOCV)
# Usar feature_importance.tsv para selecionar biomarcadores
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(caret)
  library(randomForest)
  library(glmnet)
  library(pROC)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

opt_list <- list(
  make_option("--norm_counts", type = "character", help = "Contagens normalizadas (TSV)"),
  make_option("--deseq2_sig", type = "character", help = "DEGs significativos (TSV)"),
  make_option("--metadata",   type = "character", help = "Metadados (TSV)"),
  make_option("--top_genes",  type = "integer",   default = 200),
  make_option("--outdir",     type = "character", default = "."),
  make_option("--figures_dir",type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "ml_results.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "feature_importance.tsv"))
writeLines("", file.path(opt$figures_dir, ".ml_placeholder"))

cat("═══════════════════════════════════════\n")
cat("  Machine Learning – RF + SVM + ElasticNet\n")
cat("═══════════════════════════════════════\n")

norm <- read_tsv(opt$norm_counts, show_col_types = FALSE) |>
  tibble::column_to_rownames("gene_id")
degs <- read_tsv(opt$deseq2_sig, show_col_types = FALSE)
meta <- read_tsv(opt$metadata,   show_col_types = FALSE) |>
  as.data.frame() |> tibble::column_to_rownames("sample")

# Strip versão
degs$gene_id <- gsub("\\.[0-9]+$", "", degs$gene_id)
rownames(norm) <- gsub("\\.[0-9]+$", "", rownames(norm))

# Seleciona top DEGs por |lfc|
top_genes <- degs |>
  arrange(desc(abs(log2FoldChange))) |>
  slice_head(n = opt$top_genes) |>
  pull(gene_id)
top_genes <- intersect(top_genes, rownames(norm))

if (length(top_genes) < 5) {
  message("Poucos genes DEG para ML (", length(top_genes), "); abortando.")
  quit(status = 0)
}

# Matriz: amostras × genes
X <- t(as.matrix(norm[top_genes, rownames(meta)]))
y <- factor(meta$condition)
levels(y) <- make.names(levels(y))

n_samples <- nrow(X)
cat(sprintf("Amostras: %d | Genes: %d | Classes: %s\n",
            n_samples, ncol(X), paste(levels(y), collapse = " vs ")))

# Cross-validation adaptada ao n
if (n_samples < 10) {
  ctrl <- trainControl(method = "LOOCV", classProbs = TRUE,
                       summaryFunction = twoClassSummary, savePredictions = TRUE)
  cat("Usando LOOCV (n < 10)\n")
} else {
  k <- min(5L, floor(n_samples / 2L))
  ctrl <- trainControl(method = "cv", number = k, classProbs = TRUE,
                       summaryFunction = twoClassSummary, savePredictions = TRUE)
  cat(sprintf("Usando %d-fold CV\n", k))
}

train_model <- function(method, tuneGrid = NULL) {
  tryCatch(
    train(X, y, method = method, trControl = ctrl,
          metric = "ROC", tuneGrid = tuneGrid,
          preProcess = c("center", "scale")),
    error = function(e) { message(method, " falhou: ", e$message); NULL }
  )
}

rf_mod  <- train_model("rf",  data.frame(mtry = floor(sqrt(ncol(X)))))
svm_mod <- train_model("svmRadial")
en_mod  <- train_model("glmnet")

# Importâncias
get_importance <- function(mod, label) {
  if (is.null(mod)) return(data.frame())
  imp <- tryCatch(varImp(mod)$importance, error = function(e) NULL)
  if (is.null(imp)) return(data.frame())
  imp |>
    tibble::rownames_to_column("gene") |>
    mutate(model = label) |>
    arrange(desc(Overall))
}

imp_all <- bind_rows(
  get_importance(rf_mod,  "RandomForest"),
  get_importance(svm_mod, "SVM"),
  get_importance(en_mod,  "ElasticNet")
)
write_tsv(imp_all, file.path(opt$outdir, "feature_importance.tsv"))

# Top 20 genes por RF
top20 <- imp_all |> filter(model == "RandomForest") |> slice_head(n = 20)

p_imp <- ggplot(top20, aes(x = reorder(gene, Overall), y = Overall)) +
  geom_col(fill = "#2c7bb6") + coord_flip() +
  labs(title = "Top 20 Biomarcadores – Random Forest",
       x = NULL, y = "Importância") + theme_bw()

ggsave(file.path(opt$figures_dir, "feature_importance.pdf"), p_imp, width = 8, height = 6)
ggsave(file.path(opt$figures_dir, "feature_importance.png"), p_imp, width = 8, height = 6,
       dpi = 300)

# Resultados AUC (nota: AUC=1 trivial com n pequeno)
results <- bind_rows(lapply(list(rf_mod, svm_mod, en_mod), function(m) {
  if (is.null(m)) return(data.frame())
  data.frame(model = m$method,
             AUC   = max(m$results$ROC, na.rm = TRUE))
}))
write_tsv(results, file.path(opt$outdir, "ml_results.tsv"))

if (n_samples < 10) {
  cat("AVISO: n < 10 amostras → AUC pode ser 1.0 (trivial com LOOCV).\n")
  cat("Use feature_importance.tsv para seleção de biomarcadores.\n")
}
cat("Concluído.\n")
