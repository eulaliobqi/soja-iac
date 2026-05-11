#!/usr/bin/env Rscript
# ============================================================
# install_r_bioc_packages.R
# Instala pacotes Bioconductor não disponíveis no conda
# Execute UMA VEZ após criar o env r-analysis:
#   mamba run -n r-analysis Rscript scripts/install_r_bioc_packages.R
# ============================================================

cat("Instalando pacotes Bioconductor via BiocManager...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# org.Gmax.eg.db removido do Bioconductor 3.20+; será carregado via AnnotationHub em runtime
bioc_pkgs <- c(
  "WGCNA",
  "GO.db",
  "impute",
  "preprocessCore"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Instalando %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("  ✓ %s já instalado\n", pkg))
  }
}

# Valida
cat("\nValidação:\n")
for (pkg in bioc_pkgs) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %s %s\n", ifelse(ok, "✓", "✗"), pkg))
}

cat("\nConcluído.\n")
