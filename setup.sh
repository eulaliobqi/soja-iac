#!/usr/bin/env bash
# ============================================================
# setup.sh – Instalação dos ambientes mamba e Nextflow
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "═══════════════════════════════════════════"
echo "  RNASeq Insight – Setup de Ambiente"
echo "═══════════════════════════════════════════"

# ── Verifica dependências ─────────────────────────────────────
for cmd in mamba conda; do
    if command -v "$cmd" &>/dev/null; then
        CONDA_CMD="$cmd"
        break
    fi
done

if [ -z "${CONDA_CMD:-}" ]; then
    echo "ERRO: mamba ou conda não encontrado."
    echo "Instale o Miniforge: https://github.com/conda-forge/miniforge"
    exit 1
fi

echo "Usando: $CONDA_CMD"

# ── Instala Nextflow ──────────────────────────────────────────
if ! command -v nextflow &>/dev/null; then
    echo ""
    echo "Instalando Nextflow..."
    $CONDA_CMD install -y -c bioconda -c conda-forge nextflow
fi
echo "Nextflow: $(nextflow -version 2>&1 | head -1)"

# ── Instala ambientes ─────────────────────────────────────────
install_env() {
    local yml="$1"
    local env_name
    env_name=$(grep "^name:" "$yml" | awk '{print $2}')

    if $CONDA_CMD env list | grep -q "^${env_name} "; then
        echo "  Ambiente '$env_name' já existe – pulando."
    else
        echo "  Instalando '$env_name'..."
        $CONDA_CMD env create -f "$yml"
    fi
}

echo ""
echo "Instalando ambientes Mamba..."
install_env "$SCRIPT_DIR/envs/rnaseq-tools.yml"
install_env "$SCRIPT_DIR/envs/r-analysis.yml"

# ── Cria estrutura de diretórios de resultado ─────────────────
echo ""
echo "Criando estrutura de diretórios..."
for d in results/qc results/trimmed results/aligned results/counts \
          results/deseq2/figures results/enrichment/figures \
          results/splicing results/wgcna/figures \
          results/integration/figures results/report; do
    mkdir -p "$SCRIPT_DIR/$d"
done

# ── Valida instalação ─────────────────────────────────────────
echo ""
echo "Validando instalação..."

validate_tool() {
    local env="$1"
    local tool="$2"
    if $CONDA_CMD run -n "$env" which "$tool" &>/dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool (não encontrado em $env)"
    fi
}

validate_tool rnaseq-tools fastqc
validate_tool rnaseq-tools fastp
validate_tool rnaseq-tools hisat2
validate_tool rnaseq-tools samtools
validate_tool rnaseq-tools featureCounts
validate_tool rnaseq-tools multiqc
validate_tool r-analysis   Rscript

echo ""
echo "═══════════════════════════════════════════"
echo "  Setup concluído!"
echo ""
echo "  Próximos passos:"
echo "  1. Edite params.yaml com seus dados"
echo "  2. Edite samplesheet.csv com suas amostras"
echo "  3. Execute: bash run_pipeline.sh"
echo "═══════════════════════════════════════════"
