#!/usr/bin/env bash
# ============================================================
# run_pipeline.sh – Execução do pipeline RNASeq
# Uso: bash run_pipeline.sh [local|slurm] [-resume|--resume]
# ============================================================
set -euo pipefail

PROFILE="${1:-local}"
shift || true
# Mapeia --resume → -resume (Nextflow usa traço simples)
EXTRA_ARGS=""
for arg in "$@"; do
    [[ "$arg" == "--resume" ]] && arg="-resume"
    EXTRA_ARGS="$EXTRA_ARGS $arg"
done
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$SCRIPT_DIR/logs/nextflow_${TIMESTAMP}.log"

mkdir -p "$SCRIPT_DIR/logs"

echo "═══════════════════════════════════════════"
echo "  RNASeq Insight – Execução do Pipeline"
echo "  Perfil   : $PROFILE"
echo "  Timestamp: $TIMESTAMP"
echo "═══════════════════════════════════════════"
echo ""

# Valida perfil
if [[ ! "$PROFILE" =~ ^(local|slurm|test)$ ]]; then
    echo "ERRO: perfil inválido '$PROFILE'. Use: local | slurm | test"
    exit 1
fi

# Localiza nextflow: prefere o do env rnaseq-tools (Java correto) → PATH → erro
if mamba run -n rnaseq-tools which nextflow &>/dev/null 2>&1; then
    NF_CMD="mamba run -n rnaseq-tools nextflow"
elif command -v nextflow &>/dev/null; then
    NF_CMD="nextflow"
else
    echo "ERRO: nextflow não encontrado."
    echo "Execute: bash setup.sh primeiro."
    exit 1
fi

echo "  Nextflow : $($NF_CMD -version 2>&1 | head -1)"
echo ""

# Executa pipeline
$NF_CMD run "$SCRIPT_DIR/main.nf" \
    -profile "$PROFILE" \
    -params-file "$SCRIPT_DIR/params.yaml" \
    -with-report "$SCRIPT_DIR/results/nextflow_report.html" \
    -with-trace  "$SCRIPT_DIR/results/nextflow_trace.txt" \
    -with-timeline "$SCRIPT_DIR/results/nextflow_timeline.html" \
    -with-dag    "$SCRIPT_DIR/results/nextflow_dag.html" \
    $EXTRA_ARGS \
    2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "═══════════════════════════════════════════"
    echo "  Pipeline concluído com SUCESSO!"
    echo "  Resultados: $SCRIPT_DIR/results/"
    echo ""
    echo "  Para visualizar o dashboard:"
    echo "  RESULTS_DIR=$SCRIPT_DIR/results \\"
    echo "  Rscript -e \"shiny::runApp('dashboard/app.R', port=3838)\""
    echo "═══════════════════════════════════════════"
else
    echo "═══════════════════════════════════════════"
    echo "  ERRO: pipeline falhou (código $EXIT_CODE)."
    echo "  Verifique: $LOG_FILE"
    echo "  Para retomar: bash run_pipeline.sh $PROFILE -resume"
    echo "═══════════════════════════════════════════"
    exit $EXIT_CODE
fi
