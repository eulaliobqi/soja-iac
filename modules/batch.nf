// ============================================================
// Módulo: ComBat-Seq – Correção de batch
// ============================================================

process COMBAT_SEQ {
    label 'medium_mem'
    publishDir "${params.outdir}/batch_correction", mode: 'copy'

    input:
    path(counts)
    path(metadata)

    output:
    path("counts_corrected.tsv"), emit: counts_corrected
    path("pca_before_batch.pdf"), emit: pca_before
    path("pca_after_batch.pdf"),  emit: pca_after
    path("batch_report.txt"),     emit: report

    script:
    """
    Rscript ${projectDir}/scripts/05_batch_correction.R \
        --counts   ${counts} \
        --metadata ${metadata} \
        --outdir   .
    """
}
