// ============================================================
// Módulo: Controle de Qualidade – FastQC + MultiQC
// ============================================================

process FASTQC {
    tag "${meta.sample}"
    publishDir "${params.outdir}/qc/${stage}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    val(stage)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    def prefix = meta.sample
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path(reports)
    val(stage)

    output:
    path("multiqc_${stage}.html"),   emit: report
    path("multiqc_${stage}_data/"),  emit: data

    script:
    """
    multiqc \\
        --title "${params.report_title} – QC ${stage}" \\
        --filename "multiqc_${stage}.html" \\
        --force \\
        .
    """
}
