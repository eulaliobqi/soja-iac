// ============================================================
// Módulo: Salmon – Pseudoalinhamento e quantificação
// ============================================================

process SALMON_INDEX {
    label 'medium_mem'
    publishDir "${params.outdir}/salmon_index", mode: 'copy'

    input:
    path(transcriptome_fa)

    output:
    path("salmon_index/"), emit: index

    script:
    """
    salmon index \
        --threads ${task.cpus} \
        --transcripts ${transcriptome_fa} \
        --index salmon_index \
        --gencode
    """
}

process SALMON_QUANT {
    label 'medium_mem'
    tag { meta.sample }
    publishDir "${params.outdir}/counts/salmon", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("${meta.sample}/"), emit: quant
    path("${meta.sample}/quant.sf"),           emit: quant_sf
    path("${meta.sample}/logs/"),              emit: logs

    script:
    def sample = meta.sample
    def r1 = reads[0]
    def r2 = reads.size() > 1 ? "-2 ${reads[1]}" : ""
    def lib = params.strandedness == 1 ? "SF" :
              params.strandedness == 2 ? "SR" : "A"
    """
    mkdir -p ${sample}/logs
    salmon quant \
        --index ${index} \
        --libType ${lib} \
        -1 ${r1} ${r2} \
        --threads ${task.cpus} \
        --validateMappings \
        --gcBias \
        --output ${sample} \
        2> ${sample}/logs/salmon_quant.log
    """
}

process TXIMPORT {
    label 'low_mem'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(quant_dirs)
    path(samplesheet)
    path(gtf)

    output:
    path("salmon_counts.tsv"), emit: counts
    path("salmon_tpm.tsv"),    emit: tpm

    script:
    """
    Rscript ${projectDir}/scripts/00_tximport.R \
        --salmon_dir . \
        --samplesheet ${samplesheet} \
        --gtf ${gtf} \
        --outdir .
    """
}
