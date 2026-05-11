// ============================================================
// Módulo: Alinhamento – gffread + HISAT2 + samtools
// ============================================================

process GFFREAD {
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path(gff3)
    path(genome_fasta)

    output:
    path("annotation.gtf"), emit: gtf

    script:
    """
    gffread ${gff3} -T -o annotation.gtf
    """
}

process HISAT2_BUILD {
    publishDir "${params.outdir}/genome/hisat2_index", mode: 'copy'

    input:
    path(genome_fasta)
    path(gtf)

    output:
    path("genome_index*.ht2"), emit: index
    path("genome_index"),      emit: index_dir

    script:
    """
    # Extrai splice sites e exons do GTF para melhorar o alinhamento
    hisat2_extract_splice_sites.py ${gtf} > genome_splice_sites.txt
    hisat2_extract_exons.py ${gtf}         > genome_exons.txt

    hisat2-build \\
        -p ${task.cpus} \\
        --ss genome_splice_sites.txt \\
        --exon genome_exons.txt \\
        ${genome_fasta} \\
        genome_index

    mkdir -p genome_index
    mv genome_index.*.ht2 genome_index/ 2>/dev/null || true
    """
}

process HISAT2_ALIGN {
    tag "${meta.sample}"
    publishDir "${params.outdir}/aligned", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }

    input:
    tuple val(meta), path(reads)
    path(index_dir)

    output:
    tuple val(meta), path("${meta.sample}.bam"),     emit: bam
    path("${meta.sample}_hisat2.log"),               emit: log

    script:
    def (r1, r2) = reads
    def index_prefix = "${index_dir}/genome_index"
    """
    hisat2 \\
        -p ${task.cpus} \\
        -x ${index_prefix} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --max-intronlen ${params.hisat2_max_intron} \\
        --dta \\
        --rna-strandness RF \\
        --new-summary \\
        --summary-file ${meta.sample}_hisat2.log \\
        2>> ${meta.sample}_hisat2.log \\
    | samtools view -bS -@ ${task.cpus} - > ${meta.sample}.bam

    # Valida taxa de alinhamento (mínimo 70%)
    RATE=\$(grep "Overall alignment rate" ${meta.sample}_hisat2.log | grep -oP '[\\d.]+(?=%)' | tail -1)
    echo "Alignment rate for ${meta.sample}: \${RATE}%"
    awk "BEGIN { if (\${RATE} < 70) { print \"ERROR: alignment rate below 70% for ${meta.sample}\"; exit 1 } }"
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "${meta.sample}"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}_sorted.bam"),      emit: bam
    tuple val(meta), path("${meta.sample}_sorted.bam.bai"),  emit: bai
    tuple val(meta), path("${meta.sample}_flagstat.txt"),     emit: flagstat

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${meta.sample}_sorted.bam \\
        ${bam}

    samtools index \\
        -@ ${task.cpus} \\
        ${meta.sample}_sorted.bam

    samtools flagstat \\
        -@ ${task.cpus} \\
        ${meta.sample}_sorted.bam > ${meta.sample}_flagstat.txt
    """
}

process INFER_STRANDEDNESS {
    tag "${meta.sample}"
    publishDir "${params.outdir}/aligned/strandedness", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(bed)

    output:
    tuple val(meta), path("${meta.sample}_strandedness.txt"), emit: result

    script:
    """
    infer_experiment.py \\
        -i ${bam} \\
        -r ${bed} \\
        -s 200000 \\
        > ${meta.sample}_strandedness.txt 2>&1
    """
}
