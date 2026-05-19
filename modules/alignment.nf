// ============================================================
// Módulo: Alinhamento – gffread + HISAT2 + samtools
// (STAR mantido no arquivo mas não usado: CPU sem suporte AVX2)
// ============================================================

process GFFREAD {
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path(gff3)
    path(genome_fasta)

    output:
    path("annotation.gtf"), emit: gtf

    script:
    // Lição 3: corrige chrChr gerado pelo gffread em Glycine max (Chr01→chrChr01)
    """
    gffread ${gff3} -T -o - | sed 's/^chrChr/Chr/; s/^chrchr/Chr/' > annotation.gtf
    """
}

process STAR_INDEX {
    label 'high_mem'
    publishDir "${params.outdir}/genome/star_index", mode: 'copy'

    input:
    path(genome_fasta)
    path(gtf)

    output:
    path("star_index/"), emit: index

    script:
    // STAR 2.7.10b obrigatório — 2.7.11b causa SIGSEGV em CPUs sem AVX2 (Lição 1)
    """
    mkdir -p star_index
    STAR \
        --runMode genomeGenerate \
        --runThreadN ${task.cpus} \
        --genomeDir star_index \
        --genomeFastaFiles ${genome_fasta} \
        --sjdbGTFfile ${gtf} \
        --genomeSAindexNbases 13
    """
}

process STAR_ALIGN {
    label 'high_mem'
    tag "${meta.sample}"
    publishDir "${params.outdir}/aligned", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') || fn.endsWith('.tab') ? "logs/${fn}" : fn }

    input:
    tuple val(meta), path(reads)
    path(index_dir)

    output:
    tuple val(meta), path("${meta.sample}.bam"), emit: bam
    path("${meta.sample}_Log.final.out"),         emit: log

    script:
    def (r1, r2) = reads
    def strand = params.strandedness == 1 ? "Forward" :
                 params.strandedness == 2 ? "Reverse" : "Unstranded"
    """
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${index_dir} \
        --readFilesIn ${r1} ${r2} \
        --readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NH HI AS NM MD \
        --outFilterMultimapNmax 10 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFileNamePrefix ${meta.sample}_ \
        --outStd BAM_Unsorted \
        --soloStrand ${strand} \
    > ${meta.sample}.bam

    # Valida taxa de alinhamento
    RATE=\$(grep "Uniquely mapped reads %" ${meta.sample}_Log.final.out | awk '{print \$NF}' | tr -d '%')
    echo "STAR alignment rate for ${meta.sample}: \${RATE}%"
    awk -v rate="\${RATE}" -v sample="${meta.sample}" -v min="${params.min_align_rate}" \
        'BEGIN { if (rate+0 < min+0) { print "WARN: alignment rate " rate "% below " min "% for " sample } }'
    """
}

process HISAT2_BUILD {
    publishDir "${params.outdir}/genome/hisat2_index", mode: 'copy'

    input:
    path(genome_fasta)
    path(gtf)

    output:
    path("genome_index"), emit: index_dir

    script:
    """
    hisat2_extract_splice_sites.py ${gtf} > splice_sites.txt
    hisat2_extract_exons.py ${gtf}         > exons.txt

    mkdir -p genome_index
    hisat2-build \\
        -p ${task.cpus} \\
        --ss splice_sites.txt \\
        --exon exons.txt \\
        ${genome_fasta} \\
        genome_index/genome_index
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

    # Valida taxa de alinhamento
    RATE=\$(grep "Overall alignment rate" ${meta.sample}_hisat2.log | grep -oP '[\\d.]+(?=%)' | tail -1)
    echo "Alignment rate for ${meta.sample}: \${RATE}%"
    awk -v rate="\${RATE}" -v sample="${meta.sample}" -v min="${params.min_align_rate}" \
        'BEGIN { if (rate+0 < min+0) { print "ERROR: alignment rate " rate "% below minimum " min "% for " sample; exit 1 } }'
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
