// ============================================================
// Módulo: StringTie + GFFcompare – Montagem de transcritos
// ============================================================

process STRINGTIE {
    label 'medium_mem'
    tag { sample }
    publishDir "${params.outdir}/stringtie", mode: 'copy'

    input:
    tuple val(sample), path(bam)
    path(gtf)

    output:
    tuple val(sample), path("${sample}.gtf"),              emit: gtf
    path("${sample}_abundance.tsv"),                       emit: abundance

    script:
    def strand_flag = params.strandedness == 1 ? "--fr" :
                      params.strandedness == 2 ? "--rf" : ""
    """
    stringtie ${bam} \
        -G ${gtf} \
        ${strand_flag} \
        -o ${sample}.gtf \
        -A ${sample}_abundance.tsv \
        -p ${task.cpus}
    """
}

process GFFCOMPARE {
    label 'low_mem'
    publishDir "${params.outdir}/stringtie", mode: 'copy'

    input:
    path(gtf_list)
    path(ref_gtf)

    output:
    path("gffcmp.annotated.gtf"), emit: annotated, optional: true
    path("gffcmp.tracking"),      emit: tracking
    path("gffcmp.stats"),         emit: stats
    path("gffcmp.loci"),          emit: loci
    path("merged.gtf"),           emit: merged, optional: true

    script:
    """
    # Cria lista de GTFs para merge
    ls *.gtf | grep -v ref > gtf_list.txt

    # Merge com StringTie
    stringtie --merge \
        -G ${ref_gtf} \
        -o merged.gtf \
        gtf_list.txt

    # GFFcompare: gera gffcmp.annotated.gtf (versão >=0.12)
    gffcompare \
        -r ${ref_gtf} \
        -o gffcmp \
        merged.gtf
    """
}
