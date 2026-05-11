#!/usr/bin/env nextflow
// ============================================================
// RNASeq Insight Platform – Glycine max
// Autor: Eulalio Santos
// ============================================================

nextflow.enable.dsl = 2

// Importa módulos
include { FASTQC as FASTQC_RAW      } from './modules/qc'
include { FASTQC as FASTQC_TRIMMED  } from './modules/qc'
include { MULTIQC as MULTIQC_RAW    } from './modules/qc'
include { MULTIQC as MULTIQC_TRIM   } from './modules/qc'

include { FASTP                     } from './modules/trimming'

include { GFFREAD                   } from './modules/alignment'
include { HISAT2_BUILD              } from './modules/alignment'
include { HISAT2_ALIGN              } from './modules/alignment'
include { SAMTOOLS_SORT_INDEX       } from './modules/alignment'

include { FEATURECOUNTS             } from './modules/quantification'
include { PARSE_COUNTS              } from './modules/quantification'

include { RMATS                     } from './modules/splicing'
include { PARSE_RMATS               } from './modules/splicing'

include { DESEQ2                    } from './modules/analysis'
include { ENRICHMENT                } from './modules/analysis'
include { WGCNA                     } from './modules/analysis'
include { INTEGRATION               } from './modules/analysis'
include { QUARTO_REPORT             } from './modules/analysis'

// ── Funções auxiliares ────────────────────────────────────────

def validate_params() {
    if (!params.samplesheet) {
        error "ERROR: --samplesheet é obrigatório"
    }
    if (!params.genome_fasta && !params.genome_index) {
        error "ERROR: --genome_fasta ou --genome_index é obrigatório"
    }
    if (!params.genome_gff3 && !params.genome_gtf) {
        error "ERROR: --genome_gff3 ou --genome_gtf é obrigatório"
    }
}

def parse_samplesheet(csv_path) {
    Channel.fromPath(csv_path)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample:    row.sample,
                condition: row.condition,
                replicate: row.replicate as Integer
            ]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            // Valida existência dos arquivos
            if (!reads[0].exists()) error "FASTQ não encontrado: ${reads[0]}"
            if (!reads[1].exists()) error "FASTQ não encontrado: ${reads[1]}"
            return [meta, reads]
        }
}

// ── Workflow principal ────────────────────────────────────────

workflow {

    validate_params()

    // Banner
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║   RNASeq Insight Platform – Glycine max              ║
    ║   Autor: ${params.report_author}
    ╚══════════════════════════════════════════════════════╝
    Samplesheet : ${params.samplesheet}
    Genoma      : ${params.genome_fasta ?: 'índice pré-existente'}
    Anotação    : ${params.genome_gff3 ?: params.genome_gtf}
    Resultados  : ${params.outdir}
    """.stripIndent()

    // ── 1. Leitura do samplesheet ─────────────────────────────
    reads_ch = parse_samplesheet(params.samplesheet)

    // ── 2. Preparação do genoma ───────────────────────────────
    genome_fasta = file(params.genome_fasta)

    if (params.genome_gtf) {
        gtf = file(params.genome_gtf)
    } else {
        gff3 = file(params.genome_gff3)
        GFFREAD(gff3, genome_fasta)
        gtf = GFFREAD.out.gtf
    }

    if (params.genome_index) {
        index_dir = file(params.genome_index)
    } else {
        HISAT2_BUILD(genome_fasta, gtf)
        index_dir = HISAT2_BUILD.out.index_dir
    }

    // ── 3. QC pré-trimagem ────────────────────────────────────
    FASTQC_RAW(reads_ch, 'pre_trim')
    raw_qc_reports = FASTQC_RAW.out.zip.map { meta, zip -> zip }.collect()
    MULTIQC_RAW(raw_qc_reports, 'pre_trim')

    // ── 4. Trimagem ───────────────────────────────────────────
    FASTP(reads_ch)
    trimmed_ch = FASTP.out.reads

    // ── 5. QC pós-trimagem ────────────────────────────────────
    FASTQC_TRIMMED(trimmed_ch, 'post_trim')
    trim_qc_reports = FASTQC_TRIMMED.out.zip.map { meta, zip -> zip }.collect()
        .mix(FASTP.out.json.collect())
    MULTIQC_TRIM(trim_qc_reports, 'post_trim')

    // ── 6. Alinhamento ────────────────────────────────────────
    HISAT2_ALIGN(trimmed_ch, index_dir)
    SAMTOOLS_SORT_INDEX(HISAT2_ALIGN.out.bam)
    sorted_bams_ch = SAMTOOLS_SORT_INDEX.out.bam

    // ── 7. Quantificação ──────────────────────────────────────
    all_bams = sorted_bams_ch.map { meta, bam -> bam }.collect()
    FEATURECOUNTS(all_bams, gtf, params.strandedness)
    PARSE_COUNTS(FEATURECOUNTS.out.counts, file(params.samplesheet))

    // ── 8. Splicing alternativo (rMATS) ───────────────────────
    ctrl_bams = sorted_bams_ch
        .filter { meta, bam -> meta.condition == params.control_group }
        .map    { meta, bam -> bam }
        .collect()

    treat_bams = sorted_bams_ch
        .filter { meta, bam -> meta.condition == params.treatment_group }
        .map    { meta, bam -> bam }
        .collect()

    RMATS(ctrl_bams, treat_bams, gtf)
    PARSE_RMATS(RMATS.out.results_dir)

    // ── 9. Expressão diferencial ──────────────────────────────
    DESEQ2(
        PARSE_COUNTS.out.counts,
        PARSE_COUNTS.out.metadata
    )

    // ── 10. Enriquecimento ────────────────────────────────────
    ENRICHMENT(
        DESEQ2.out.results_all,
        DESEQ2.out.norm_counts
    )

    // ── 11. Co-expressão (WGCNA) ──────────────────────────────
    WGCNA(
        DESEQ2.out.norm_counts,
        PARSE_COUNTS.out.metadata
    )

    // ── 12. Integração multi-ômica ────────────────────────────
    INTEGRATION(
        DESEQ2.out.results_sig,
        PARSE_RMATS.out.significant,
        ENRICHMENT.out.go_bp,
        ENRICHMENT.out.kegg,
        WGCNA.out.modules,
        WGCNA.out.hub_genes
    )

    // ── 13. Relatório automático ──────────────────────────────
    QUARTO_REPORT(
        DESEQ2.out.figures,
        ENRICHMENT.out.figures,
        PARSE_RMATS.out.significant.map { it.parent },
        WGCNA.out.figures,
        INTEGRATION.out.figures
    )

    // ── Resumo final ──────────────────────────────────────────
    workflow.onComplete {
        log.info """
        ═══════════════════════════════════════
         Pipeline concluído com sucesso!
         Duração  : ${workflow.duration}
         Status   : ${workflow.success ? 'OK' : 'FALHOU'}
         Resultados: ${params.outdir}/
        ═══════════════════════════════════════
        """.stripIndent()
    }
}
