#!/usr/bin/env nextflow
// ============================================================
// RNASeq Insight Platform – Glycine max Wm82.a4.v1
// Autor: Eulalio Santos | UFV
// ============================================================

nextflow.enable.dsl = 2

// ── Importa módulos ───────────────────────────────────────────
include { FASTQC as FASTQC_RAW     } from './modules/qc'
include { FASTQC as FASTQC_TRIMMED } from './modules/qc'
include { MULTIQC as MULTIQC_RAW   } from './modules/qc'
include { MULTIQC as MULTIQC_TRIM  } from './modules/qc'

include { FASTP                    } from './modules/trimming'

include { GFFREAD                  } from './modules/alignment'
include { HISAT2_BUILD             } from './modules/alignment'
include { HISAT2_ALIGN             } from './modules/alignment'
include { SAMTOOLS_SORT_INDEX      } from './modules/alignment'

include { FEATURECOUNTS            } from './modules/quantification'
include { PARSE_COUNTS             } from './modules/quantification'

include { SALMON_INDEX             } from './modules/salmon'
include { SALMON_QUANT             } from './modules/salmon'
include { TXIMPORT                 } from './modules/salmon'

include { COMBAT_SEQ               } from './modules/batch'

include { RMATS                    } from './modules/splicing'
include { PARSE_RMATS              } from './modules/splicing'

include { STRINGTIE                } from './modules/assembly'
include { GFFCOMPARE               } from './modules/assembly'

include { LNCRNA_PRED              } from './modules/lncrna'

include { DESEQ2                   } from './modules/analysis'
include { ENRICHMENT               } from './modules/analysis'
include { GOSEQ                    } from './modules/analysis'
include { WGCNA                    } from './modules/analysis'
include { INTEGRATION              } from './modules/analysis'
include { MACHINE_LEARNING         } from './modules/analysis'
include { PPI_NETWORK              } from './modules/analysis'
include { PLANTFDB                 } from './modules/analysis'
include { GENIE3                   } from './modules/analysis'
include { METANALYSIS              } from './modules/analysis'
include { QUARTO_REPORT            } from './modules/analysis'

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
            if (!reads[0].exists()) error "FASTQ não encontrado: ${reads[0]}"
            if (!reads[1].exists()) error "FASTQ não encontrado: ${reads[1]}"
            return [meta, reads]
        }
}

// ── Workflow principal ────────────────────────────────────────

workflow {

    validate_params()

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║   RNASeq Insight Platform – Glycine max Wm82.a4.v1  ║
    ║   Autor: ${params.report_author}
    ╚══════════════════════════════════════════════════════╝
    Samplesheet : ${params.samplesheet}
    Genoma      : ${params.genome_fasta ?: 'índice pré-existente'}
    Anotação    : ${params.genome_gff3 ?: params.genome_gtf}
    KEGG org    : ${params.kegg_organism}
    Resultados  : ${params.outdir}
    StringTie   : ${params.run_stringtie ?: false}
    ML/PPI      : ${params.run_ml ?: false} / ${params.run_ppi ?: false}
    """.stripIndent()

    // ── 1. Samplesheet ───────────────────────────────────────
    reads_ch = parse_samplesheet(params.samplesheet)

    // ── 2. Genoma / GTF ──────────────────────────────────────
    genome_fasta = file(params.genome_fasta)

    if (params.genome_gtf) {
        gtf = file(params.genome_gtf)
    } else {
        GFFREAD(file(params.genome_gff3), genome_fasta)
        gtf = GFFREAD.out.gtf
    }

    // ── 3. Índice alinhador ───────────────────────────────────
    if (params.genome_index) {
        index_dir = file(params.genome_index)
    } else {
        HISAT2_BUILD(genome_fasta, gtf)
        index_dir = HISAT2_BUILD.out.index_dir
    }

    // ── 4. QC pré-trimagem ────────────────────────────────────
    FASTQC_RAW(reads_ch, 'pre_trim')
    MULTIQC_RAW(FASTQC_RAW.out.zip.map { m, z -> z }.collect(), 'pre_trim')

    // ── 5. Trimagem ───────────────────────────────────────────
    FASTP(reads_ch)
    trimmed_ch = FASTP.out.reads

    // ── 6. QC pós-trimagem ────────────────────────────────────
    FASTQC_TRIMMED(trimmed_ch, 'post_trim')
    MULTIQC_TRIM(
        FASTQC_TRIMMED.out.zip.map { m, z -> z }.collect()
            .mix(FASTP.out.json.collect()),
        'post_trim'
    )

    // ── 7. Alinhamento ────────────────────────────────────────
    HISAT2_ALIGN(trimmed_ch, index_dir)
    SAMTOOLS_SORT_INDEX(HISAT2_ALIGN.out.bam)
    sorted_bams_ch = SAMTOOLS_SORT_INDEX.out.bam

    // ── 8. Salmon (pseudoalinhamento paralelo) ────────────────
    if (params.run_salmon != false) {
        SALMON_INDEX(genome_fasta)
        SALMON_QUANT(trimmed_ch, SALMON_INDEX.out.index)
        quant_dirs = SALMON_QUANT.out.quant.map { s, d -> d }.collect()
        TXIMPORT(quant_dirs, file(params.samplesheet), gtf)
    }

    // ── 9. featureCounts ──────────────────────────────────────
    all_bams = sorted_bams_ch.map { meta, bam -> bam }.collect()
    FEATURECOUNTS(all_bams, gtf, params.strandedness)
    PARSE_COUNTS(FEATURECOUNTS.out.counts, file(params.samplesheet))

    // ── 10. Correção de batch ─────────────────────────────────
    if (params.run_combat_seq) {
        COMBAT_SEQ(PARSE_COUNTS.out.counts, PARSE_COUNTS.out.metadata)
        counts_for_de = COMBAT_SEQ.out.counts_corrected
    } else {
        counts_for_de = PARSE_COUNTS.out.counts
    }

    // ── 11. Splicing alternativo ──────────────────────────────
    ctrl_bams = sorted_bams_ch
        .filter { meta, bam -> meta.condition == params.control_group }
        .map    { meta, bam -> bam }.collect()
    treat_bams = sorted_bams_ch
        .filter { meta, bam -> meta.condition == params.treatment_group }
        .map    { meta, bam -> bam }.collect()

    RMATS(ctrl_bams, treat_bams, gtf)
    PARSE_RMATS(RMATS.out.results_dir)

    // ── 12. Expressão diferencial ─────────────────────────────
    DESEQ2(counts_for_de, PARSE_COUNTS.out.metadata)

    // ── 13. Enriquecimento funcional ──────────────────────────
    ENRICHMENT(DESEQ2.out.results_all, DESEQ2.out.norm_counts)
    GOSEQ(DESEQ2.out.results_all, gtf)

    // ── 14. Co-expressão ──────────────────────────────────────
    WGCNA(DESEQ2.out.norm_counts, PARSE_COUNTS.out.metadata)

    // ── 15. Integração ────────────────────────────────────────
    INTEGRATION(
        DESEQ2.out.results_sig,
        PARSE_RMATS.out.significant,
        ENRICHMENT.out.go_bp,
        ENRICHMENT.out.kegg,
        WGCNA.out.modules,
        WGCNA.out.hub_genes
    )

    // ── 16. StringTie + lncRNA (run_stringtie=true) ───────────
    if (params.run_stringtie) {
        STRINGTIE(sorted_bams_ch, gtf)
        gtf_list = STRINGTIE.out.gtf.map { s, g -> g }.collect()
        GFFCOMPARE(gtf_list, gtf)

        annotated_gtf = GFFCOMPARE.out.annotated.ifEmpty(file('NO_FILE'))
        LNCRNA_PRED(
            genome_fasta,
            annotated_gtf,
            GFFCOMPARE.out.merged,
            GFFCOMPARE.out.tracking
        )
    }

    // ── 17. Fase 3 – Alto Impacto ─────────────────────────────
    if (params.run_ml) {
        MACHINE_LEARNING(
            DESEQ2.out.norm_counts,
            DESEQ2.out.results_sig,
            PARSE_COUNTS.out.metadata
        )
    }

    if (params.run_ppi) {
        PPI_NETWORK(DESEQ2.out.results_sig)
    }

    PLANTFDB(DESEQ2.out.results_sig, DESEQ2.out.results_all)
    GENIE3(DESEQ2.out.norm_counts, PLANTFDB.out.tf_classified)

    if (params.geo_accessions) {
        METANALYSIS(DESEQ2.out.results_sig)
    }

    // ── 18. Relatório automático ──────────────────────────────
    QUARTO_REPORT(
        DESEQ2.out.results_all,
        DESEQ2.out.results_sig,
        DESEQ2.out.norm_counts,
        ENRICHMENT.out.go_bp,
        ENRICHMENT.out.go_mf,
        ENRICHMENT.out.go_cc,
        ENRICHMENT.out.kegg,
        ENRICHMENT.out.gsea_go,
        ENRICHMENT.out.gsea_kegg,
        PARSE_RMATS.out.significant,
        WGCNA.out.modules,
        WGCNA.out.hub_genes,
        INTEGRATION.out.candidates
    )
}
