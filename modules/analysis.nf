// ============================================================
// Módulo: Análises R – DESeq2, Enriquecimento, WGCNA, Integração
// ============================================================

process DESEQ2 {
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path(counts)
    path(metadata)

    output:
    path("deseq2_results_all.tsv"),      emit: results_all
    path("deseq2_results_sig.tsv"),      emit: results_sig
    path("normalized_counts.tsv"),       emit: norm_counts
    path("figures/"),                    emit: figures
    path("deseq2_summary.txt"),          emit: summary

    script:
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/01_deseq2.R \\
        --counts         ${counts} \\
        --metadata       ${metadata} \\
        --control        ${params.control_group} \\
        --treatment      ${params.treatment_group} \\
        --padj           ${params.padj_cutoff} \\
        --lfc            ${params.lfc_cutoff} \\
        --outdir         . \\
        --figures_dir    figures
    """
}

process ENRICHMENT {
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path(deseq2_results)
    path(norm_counts)

    output:
    path("go_bp_results.tsv"),   emit: go_bp
    path("go_mf_results.tsv"),   emit: go_mf
    path("go_cc_results.tsv"),   emit: go_cc
    path("kegg_results.tsv"),    emit: kegg
    path("gsea_go_results.tsv"), emit: gsea_go
    path("gsea_kegg_results.tsv"), emit: gsea_kegg
    path("figures/"),            emit: figures

    script:
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/02_enrichment.R \\
        --deseq2      ${deseq2_results} \\
        --norm_counts ${norm_counts} \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --organism    ${params.kegg_organism} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process WGCNA {
    publishDir "${params.outdir}/wgcna", mode: 'copy'

    input:
    path(norm_counts)
    path(metadata)

    output:
    path("wgcna_modules.tsv"),     emit: modules
    path("wgcna_hub_genes.tsv"),   emit: hub_genes
    path("wgcna_eigengenes.tsv"),  emit: eigengenes
    path("figures/"),              emit: figures

    script:
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/03_wgcna.R \\
        --norm_counts ${norm_counts} \\
        --metadata    ${metadata} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process INTEGRATION {
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path(deseq2_sig)
    path(splicing_sig)
    path(go_bp)
    path(kegg)
    path(wgcna_modules)
    path(hub_genes)

    output:
    path("integrated_genes.tsv"),  emit: integrated
    path("gene_ranking.tsv"),      emit: ranking
    path("key_candidates.tsv"),    emit: candidates
    path("figures/"),              emit: figures

    script:
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/04_integration.R \\
        --deseq2    ${deseq2_sig} \\
        --splicing  ${splicing_sig} \\
        --go        ${go_bp} \\
        --kegg      ${kegg} \\
        --wgcna     ${wgcna_modules} \\
        --hub_genes ${hub_genes} \\
        --outdir    . \\
        --figures_dir figures
    """
}

process QUARTO_REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path(deseq2_dir,      stageAs: 'deseq2_figures')
    path(enrichment_dir,  stageAs: 'enrichment_figures')
    path(splicing_dir,    stageAs: 'splicing_dir')
    path(wgcna_dir,       stageAs: 'wgcna_figures')
    path(integration_dir, stageAs: 'integration_figures')

    output:
    path("rnaseq_report.html"), emit: report

    script:
    """
    cp -r ${projectDir}/report/rnaseq_report.qmd .

    quarto render rnaseq_report.qmd \\
        -P deseq2_dir:deseq2_figures \\
        -P enrichment_dir:enrichment_figures \\
        -P splicing_dir:splicing_dir \\
        -P wgcna_dir:wgcna_figures \\
        -P integration_dir:integration_figures \\
        -P report_title:"${params.report_title}" \\
        -P report_author:"${params.report_author}" \\
        --output rnaseq_report.html
    """
}
