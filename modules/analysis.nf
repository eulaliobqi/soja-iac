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
    def go_annot_arg = params.go_annot ? "--go_annot ${params.go_annot}" : ""
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/02_enrichment.R \\
        --deseq2      ${deseq2_results} \\
        --norm_counts ${norm_counts} \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --organism    ${params.kegg_organism} \\
        --outdir      . \\
        --figures_dir figures \\
        ${go_annot_arg}
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
    path deseq2_all
    path deseq2_sig
    path norm_counts
    path go_bp
    path go_mf
    path go_cc
    path kegg
    path gsea_go
    path gsea_kegg
    path splicing_sig
    path wgcna_modules
    path wgcna_hub_genes
    path integration_candidates

    output:
    path("rnaseq_report.html"), emit: report

    script:
    """
    mkdir -p deseq2_data enrichment_data splicing_data wgcna_data integration_data

    cp ${deseq2_all}             deseq2_data/deseq2_results_all.tsv
    cp ${deseq2_sig}             deseq2_data/deseq2_results_sig.tsv
    cp ${norm_counts}            deseq2_data/normalized_counts.tsv
    cp ${go_bp}                  enrichment_data/go_bp_results.tsv
    cp ${go_mf}                  enrichment_data/go_mf_results.tsv
    cp ${go_cc}                  enrichment_data/go_cc_results.tsv
    cp ${kegg}                   enrichment_data/kegg_results.tsv
    cp ${gsea_go}                enrichment_data/gsea_go_results.tsv
    cp ${gsea_kegg}              enrichment_data/gsea_kegg_results.tsv
    cp ${splicing_sig}           splicing_data/splicing_significant.tsv
    cp ${wgcna_modules}          wgcna_data/wgcna_modules.tsv
    cp ${wgcna_hub_genes}        wgcna_data/wgcna_hub_genes.tsv
    cp ${integration_candidates} integration_data/key_candidates.tsv

    cp -r ${projectDir}/report/rnaseq_report.qmd .

    quarto render rnaseq_report.qmd \\
        -P deseq2_dir:deseq2_data \\
        -P enrichment_dir:enrichment_data \\
        -P splicing_dir:splicing_data \\
        -P wgcna_dir:wgcna_data \\
        -P integration_dir:integration_data \\
        -P report_title:"${params.report_title}" \\
        -P report_author:"${params.report_author}" \\
        --output rnaseq_report.html
    """
}
