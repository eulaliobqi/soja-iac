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

process GOSEQ {
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path(deseq2_all)
    path(gtf)

    output:
    path("goseq_bp_results.tsv"), emit: go_bp
    path("goseq_mf_results.tsv"), emit: go_mf
    path("goseq_cc_results.tsv"), emit: go_cc

    script:
    def go_annot_arg = params.go_annot ? "--go_annot ${params.go_annot}" : ""
    """
    Rscript ${projectDir}/scripts/02b_goseq.R \
        --deseq2 ${deseq2_all} \
        --gtf    ${gtf} \
        --padj   ${params.padj_cutoff} \
        --lfc    ${params.lfc_cutoff} \
        --outdir . \
        ${go_annot_arg}
    """
}

process MACHINE_LEARNING {
    label 'medium_mem'
    publishDir "${params.outdir}/ml", mode: 'copy'

    input:
    path(norm_counts)
    path(deseq2_sig)
    path(metadata)

    output:
    path("ml_results.tsv"),        emit: results
    path("feature_importance.tsv"),emit: features
    path("figures/"),              emit: figures

    script:
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/06_machine_learning.R \
        --norm_counts ${norm_counts} \
        --deseq2_sig  ${deseq2_sig} \
        --metadata    ${metadata} \
        --outdir      . \
        --figures_dir figures
    """
}

process PPI_NETWORK {
    label 'medium_mem'
    publishDir "${params.outdir}/network", mode: 'copy'

    input:
    path(deseq2_sig)

    output:
    path("ppi_edges.tsv"),      emit: edges
    path("ppi_nodes.tsv"),      emit: nodes
    path("hub_genes.tsv"),      emit: hub_genes
    path("network_summary.txt"),emit: summary
    path("figures/"),           emit: figures

    script:
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/07_ppi_network.R \
        --deseq2_sig  ${deseq2_sig} \
        --outdir      . \
        --figures_dir figures
    """
}

process PLANTFDB {
    label 'low_mem'
    publishDir "${params.outdir}/plantfdb", mode: 'copy'

    input:
    path(deseq2_sig)
    path(deseq2_all)

    output:
    path("tf_deg_classified.tsv"),   emit: tf_classified
    path("tf_family_summary.tsv"),   emit: family_summary
    path("tf_family_enrichment.tsv"),emit: family_enrichment
    path("plantfdb_summary.txt"),    emit: summary
    path("figures/"),                emit: figures

    script:
    def pf_arg = params.plantfdb_file ? "--plantfdb_file ${params.plantfdb_file}" : ""
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/08_plantfdb.R \
        --deseq2_sig ${deseq2_sig} \
        --deseq2_all ${deseq2_all} \
        --outdir      . \
        --figures_dir figures \
        ${pf_arg}
    """
}

process GENIE3 {
    label 'high_mem'
    publishDir "${params.outdir}/genie3", mode: 'copy'

    input:
    path(norm_counts)
    path(tf_classified)

    output:
    path("genie3_network.tsv"),  emit: network
    path("genie3_hub_tfs.tsv"),  emit: hub_tfs
    path("genie3_summary.txt"),  emit: summary
    path("figures/"),            emit: figures

    script:
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/09_genie3.R \
        --norm_counts   ${norm_counts} \
        --tf_classified ${tf_classified} \
        --n_trees       ${params.genie3_trees ?: 500} \
        --n_links       ${params.genie3_links ?: 5000} \
        --ncores        ${task.cpus} \
        --outdir        . \
        --figures_dir   figures
    """
}

process METANALYSIS {
    label 'medium_mem'
    publishDir "${params.outdir}/metanalysis", mode: 'copy'

    input:
    path(deseq2_sig)

    output:
    path("metanalysis_validated_genes.tsv"), emit: validated
    path("metanalysis_overlap.tsv"),         emit: overlap
    path("metanalysis_summary.tsv"),         emit: summary
    path("metanalysis_report.txt"),          emit: report
    path("figures/"),                        emit: figures

    script:
    def geo_arg = params.geo_accessions ? "--geo_accessions '${params.geo_accessions}'" : ""
    """
    mkdir -p figures
    Rscript ${projectDir}/scripts/10_metanalysis.R \
        --deseq2_sig ${deseq2_sig} \
        --outdir      . \
        --figures_dir figures \
        ${geo_arg}
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
