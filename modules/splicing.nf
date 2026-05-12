// ============================================================
// Módulo: Splicing Alternativo – rMATS
// ============================================================

process RMATS {
    publishDir "${params.outdir}/splicing", mode: 'copy'

    input:
    path(ctrl_bams)     // lista de BAMs do grupo controle
    path(treat_bams)    // lista de BAMs do grupo tratamento
    path(gtf)

    output:
    path("rmats_output/"),      emit: results_dir
    path("rmats_summary.tsv"),  emit: summary

    script:
    def b1 = ctrl_bams.collect  { it.toString() }.join(',')
    def b2 = treat_bams.collect { it.toString() }.join(',')
    def lib_type = params.rmats_type == "paired" ? "fr-firststrand" : "fr-unstranded"
    """
    # Cria arquivos de lista de BAMs (rMATS exige CSV em linha única)
    echo "${b1}" > b1_bams.txt
    echo "${b2}" > b2_bams.txt

    rmats.py \\
        --b1 b1_bams.txt \\
        --b2 b2_bams.txt \\
        --gtf ${gtf} \\
        --od rmats_output \\
        --tmp rmats_tmp \\
        -t ${params.rmats_type} \\
        --libType ${lib_type} \\
        --readLength ${params.read_length} \\
        --nthread ${task.cpus} \\
        --tstat ${task.cpus} \\
        --statoff

    # Gera resumo dos eventos significativos
    python3 - <<'PYEOF'
import os, glob, pandas as pd

event_types = ['SE','A5SS','A3SS','MXE','RI']
summary = []

for evt in event_types:
    f = f"rmats_output/{evt}.MATS.JC.txt"
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f, sep="\\t")
    if df.empty:
        continue
    df['FDR'] = pd.to_numeric(df.get('FDR', df.get('PValue', 1)), errors='coerce').fillna(1)
    df['IncLevelDifference'] = pd.to_numeric(df.get('IncLevelDifference', 0), errors='coerce').fillna(0)
    sig = df[(df['FDR'] < 0.05) & (df['IncLevelDifference'].abs() > 0.1)]
    summary.append({'event_type': evt, 'total': len(df), 'significant': len(sig)})

pd.DataFrame(summary).to_csv("rmats_summary.tsv", sep="\\t", index=False)
print("rMATS summary:")
print(pd.DataFrame(summary).to_string(index=False))
PYEOF
    """
}

process PARSE_RMATS {
    publishDir "${params.outdir}/splicing", mode: 'copy'

    input:
    path(rmats_dir)

    output:
    path("splicing_significant.tsv"), emit: significant
    path("splicing_all.tsv"),         emit: all_events

    script:
    """
    #!/usr/bin/env python3
    import os, glob, pandas as pd

    event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']
    all_dfs, sig_dfs = [], []

    for evt in event_types:
        fpath = os.path.join("${rmats_dir}", f"{evt}.MATS.JC.txt")
        if not os.path.exists(fpath):
            continue
        df = pd.read_csv(fpath, sep="\\t")
        if df.empty:
            continue
        df['event_type'] = evt
        df['FDR'] = pd.to_numeric(df.get('FDR', df.get('PValue', 1)), errors='coerce').fillna(1)
        df['IncLevelDifference'] = pd.to_numeric(
            df.get('IncLevelDifference', 0), errors='coerce').fillna(0)
        all_dfs.append(df)
        sig = df[(df['FDR'] < 0.05) & (df['IncLevelDifference'].abs() > 0.1)]
        sig_dfs.append(sig)

    if all_dfs:
        pd.concat(all_dfs, ignore_index=True).to_csv("splicing_all.tsv", sep="\\t", index=False)
    if sig_dfs:
        pd.concat(sig_dfs, ignore_index=True).to_csv("splicing_significant.tsv", sep="\\t", index=False)
        print(f"Eventos significativos: {sum(len(d) for d in sig_dfs)}")
    else:
        pd.DataFrame().to_csv("splicing_significant.tsv", sep="\\t", index=False)
        print("Nenhum evento de splicing significativo encontrado.")
    """
}
