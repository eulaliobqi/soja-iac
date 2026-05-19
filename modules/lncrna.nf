// ============================================================
// Módulo: LNCRNA_PRED – Predição de lncRNAs novos
// Fallback: annotated.gtf → merged.gtf → tracking → merged completo
// ============================================================

process LNCRNA_PRED {
    label 'medium_mem'
    publishDir "${params.outdir}/lncrna", mode: 'copy'

    input:
    path(genome_fa)
    path(annotated_gtf)    // gffcmp.annotated.gtf (pode ser NO_FILE)
    path(merged_gtf)       // merged.gtf
    path(tracking)         // gffcmp.tracking

    output:
    path("novel_transcripts.fa"),    emit: novel_fa, optional: true
    path("novel_transcripts.gtf"),   emit: novel_gtf, optional: true
    path("lncrna_candidates.tsv"),   emit: candidates
    path("lncrna_all.tsv"),          emit: all
    path("lncrna_summary.txt"),      emit: summary
    path("figures/"),                emit: figures

    script:
    // Nunca usar \n em strings Python dentro de Nextflow (Lição 11)
    // Usar rstrip() sem argumento; sep = chr(9) para tab
    """
    mkdir -p figures

    python3 - << 'PYEOF'
import os, re, sys

novel_codes = {'u', 'i', 'x', 's', 'o'}
sep = chr(9)

out_lines   = []
novel_ids   = set()

# Cadeia de fallback (Lição 10):
# 1. gffcmp.annotated.gtf (gffcompare >=0.12)
# 2. gffcmp.merged.gtf (versão antiga)
# 3. Parse tracking + merged.gtf por TCONS novel
# 4. merged.gtf completo

annotated_candidates = ("${annotated_gtf}", "gffcmp.merged.gtf")
for gtf_path in annotated_candidates:
    if os.path.exists(gtf_path) and os.path.getsize(gtf_path) > 0:
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                line = line.rstrip()
                if not line:
                    continue
                parts = line.split(sep)
                if len(parts) < 9:
                    continue
                attr = parts[8]
                # class_code "u","i","x","s","o" = novel
                m = re.search(r'class_code "([^"]+)"', attr)
                if m and m.group(1) in novel_codes:
                    out_lines.append(line)
                    m2 = re.search(r'transcript_id "([^"]+)"', attr)
                    if m2:
                        novel_ids.add(m2.group(1))
        if out_lines:
            print(f"Novel transcripts (from {gtf_path}): {len(novel_ids)}")
            break

if not out_lines and os.path.exists("${tracking}"):
    print("Fallback: parsing tracking file...")
    with open("${tracking}") as fh:
        for line in fh:
            parts = line.rstrip().split(sep)
            if len(parts) >= 4 and parts[3] in novel_codes:
                novel_ids.add(parts[0])
    if novel_ids and os.path.exists("${merged_gtf}"):
        with open("${merged_gtf}") as fh:
            for line in fh:
                line = line.rstrip()
                if any(nid in line for nid in novel_ids):
                    out_lines.append(line)
        print(f"Novel transcripts (from tracking): {len(novel_ids)}")

if not out_lines and os.path.exists("${merged_gtf}"):
    print("Ultimo fallback: usando merged.gtf completo")
    with open("${merged_gtf}") as fh:
        out_lines = [l.rstrip() for l in fh if not l.startswith('#') and l.strip()]

if out_lines:
    with open("novel_transcripts.gtf", "w") as fh:
        fh.write("\\n".join(out_lines) + "\\n")
    print(f"GTF novel escrito: {len(out_lines)} linhas")
else:
    open("novel_transcripts.gtf", "w").close()
    print("Nenhum transcrito novel encontrado")
PYEOF

    # Extrai sequências se GTF não vazio
    if [ -s novel_transcripts.gtf ]; then
        gffread novel_transcripts.gtf -g ${genome_fa} -w novel_transcripts.fa || \
            touch novel_transcripts.fa
    else
        touch novel_transcripts.fa
    fi

    # Predição lncRNA via R
    Rscript ${projectDir}/scripts/11_lncrna.R \
        --fasta       novel_transcripts.fa \
        --min_len_nt  200 \
        --max_orf_aa  100 \
        --outdir      . \
        --figures_dir figures
    """
}
