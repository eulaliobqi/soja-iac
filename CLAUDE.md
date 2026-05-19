# RNASeq Insight Platform – *Glycine max* Wm82.a4.v1

Pipeline Nextflow DSL2 para análise de RNA-Seq em soja (*Glycine max* Wm82.a4.v1).
Adaptado do pipeline Arabidopsis TAIR10 com todas as correções de bugs documentadas.

## Parâmetros organismo-específicos

```yaml
kegg_organism:   "gmx"          # código KEGG correto para Glycine max
string_species:  3847           # taxon ID STRING
plantfdb_prefix: "Gma"
biomart_dataset: "gmax_eg_gene" # Ensembl Plants
```

**IDs de genes:** formato `Glyma.01G000100` — sempre strip do sufixo numérico:
```r
gene_id <- gsub("\\.[0-9]+$", "", gene_id)
```

**OrgDb:** `org.Gmax.eg.db` removido do Bioconductor 3.20+. Usar `biomaRt` + `enricher()` com TERM2GENE manual.

## Agentes especializados

- `/validate-qc`       — Valida QC (FastQC/MultiQC + STAR + Salmon)
- `/interpret-results` — Interpretação biológica em contexto de Glycine max
- `/debug-pipeline`    — Diagnóstico de falhas Nextflow (logs + work dir)
- `/review-code`       — Revisão de código R/Nextflow antes de commit
- `/check-ids`         — Validação de consistência dos gene IDs Glyma

### `/validate-qc`

Verificar após MULTIQC_POST:
- Per-base Q ≥ Q20 em >90% das posições
- %GC ±5% do esperado (~44% para Glycine max)
- Duplicação < 50%
- Reads pós-trimagem ≥ 36 bp
- STAR alignment rate ≥ 70%
- Salmon mapping rate ≥ 60%
- Correlação inter-réplicas ≥ 0.95

```bash
grep -i "overall alignment rate" logs/star_align_*.log
grep "Mapping rate" results/counts/salmon/*/logs/salmon_quant.log
```

### `/interpret-results`

Trigger: Após INTEGRATION + ENRICHMENT.

1. DEGs: up vs down; genes |lfc| > 2 mais responsivos; GmWRKY, GmERF, GmNAC
2. GO-BP: resposta a estresse, desenvolvimento, isoflavonoides
3. KEGG: nitrogênio (gmx00910), flavonoides (gmx00941), fotossíntese
4. WGCNA: módulos |r| > 0.7, p < 0.05; hub genes
5. Integração: score > 4.0, ≥3 camadas → candidatos fortes
6. TFs: WRKY, ERF/AP2, MYB, NAC, bHLH, ARF

### `/debug-pipeline`

Trigger: Qualquer processo com exit status ≠ 0.

```bash
cat .nextflow.log | grep "ERROR\|failed\|exit"
cd work/XX/YYYYYY*
cat .command.sh; cat .command.out; cat .command.err
```

| Erro | Causa | Solução |
|------|-------|---------|
| exit 139 STAR | STAR-avx2 sem AVX2 | `star=2.7.10b` em rnaseq-tools.yml |
| Missing output file | Script R falhou | Ler `.command.err` |
| cannot rescale constant column | prcomp variância zero | `vst_mat[apply(vst_mat,1,var)>0,]` |
| makeTxDbFromGFF deprecated | GenomicFeatures ≥1.61 | helper `make_txdb()` com txdbmaker |
| FileNotFoundError: gffcmp.merged.gtf | gffcompare ≥0.12 | fallback para `.annotated.gtf` |

### `/review-code`

Checklist obrigatório antes de `git push`:

**R:**
- [ ] Outputs inicializados vazios no início?
- [ ] gene_id strip de versão (.1, .2)?
- [ ] prcomp precedido de filtro var > 0?
- [ ] makeTxDbFromGFF com helper make_txdb()?
- [ ] GOseq com gene2cat customizado (NÃO genome/id)?
- [ ] clusterProfiler: enricher() com TERM2GENE (SEM OrgDb)?
- [ ] GENIE3 com run_genie3() e fallback nCores=1?
- [ ] STRINGdb com species=3847?
- [ ] PlantTFDB prefixo Gma_ ?
- [ ] KEGG organism="gmx"?
- [ ] tryCatch em chamadas biomaRt?

**Nextflow:**
- [ ] Nenhum `\n` em strings Python dentro de heredoc?
- [ ] rstrip() sem argumento?
- [ ] b1.txt/b2.txt: BAMs em uma linha com vírgulas?
- [ ] gffread → `sed 's/^chrChr/Chr/'` após conversão?
- [ ] LNCRNA_PRED com fallback para gffcmp.annotated.gtf?

### `/check-ids`

```r
head(gene_ids)  # Esperado: "Glyma.01G000100"
test_map <- getBM(attributes="ensembl_gene_id", filters="ensembl_gene_id",
                  values=gene_ids[1:10], mart=mart)
cat("IDs reconhecidos:", nrow(test_map), "/ 10\n")
go_coverage <- sum(gene_ids %in% go_map$ensembl_gene_id) / length(gene_ids)
cat("Cobertura GO:", round(go_coverage*100, 1), "%\n")  # Esperado: >60%
```

## Critérios de qualidade

```
padj_cutoff:  0.05
lfc_cutoff:   1.0
splicing_fdr: 0.05
splicing_dpsi: 0.1
wgcna_r2:     0.85
key_candidates: ≥2 camadas de evidência
```

## Lições críticas (resumo)

1. **STAR:** usar 2.7.10b (não 2.7.11b — SIGSEGV em CPUs sem AVX2)
2. **Gene IDs:** strip `gsub("\\.[0-9]+$", "")` em todos os scripts
3. **Cromossomos:** gffread → `sed 's/^chrChr/Chr/'`
4. **OrgDb:** não usar `org.Gmax.eg.db` — usar biomaRt + enricher()
5. **ComBat-Seq:** filtrar var > 0 antes de prcomp
6. **GOseq:** nunca usar `genome="gma"` — sempre `gene2cat` customizado
7. **txdbmaker:** helper `make_txdb()` compatível com GenomicFeatures ≥1.61
8. **GOseq nullp:** fallback `method="Hypergeometric"` se < 6 comprimentos únicos
9. **rMATS:** b1.txt/b2.txt uma linha, vírgulas
10. **lncRNA:** fallback gffcmp.annotated.gtf → merged.gtf → tracking
11. **Heredoc:** nunca `\n` em Python dentro de Nextflow; usar `rstrip()`
12. **GENIE3:** sempre `run_genie3()` com fallback nCores=1
13. **Outputs:** inicializar vazios no início de todo script R
14. **STRINGdb:** species=3847; inicializar outputs antes de tentar conexão
