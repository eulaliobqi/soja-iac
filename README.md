# RNASeq Insight Platform – *Glycine max*

Pipeline automatizado e reprodutível para análise completa de RNA-Seq em soja, incluindo expressão diferencial, splicing alternativo, enriquecimento funcional, co-expressão gênica e dashboard interativo.

---

## Visão Geral

```
FASTQ
  │
  ├─ FastQC (pré) ──────────────────────────────────── MultiQC
  │
  ├─ fastp (trimagem)
  │
  ├─ FastQC (pós) ─────────────────────────────────── MultiQC
  │
  ├─ HISAT2 (alinhamento splice-aware)
  │
  ├─ samtools (sort + index + flagstat)
  │
  ├─ featureCounts ──── DESeq2 ──── GO/KEGG/GSEA
  │                        │
  │                      WGCNA
  │
  └─ rMATS (splicing) ────────── Integração Multi-ômica
                                        │
                               Dashboard Shiny + Relatório Quarto
```

---

## Funcionalidades

| Módulo | Ferramenta | Output |
|---|---|---|
| QC | FastQC + MultiQC | Relatórios HTML |
| Trimagem | fastp | Reads limpos + métricas |
| Alinhamento | HISAT2 | BAM ordenados e indexados |
| Quantificação | featureCounts | Matriz de contagens |
| Expressão diferencial | DESeq2 | Tabelas + PCA, volcano, heatmap, MA |
| Enriquecimento | clusterProfiler | GO-BP/MF/CC, KEGG, GSEA |
| Splicing alternativo | rMATS | SE, A5SS, A3SS, MXE, RI |
| Co-expressão | WGCNA | Módulos + hub genes |
| Integração | Script R | Gene ranking + candidatos |
| Dashboard | Shiny + Plotly | App interativo |
| Relatório | Quarto | HTML completo |

---

## Requisitos

- Linux (Ubuntu 20.04+ recomendado)
- [Miniforge / Mamba](https://github.com/conda-forge/miniforge)
- [Nextflow](https://www.nextflow.io/) ≥ 24.04
- Java ≥ 11

---

## Instalação

```bash
# Clone o repositório
git clone https://github.com/eulaliobqi/soja-iac.git
cd soja-iac

# Instala ambientes e valida ferramentas
bash setup.sh
```

---

## Uso

### 1. Configure suas amostras

Edite [samplesheet.csv](samplesheet.csv):

```csv
sample,fastq_1,fastq_2,condition,replicate
ctrl_rep1,data/ctrl_rep1_R1.fastq.gz,data/ctrl_rep1_R2.fastq.gz,control,1
ctrl_rep2,data/ctrl_rep2_R1.fastq.gz,data/ctrl_rep2_R2.fastq.gz,control,2
ctrl_rep3,data/ctrl_rep3_R1.fastq.gz,data/ctrl_rep3_R2.fastq.gz,control,3
treat_rep1,data/treat_rep1_R1.fastq.gz,data/treat_rep1_R2.fastq.gz,treatment,1
treat_rep2,data/treat_rep2_R1.fastq.gz,data/treat_rep2_R2.fastq.gz,treatment,2
treat_rep3,data/treat_rep3_R1.fastq.gz,data/treat_rep3_R2.fastq.gz,treatment,3
```

### 2. Configure os parâmetros

Edite [params.yaml](params.yaml):

```yaml
genome_fasta: "/path/to/Gmax_508_v4.0.fa"
genome_gff3:  "/path/to/Gmax_508_Wm82.a4.v1.gene_exons.gff3"
read_length:  150
```

### 3. Execute o pipeline

```bash
# Execução local
bash run_pipeline.sh local

# Execução em cluster SLURM
bash run_pipeline.sh slurm

# Retomar execução interrompida
bash run_pipeline.sh local --resume
```

Ou diretamente com Nextflow:

```bash
nextflow run main.nf \
  -profile local \
  -params-file params.yaml \
  -with-report results/nextflow_report.html \
  -with-trace  results/nextflow_trace.txt \
  -with-timeline results/nextflow_timeline.html
```

---

## Estrutura de Resultados

```
results/
├── qc/
│   ├── pre_trim/          FastQC antes da trimagem
│   ├── post_trim/         FastQC após trimagem
│   └── multiqc/           Relatórios MultiQC consolidados
├── trimmed/               Reads trimados + relatórios fastp
├── aligned/               BAMs ordenados, indexados, flagstat
├── counts/                Matriz de contagens (featureCounts)
├── deseq2/
│   ├── deseq2_results_all.tsv
│   ├── deseq2_results_sig.tsv
│   ├── normalized_counts.tsv
│   └── figures/           PCA, volcano, heatmap, MA plot
├── enrichment/
│   ├── go_bp_results.tsv
│   ├── kegg_results.tsv
│   ├── gsea_*.tsv
│   └── figures/
├── splicing/
│   ├── splicing_significant.tsv
│   ├── splicing_all.tsv
│   └── rmats_output/
├── wgcna/
│   ├── wgcna_modules.tsv
│   ├── wgcna_hub_genes.tsv
│   └── figures/
├── integration/
│   ├── gene_ranking.tsv
│   ├── key_candidates.tsv
│   └── figures/
├── report/
│   └── rnaseq_report.html  Relatório completo
└── nextflow_*.html         Logs e timeline Nextflow
```

---

## Dashboard Interativo

Após a análise, visualize os resultados no dashboard:

```bash
RESULTS_DIR=results Rscript -e "shiny::runApp('dashboard/app.R', port=3838)"
```

Acesse em `http://localhost:3838`

**Funcionalidades do dashboard:**
- PCA interativo
- Volcano plot com filtros dinâmicos
- MA plot
- Tabelas filtráveis e exportáveis (CSV/Excel)
- Dotplots GO e KEGG
- Eventos de splicing
- Módulos WGCNA e hub genes
- Ranking de genes candidatos

---

## Genoma de Referência

Para *Glycine max* Wm82.a4.v1:

```bash
# Download do Phytozome (requer conta gratuita)
# https://phytozome-next.jgi.doe.gov/

# Ou NCBI RefSeq:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.gff.gz
```

---

## Critérios de Análise

| Parâmetro | Valor padrão | Ajustável em |
|---|---|---|
| FDR (padj) | < 0.05 | `params.yaml` |
| \|log2FC\| | > 1.0 | `params.yaml` |
| Taxa alinhamento mínima | 70% | `modules/alignment.nf` |
| rMATS FDR | < 0.05 | `scripts/04_integration.R` |
| rMATS \|ΔPSI\| | > 0.10 | `scripts/04_integration.R` |

---

## Citações

Se utilizar este pipeline, cite:

- **DESeq2**: Love MI, Huber W, Anders S (2014). *Genome Biology*, 15:550.
- **HISAT2**: Kim D et al. (2019). *Nature Methods*, 16:3.
- **clusterProfiler**: Wu T et al. (2021). *The Innovation*, 2(3):100141.
- **rMATS**: Shen S et al. (2014). *PNAS*, 111:E5593.
- **WGCNA**: Langfelder P, Horvath S (2008). *BMC Bioinformatics*, 9:559.
- **fastp**: Chen S et al. (2018). *Bioinformatics*, 34:i884.

---

## Autor

**Eulalio Santos** | Universidade Federal de Viçosa  
GitHub: [@eulaliobqi](https://github.com/eulaliobqi)  
Email: eulalio.santos@ufv.br

---

*Desenvolvido com [Claude Code](https://claude.ai/claude-code) – Anthropic*
