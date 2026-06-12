# Preeclampsia count-matrix statistical analysis

R/DESeq2 workflow for placental RNA-seq count matrices from a preeclampsia study.

This is **Part II** of the project and starts from count matrices.

## Design

- `placenta`: Fetal / Materna
- `grupo`: Control / Temprana / Tardia
- `parto`: Simple / Gemelar
- `sexo_rn`: 0 varón / 1 mujer

Main model within each placental compartment:

```r
~ sexo_rn + parto + grupo
```

Global interaction model:

```r
~ sexo_rn + parto + placenta + grupo + placenta:grupo
```

## Quick start

```bash
mamba env create -f environment.yml
mamba activate preeclampsia-stats
Rscript scripts/run_all.R
```

## Expected inputs

```text
data/raw_counts/genes_counts.csv
data/raw_counts/transcripts_counts.csv   # optional
data/metadata/sample_sheet.csv
```

The count matrix must have genes/transcripts as rows and samples as columns:

```text
gene_id,SAMPLE01,SAMPLE02,...
```

## Analyses

1. Metadata/count matching
2. PCA and sample-distance QC
3. DESeq2 by placental compartment
4. Maternal vs fetal contrasts within clinical groups
5. Placenta × group interaction
6. ncRNA filtering module
7. Isoform/transcript-level module
8. Hallmark enrichment
9. Ordinal gradient: Control → Temprana → Tardia
10. Quarto report

## Outputs

```text
results/deseq2/
results/enrichment/
figures/qc/
figures/volcano/
reports/preeclampsia_report.html
```
