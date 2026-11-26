# Breast Cancer Gene Expression Analysis

## Project overview
A reproducible bioinformatics project focused on apoptosis-associated genes in breast cancer, especially **BAX, BCL2, CASP3, and TP53**. The work is motivated by experimental research involving MCF-7 cells, beta-cyclocitral, and nanoparticle-based treatments.

## Selected public dataset: GSE103001
We selected **GSE103001**, a human RNA-seq study of 22 estrogen-receptor-positive (ER+) primary breast cancers and their 22 paired adjacent non-tumour breast tissues (44 samples total). GEO provides per-sample HTSeq gene-count files, making this dataset suitable for a count-based workflow.

- GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103001
- Publication: https://www.nature.com/articles/s41598-017-17811-2
- Genome annotation reported by GEO: GRCh37.75 / hg19
- Study design: paired tumour and adjacent normal tissue from the same patients.

**Scope limitation:** this cohort represents ER-positive primary breast cancer, not all breast-cancer subtypes. It is patient-tissue data, not an MCF-7 cell experiment and not a beta-cyclocitral treatment experiment.

## Repository contents
- `data/GSE103001_sample_metadata.csv`: 44 sample identifiers with paired patient IDs and tissue groups.
- `data/sample_metadata_template.csv`: generic template for future datasets.
- `scripts/analyze_GSE103001.R`: downloads the GEO supplementary count archive, combines per-sample HTSeq counts, performs paired differential-expression analysis with DESeq2, and generates a volcano plot and heatmap.
- `scripts/differential_expression.R`: generic starter workflow for a separately supplied count matrix.

## Run the GSE103001 analysis
Install R, then install the required packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GEOquery"))
install.packages(c("ggplot2", "pheatmap", "matrixStats"))
```

From the repository root, run:

```bash
Rscript scripts/analyze_GSE103001.R
```

The script downloads the GEO supplementary archive and writes:
- `data/GSE103001_counts_matrix.csv`
- `results/GSE103001/differential_expression.csv`
- `results/GSE103001/apoptosis_gene_results.csv`
- `results/GSE103001/volcano_plot.png`
- `results/GSE103001/top_variable_genes_heatmap.png`

The DESeq2 model includes patient identity (`~ patient_id + group`) to account for the paired study design. The reported contrast is Tumour versus Normal.

## Scientific safeguards
- Use the raw integer HTSeq counts for DESeq2; do not substitute TPM, FPKM, or other normalized expression values.
- Review sample pairing, file completeness, count format, and quality-control metrics before interpreting results.
- Gene identifiers in the count files may be Ensembl IDs. Gene-symbol mapping should be performed with an annotation release compatible with the study's GRCh37.75 reference before extracting BAX, BCL2, CASP3, and TP53. The current script's direct symbol filter is a placeholder and may return no rows until this mapping is added.
- A tumour-versus-adjacent-normal association does not show that beta-cyclocitral or nanoparticles caused an expression change. Those effects require appropriate treated-versus-control experimental data.
- This is a research workflow, not a clinical diagnostic tool.

## Genes of interest
- **BAX** — pro-apoptotic BCL-2 family member
- **BCL2** — anti-apoptotic BCL-2 family member
- **CASP3** — executioner caspase
- **TP53** — tumour-suppressor transcription factor

## Author
Alireza Shabani  
PhD in Biotechnology
