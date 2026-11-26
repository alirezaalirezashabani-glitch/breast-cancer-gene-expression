# Generate the GSE103001 figures

The repository contains `scripts/make_GSE103001_figures.R`, a reproducible paired RNA-seq workflow. It downloads the public GEO count archive, validates sample pairing and raw counts, runs DESeq2, maps Ensembl IDs to gene symbols, and generates analysis tables and figures.

## 1. Install R packages

Run in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2", "AnnotationDbi", "org.Hs.eg.db"))
install.packages(c("ggplot2", "pheatmap", "matrixStats"))
```

## 2. Run from the repository folder

```bash
Rscript scripts/make_GSE103001_figures.R
```

The script generates:

- `data/GSE103001_counts_matrix.csv`
- `results/GSE103001/differential_expression.csv`
- `results/GSE103001/apoptosis_gene_results.csv`
- `results/GSE103001/volcano_plot.png`
- `results/GSE103001/top_variable_genes_heatmap.png`
- `results/GSE103001/pca_plot.png`
- `results/GSE103001/apoptosis_gene_expression.png` (when target genes map)
- `results/GSE103001/analysis_sample_metadata.csv`
- `results/GSE103001/sessionInfo.txt`

The DESeq2 model is paired (`~ patient_id + group)), comparing tumour with adjacent normal tissue. The script checks that each patient has exactly one sample in each group, uses raw integer counts for DESeq2, and uses normalized counts only for visualization.

## 3. Download outputs from GitHub Actions

The workflow automatically uploads generated outputs as the artifact `GSE103001-RNA-seq-results`:

1. Open the repository's **Actions** tab.
2. Select the latest successful **Run GSE103001 RNA-seq analysis**.
3. Download the artifact from the **Artifacts** section of the run summary.

Artifacts are retained for 30 days. The workflow does not try to commit results from inside the Bioconductor container; this avoids the earlier Git working-directory error.

**Interpretation boundary:** GSE103001 is a paired patient-tissue cohort of ER-positive breast cancers. It is not an MCF-7 beta-cyclocitral/nanoparticle treatment experiment, so its results cannot establish treatment effects.
