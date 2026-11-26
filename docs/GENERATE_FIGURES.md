# Generate the GSE103001 figures

The repository now includes `scripts/make_GSE103001_figures.R`. It downloads the public GEO count archive, uses the paired-patient design, maps Ensembl IDs to gene symbols, and generates the analysis tables and plots.

## 1. Install R packages

Run in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2", "org.Hs.eg.db"))
install.packages(c("ggplot2", "pheatmap", "matrixStats"))
```

## 2. Run from the repository folder

```bash
Rscript scripts/make_GSE103001_figures.R
```

The script downloads the public supplementary archive and creates:

- `results/GSE103001/differential_expression.csv`
- `results/GSE103001/apoptosis_gene_results.csv`
- `results/GSE103001/volcano_plot.png`
- `results/GSE103001/top_variable_genes_heatmap.png`
- `results/GSE103001/apoptosis_gene_expression.png` (if the four symbols map successfully)
- `data/GSE103001_counts_matrix.csv`

The DESeq2 model is paired (`~ patient_id + group`), comparing tumour with adjacent normal tissue. Plots are only created after the R script has been run with the real GEO count files; no simulated expression values or placeholder results are presented as experimental findings.

**Interpretation boundary:** GSE103001 is a patient-tissue cohort of ER-positive breast cancers. It is not an MCF-7 beta-cyclocitral/nanoparticle treatment experiment, so its results cannot establish treatment effects.
