# Breast Cancer Gene Expression Analysis

## Project overview
A reproducible bioinformatics project for exploring apoptosis-associated gene expression in breast cancer, with a biological focus on **BAX, BCL2, CASP3, and TP53**. The project is motivated by experimental work involving MCF-7 cells, beta-cyclocitral, and nanoparticle-based treatments.

## Current repository contents
- `data/sample_metadata_template.csv`: example sample sheet. Replace the example IDs and groups with the real sample information.
- `scripts/differential_expression.R`: starter RNA-seq differential-expression workflow using DESeq2, with volcano-plot and heatmap outputs.

## Planned workflow
1. Obtain an appropriate public breast-cancer RNA-seq dataset (or use a properly prepared count matrix).
2. Document sample identifiers, groups, and relevant clinical/technical covariates.
3. Perform quality control and exploratory analysis.
4. Test differential expression with an appropriate statistical model.
5. Visualize results (PCA, volcano plot, heatmap, and selected apoptosis genes).
6. Interpret findings in the context of apoptosis biology and the experimental MCF-7 research.

## Input requirements
The starter R script expects:
- `data/counts_matrix.csv`: a gene-by-sample matrix; first column is gene ID and remaining columns are sample IDs.
- `data/sample_metadata_template.csv`: columns named `sample_id` and `group`; sample IDs must exactly match count-matrix column names.

**Important:** DESeq2 requires raw, un-normalized RNA-seq counts (or supported estimated counts), not TPM, FPKM, or already library-size-normalized values. Keep the sample metadata rows aligned to the count-matrix columns. See the official [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

## Run the starter analysis
Install R and the required packages, then run:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages(c("ggplot2", "pheatmap", "matrixStats"))
```

From the repository root:

```bash
Rscript scripts/differential_expression.R
```

The script writes results to `results/`:
- `differential_expression_results.csv`
- `volcano_plot.png`
- `top_variable_genes_heatmap.png`

## Scientific safeguards
- No biological conclusions are drawn from the example metadata template.
- Do not run statistical testing until real, correctly matched sample-level data are supplied.
- Confirm the dataset platform, count-generation method, sample groups, biological replicates, and covariates before choosing the final model.
- Public cohort associations do not establish that beta-cyclocitral or nanoparticles caused an expression change. Such causal claims require suitable experimental data.

## Genes of interest
- **BAX** — pro-apoptotic BCL-2 family member
- **BCL2** — anti-apoptotic BCL-2 family member
- **CASP3** — executioner caspase
- **TP53** — tumor-suppressor transcription factor

## Author
Alireza Shabani  
PhD in Biotechnology
