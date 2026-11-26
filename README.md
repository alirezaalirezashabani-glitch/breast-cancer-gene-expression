# Breast Cancer Gene Expression Analysis

**Author:** Alireza Shabani, PhD in Agricultural Biotechnology

A reproducible RNA-seq workflow examining tumour-versus-adjacent-normal gene expression in an ER-positive breast cancer cohort, with a focused look at apoptosis-associated genes (**BAX, BCL2, CASP3, TP53**).

## Dataset

- **GEO accession:** [GSE103001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103001)
- **Study publication:** [Scientific Reports](https://www.nature.com/articles/s41598-017-17811-2)
- **Samples:** 22 ER-positive primary breast tumours and 22 paired adjacent non-tumour tissues (44 samples; 22 patient pairs).
- **Data:** per-sample HTSeq raw gene counts.
- **Reference annotation reported for the study:** GRCh37.75 / hg19.

The sample-to-patient mapping is recorded in `data/GSE103001_sample_metadata.csv`.

## Analysis design

The main pipeline is `scripts/make_GSE103001_figures.R`. It downloads the GEO supplementary archive, validates the expected sample files, constructs a raw integer count matrix, performs paired differential expression using DESeq2, maps Ensembl identifiers to gene symbols, and generates tables and figures.

The model is:

```r
design = ~ patient_id + group
```

The reported contrast is **Tumor versus Normal**. The volcano plot uses adjusted p-value and absolute log2 fold-change thresholds for visual highlighting; the full result table retains all tested genes and statistics.

## Run locally

Use an R installation with Bioconductor available. The GitHub Actions workflow uses the Bioconductor Docker image and installs:

- Bioconductor: `GEOquery`, `DESeq2`, `AnnotationDbi`, `org.Hs.eg.db`
- CRAN: `ggplot2`, `pheatmap`, `matrixStats`

Then run from the repository root:

```bash
Rscript scripts/make_GSE103001_figures.R
```

Or open the **Actions** tab and run **Run GSE103001 RNA-seq analysis** manually with `workflow_dispatch`.

## Expected outputs

After a successful run, the pipeline writes:

- `data/GSE103001_counts_matrix.csv` — raw integer count matrix.
- `results/GSE103001/differential_expression.csv` — genome-wide DESeq2 results.
- `results/GSE103001/apoptosis_gene_results.csv` — results for BAX, BCL2, CASP3 and TP53 when symbol mapping succeeds.
- `results/GSE103001/volcano_plot.png` — differential-expression overview.
- `results/GSE103001/top_variable_genes_heatmap.png` — top variable genes across samples.
- `results/GSE103001/apoptosis_gene_expression.png` — normalized counts for apoptosis genes when mapped.
- `results/GSE103001/sessionInfo.txt` — R and package versions for reproducibility.

**Important:** output files are not considered available until a successful Actions run has generated and committed them. Check the Actions run status before citing or interpreting numerical results.

## Scientific interpretation and limitations

1. This is a comparison of patient tumour tissue with paired adjacent non-tumour tissue. It is **not** an MCF-7 cell-line treatment experiment.
2. The cohort is ER-positive; findings should not be generalized automatically to every breast-cancer subtype.
3. This dataset cannot establish that beta-cyclocitral, iron-oxide nanoparticles, or any other treatment caused expression changes. Treatment effects require a suitable treated-versus-control experiment.
4. Raw integer counts are used for DESeq2. Normalized counts are used only for visualization.
5. Gene-symbol annotation should be checked for compatibility with the study's genome build before biological interpretation.
6. This workflow is a research analysis, not a clinical diagnostic test.

## Genes of interest

- **BAX:** pro-apoptotic BCL-2 family member.
- **BCL2:** anti-apoptotic BCL-2 family member.
- **CASP3:** executioner caspase.
- **TP53:** tumour-suppressor transcription factor.

## Repository structure

- `data/`: sample metadata and generated count matrix.
- `scripts/`: analysis scripts.
- `docs/GENERATE_FIGURES.md`: execution and output guide.
- `.github/workflows/`: automated GitHub Actions workflow.
- `results/GSE103001/`: generated result tables and figures (created after a successful run).
