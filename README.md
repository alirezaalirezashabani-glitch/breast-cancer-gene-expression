# Breast Cancer Gene Expression Analysis

**Author:** Alireza Shabani, PhD in Agricultural Biotechnology

A reproducible paired RNA-seq workflow for tumour-versus-adjacent-normal analysis in an ER-positive breast cancer cohort, with a focused analysis of apoptosis-associated genes (**BAX, BCL2, CASP3, TP53**).

## Dataset

- **GEO accession:** [GSE103001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103001)
- **Study publication:** [Scientific Reports](https://www.nature.com/articles/s41598-017-17811-2)
- **Samples:** 22 ER-positive primary breast tumours and 22 paired adjacent non-tumour tissues (44 samples; 22 patient pairs).
- **Data:** per-sample HTSeq raw gene counts.
- **Reference annotation reported for the study:** GRCh37.75 / hg19.
- **Sample-to-patient mapping:** `data/GSE103001_sample_metadata.csv`.

## Reproducible analysis

The main pipeline is `scripts/make_GSE103001_figures.R`. It downloads the GEO supplementary archive, validates the sample files and paired design, builds an integer count matrix, performs paired differential expression with DESeq2, maps Ensembl IDs to gene symbols, and generates tables and figures.

The model is:

```r
design = ~ patient_id + group
```

The contrast is **Tumor versus Normal**. The volcano plot highlights genes using adjusted p-value and absolute log2 fold-change thresholds; the full result table retains all tested genes and statistics.

## Run locally

Use R with a compatible Bioconductor installation. The workflow uses the Bioconductor Docker image `RELEASE_3_20` and installs:

- Bioconductor: `GEOquery`, `DESeq2`, `AnnotationDbi`, `org.Hs.eg.db`
- CRAN: `ggplot2`, `pheatmap`, `matrixStats`

From the repository root:

```bash
Rscript scripts/make_GSE103001_figures.R
```

Or open the **Actions** tab and run **Run GSE103001 RNA-seq analysis** manually with `workflow_dispatch`.

## Download generated outputs

The workflow uploads the real analysis outputs as a GitHub Actions artifact named **GSE103001-RNA-seq-results**:

1. Open the repository's **Actions** tab.
2. Select the latest successful **Run GSE103001 RNA-seq analysis** run.
3. In the run summary, find **Artifacts** and download `GSE103001-RNA-seq-results`.

Artifacts are retained for 30 days. Generated results are not automatically committed to the source branch; this keeps the workflow independent of write permissions and avoids the previous `fatal: not in a git directory` failure in the container job.

## Expected outputs

- `data/GSE103001_counts_matrix.csv` — raw integer count matrix.
- `results/GSE103001/differential_expression.csv` — genome-wide DESeq2 results.
- `results/GSE103001/apoptosis_gene_results.csv` — results for BAX, BCL2, CASP3 and TP53 when mapped.
- `results/GSE103001/volcano_plot.png` — differential-expression overview.
- `results/GSE103001/top_variable_genes_heatmap.png` — top variable genes.
- `results/GSE103001/pca_plot.png` — PCA of variance-stabilized counts.
- `results/GSE103001/apoptosis_gene_expression.png` — paired normalized expression plot for apoptosis genes when mapped.
- `results/GSE103001/analysis_sample_metadata.csv` — sample/patient/group mapping used in analysis.
- `results/GSE103001/sessionInfo.txt` — R and package versions.

**Do not cite numerical findings until the latest workflow run succeeds and the downloaded result tables have been checked.**

## Scientific interpretation and limitations

1. This is a comparison of patient tumour tissue with paired adjacent non-tumour tissue. It is **not** an MCF-7 cell-line treatment experiment.
2. The cohort is ER-positive; findings should not automatically be generalized to every breast-cancer subtype.
3. This dataset cannot establish that beta-cyclocitral, iron-oxide nanoparticles, or another treatment caused expression changes. Treatment effects require an appropriate treated-versus-control experiment.
4. DESeq2 is supplied raw integer counts. Normalized counts are used only for visualization.
5. Gene-symbol annotation and the study genome build should be reviewed before biological interpretation.
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
- `results/GSE103001/`: generated tables and figures, included in the workflow artifact after a successful run.
