# Preliminary results: GSE103001 paired RNA-seq

**Analysis run:** GitHub Actions run #5 (26 November 2026), successful  
**Dataset:** 22 ER-positive breast tumour samples and 22 matched adjacent non-tumour samples  
**Model:** DESeq2, paired design `~ patient_id + group`; contrast Tumour vs Normal

## Genome-wide overview

- 48,517 genes were retained in the differential-expression result table.
- 31,585 genes had a non-missing adjusted p-value.
- 4,822 genes met both thresholds: adjusted p-value (Benjamini–Hochberg) < 0.05 and absolute log2 fold change ≥ 1.

The downloadable workflow artifact contains the complete differential-expression table, raw count matrix, sample metadata, and plots. These numerical summaries are reported from the successful run artifact; they are not simulated.

## Apoptosis-associated genes

| Gene | log2 fold change (Tumour / Normal) | Approx. fold change | Adjusted p-value | Interpretation at FDR < 0.05 |
| --- | ---: | ---: | ---: | --- |
| CASP3 | +0.484 | 1.40× | 2.96 × 10⁻⁷ | Higher in tumour |
| BAX | +0.526 | 1.44× | 0.00155 | Higher in tumour |
| TP53 | +0.385 | 1.31× | 0.05003 | Not significant at strict 0.05 threshold |
| BCL2 | −0.140 | 0.91× | 0.573 | Not significant |

## Interpretation

BAX and CASP3 show higher expression in the tumour tissue cohort relative to paired adjacent tissue. TP53 is borderline but does not pass the prespecified adjusted-p-value threshold of 0.05; BCL2 is not significantly different in this analysis.

These are **observational tissue-level differences**, not evidence that apoptosis was experimentally induced. The dataset does not contain beta-cyclocitral or iron-oxide nanoparticle treatment groups and cannot test the user's MCF-7 treatment hypothesis. Differences in cell composition, tumour biology, and adjacent-tissue characteristics may contribute to the observed signals.

## Reproducibility and limitations

- Raw HTSeq integer counts were used as input to DESeq2.
- Normalized counts were used for visualization.
- The design accounts for patient pairing.
- Gene-symbol annotation and compatibility with the reported GRCh37.75 / hg19 reference should be reviewed before publication.
- This analysis is exploratory and is not a clinical diagnostic or treatment-response test.

## Download

The full result bundle is available from the successful GitHub Actions run under the artifact name **GSE103001-RNA-seq-results**. Artifacts are retained for 30 days.
