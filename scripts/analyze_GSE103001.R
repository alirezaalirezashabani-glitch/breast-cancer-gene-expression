#!/usr/bin/env Rscript

# Download and analyse GSE103001: 22 paired ER-positive breast tumours
# and adjacent non-tumour breast tissues.
# GEO source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103001

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
})

dir.create("data/GSE103001", recursive = TRUE, showWarnings = FALSE)
dir.create("results/GSE103001", recursive = TRUE, showWarnings = FALSE)

meta <- read.csv("data/GSE103001_sample_metadata.csv",
                 stringsAsFactors = FALSE)
meta$group <- factor(meta$group, levels = c("Normal", "Tumor"))
meta$patient_id <- factor(meta$patient_id)

# Download the GEO supplementary archive (approximately 9.2 MB at the time
# this dataset was checked) and extract the per-sample HTSeq count files.
supp <- getGEOSuppFiles("GSE103001", makeDirectory = FALSE,
                        baseDir = "data/GSE103001")
archive <- rownames(supp)[grepl("GSE103001_RAW.tar$", rownames(supp))]
if (length(archive) == 0) {
  candidates <- list.files("data/GSE103001", pattern = "RAW.tar$",
                           recursive = TRUE, full.names = TRUE)
  if (length(candidates) == 0) stop("Could not locate GSE103001_RAW.tar.")
  archive <- candidates[1]
}
untar(archive[1], exdir = "data/GSE103001/extracted")

files <- list.files("data/GSE103001/extracted",
                    pattern = "\\.(txt|txt.gz)$",
                    recursive = TRUE, full.names = TRUE)
count_files <- setNames(files[
  vapply(files, function(f) {
    id <- regmatches(basename(f), regexpr("GSM[0-9]+", basename(f)))
    length(id) > 0 && id %in% meta$sample_id
  }, logical(1))
], vapply(files[
  vapply(files, function(f) {
    id <- regmatches(basename(f), regexpr("GSM[0-9]+", basename(f)))
    length(id) > 0 && id %in% meta$sample_id
  }, logical(1))
], function(f) regmatches(basename(f), regexpr("GSM[0-9]+", basename(f))),
character(1)))

if (!all(meta$sample_id %in% names(count_files))) {
  missing <- setdiff(meta$sample_id, names(count_files))
  stop("Missing count files for: ", paste(missing, collapse = ", "),
       ". Check the extracted GEO archive and file names.")
}

read_htseq <- function(path) {
  x <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  x <- x[!grepl("^__", x[[1]]), 1:2]
  colnames(x) <- c("gene_id", "count")
  x
}

tables <- lapply(meta$sample_id, function(id) read_htseq(count_files[[id]]))
names(tables) <- meta$sample_id
gene_ids <- Reduce(intersect, lapply(tables, `[[`, "gene_id"))
count_mat <- do.call(cbind, lapply(tables, function(x) {
  x[match(gene_ids, x$gene_id), "count"]
}))
rownames(count_mat) <- gene_ids
colnames(count_mat) <- meta$sample_id
storage.mode(count_mat) <- "integer"

write.csv(data.frame(gene_id = rownames(count_mat), count_mat,
                     check.names = FALSE),
          "data/GSE103001_counts_matrix.csv", row.names = FALSE)

rownames(meta) <- meta$sample_id
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ patient_id + group
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("group", "Tumor", "Normal"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, "results/GSE103001/differential_expression.csv",
          row.names = FALSE)

# Highlight the apoptosis-related genes of interest where present.
target_genes <- c("BAX", "BCL2", "CASP3", "TP53")
write.csv(res_df[res_df$gene_id %in% target_genes, ],
          "results/GSE103001/apoptosis_gene_results.csv", row.names = FALSE)

volcano_df <- res_df
volcano_df$significance <- ifelse(
  !is.na(volcano_df$padj) & volcano_df$padj < 0.05 &
    abs(volcano_df$log2FoldChange) >= 1,
  "Significant", "Not significant"
)
p <- ggplot(volcano_df, aes(log2FoldChange,
                            -log10(pmax(padj, 1e-300)),
                            color = significance)) +
  geom_point(alpha = 0.65, size = 1.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(title = "GSE103001: Tumour vs adjacent normal",
       x = "log2 fold change (Tumour / Normal)",
       y = "-log10 adjusted p-value", color = NULL)
ggsave("results/GSE103001/volcano_plot.png", p,
       width = 7, height = 5, dpi = 300)

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
top <- head(order(rowVars(mat), decreasing = TRUE), min(50, nrow(mat)))
ann <- as.data.frame(colData(dds)[, c("group"), drop = FALSE])
png("results/GSE103001/top_variable_genes_heatmap.png",
    width = 1800, height = 1600, res = 220)
pheatmap(mat[top, , drop = FALSE], scale = "row",
         annotation_col = ann, show_rownames = FALSE,
         main = "GSE103001: Top variable genes")
dev.off()

message("Finished. See results/GSE103001/.")
