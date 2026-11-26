#!/usr/bin/env Rscript

# Differential expression starter workflow for breast cancer RNA-seq counts.
# Input 1: data/counts_matrix.csv (genes in rows, sample IDs in columns;
#          first column must be gene_id; values must be raw integer counts)
# Input 2: data/sample_metadata_template.csv (sample_id, group)
# Replace the template metadata with the real sample IDs and group labels.
# Do not use normalized/TPM/FPKM values as input to DESeq2.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

counts_file <- "data/counts_matrix.csv"
metadata_file <- "data/sample_metadata_template.csv"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(counts_file)) {
  stop("Missing data/counts_matrix.csv. Add a raw-count matrix before running.")
}

counts_df <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, stringsAsFactors = TRUE)

required_cols <- c("sample_id", "group")
if (!all(required_cols %in% colnames(metadata))) {
  stop("Metadata must contain columns: sample_id and group.")
}
if (!all(metadata$sample_id %in% colnames(counts_df))) {
  stop("Every sample_id in metadata must match a column in the count matrix.")
}
if (length(unique(metadata$group)) != 2) {
  stop("This starter script expects exactly two groups.")
}

counts <- as.matrix(counts_df[, metadata$sample_id, drop = FALSE])
storage.mode(counts) <- "integer"
if (anyNA(counts) || any(counts < 0)) {
  stop("Counts must be non-negative integers without missing values.")
}

metadata$group <- relevel(factor(metadata$group), ref = levels(factor(metadata$group))[1])
rownames(metadata) <- metadata$sample_id

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ group
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(out_dir, "differential_expression_results.csv"), row.names = FALSE)

# Volcano plot
volcano_df <- res_df
volcano_df$significance <- ifelse(
  !is.na(volcano_df$padj) & volcano_df$padj < 0.05 &
    abs(volcano_df$log2FoldChange) >= 1,
  "Significant", "Not significant"
)
p_volcano <- ggplot(volcano_df, aes(log2FoldChange, -log10(pmax(padj, 1e-300)),
                                    color = significance)) +
  geom_point(alpha = 0.65, size = 1.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Differential expression", x = "log2 fold change",
       y = "-log10 adjusted p-value", color = NULL) +
  theme_classic()
ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width = 7, height = 5, dpi = 300)

# Variance-stabilized expression heatmap for the most variable genes
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
top_n <- min(50, nrow(mat))
top_genes <- names(sort(matrixStats::rowVars(mat), decreasing = TRUE))[seq_len(top_n)]
png(file.path(out_dir, "top_variable_genes_heatmap.png"), width = 1800, height = 1600, res = 220)
pheatmap(mat[top_genes, , drop = FALSE], scale = "row",
         annotation_col = as.data.frame(colData(dds)[, "group", drop = FALSE]),
         main = "Top variable genes")
dev.off()

message("Analysis complete. Results are in: ", out_dir)
