#!/usr/bin/env Rscript
# Reproducible paired RNA-seq analysis for GSE103001.
# Run from repository root: Rscript scripts/make_GSE103001_figures.R

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

dir.create("data/GSE103001/extracted", recursive=TRUE, showWarnings=FALSE)
dir.create("results/GSE103001", recursive=TRUE, showWarnings=FALSE)

meta <- read.csv("data/GSE103001_sample_metadata.csv", stringsAsFactors=FALSE)
required <- c("sample_id", "patient_id", "group")
if (!all(required %in% names(meta))) stop("Metadata must contain: ", paste(required, collapse=", "))
if (anyDuplicated(meta$sample_id)) stop("Duplicate sample IDs in metadata.")
if (!all(meta$group %in% c("Normal", "Tumor"))) stop("Unexpected group labels.")
pair_check <- table(meta$patient_id, meta$group)
if (!all(c("Normal", "Tumor") %in% colnames(pair_check)) ||
    any(pair_check[, "Normal"] != 1L) || any(pair_check[, "Tumor"] != 1L)) {
  stop("Each patient must have exactly one Normal and one Tumor sample.")
}
meta$group <- factor(meta$group, levels=c("Normal", "Tumor"))
meta$patient_id_original <- meta$patient_id
# Safe factor levels prevent formula/column-name warnings from hyphenated patient IDs.
meta$patient_id <- factor(make.names(meta$patient_id))
rownames(meta) <- meta$sample_id

supp <- GEOquery::getGEOSuppFiles("GSE103001", makeDirectory=FALSE,
                                  baseDir="data/GSE103001")
arc <- rownames(supp)[grepl("GSE103001_RAW.tar$", rownames(supp))]
if (!length(arc)) arc <- list.files("data/GSE103001", "GSE103001_RAW.tar$",
                                    recursive=TRUE, full.names=TRUE)
if (!length(arc)) stop("GEO RAW.tar archive was not found.")
untar(arc[1], exdir="data/GSE103001/extracted")

ff <- list.files("data/GSE103001/extracted", "\\.(txt|txt.gz)$",
                 recursive=TRUE, full.names=TRUE)
ids <- vapply(ff, function(f) {
  z <- regmatches(basename(f), regexpr("GSM[0-9]+", basename(f)))
  if (length(z)) z else NA_character_
}, character(1))
ok <- !is.na(ids) & ids %in% meta$sample_id
ff <- ff[ok]
ids <- ids[ok]
if (anyDuplicated(ids)) stop("Duplicate GSM count files detected.")
files <- setNames(ff, ids)
missing_ids <- setdiff(meta$sample_id, names(files))
if (length(missing_ids)) stop("Missing sample files: ", paste(missing_ids, collapse=", "))

read_htseq <- function(f) {
  con <- if (grepl("\\.gz$", f)) gzfile(f, "rt") else f
  x <- read.delim(con, header=FALSE, stringsAsFactors=FALSE)
  if (inherits(con, "connection")) close(con)
  x <- x[!grepl("^__", x[[1]]), 1:2]
  names(x) <- c("gene_id", "count")
  x$count <- as.numeric(x$count)
  if (any(!is.finite(x$count)) || any(x$count < 0) ||
      any(abs(x$count - round(x$count)) > 1e-8)) {
    stop("HTSeq file contains invalid or non-integer raw counts: ", f)
  }
  x
}

tabs <- lapply(meta$sample_id, function(id) read_htseq(files[[id]]))
names(tabs) <- meta$sample_id
genes <- Reduce(intersect, lapply(tabs, function(x) x$gene_id))
if (!length(genes)) stop("No shared gene IDs found across sample count files.")
cm <- do.call(cbind, lapply(tabs, function(x)
  as.integer(x$count[match(genes, x$gene_id)])))
rownames(cm) <- genes
colnames(cm) <- meta$sample_id
storage.mode(cm) <- "integer"
write.csv(data.frame(gene_id=rownames(cm), cm, check.names=FALSE),
          "data/GSE103001_counts_matrix.csv", row.names=FALSE)

dds <- DESeqDataSetFromMatrix(countData=cm, colData=meta,
                               design=~ patient_id + group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
rr <- as.data.frame(results(dds, contrast=c("group", "Tumor", "Normal")))
rr$ensembl_id <- sub("\\.[0-9]+$", "", rownames(rr))

# Resolve one-to-many annotation explicitly instead of silently duplicating genes.
symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db, keys=unique(rr$ensembl_id),
  keytype="ENSEMBL", column="SYMBOL", multiVals="first"
)
rr$gene_symbol <- unname(symbol_map[rr$ensembl_id])
rr <- rr[order(rr$padj, na.last=TRUE), ]
write.csv(rr, "results/GSE103001/differential_expression.csv", row.names=FALSE)

targets <- c("BAX", "BCL2", "CASP3", "TP53")
target_results <- rr[!is.na(rr$gene_symbol) & rr$gene_symbol %in% targets, ]
write.csv(target_results, "results/GSE103001/apoptosis_gene_results.csv", row.names=FALSE)

# Volcano plot: exclude genes without valid test statistics from the plotting layer.
vd <- rr[is.finite(rr$log2FoldChange) & is.finite(rr$padj), ]
vd$call <- ifelse(vd$padj < 0.05 & abs(vd$log2FoldChange) >= 1,
                  "Significant", "Not significant")
p1 <- ggplot(vd, aes(log2FoldChange, -log10(pmax(padj, 1e-300)), color=call)) +
  geom_point(alpha=.65, size=1.25) +
  geom_vline(xintercept=c(-1, 1), linetype=2) +
  geom_hline(yintercept=-log10(.05), linetype=2) +
  theme_classic(base_size=12) +
  labs(title="GSE103001 | Tumour vs paired adjacent normal",
       x="log2 fold change (Tumour / Normal)",
       y="-log10 adjusted p-value", color=NULL)
ggsave("results/GSE103001/volcano_plot.png", p1, width=7.2, height=5.4, dpi=300)

vs <- vst(dds, blind=FALSE)
vm <- assay(vs)
ti <- head(order(rowVars(vm), decreasing=TRUE), min(50, nrow(vm)))
png("results/GSE103001/top_variable_genes_heatmap.png", 2000, 1700, res=220)
pheatmap(vm[ti, , drop=FALSE], scale="row",
         annotation_col=data.frame(group=meta$group, row.names=meta$sample_id),
         show_rownames=FALSE, main="GSE103001 | Top variable genes")
dev.off()

# PCA overview of sample-level variation.
pca_data <- plotPCA(vs, intgroup="group", returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
p3 <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
  geom_point(size=3, alpha=.85) +
  theme_classic(base_size=12) +
  labs(title="GSE103001 | PCA of variance-stabilized counts",
       x=paste0("PC1 (", percent_var[1], "%)"),
       y=paste0("PC2 (", percent_var[2], "%)"), color=NULL)
ggsave("results/GSE103001/pca_plot.png", p3, width=7.2, height=5.6, dpi=300)

# Paired normalized expression for apoptosis genes, retaining patient links.
valid_targets <- targets[targets %in% rr$gene_symbol]
if (length(valid_targets)) {
  target_ids <- rr$ensembl_id[match(valid_targets, rr$gene_symbol)]
  nc <- counts(dds, normalized=TRUE)
  nc_ids <- sub("\\.[0-9]+$", "", rownames(nc))
  row_idx <- match(target_ids, nc_ids)
  keep <- !is.na(row_idx)
  valid_targets <- valid_targets[keep]
  row_idx <- row_idx[keep]
  if (length(valid_targets)) {
    dat <- do.call(rbind, lapply(seq_along(valid_targets), function(i) {
      data.frame(gene=valid_targets[i], sample=colnames(nc),
                 patient_id=meta$patient_id[match(colnames(nc), meta$sample_id)],
                 count=as.numeric(nc[row_idx[i], ]),
                 group=meta$group[match(colnames(nc), meta$sample_id)])
    }))
    p2 <- ggplot(dat, aes(group, count)) +
      geom_line(aes(group=patient_id), color="grey65", alpha=.55, linewidth=.45) +
      geom_boxplot(aes(fill=group), outlier.shape=NA, alpha=.35, width=.55) +
      geom_point(aes(color=group), position=position_jitter(width=.08), size=1.5) +
      facet_wrap(~gene, scales="free_y") + scale_y_log10() +
      theme_classic(base_size=12) +
      labs(title="Apoptosis-related genes | paired patient samples",
           x=NULL, y="Normalized counts (log10 scale)", color=NULL, fill=NULL)
    ggsave("results/GSE103001/apoptosis_gene_expression.png", p2,
           width=8, height=5.5, dpi=300)
  }
} else {
  warning("No target symbols mapped; inspect Ensembl annotation.")
}

write.csv(data.frame(sample_id=meta$sample_id,
                     patient_id=meta$patient_id_original,
                     group=as.character(meta$group)),
          "results/GSE103001/analysis_sample_metadata.csv", row.names=FALSE)
writeLines(capture.output(sessionInfo()), "results/GSE103001/sessionInfo.txt")
message("Analysis completed successfully. Outputs are in results/GSE103001/")
