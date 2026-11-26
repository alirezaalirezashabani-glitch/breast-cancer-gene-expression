#!/usr/bin/env Rscript
# GSE103001 paired RNA-seq analysis and publication-quality plots.
# Run at repository root: Rscript scripts/make_GSE103001_figures.R
suppressPackageStartupMessages({
  library(GEOquery); library(DESeq2); library(ggplot2)
  library(pheatmap); library(matrixStats)
  library(AnnotationDbi); library(org.Hs.eg.db)
})
dir.create("data/GSE103001/extracted", recursive=TRUE, showWarnings=FALSE)
dir.create("results/GSE103001", recursive=TRUE, showWarnings=FALSE)
meta <- read.csv("data/GSE103001_sample_metadata.csv", stringsAsFactors=FALSE)
meta$group <- factor(meta$group, levels=c("Normal","Tumor"))
meta$patient_id <- factor(meta$patient_id)
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
  if(length(z)) z else NA_character_
}, character(1))
ok <- !is.na(ids) & ids %in% meta$sample_id
ff <- ff[ok]; ids <- ids[ok]
if(anyDuplicated(ids)) stop("Duplicate GSM count files detected.")
files <- setNames(ff, ids)
miss <- setdiff(meta$sample_id, names(files))
if(length(miss)) stop("Missing sample files: ", paste(miss, collapse=", "))
read_htseq <- function(f) {
  con <- if(grepl("\\.gz$", f)) gzfile(f, "rt") else f
  x <- read.delim(con, header=FALSE, stringsAsFactors=FALSE)
  if(inherits(con,"connection")) close(con)
  x <- x[!grepl("^__", x[[1]]), 1:2]
  names(x) <- c("gene_id","count"); x
}
tabs <- lapply(meta$sample_id, function(id) read_htseq(files[[id]]))
names(tabs) <- meta$sample_id
genes <- Reduce(intersect, lapply(tabs, function(x) x$gene_id))
cm <- do.call(cbind, lapply(tabs, function(x)
  as.numeric(x$count[match(genes,x$gene_id)])))
rownames(cm) <- genes; colnames(cm) <- meta$sample_id
storage.mode(cm) <- "integer"
write.csv(data.frame(gene_id=rownames(cm),cm,check.names=FALSE),
          "data/GSE103001_counts_matrix.csv",row.names=FALSE)
rownames(meta) <- meta$sample_id
dds <- DESeqDataSetFromMatrix(cm, meta, design=~patient_id+group)
dds <- dds[rowSums(counts(dds))>=10,]; dds <- DESeq(dds)
rr <- as.data.frame(results(dds,contrast=c("group","Tumor","Normal")))
rr$ensembl_id <- sub("\\.[0-9]+$","",rownames(rr))
mp <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(rr$ensembl_id),
  keytype="ENSEMBL", columns="SYMBOL")
mp <- mp[!duplicated(mp$ENSEMBL),]
rr$gene_symbol <- mp$SYMBOL[match(rr$ensembl_id,mp$ENSEMBL)]
rr <- rr[order(rr$padj),]
write.csv(rr,"results/GSE103001/differential_expression.csv",row.names=FALSE)
targets <- c("BAX","BCL2","CASP3","TP53")
tr <- rr[rr$gene_symbol %in% targets,]
write.csv(tr,"results/GSE103001/apoptosis_gene_results.csv",row.names=FALSE)
vd <- rr
vd$call <- ifelse(!is.na(vd$padj)&vd$padj<0.05&abs(vd$log2FoldChange)>=1,
                  "Significant","Not significant")
p1 <- ggplot(vd,aes(log2FoldChange,-log10(pmax(padj,1e-300)),color=call))+
 geom_point(alpha=.65,size=1.25)+geom_vline(xintercept=c(-1,1),linetype=2)+
 geom_hline(yintercept=-log10(.05),linetype=2)+theme_classic(base_size=12)+
 labs(title="GSE103001 | Tumour vs adjacent normal",
 x="log2 fold change (Tumour / Normal)",y="-log10 adjusted p-value",color=NULL)
ggsave("results/GSE103001/volcano_plot.png",p1,width=7.2,height=5.4,dpi=300)
vs <- vst(dds,blind=FALSE); vm <- assay(vs)
ti <- head(order(rowVars(vm),decreasing=TRUE),min(50,nrow(vm)))
png("results/GSE103001/top_variable_genes_heatmap.png",2000,1700,res=220)
pheatmap(vm[ti,,drop=FALSE],scale="row",
 annotation_col=data.frame(group=meta$group,row.names=meta$sample_id),
 show_rownames=FALSE,main="GSE103001 | Top variable genes")
dev.off()
# Plot normalized expression of the four apoptosis genes, if mapped.
valid <- targets %in% rr$gene_symbol
if(any(valid)) {
  sel <- rr$ensembl_id[match(targets[valid],rr$gene_symbol)]
  nc <- counts(dds,normalized=TRUE)
  dat <- do.call(rbind,lapply(seq_along(sel),function(i)
    data.frame(gene=targets[valid][i],sample=colnames(nc),
      count=as.numeric(nc[sel[i],]),
      group=meta$group[match(colnames(nc),meta$sample_id)])))
  p2 <- ggplot(dat,aes(group,count,fill=group))+
    geom_boxplot(outlier.shape=NA,alpha=.5)+
    geom_point(position=position_jitter(width=.12),size=1.2)+
    facet_wrap(~gene,scales="free_y")+scale_y_log10()+
    theme_classic(base_size=12)+
    labs(title="Apoptosis-related genes",x=NULL,
         y="Normalized counts (log10 scale)",fill=NULL)
  ggsave("results/GSE103001/apoptosis_gene_expression.png",p2,
         width=8,height=5.5,dpi=300)
} else warning("No target symbols mapped; inspect Ensembl annotation.")
writeLines(capture.output(sessionInfo()),"results/GSE103001/sessionInfo.txt")
message("Done: see results/GSE103001/")
