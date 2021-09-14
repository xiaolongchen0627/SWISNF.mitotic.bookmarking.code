load('CnR.bookmarkers.Rdata') ##### saved image from CnR.boosmarkers.SMARCE1MD.cells.r
library(data.table)
library(ggplot2)
library(ggsci)
library(edgeR)
library(limma)


fpkm <- all.ESRRB;fc <- 2;
fpkm[,6:13] <- log2(fpkm[,6:13]);fpkm <- fpkm[order(fpkm$MDmean),];fpkm[,6:13][fpkm[,6:13]>=5] <- 5;fpkm[,6:13][fpkm[,6:13]<1] <- 0
de <- fpkm[fpkm$FDR<0.01&fpkm$logFC >= fc&fpkm$pvalue<0.01,]
print(dim(de));print(dim(de[de$MDmean<3,]))
#ph <- pheatmap(de[,6:9],cluster_cols=T,cluster_rows=T,show_rownames=F,show_colnames=T,silent = TRUE,kmeans=2)
#ct <- cutree(ph$tree_row,k=2)
#pdf('ESRRB.row.scaled.pdf',8,10);
#pheatmap(de[,6:13],cluster_cols=T,cluster_rows=F,show_rownames=F,show_colnames=T,scale='row')
#dev.off()
pdf('ESRRB.log2.fpkm.pdf',8,10);
pheatmap(de[,6:13],cluster_cols=T,cluster_rows=F,show_rownames=F,show_colnames=T)
dev.off()

fpkm <- all.EZH2;fc <- 1;
fpkm[,6:13] <- log2(fpkm[,6:13]);fpkm <- fpkm[order(fpkm$MDmean),];fpkm[,6:13][fpkm[,6:13]>=5] <- 5;fpkm[,6:13][fpkm[,6:13]<1] <- 0
de <- fpkm[fpkm$FDR<0.01&fpkm$logFC <= -1* fc&fpkm$pvalue<0.01,]
print(dim(de));print(dim(de[de$R42A<3,]))
#ph=pheatmap(de[,6:13],cluster_cols=T,cluster_rows=T,show_rownames=F,show_colnames=T,silent = TRUE,kmeans=2)
#ct=cutree(ph$tree_row,k=2)
#pdf('EZH2.row.scaled.pdf',8,10);
#pheatmap(de[,6:13],cluster_cols=T,cluster_rows=F,show_rownames=F,show_colnames=T,scale='row')
#dev.off()
pdf('EZH2.log2.fpkm1.pdf',8,10);
pheatmap(de[,6:13],cluster_cols=T,cluster_rows=F,show_rownames=F,show_colnames=T)
dev.off()

fpkm <- all.SOX2;fc <- 1;
fpkm[,6:13] <- log2(fpkm[,6:13]);fpkm <- fpkm[order(fpkm$MDmean),];fpkm[,6:13][fpkm[,6:13]>=5] <- 5;fpkm[,6:13][fpkm[,6:13]<1] <- 0
de <- fpkm[fpkm$FDR<0.01&fpkm$logFC >= fc&fpkm$pvalue<0.01,];
print(dim(de));print(dim(de[de$MDmean<3,]))
#ph=pheatmap(de[,6:9],cluster_cols=T,cluster_rows=T,show_rownames=F,show_colnames=T,silent = TRUE,kmeans=2)
#ct=cutree(ph$tree_row,k=2)
#pdf('SOX2.row.scaled.pdf',8,10);
#pheatmap(de[,6:13],cluster_cols=T,cluster_rows=F,show_rownames=F,show_colnames=T,scale='row')
#dev.off()
pdf('SOX2.log2.fpkm1.pdf',8,10);
pheatmap(de[,6:13],cluster_cols=T,cluster_rows=F,show_rownames=F,show_colnames=T)
dev.off()


