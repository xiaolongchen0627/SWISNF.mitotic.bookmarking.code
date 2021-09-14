library(DESeq2)
library(data.table)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

ns <- fread('ROBERTS-222993-STRANDED_RSEM_count.txt')
ns.cts <- as.matrix(ns[,5:ncol(ns)])
ns.cts <- ns.cts[complete.cases(ns.cts),]
ns.cts <- round(ns.cts,0)
ns.cts <- data.frame(ns.cts)
ns.cts <- apply(ns.cts,2,as.numeric)
rownames(ns.cts) <- ns$geneSymbol
coldatans <- data.frame(
	name=colnames(ns.cts),
    condition=c('Mut1','Mut1','Mut2','Mut2','WT1','WT1','WT2','WT2'),
    rep=c('rep1','rep2','rep1','rep2','rep1','rep2','rep1','rep2'),
    MD=c('R42A','R42A','R42A','R42A','WT','WT','WT','WT')
    )

rownames(coldatans) <- coldatans[,1]
ddsns <- DESeqDataSetFromMatrix(countData = ns.cts,colData=coldatans,design=formula(~MD))
ddsns <- DESeq(ddsns)
resns <- results(ddsns)
resns <- resns[complete.cases(resns),]

# resultsNames(ddsns)
pc <- fread('mm10.proteincoding.gene.bed') ### Protein_coding gene from gencode GRCm38 M25
# chr start end geneName ensembl_gene_ID strand

resns.pc <- resns[rownames(resns)%in%pc$symbol,]
rld.ns <- assay(rlog(ddsns))
dim(resns.pc[resns.pc$pvalue<0.05&resns.pc$log2FoldChange>1,]);dim(resns.pc[resns.pc$pvalue<0.05&resns.pc$log2FoldChange<(-1),]);

dens <- merge(data.frame(gene=rownames(resns.pc[resns.pc$pvalue<0.05,]),resns.pc[resns.pc$pvalue<0.05,]),data.frame(gene=rownames(rld.ns),rld.ns),by='gene')

dens <- dens[abs(dens$log2FoldChange)>1,];dens=dens[order(dens$log2FoldChange),]

pdf('BulkRNA.ns.heatmap.pdf',10,8,colormodel = 'cmyk'); 
pheatmap(dens[,8:15],cluster_rows = F,scale = 'row',show_rownames = F,
	color = colorRampPalette(c('blue','white','red'))(100));
dev.off()




resns.pc <- as.data.frame(resns.pc);
resns.pc$pvalue[resns.pc$pvalue<1e-100] <- 1e-100
resns.pc$log2FoldChange[resns.pc$log2FoldChange>6] <- 6
resns.pc$log2FoldChange[resns.pc$log2FoldChange< -6] <- (-6)

gp <- ggplot(data.frame(resns.pc),aes(x=log2FoldChange,y=-log10(pvalue)))+geom_point(color='grey')+
  geom_point(data=resns.pc[resns.pc$log2FoldChange>1&resns.pc$pvalue <0.05,],color='red')+
  geom_point(data=resns.pc[resns.pc$log2FoldChange<(-1)&resns.pc$pvalue<0.05,],color='blue')+
  theme_bw()+geom_hline(yintercept = -log10(0.05))+  xlim(-6,6)+
  annotate('text',x=5,y=60,label='MD high reg n = 1293',col='red')+
  annotate('text', x= -5 ,y=60 ,label='MD down reg n = 976',col='blue')

pdf('BulkRNA.ns.volcanoplot.pdf',8,8,colormodel = 'cmyk')
plot(gp)
dev.off()

MDupgons <- enrichGO(rownames(resns.pc[resns.pc$pvalue<0.05&resns.pc$log2FoldChange>1,]),OrgDb = org.Mm.eg.db,keyType = 'SYMBOL',ont = 'BP')
MDdowngons <- enrichGO(rownames(resns.pc[resns.pc$pvalue<0.05&resns.pc$log2FoldChange<(-1),]),OrgDb = org.Mm.eg.db,keyType = 'SYMBOL',ont = 'BP')

pdf('BulkRNA.ns.go.bar.mdup.pdf',height = 5,8,colormodel = 'cmyk')
barplot(MDupgons)
dev.off()

pdf('BulkRNA.ns.go.bar.mddown.pdf',height = 5,8,colormodel = 'cmyk')
barplot(MDdowngons)
dev.off()





