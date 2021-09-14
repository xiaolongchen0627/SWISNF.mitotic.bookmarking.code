library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(data.table)


setwd('BETA/')

beta <- fread('minus.ctcf.100000_targets_associated_peaks.txt') #results from BETA.sh 

dir <- '/Volumes/Xiaolong/Nascent_RNA/STAR' ### where to put the counts file from HT-seq 
sf <- read.delim('/Volumes/Xiaolong/Nascent_RNA/STAR/NascentRNA.scalar.txt' ,header = T, stringsAsFactor =F ) ### information about the samples, sampleNames and reads staistics


pc <- fread('mm10.proteincoding.gene.bed') ### Protein_coding gene from gencode GRCm38 M25
# chr start end geneName ensembl_gene_ID strand

time <- 90
#time <- 45 ####### same codes for 45 min
sf.min90 <- subset(sf,time == time)
counts.min90 <- DESeqDataSetFromHTSeqCount(sample=sf.min90,directory = dir , design = formula(~ group))
sizeFactors(counts.min90) <- sf.min90$sf
dds.min90 <- DESeq(counts.min90)
vst.min90 <- assay(vst(dds.min90))
dds.min90 <- results(dds.min90)

dds.min90.pc <- merge(data.frame(symbol=rownames(dds.min90),data.frame(dds.min90)),pc,by='symbol')
# print(dim(dds.min90[complete.cases(dds.min90) & dds.min90$pvalue<0.05 & dds.min90$log2FoldChange >=1 ,]))
# print(dim(dds.min90[complete.cases(dds.min90) & dds.min90$pvalue<0.05 & dds.min90$log2FoldChange <= -1 ,]))
# print(dim(dds.min90.pc[complete.cases(dds.min90.pc) & dds.min90.pc$pvalue<0.05 & dds.min90.pc$log2FoldChange >=1 ,]))
# print(dim(dds.min90.pc[complete.cases(dds.min90.pc) & dds.min90.pc$pvalue<0.05 & dds.min90.pc$log2FoldChange <= -1 ,]))
dds.min90.pc.refseq <- bitr(dds.min90.pc$symbol,fromType = 'SYMBOL',toType = 'REFSEQ',OrgDb = org.Mm.eg.db)
dds.min90.pc.refseq <- merge(dds.min90.pc.refseq[grep(dds.min90.pc.refseq$REFSEQ,pa='NM_'),],dds.min90.pc,by.x='SYMBOL',by.y='symbol')
dds.min90.pc.refseq <- dds.min90.pc.refseq[order(dds.min90.pc.refseq$SYMBOL,dds.min90.pc.refseq$REFSEQ),]

beta <- merge(beta,dds.min90.pc.refseq,by.x='refseqID',by.y='REFSEQ')
beta[,c('peaks','tmp'):=tstrsplit(peak,'_ID\\.')]
beta.up <- unique(beta[padj<0.01&log2FoldChange>1,])
beta.up <- beta.up[order(Gene_Symbol,peak,score,decreasing = T),];beta.up[,sum.score:=mean(score),by=list(peaks,Gene_Symbol)]
pdf('BETA.90min.potentialscore.pdf',10,10,colormode = 'cmyk')
ggplot(beta.up,aes(x=peaks,y=sum.score,fill=peaks))+ geom_boxplot(outlier.shape = NA)+theme_bw()+scale_color_aaas() #+facet_grid(.~peaks,scales = 'free_x')
dev.off()