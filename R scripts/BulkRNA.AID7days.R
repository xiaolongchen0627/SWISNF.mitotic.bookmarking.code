library(DESeq2)
library(data.table)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

counts1=fread('ROBERTS-269902-STRANDED_RSEM_gene_count.2022-06-30_23-11-17.txt')
counts2=fread('ROBERTS-271135-STRANDED_RSEM_gene_count.2022-07-06_08-33-17.txt')
counts=merge(counts1,counts2,by=colnames(counts1)[1:4])

coldata=data.table(data.frame(matrix(c('aid6_0min','2431425',
                                       'aid6_30min','2431426',
                                       'aid6_60min','2431427',
                                       'aid6_240min','2431428',
                                       'aid6_480min','2431429',
                                       'aid6_1day','2431430',
                                       'aid6_3day','2431431',
                                       'aid23_0min','2431432',
                                       'aid23_30min','2431433',
                                       'aid23_60min','2431434',
                                       'aid23_240min','2431435',
                                       'aid23_480min','2431436',
                                       'aid23_1day','2431437',
                                       'aid23_3day','2431438',
                                       'aid6_7day','2435738',
                                       'aid23_7day','2435739'),ncol = 2,byrow = T)))
colnames(coldata)=c('V1','V2')
coldata[,rep:=paste(V1,'_rep',1:.N,sep=''),by='V1']

colnames(counts)[5:ncol(counts)]=gsub(colnames(counts)[5:ncol(counts)],pa='_.*',rep='')
colnames(counts)[5:ncol(counts)]=unlist(lapply(colnames(counts)[5:ncol(counts)],function(x)coldata$V1[match(x,coldata$V2)]))

tmp=table(counts$geneSymbol)
counts=counts[counts$geneSymbol%in%names(tmp[tmp==1]),]
counts.cts=as.matrix(counts[,5:ncol(counts)])
counts.cts=counts.cts[complete.cases(counts.cts),]
counts.cts=round(counts.cts,0)
counts.cts=data.frame(counts.cts)
counts.cts=apply(counts.cts,2,as.numeric)

rownames(counts.cts)=counts$geneSymbol
counts.cts=counts.cts[,c(1:7,15,8:14,16)]
coldata.iaa.timecourse=data.frame(name=colnames(counts.cts),
                                  cells=rep(c('aid6','aid23'),each=8),
                                  time=rep(c('0m','30m','60','240m','480m','1d','3d','7d'),2))

dds.iaa.tc=DESeqDataSetFromMatrix(countData = counts.cts,colData=coldata.iaa.timecourse,design=formula(~ time))
dds.iaa.tc=DESeq(dds.iaa.tc,test='LRT',reduced=~1)
rlog.iaatc=assay(vst(dds.iaa.tc))
res.iaa.tc=results(dds.iaa.tc)
res.iaa.tc=res.iaa.tc[complete.cases(res.iaa.tc),]
res.iaa.tc=data.table(data.frame(gene_name=rownames(res.iaa.tc),data.frame(res.iaa.tc)))

de.iaa.tc = res.iaa.tc[complete.cases(res.iaa.tc)&pvalue<0.05&abs(log2FoldChange)>1,]
km=kmeans(t(apply(rlog.iaatc[de.iaa.tc$gene_name,],1,scale)),centers = 3)

pdf('phetmap.iaa.bulkRNAseq.kmeans.3.pdf',5,height=3)
pheatmap(rlog.iaatc[de.iaa.tc$gene_name,][order(km$cluster),c(rbind(1:8,9:16))],cluster_rows = F,cluster_cols = F,
         scale = 'row',show_rownames = F,show_colnames = T,legend =T)
dev.off()





