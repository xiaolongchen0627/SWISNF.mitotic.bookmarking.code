library(DESeq2)
library(data.table)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
ns=fread('ROBERTS-222993-STRANDED_RSEM_count.txt')
bmp4.re1=fread('ROBERTS-260503-STRANDED_RSEM_gene_count.2022-04-08_22-39-33.txt')
bmp4.re2=fread('ROBERTS-260917-STRANDED_RSEM_gene_count.2022-04-14_19-31-35.txt')
bmp4.re=merge(bmp4.re1,bmp4.re2,by=colnames(bmp4.re1)[1:4])

tmp=table(ns$geneSymbol)
ns=ns[ns$geneSymbol%in%names(tmp[tmp==1]),]
ns.cts=as.matrix(ns[,5:ncol(ns)])
ns.cts=ns.cts[complete.cases(ns.cts),]
ns.cts=round(ns.cts,0)
ns.cts=data.frame(ns.cts)
ns.cts=apply(ns.cts,2,as.numeric)
rownames(ns.cts)=ns$geneSymbol
coldatans=data.frame(name=colnames(ns.cts),
                     condition=c('R42A','R42A','R42A','R42A','WT','WT','WT','WT'),
                     rep=c('rep1','rep2','rep1','rep2','rep1','rep2','rep1','rep2'),
                     cell=c('A04','A04','A10','A10','MD09','MD09','MD30','MD30'),
                     bmp4='ES')


tmp=table(bmp4.re$geneSymbol)
bmp4.re=bmp4.re[bmp4.re$geneSymbol%in%names(tmp[tmp==1]),]
bmp4.re.cts=as.matrix(bmp4.re[,5:ncol(bmp4.re)])
bmp4.re.cts=bmp4.re.cts[complete.cases(bmp4.re.cts),]
bmp4.re.cts=round(bmp4.re.cts,0)
bmp4.re.cts=data.frame(bmp4.re.cts)
bmp4.re.cts=apply(bmp4.re.cts,2,as.numeric)

coldatabmp4.re=data.frame(
  name=c('A04_0ng_rep1','A04_0.1ng_rep1','A04_0.2ng_rep1','A10_0ng_rep1','A10_0.1ng_rep1','A10_0.2ng_rep1',
         'MD09_0ng_rep1','MD09_0.1ng_rep1','MD09_0.2ng_rep1','MD30_0ng_rep1','MD30_0.1ng_rep1','MD30_0.2ng_rep1',
         'A04_0ng_rep2','A04_0.1ng_rep2','A04_0.2ng_rep2','A10_0ng_rep2','A10_0.1ng_rep2','A10_0.2ng_rep2',
         'MD09_0ng_rep2','MD09_0.1ng_rep2','MD09_0.2ng_rep2','MD30_0ng_rep2','MD30_0.1ng_rep2','MD30_0.2ng_rep2'),
  condition=c('R42A','R42A','R42A','R42A','R42A','R42A','WT','WT','WT','WT','WT','WT','R42A','R42A','R42A','R42A','R42A','R42A','WT','WT','WT','WT','WT','WT'),
  rep=c(rep('rep1',12),rep('rep2',12)),
  id=gsub(colnames(bmp4.re.cts),pa='X',rep=''),
  cell=c('A04','A04','A04','A10','A10','A10','MD09','MD09','MD09','MD30','MD30','MD30','A04','A04','A04','A10','A10','A10','MD09','MD09','MD09','MD30','MD30','MD30'),
  bmp4=c('0ng','0.1ng','0.2ng','0ng','0.1ng','0.2ng','0ng','0.1ng','0.2ng','0ng','0.1ng','0.2ng','0ng','0.1ng','0.2ng','0ng','0.1ng','0.2ng','0ng','0.1ng','0.2ng','0ng','0.1ng','0.2ng'))
colnames(bmp4.re.cts)=c('A04_0ng_rep1','A04_0.1ng_rep1','A04_0.2ng_rep1','A10_0ng_rep1','A10_0.1ng_rep1','A10_0.2ng_rep1',
                        'MD09_0ng_rep1','MD09_0.1ng_rep1','MD09_0.2ng_rep1','MD30_0ng_rep1','MD30_0.1ng_rep1','MD30_0.2ng_rep1',
                        'A04_0ng_rep2','A04_0.1ng_rep2','A04_0.2ng_rep2','A10_0ng_rep2','A10_0.1ng_rep2','A10_0.2ng_rep2',
                        'MD09_0ng_rep2','MD09_0.1ng_rep2','MD09_0.2ng_rep2','MD30_0ng_rep2','MD30_0.1ng_rep2','MD30_0.2ng_rep2')
rownames(bmp4.re.cts)=bmp4.re$geneSymbol

ddsns=DESeqDataSetFromMatrix(countData = ns.cts,colData=coldatans,design=formula(~condition))
ddsns=DESeq(ddsns)
resns=results(ddsns)
resns=resns[complete.cases(resns),]
# resultsNames(ddsns)
resns.pc=resns[rownames(resns)%in%pc$symbol,]
rld.ns=assay(rlog(ddsns))

ddsbmp4.re=DESeqDataSetFromMatrix(countData = bmp4.re.cts,colData=coldatabmp4.re,design=formula(~condition+bmp4+bmp4:condition))
ddsbmp4.re=DESeq(ddsbmp4.re)
resbmp4.re=results(ddsbmp4.re)
resbmp4.re=resbmp4.re[complete.cases(resbmp4.re),]
# resultsNames(ddsbmp4.re)
resbmp4.re.pc=resbmp4.re[rownames(resbmp4.re)%in%pc$symbol,]
rld.bmp4.re=assay(rlog(ddsbmp4.re))

bmp4.0.1=DESeqDataSetFromMatrix(countData = bmp4.re.cts[,c(7,8,10,11,19,20,22,23)],colData = coldatabmp4.re[c(7,8,10,11,19,20,22,23),],design=formula(~bmp4))
bmp4.0.1=DESeq(bmp4.0.1)
bmp4.0.1=results(bmp4.0.1)
bmp4.0.1=bmp4.0.1[complete.cases(bmp4.0.1),]

rld.bmp4.re1=rld.bmp4.re[rownames(rld.bmp4.re)%in%
                           intersect(
                             union(rownames(bmp4.0.1[bmp4.0.1$padj<0.01&abs(bmp4.0.1$log2FoldChange)>1,]),
                                   rownames(bmp4.0.2[bmp4.0.2$padj<0.01&abs(bmp4.0.2$log2FoldChange)>1,])),
                             rownames(bmp4.mdr42a[bmp4.mdr42a$padj<0.01&abs(bmp4.mdr42a$log2FoldChange)>1,]))
                         ,]
# rld.bmp4.re1=rld.bmp4.re1[apply(rld.bmp4.re1,1,function(x)length(x[x>0]))>3,]
pca=prcomp(t(rld.bmp4.re1[,c(1,4,7:13,16,19:24)]),scale. = T)
pca1=data.frame(pca$x,rbind(coldatabmp4.re[c(1,4,7:13,16,19:24),]))

ggplot(pca1,aes(x=PC1,y=PC2,col=condition,shape=bmp4))+geom_point()+geom_text_repel(aes(label=name))+theme_bw()+scale_color_aaas()
