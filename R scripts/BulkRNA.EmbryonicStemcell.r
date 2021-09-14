library(data.table)
library(DESeq2)
library(sva)
library(data.table)

batch1 <- fread('ROBERTS-208124-STRANDED_RSEM_count.txt')
batch2 <- fread('ROBERTS-214059-STRANDED_RSEM_count.txt')
batch2 <- batch2[order(geneID),]
batch1 <- batch1[order(geneID),]

cts <- merge(data.frame(batch1),data.frame(batch2),by=colnames(batch1)[1:4])
rownames(cts) <- cts[,2]
cts1 <- as.matrix(cts[,5:ncol(cts)])
cts2 <- cts1
cts2 <- apply(cts2,2,as.numeric)
cts2 <- cts2[complete.cases(cts2),]
cts2 <- round(cts2,0)
rownames(cts2) <- rownames(cts1)
colnames(cts1) <- gsub(colnames(cts1),pa='.*BAF',rep='BAF')
colnames(cts2) <- gsub(colnames(cts2),pa='.*BAF',rep='BAF')


coldata <- read.table('BulkRNA.ES.coldata',stringsAsFactors = F,header=F)
colnames(coldata) <- c('sample','condition','rep','cells','batch')
rownames(coldata) <- coldata[,1]

ddses.4cells <- DESeqDataSetFromMatrix(countData = cts2,colData=coldata,design=formula(~batch+condition))
ddses.4cells <- DESeq(ddses.4cells)
reses.4cells <- results(ddses.4cells)
reses.4cells <- reses.4cells[complete.cases((reses.4cells)),]
reses.4cells.pc <- reses.4cells[rownames(reses.4cells)%in%pc$symbol,]
rldes.4cells <- assay(vst(ddses.4cells))
rldes.4cells <- rldes.4cells[rowSums(rldes.4cells)>0,];
rldes.4cells <- rldes.4cells[rowSums(rldes.4cells[,1:12])>0,]

cts3 <- cts2[,colnames(cts2)%in%rownames(coldata[coldata$batch=='batch1',])]
###### De gene was calculated using the original two cells. The extra two cells was used to eleminate possible variation.
ddses.2cells <- DESeqDataSetFromMatrix(countData = cts3,colData=coldata[coldata$batch=='batch1',],design=formula(~condition))
ddses.2cells <- DESeq(ddses.2cells )
reses.2cells <- results(ddses.2cells)
reses.2cells <- reses.2cells[complete.cases((reses.2cells)),]
reses.2cells.pc <- reses.2cells[rownames(reses.2cells)%in%pc$symbol,]
rldes.2cells <- assay(vst(ddses.2cells))
rldes.2cells <- rldes.2cells[rowSums(rldes.2cells)>0,];
rldes.2cells <- rldes.2cells[rowSums(rldes.2cells[,1:12])>0,]
mean.2es.rlog <- rowGrpMeans(rldes.2cells[,c(1:12)],gr=factor(rep(c(1,2,3,4),3)));
colnames(mean.2es.rlog) <- c('R42A_A04','R42A_A10','MD_09','MD_30')

########## remove batch effect to do pca analysis

modddses <- model.matrix(~as.factor(condition),coldata)
mod0ddses <- model.matrix(~1,coldata)
batch <- c(rep(1,12),rep(2,12))

rldes.combat <- ComBat(rldes.4cells,batch = batch,mod = mod0ddses,par.prior = T,prior.plots = T)
rldes.combat <- rldes.combat[rowSums(rldes.combat)>0,]
sd1 <- apply(rldes.combat,1,sd) #remove genes with no changes
rldes.combat <- rldes.combat[rownames(rldes.combat)%in%names(sd1[sd1>0]),]

pca=prcomp(t(rldes.combat),scale. = T)
pca1=data.frame(pca$x,coldata)
pca1=data.table(pca1)
pca1$cellrep=NULL
pca1[cells=='A04',cellrep:='rep1']
pca1[cells=='A10',cellrep:='rep2']
pca1[cells=='A01',cellrep:='rep3']
pca1[cells=='A02',cellrep:='rep4']

pca1[cells=='09',cellrep:='rep1']
pca1[cells=='30',cellrep:='rep2']
pca1[cells=='45',cellrep:='rep3']
pca1[cells=='47',cellrep:='rep4']


gp <- ggplot(pca1,aes(x=PC1,y=PC2,col=condition,shape=cellrep))+geom_point()+
	geom_text_repel(aes(label=sample))+theme_bw()+scale_color_aaas()
pdf('BulkRNA.ES.pca.pdf',8,8,colormodel = 'cmyk')
plot(gp)
dev.off()

reses.2cells.pc <- data.frame(reses.pc)
reses.2cells.pc$pvalue[reses.pc$pvalue<1e-80] <- 1e-80
reses.2cells.pc$log2FoldChange[reses.pc$log2FoldChange>6] <- 6
reses.2cells.pc$log2FoldChange[reses.pc$log2FoldChange< -6] <- (-6)

gp <- gplot(data.frame(reses.pc),aes(x=-log2FoldChange,y=-log10(pvalue)))+
  geom_point(color='grey')+
  geom_point(data=reses.pc[reses.pc$log2FoldChange>1&reses.pc$pvalue   <0.05,],color='blue')+
  geom_point(data=reses.pc[reses.pc$log2FoldChange<(-1)&reses.pc$pvalue<0.05,],color='red')+
  geom_point(data=reses.pc['Bmp4',],color='black')+
  theme_bw()+geom_hline(yintercept = -log10(0.05))+  xlim(-6,6)+
  annotate('text',x=-5,y=60,label='MD down reg n = 704',col='blue')+
  annotate('text', x= 5 ,y=60 ,label='MD up reg n = 1167',col='red')

pdf('BulkRNA.ES.volcanoplot.pdf',8,8,colormodel = 'cmyk')
plot(gp)
dev.off()

###using two as DE detection and plot all the 4 
dees <- merge(data.frame(gene=rownames(reses.2cells.pc[reses.2cells.pc$pvalue<0.05,]),
	                     reses.2cells.pc[reses.2cells.pc$pvalue<0.05,]),
              data.frame(gene=rownames(rldes.2cells),rldes.2cells),
              by='gene')
dees <- dees[abs(dees$log2FoldChange)>1,];
dees <- dees[order(dees$log2FoldChange),]

pdf('BulkRNA.ES.heatmap.pdf',10,8,colormodel = 'cmyk');
pheatmap(dees[,8:31],cluster_rows = F,scale = 'row',show_rownames = F,color = colorRampPalette(c('blue','white','red'))(100));
dev.off()


MDdowngoes <- enrichGO(rownames(reses.pc[reses.pc$pvalue<0.05&reses.pc$log2FoldChange>1,]),OrgDb = org.Mm.eg.db,keyType = 'SYMBOL',ont = 'BP')
MDupgoes <- enrichGO(rownames(reses.pc[reses.pc$pvalue<0.05&reses.pc$log2FoldChange<(-1),]),OrgDb = org.Mm.eg.db,keyType = 'SYMBOL',ont = 'BP')

tmp <- MDupgoes[c(1,4,6,7,9,13,15,16,17)]
tmp <- tmp[order(tmp$p.adjust,decreasing = T),]
tmp$Description <- factor(tmp$Description,levels = tmp$Description)

pdf('BulkRNA.ES.go.bar.esmdup.selected.pdf',height=5,8,colormodel = 'cmyk')
ggplot(tmp,aes(y=Count,x=Description,fill=p.adjust))+geom_bar(stat='identity')+scale_fill_gradient(high='blue',low='red')+coord_flip()+theme(axis.text.x = element_text(angle=90))+theme_bw()
dev.off()

tmp <- MDdowngoes[c(1:8)]
tmp <- tmp[order(tmp$p.adjust,decreasing = T),]
tmp$Description <- factor(tmp$Description,levels = tmp$Description)

pdf('BulkRNA.ES.go.bar.esmddown.selected.pdf',height=5,8,colormodel = 'cmyk')
ggplot(tmp,aes(y=Count,x=Description,fill=p.adjust))+geom_bar(stat='identity')+scale_fill_gradient(high='blue',low='red')+coord_flip()+theme(axis.text.x = element_text(angle=90))+theme_bw()
dev.off()
