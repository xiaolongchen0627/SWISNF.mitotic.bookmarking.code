library(ggplot2)
library(data.table)
library(ggsci)
library(DEseq2)
library(ggrepel)

dir <- '/Volumes/Xiaolong/Nascent_RNA/STAR' ### where to put the counts file from HT-seq 
sf <- read.delim('/Volumes/Xiaolong/Nascent_RNA/STAR/NascentRNA.scalar.rev.txt' ,header = T, stringsAsFactor =F ) ### information about the samples, sampleNames and reads staistics

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sf,
                                       directory = dir,
                                       design= formula(~ group + time + group : time))

sizeFactors(ddsHTSeq) <- as.numeric(all$sf)
ddsTC <- DESeq(ddsHTSeq,test="LRT",reduced =  ~group + time)
ddsTC <- ddsHTSeq[rowSums(counts(ddsTC)) >= 10,]
rld <- vst(ddsTC)
rlog <- assay(rld)

mean.rlog.2types <- rowGrpMeans(rlog,gr=gl(12,4)); colnames(mean.rlog.2types) <- unique(paste(sapply(c('R42A','MD'),rep,4),sapply(c(0,45,90,180,240,'asyn'),rep,8),sep='.'))
mean.rlog.4cells <- rowGrpMeans(rlog,gr=gl(24,2)); colnames(mean.rlog.4cells) <- unique(paste(sapply(c('R42AA04','R42AA10','MD09','MD30'),rep,2),sapply(c(0,45,90,180,240,'asyn'),rep,8),sep='.'))

m.rlog.2types <- melt(mean.rlog.2types);
m.rlog.2types <- data.table(m.rlog.2types);
m.rlog.2types[,c('type','time'):=tstrsplit(Var2,'\\.')]
m.rlog.2types$time <- factor(m.rlog.2types$time,levels=c(0,45,90,180,240,'asyn'))
m.rlog.2types$type <- factor(m.rlog.2types$type,levels=c('R42A','MD'))

pdf('overall.average.rlog.pdf',8,8,colormodel = 'cmyk')
ggplot(m.rlog.2types,aes(x=time,y=value,fill=type))+geom_boxplot(outlier.shape = NA)+scale_fill_aaas()
dev.off()

#### pca analysis 

pca <- prcomp(t(rlog),scale. = T)
# autoplot(pca)
pca1 <- merge(data.frame(Name=rownames(pca$x),data.frame(pca$x)),sf[,c(1,8,9,10,11)],by='Name')
pca1$time <-factor(pca1$time,level=c('asyn','0','45','90','180','240'))
pdf('NascentRNA.pca.pdf',width=10,8)
ggplot(pca1[pca1$batch!='1',],aes(x=PC1,y=PC2,col=time, shape=group))+geom_point() +
  geom_text_repel(aes(label=gsub(gsub(Name,pa='_rep.*',rep=''),pa='^(\\d*)_',rep='')))+
  theme_bw()+scale_color_aaas()
dev.off()