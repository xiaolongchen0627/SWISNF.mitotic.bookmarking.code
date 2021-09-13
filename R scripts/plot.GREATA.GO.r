library(ggExtra)
library(ggplot2)
library(data.table)

go <- fread('smarcb1.peaks.GO.GREAT.tsv')
##################################################
# summarized Gene ontology table was saved manually from GREAT results by adding annotation columns,
# go ; 'cc' 'bp' or 'mf'
# time ; 'asyn' , 'overlap' or 'mit'
##############
go <- go[go=='bp',1:15]

asyn <- go[time=='asyn',];
asyn <- asyn[order(asyn$`binom pvalue`),];
asyn$Term <- factor(asyn$Term,levels = asyn$Term)

mit <- go[time%in%c('mit','Mit'),];
mit <- mit[order(mit$`binom pvalue`),];
mit$Term <- factor(mit$Term,levels = mit$Term)

overlap <- go[time=='overlap',];
overlap <- overlap[order(overlap$`binom pvalue`),];
overlap$Term <- factor(overlap$Term,levels = overlap$Term)

gp1 <- ggplot(asyn,aes(y=`binom hist`,x=Term,fill=log10(`binom pvalue`)))+
  geom_bar(stat='identity')+
  scale_fill_gradient(high='blue',low='red')+
  coord_flip()+theme(axis.text.x = element_text(angle=90))+theme_bw()+
  facet_wrap(~time)

gp2 <- ggplot(mit,aes(y=`binom hist`,x=Term,fill=log10(`binom pvalue`)))+ 
geom_bar(stat='identity')+
  scale_fill_gradient(high='blue',low='red')+
  coord_flip()+theme(axis.text.x = element_text(angle=90))+theme_bw()+
  facet_wrap(~time)

gp3 <- ggplot(overlap,aes(y=`binom hist`,x=Term,fill=log10(`binom pvalue`)))+
geom_bar(stat='identity')+
  scale_fill_gradient(high='blue',low='red')+
  coord_flip()+theme(axis.text.x = element_text(angle=90))+theme_bw()+
  facet_wrap(~time)

pdf('SMARCB1.peaks.GREAT.GO.pdf',10,10, colormode='cmyk')
grid.arrange(gp1,gp2,gp3)
dev.off()
