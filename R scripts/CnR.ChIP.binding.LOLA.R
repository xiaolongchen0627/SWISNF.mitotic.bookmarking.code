library(LOLA)
background=readBed('ENCFF950SQZ.bed') #########ATAC as background 
database=loadRegionDB('all/') ############put public files as bed in the folder.
files=list.files('usrset/')
files=files[grep(files,pa='SWISNF')]
names=gsub(files,pa='.SWISNF.peak.bed',rep='')
for( i in files){
  name=gsub(i,pa='.SWISNF.peak.bed',rep='')
  tmp=readBed(paste('usrset/',i,sep=''))
  assign(name,tmp)
}
userSets=GRangesList(name.ESRRB.Asyn.only = ESRRB.Asyn.only,
                     name.ESRRB.Mit.retain = ESRRB.Mit.retain,
                     name.SMARCB1.Asyn.only = SMARCB1.Asyn.only,name.SMARCB1.Mit.only = SMARCB1.Mit.only,
                     name.SMARCB1.Mit.retain = SMARCB1.Mit.retain,name.SMARCB1.overlap = SMARCB1.overlap,
                     name.SMARCE1.Asyn.only = SMARCE1.Asyn.only,name.SMARCE1.Mit.only = SMARCE1.Mit.only,
                     name.SMARCE1.Mit.retain = SMARCE1.Mit.retain,name.SMARCE1.overlap = SMARCE1.overlap,
                     name.SOX2.Asyn.only = SOX2.Asyn.only,
                     name.SOX2.Mit.retain = SOX2.Mit.retain)

locResults.enrich = runLOLA(userSets, background, database, cores=1,direction = 'two.sided')

gp=ggplot(locResults.enrich,aes(x=userSet,y=anno,fill=oddsRatio))+
  geom_tile()+
  theme(axis.text.x = element_text(angle=90))+
  scale_fill_gradientn(colours = color,breaks=c(0,1,5),limit=c(0,5))+
  geom_text(aes(label=paste(ifelse(support<100,'#',' '),ifelse(qValue>0.05,'n.s',format(qValue, digits=2)),sep='')))
gp