library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyr)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
smarce1.overlap.peaks <- readPeakFile('SMARCE1/overlap.peaks.bed')
smarce1.asyn_only.peaks <- readPeakFile('SMARCE1/Asyn_only.peaks.bed')
smarce1.mit_only.peaks <- readPeakFile('SMARCE1/Mit_only.peaks.bed')

# peakAnnoList <- lapply(list(Asyn_only.SMARCE1=smarce1.asyn_only.peaks,
# overlap.SMARCE1=smarce1.overlap.peaks,
# Mit_only.SMARCE1=smarce1.mit_only.peaks),annotatePeak,TxDb=txdb,tssRegion=c(-2000,2000),verbose=F)

smarcb1.overlap.peaks <- readPeakFile('SMARCB1/overlap.peaks.bed')
smarcb1.asyn_only.peaks <- readPeakFile('SMARCB1/Asyn_only.peaks.bed')
smarcb1.mit_only.peaks <- readPeakFile('SMARCB1/Mit_only.peaks.bed')

peakAnnoList.smarc <- lapply(list(Asyn_only.SMARCE1=smarce1.asyn_only.peaks,
                            overlap.SMARCE1=smarce1.overlap.peaks,
                            Mit_only.SMARCE1=smarce1.mit_only.peaks,
                            Asyn_only.SMARCB1=smarcb1.asyn_only.peaks,
                            overlap.SMARCB1=smarcb1.overlap.peaks,
                            Mit_only.SMARCB1=smarcb1.mit_only.peaks),
                            annotatePeak,
                            TxDb=txdb,tssRegion=c(-2000,2000),verbose=F
                            )

pdf('SMARCB1.SMARCE1.annotated.ChIPseeker.pdf',colormodel = 'cmyk',width=10,height = 8)
plotAnnoBar(peakAnnoList.smarc)
dev.off()

setwd('SOX2.peaks')

chipsox2 <- readPeakFile('Asyn.SOX2.ChIP.peaks.bed')
cnrsox2 <- readPeakFile('Asyn.SOX2.CnR.peaks.bed')
cnrsox2.mit <- readPeakFile('Mit.SOX2.CnR.peaks.bed')

peakAnnoList.sox2 <- lapply(list(SOX2.ChIP=chipsox2,
                            SOX2.CnR=cnrsox2,
                            SOX2.CnR.mit=cnrsox2.mit),
                            annotatePeak,
                            TxDb=txdb,tssRegion=c(-2000,2000),verbose=F
                            )
pdf('SMARCB1.SMARCE1.annotated.ChIPseeker.pdf',colormodel = 'cmyk',width=10,height = 8)
plotAnnoBar(peakAnnoList.sox2)
dev.off()                            
