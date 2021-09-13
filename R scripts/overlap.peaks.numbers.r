library(ggsci)
library(eulerr)
library(ggplot2)
library(gridExtra)

aaas.color_pallete <- pal_aaas(alpha=0.5)(10)

smarca4.overlap <- plot(euler(c(Mit=382,Asyn=14362,"Mit&Asyn"=1890)), 
    fill=aaas.color_pallete[1:2],
    legend=NULL,
    edges = list(lty = 1,lwd=0.5),
    quantities = list(type = c("counts", "percent"), 
    font=1, round=2, cex=0.8))

smarcb1.overlap <- plot(euler(c(Mit=833,Asyn=7635,"Mit&Asyn"=1466)),
    fill=aaas.color_pallete[1:2],
    legend=NULL,
    edges = list(lty = 1,lwd=0.5),
    quantities = list(type = c("counts", "percent"), 
    font=1, round=2, cex=0.8))

smarce1.overlap <- plot(euler(c(Mit=3176,Asyn=4864,"Mit&Asyn"=3501)),
    fill=aaas.color_pallete[1:2],
    legend=NULL,
    edges = list(lty = 1,lwd=0.5),
    quantities = list(type = c("counts", "percent"), 
    font=1, round=2, cex=0.8))

sox2.overlap <- plot(euler(c(Mit=1913,Asyn=11805,"Mit&Asyn"=3382)),
    fill=aaas.color_pallete[1:2],
    legend=NULL,
    edges = list(lty = 1,lwd=0.5),
    quantities = list(type = c("counts", "percent"), 
    font=1, round=2, cex=0.8))

pdf('CnR.Mit.Asny.overlap.pdf',colormodel = 'cmyk',width = 10,height = 8)
grid.arrange(smarca4.overlap,smarcb1.overlap,smarce1.overlap,nrow=1)
dev.off()


pdf('supple.SMARCA4.chip.overlap.pdf',colormodel='cmyk',width=6,height = 6)
smarca4.chip.overlap <- plot(euler(c(Mit=196,Asyn=18359,"Mit&Asyn"=339)),
    fill=aaas.color_pallete[1:2],
    legend=NULL,
    edges = list(lty = 1,lwd=0.5),
    quantities = list(type = c("counts", "percent"), 
    font=1, round=2, cex=0.8))
smarca4.chip.overlap
dev.off()

pdf('supple.sox2.chip.overlap.pdf',colormodel='cmyk',width=6,height = 6)
sox2.chip.overlap <- plot(euler(c(Mit=13,Asyn=32599,"Mit&Asyn"=383)),
    fill=aaas.color_pallete[1:2],
    legend=NULL,
    edges = list(lty = 1,lwd=0.5),
    quantities = list(type = c("counts", "percent"), 
    font=1, round=2, cex=0.8))
sox2.chip.overlap
dev.off()

pdf('supple.sox2.cnr.overlap.pdf',colormodel='cmyk',width=6,height = 6)
sox2.cnr.overlap <- plot(euler(c(Mit=1913,Asyn=11805,"Mit&Asyn"=3382)),
    fill=aaas.color_pallete[1:2],
    legend=NULL,
    edges = list(lty = 1,lwd=0.5),
    quantities = list(type = c("counts", "percent"), 
    font=1, round=2, cex=0.8))
sox2.cnr.overlap
dev.off()
