library(data.table)
library(ggplot2)
library(ggsci)
library(edgeR)
library(limma)

setwd('bookmarkers.SMARCE1.MD')   
# there is predicted peaks and results of multibigwigsummary results for peaks
# as well as reads counts and normalized counts matrix to the peaks 


 dataname <- c('ESRRB','EZH2','SOX2')

 for(i in dataname){
 	tmp <- fread(paste('../',i,".3types.peaks.bed",sep=''));
 	assign(paste('ann',i,sep='.'),tmp)}
 
 for( i in dataname){
 	tmp <- fread(paste(i,'.fpkm',sep=''));
 	colnames(tmp) <- gsub(gsub(colnames(tmp),pa='\'',rep=''),pa='#',rep='');
 	tmp <- tmp[complete.cases(tmp),];
 	assign(i,tmp)
 }
#### fpkm to peaks were pre-calculated 
 for(i in dataname){
 	tmp1 <- get(paste('ann',i,sep='.'));
 	tmp2 <- get(i);
 	tmp3 <- merge(tmp1,tmp2,by.x=colnames(tmp1)[1:3],by.y=colnames(tmp2)[1:3]);
 	assign(paste('fpkm.',i,sep=''),tmp3)
 }


 for(i in dataname){
 	tmp <- get(paste('fpkm.',i,sep=''));
 	MDcol <- colnames(tmp)[grep(colnames(tmp),pa='MD')];
 	R42Acol <- colnames(tmp)[grep(colnames(tmp),pa='R42A')];
 	tmp$MDmean <- rowMeans(tmp[,..MDcol]);
 	tmp$R42Amean <- rowMeans(tmp[,..R42Acol]);
 	assign(paste('fpkm.',i,sep=''),tmp)
 }

#### all normalized counts and location in fpkm.ESRRB ,FPKM.SOX2, FPKM.EZH2

# for(i in dataname){
# 	tmp <- get(paste('fpkm.',i,sep=''));
# 	tmp$pvalue <- apply(tmp[,c(6:13)],1,function(x) t.test(x[1:4],x[5:8])$p.value);
# 	assign(paste('fpkm.',i,sep=''),tmp)
# 	}
##### t.test pvalue for normalized counts 

# for(i in dataname){
# 	tmp <- get(paste('fpkm.',i,sep=''));
# 	tmp[,rsq:=summary(lm(MDmean~R42Amean))$r.squared]
# 	assign(paste('fpkm.',i,sep=''),tmp)
# 	}
##### R squared 


dataname <- c('ESRRB','EZH2','SOX2')
sf <- fread('CnR.bookmark.SMARCE1MD.scale.factor');
sf <- data.frame(sf);
rownames(sf) <- sf[,1]

for( i in dataname){
	tmp <- fread(paste(i,'.counts',sep=''));
	colnames(tmp) <- gsub(gsub(colnames(tmp),pa='\'',rep=''),pa='#',rep='');
	tmp <- tmp[complete.cases(tmp),];
	assign(i,tmp)
	}
### counts for peaks was to calculate pvalue in edgeR
for(i in dataname){
	tmp1 <- get(paste('ann',i,sep='.'));
	tmp2 <- get(i);tmp3=merge(tmp1,tmp2,by.x=colnames(tmp1)[1:3],by.y=colnames(tmp2)[1:3]);
	assign(paste('m.',i,sep=''),tmp3)
	}

for(i in dataname){
	tmp <- get(paste('m.',i,sep=''));
	MDcol <- colnames(tmp)[grep(colnames(tmp),pa='MD')];
	R42Acol <- colnames(tmp)[grep(colnames(tmp),pa='R42A')];
	tmp$MDmean <- rowMeans(tmp[,..MDcol]);
	tmp$R42Amean <- rowMeans(tmp[,..R42Acol]);
	assign(paste('m.',i,sep=''),tmp)
}


for(i in dataname){
	tmp <- get(paste('m.',i,sep=''));
	tmp$pvalue <- apply(tmp[,c(6:13)],1,function(x) t.test(x[1:4],x[5:8])$p.value);
	assign(paste('m.',i,sep=''),tmp)
}

for(i in dataname){
	tmp <- get(paste('m.',i,sep=''));
	tmp[,rsq:=summary(lm(MDmean~R42Amean))$r.squared]
	assign(paste('m.',i,sep=''),tmp)
}

coldata <- factor(rep(c('MD','R42A'),each=4));
mod <- model.matrix(~0+coldata);
contrasts <- makeContrasts(coldataR42A-coldataMD,levels=mod)

for(i in dataname){
	tmp <- get(paste('m.',i,sep=''));
	tmp <- data.frame(tmp);rownames(tmp)<- apply(tmp[,1:3],1,paste,collapse='',sep='');
	edger.obj <- DGEList(as.matrix(tmp[,6:13]),group=rep(c('MD','R42A'),each=4))
	edger.obj$samples$lib.size <- sf[rownames(edger.obj$samples),2];
	edger.obj$samples$norm.factors <- sf[rownames(edger.obj$samples),9]
	#edger.obj=calcNormFactors(edger.obj)
	design <- model.matrix(~0+group, data=edger.obj$samples) ; 
	contrasts <- makeContrasts(groupR42A-groupMD,levels=design);
	edger.obj <- estimateGLMCommonDisp(edger.obj,design);
	edger.obj <- estimateGLMTrendedDisp(edger.obj, design); 
	edger.obj <- estimateGLMTagwiseDisp(edger.obj, design);
	et <- exactTest(edger.obj,pair=c('MD','R42A'));
	et$table$BH <- p.adjust(et$table$PValue,method='BH')
	fit <- glmQLFit(edger.obj,design);
	fit <- glmQLFTest(fit,contrast=contrasts);
	tt <- topTags(fit,sort.by='logFC',n=Inf)  
#logcpm=data.frame(cpm(edger.obj))
#MDcol=colnames(logcpm)[grep(colnames(logcpm),pa='MD')];
#R42Acol=colnames(logcpm)[grep(colnames(logcpm),pa='R42A')];
#logcpm$MDmean=rowMeans(logcpm[,colnames(logcpm)%in%MDcol]);
#logcpm$R42Amean=rowMeans(logcpm[,colnames(logcpm)%in%R42Acol]);
#edger.obj1=as.matrix(edger.obj[,6:13]);rownames(edger.obj1)=rownames(edger.obj)
#fit=topTable(eBayes(contrasts.fit(lmFit(edger.obj[,6:13],mod),makeContrasts(coldataR42A-coldataMD,levels=c('coldataMD','coldataR42A')))),n=Inf)
#edger.obj=cbind(edger.obj,fit[rownames(edger.obj),]);
#assign(paste('d.',i,sep=''),edger.obj)
assign(paste('tt.',i,sep=''),tt)
}

for(i in dataname){
  fpkm <- data.frame(get(paste('fpkm.',i,sep='')));rownames(fpkm)=apply(fpkm[,1:3],1,paste,collapse='',sep='');test=get(paste('tt.',i,sep=''))
  fpkm <- cbind(fpkm,test$table[rownames(fpkm),])
  #write.table(paste('Add.edgeR.all.CnR.',i,'.figure5.txt',sep=''),x=fpkm,col.names=T,row.names=F,quote=F,sep='\t')
  fpkm$p.adjust <- p.adjust(fpkm$pvalue,method='BH');
  assign(paste('all.',i,sep=''),fpkm)
  fc <- ifelse(i=='ESRRB',yes=2,no=1)   #fold change cut off 
  gp <- ggplot(fpkm,aes(y=log2(MDmean),x=log2(R42Amean)))+
    geom_point(data=fpkm[!(fpkm$FDR<0.01&abs(fpkm$logFC) >= fc &fpkm$pvalue<0.01),],aes(y=log2(MDmean),x=log2(R42Amean)),color='grey',shape=16,stroke=0,alpha=0.3,size=2)+
    geom_smooth(data=fpkm,method='lm')+
    geom_point(data=fpkm[fpkm$FDR<0.01&fpkm$logFC <= -1*fc &fpkm$pvalue<0.01,],aes(y=log2(MDmean),x=log2(R42Amean)),color='blue',shape=16,stroke=0,size=2,alpha=0.5)+
    geom_point(data=fpkm[fpkm$FDR<0.01&fpkm$logFC >= fc&fpkm$pvalue<0.01,],aes(y=log2(MDmean),x=log2(R42Amean)),color='red',shape=16,stroke=0,size=2,alpha=0.5)+
    geom_abline(intercept = 0, slope = 1,linetype= 'dashed',colour='black')+
    scale_color_aaas()+
    theme_bw()+theme(legend.position = "none")+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
    	panel.border = element_rect(colour = "black", fill=NA, size=2))+
    #+coord_trans(x="log2",y="log2")+
    annotate('text',y=quantile(log2(fpkm$MDmean),0.95),x=quantile(log2(fpkm$R42Amean),0.999),label=paste('R^2 = ', round(max(fpkm$rsq),2),sep=''),color='blue');
  
  pdf(paste('CnR.',i,'.scatterplot.pdf',sep=''),width=5,height=5,colormodel='cmyk')
  plot(gp)
  dev.off()

}
