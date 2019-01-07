#PheWas=read.csv2('/Users/mrotival/WORK/06_papers/PapierSplicing_20180126/tables/new/Table_PheWAS_CYP3A5.csv',stringsAsFactor=F)
#
#PheWas$type=apply(PheWas[,c(6:10,12:21)]=='x',1,function(x){paste((colnames(PheWas)[c(6:10,12:21)])[x],collapse='//')})
#tab=sort(table(PheWas$type))
#PheWas$type[PheWas$type%in%names(tab)[tab<5]]=''
#PheWas$type[grep('Immune.disorders',PheWas$type)]='Immune.disorders'
#tab=sort(table(PheWas$type))
#
#PheWas$type[PheWas$type=='Heamatological.measurements//Immune.measurements']='Heamatological.measurements'
#
#by(as.numeric(PheWas$p.value),PheWas$type,min)
#
#library(RcolorBrewer)
#colors=rev(c(brewer.pal(8, 'Set1'),brewer.pal(8, 'Set2')))
#PheWas$color=as.numeric(as.factor(PheWas$type))
#PheWas$p.value=as.numeric(PheWas$p.value)
#
#oo=order(-PheWas$color,runif(length(PheWas$p.value)))
#
#oo=order(-PheWas$color,PheWas$p.value)
#PheWas$fdr=p.adjust(PheWas$p.value)
#plot(-log10(PheWas[oo,'p.value']),col=ifelse(PheWas$fdr[oo]<0.05,'black',colors[PheWas[oo,'color']]),pch=ifelse(PheWas$fdr[oo]<0.05,21,16),cex=ifelse(PheWas$fdr[oo]<0.05,0.6,0.4),bg=colors[PheWas[oo,'color']],axes=F,ylab=expression(-log[10](p)),xlab='')
#axis(2,las=1) 
#axis(1,at=c(0,cumsum(rev(table(PheWas$type)))),labels=FALSE) 
#abline(h=-log10(0.05/779),lty=2,col='grey')
PheWas=read.csv2(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/new/Table_PheWAS_CYP3A5.csv',HOME),stringsAsFactor=F)
PheWas$type=apply(PheWas[,c(23:33)]=='x',1,function(x){paste((gsub('.1','',colnames(PheWas)[c(23:33)]))[x],collapse='//')})
PheWas$type=as.factor(PheWas$type)
PheWas$type=relevel(PheWas$type,'Other')
levels(PheWas$type)=c("Other", "Anthropometric measurements", "Bone and joints disorders", 
"Cancer", "Cardiovascular disorders", "Digestive and metabolic disorders", 
"Heamatological measurements", "Immune disorders", "Lifestyle", 
"Nervous and psychiatric disorders", "Pulmonary disease")

colors=rev(c(brewer.pal(8, 'Set1'),brewer.pal(8, 'Set2')))
colors=colors[c(1:8,10,13,16)]
PheWas$color=as.numeric(as.factor(PheWas$type))
PheWas$p.value=as.numeric(PheWas$p.value)

oo=order(-PheWas$color,runif(length(PheWas$p.value)))

oo=order(-PheWas$color,PheWas$p.value)
PheWas$fdr=p.adjust(PheWas$p.value)
#4 x 7.75 inches
par(mar=c(5.1,4.1,0.5,2.1))
plot(-log10(PheWas[oo,'p.value']),col=ifelse(PheWas$fdr[oo]<0.05,'black',colors[PheWas[oo,'color']]),pch=ifelse(PheWas$fdr[oo]<0.05,ifelse(PheWas$Beta[oo]<0,24,25),16),cex=ifelse(PheWas$fdr[oo]<0.05,0.8,0.4),bg=colors[PheWas[oo,'color']],axes=F,ylab=expression(-log[10](p)),xlab='')
axis(2,las=1) 
axis(1,at=c(0,cumsum(rev(table(PheWas$type)))),labels=FALSE) 
abline(h=-log10(0.05/779),lty=2,col='grey')
par(xpd=T)
legend(-120,-1.5,ncol=2,bty='n',fill=rev(colors),legend=rev(levels(PheWas$type)),cex=0.7)
par(xpd=F)

