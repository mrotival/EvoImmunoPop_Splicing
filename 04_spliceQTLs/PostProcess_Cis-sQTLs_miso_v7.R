#######################################################
###				Aggregate psiQTL Pop Combined		###
#######################################################


# set the path to files that will be needed
Source='Ens70_HISAT'
allPSI_events=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep='')
Adjusted_PSI_values=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME)
# allGOterms

Require('Map_imputed')
Require('allGOterms')
# load data 
load(allPSI_events)
load(Adjusted_PSI_values)

GoToTest=list(immuneResponse='GO:0006955',
	InnateImmuneResponse='GO:0045087',
	AdaptiveImmuneResponse='GO:0002250',
	AntiviralResponse='GO:0051607',
	AntibacterialResponse='GO:0042742',
	TF_DNAbinding='GO:0003700',
	development='GO:0032502',
	olfactory_receptors='GO:0004984')
GoToTest_genes=lapply(GoToTest,function(x){allGOterms$gene[allGOterms$go==x]})
GoToTest_genes$all=unique(allGOterms$gene)

PSI_Annot$isImmune=PSI_Annot$gene_id%in%GoToTest_genes$immuneResponse

perm=0
thresholds=10^-c(seq(3,13),15,20,30,50)
RES=list()
for(pop in c('ALL')){
	for(cond in 1:5){
		for( CHR in 1:22){
		cat(pop, cond,CHR,'\n')
			load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/sQTL/MatrixEQTL-cis/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_KW_misoPSI_',Source,'_V7.Rdata',sep=''))
			RES0=RESCIS
			RES0$R2=RES0$statistic^2/(RES0$statistic^2+sum(SampleAnnot$cond==cond)-2)
			RES[[paste(cond,pop,CHR)]]=RES0
		}
	}
}
RES=do.call(rbind,RES)
colnames(RES)[colnames(RES)=='gene']='event_id'
RESobs=RES
mm=match(RESobs$event_id,PSI_Annot$event_id)
RESobs$gene=PSI_Annot$gene_id[mm]
RESobs$symbol=PSI_Annot$symbol[mm]
RESobs$event_type=PSI_Annot$event_type[mm]
RESobs$toKeep=RESobs$event_id%in% PSI_Annot$event_id[toKeep]
rm(RES)
gc()

TESTABLE=as.matrix(PSI_Annot[toKeep,paste('Testable',condIndex,sep='_')])
SUPPORT=as.matrix(PSI_Annot[toKeep,paste('Support',condIndex,sep='_')])
MEANPSI=as.matrix(PSI_Annot[toKeep,paste('MeanPSI',condIndex,sep='_')])
PCTNA=as.matrix(PSI_Annot[toKeep,paste('PctNA',condIndex,sep='_')])
JUNC=as.matrix(PSI_Annot[toKeep,paste('JuncCovered',condIndex,sep='_')])

PSI_Annot_cond=PSI_Annot[toKeep,]
testable_list=sapply(1:5,function(cond){unique(PSI_Annot_cond[TESTABLE[,cond],'gene_id'])})
luq(unlist(testable_list))
#[1] 4739
luq(setdiff(unlist(testable_list[-1]),testable_list[[1]]))
# 630 genes
GeneSymbols=PSI_Annot$symbol[toKeep]
MeanExpr=log2(1+GeneAnnot[,10:14])


lFC=MeanExpr[,-1]-MeanExpr[,1]%o%rep(1,4)
lFC_event=cbind(0,lFC[match(PSI_Annot_cond[,'gene_id'],rn(lFC)),])
TESTABLE_1=t(sapply(TESTABLE[,1],rep,5))
SUPPORT_1=t(sapply(SUPPORT[,1],rep,5))
MEANPSI_1=t(sapply(MEANPSI[,1],rep,5))
PCTNA_1=t(sapply(PCTNA[,1],rep,5))
JUNC_1=t(sapply(JUNC[,1],rep,5))

TESTABLE_DIFF=(PCTNA_1<0.05 & PCTNA<0.05) & (SUPPORT_1>10 & SUPPORT>10) & (TESTABLE | TESTABLE_1)
colnames(TESTABLE_DIFF)=paste('tested_DiffSplice',condIndex,sep='_')

TESTABLE_all=(apply(PCTNA<0.05 & SUPPORT>10, 1 ,all) & apply(TESTABLE,1,any)) %o%rep(1,5)==1

mm=match(RESobs$snps,Map_imputed$snp.name)
RESobs$pos=Map_imputed$position[mm]
RESobs$chrom=Map_imputed$chrom[mm]
RESobs$locus=Map_imputed$locus[mm]
RESobs$daf_AFB=Map_imputed$daf_char_AFB[mm]
RESobs$daf_EUB=Map_imputed$daf_char_EUB[mm]
RESobs$SNPfreq=Map_imputed$SNPfreq[mm]
mm=match(RESobs$event_id,PSI_Annot$event_id)
RESobs$event_start=as.numeric(PSI_Annot$start[mm])
RESobs$event_end=as.numeric(PSI_Annot$end[mm])
RESobs$CisDist=ifelse(RESobs$pos<RESobs$event_start,RESobs$pos-RESobs$event_start,ifelse(RESobs$pos<RESobs$event_end,0,RESobs$pos-RESobs$event_end))
RESobs$isTestable=PSI_Annot[toKeep,grep('Testable_',cn(PSI_Annot))][cbind(match(RESobs$event_id,PSI_Annot$event_id[toKeep]),RESobs$cond)]

RESobs$isExpressed=(PCTNA<0.05 & SUPPORT>10)[cbind(match(RESobs$event_id,PSI_Annot$event_id[toKeep]),RESobs$cond)]

RESobs$isTestable_DIFF=TESTABLE_DIFF[cbind(match(RESobs$event_id,PSI_Annot$event_id[toKeep]),RESobs$cond)]
RESobs$isTestable_all=TESTABLE_all[cbind(match(RESobs$event_id,PSI_Annot$event_id[toKeep]),RESobs$cond)]

RESobs=RESobs[which(!is.na(RESobs$pos)),]

max(abs(RESobs$CisDist))
mm=match(RESobs$event_id,PSI_Annot$event_id)
RESobs$CisDist=ifelse(PSI_Annot$strand[mm]=='+',1,-1)*RESobs$CisDist

NbAssocEvent=list()
NbAssocGene=list()

for (myDist in c(5e3,1e4,5e4,1e5,5e5,1e6)){
	NbAssocGene[[paste(as.character(myDist/1e3),'kb')]]=sapply(thresholds,function(th){luq(RESobs$gene[RESobs$pval<th & RESobs$toKeep & abs(RESobs$CisDist)<myDist & RESobs$isTestable])})
#	NbAssocEvent[[paste(as.character(myDist/1e3),'kb')]]=sapply(thresholds,function(th){luq(RESobs$event_id[RESobs$pval<th & RESobs$toKeep & abs(RESobs$CisDist)<myDist & RESobs$isTestable])})
}

NbAssocEventPerm=NbAssocGenePerm=list()
for(perm in 1:100){
cat(perm,'')
	RES=list()
	for( pop in c('ALL')){
		for(cond in 1:5){
			for( CHR in 1:22){
				cat(perm,pop, cond,CHR,'\n')
				load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/sQTL/MatrixEQTL-cis/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_KW_misoPSI_',Source,'_V7.Rdata',sep=''))
				RES0=RESCIS[which(RESCIS$pval<1e-3),]
				RES0$R2=RES0$statistic^2/(RES0$statistic^2+sum(SampleAnnot$cond==cond)-2)
				RES0$perm=perm
				RES[[paste(cond,pop,CHR)]]=RES0
			}
		}
	}
}
RES=do.call(rbind,RES)
colnames(RES)[colnames(RES)=='gene']='event_id'
mm=match(RES$event_id,PSI_Annot$event_id)
RES$gene=PSI_Annot$gene_id[mm]	
RES$symbol=PSI_Annot$symbol[mm]
RES$event_type=PSI_Annot$event_type[mm]
RES$toKeep=RES$event_id%in% PSI_Annot$event_id[toKeep]
mm=match(RES$snps,Map_imputed$snp.name)
RES$pos=Map_imputed$position[mm]
RES$chrom=Map_imputed$chrom[mm]
RES$locus=Map_imputed$locus[mm]
RES$daf_AFB=Map_imputed$daf_char_AFB[mm]
RES$daf_EUB=Map_imputed$daf_char_EUB[mm]
RES$SNPfreq=Map_imputed$SNPfreq[mm]
mm=match(RES$event_id,PSI_Annot$event_id)
RES$event_start=as.numeric(PSI_Annot$start[mm])
RES$event_end=as.numeric(PSI_Annot$end[mm])
RES$CisDist=ifelse(RES$pos<RES$event_start,RES$pos-RES$event_start,ifelse(RES$pos<RES$event_end,0,RES$pos-RES$event_end))
RES$CisDist=ifelse(PSI_Annot$strand[mm]=='+',1,-1)*RES$CisDist

RES$isTestable=PSI_Annot[toKeep,grep('Testable_',cn(PSI_Annot))][cbind(match(RES$event_id,PSI_Annot$event_id[toKeep]),RES$cond)]

RES=RES[which(!is.na(RES$pos)),]

for (myDist in c(5e3,1e4,5e4,1e5,5e5,1e6)){
	if(perm==1){
	NbAssocGenePerm[[paste(as.character(myDist/1e3),'kb')]]=list()
#	NbAssocEventPerm[[paste(as.character(myDist/1e3),'kb')]]=list()
	}
	NbAssocGenePerm[[paste(as.character(myDist/1e3),'kb')]][[perm]]=sapply(thresholds,function(th){luq(RES$gene[RES$pvalue<th & RES$toKeep & abs(RES$CisDist)<myDist & RES$isTestable])})
#	NbAssocEventPerm[[paste(as.character(myDist/1e3),'kb')]][[perm]]=sapply(thresholds,function(th){luq(RES$event_id[RES$pval<th & RES$toKeep & abs(RES$CisDist)<myDist & RES$isTestable])})
	}
gc()
}

NbFP_Gene=NbFP_Event=list()
for (myDist in c(5e3,1e4,5e4,1e5,5e5,1e6)){
	NbFP_Gene[[paste(as.character(myDist/1e3),'kb')]]=apply(do.call(cbind,NbAssocGenePerm[[paste(as.character(myDist/1e3),'kb')]]),1,mean)
	NbFP_Event[[paste(as.character(myDist/1e3),'kb')]]=apply(do.call(cbind,NbAssocEventPerm[[paste(as.character(myDist/1e3),'kb')]]),1,mean)
}

addcolprefix=function(data,prefix){
	colnames(data)=paste(prefix,colnames(data),sep='_')
	data
}

FDR_threshold=cbind(thresholds,addcolprefix(do.call(cbind,NbAssocEvent),'NbAssoc_Event'),addcolprefix(do.call(cbind,NbFP_Event)/do.call(cbind,NbAssocEvent),'FDR_Event'),addcolprefix(do.call(cbind,NbAssocGene),'NbAssoc_Gene'),addcolprefix(do.call(cbind,NbFP_Gene)/do.call(cbind,NbAssocGene),'FDR_Gene'))
save(RESobs,FDR_threshold,NbAssocGenePerm,NbAssocEventPerm,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V7.Rdata',EVO_IMMUNO_POP))
#load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5.Rdata',EVO_IMMUNO_POP))
#NbAssocGenePerm_1=sapply(thresholds,function(th){luq(RES$gene[ RES$pvalue<th])})
#NbAssocEventPerm_1=sapply(thresholds,function(th){luq(RES$event_id[RES$pval<th ])})

#NbFP_Gene/NbAssocGene
#
#NbFP_Event=apply(do.call(cbind,NbAssocEventPerm),1,mean)
#NbFP_Event/NbAssocEvent

##### 1e-7 is a good threshold for 5% False positive psiQTL
FDR_compute=function(pval,FDR_th=FDR_threshold[['FDR_Gene_1000 kb']]){y=approxfun(c(0,-log10(FDR_threshold$thresholds),500),c(pmin(-log10(c(1,FDR_th)),6),500))(-log10(pval)); 10^-y}
FDR_threshold=as.data.frame(FDR_threshold)
RESobs$FDR_1Mb=FDR_compute(RESobs$pvalue,FDR_threshold[['FDR_Gene_1000 kb']])

#RESobs$FDR_1Mb=FDR_compute(RESobs$pvalue,FDR_threshold[['FDR_Gene_1000 kb']])
#RESobs$FDR_5kb=FDR_compute(RESobs$pvalue,FDR_threshold[['FDR_Gene_5 kb']])

library(RColorBrewer)
library(DescTools)
acol= c(rev(brewer.pal(4,'YlOrRd')),brewer.pal(4,'YlGnBu'))
acol2=acol[c(1:3,6:8)]

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/other/Nb_detected_sQTLs_gene_MISO_V5.pdf',HOME))
thresholds=10^-c(seq(3,13,1),15,20,30,50)
thresholdsDetail=10^-c(seq(3,13,0.1),15,20,30,50)
extrapolate=function(x){approx(-log10(thresholds),x,-log10(thresholdsDetail))$y}
plot(-log10(thresholds),FDR_threshold[,'NbAssoc_Gene_1000 kb'],col='#00000000',axes=F,xlim=c(0,50),ylab='Nb Genes with sQTL')
for (i in 1:6){
	lines(-log10(thresholds),FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]],col=acol2[i],lwd=2)
	lines(-log10(thresholds),FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]]*FDR_threshold[,grep('FDR_Gene_',cn(FDR_threshold))[i]],col=acol2[i],lwd=2,lty=2)
}
axis(2,las=2)
axis(1,at=c(3,5,7,10,20,30,40,50),labels=c(3,5,7,10,20,30,40,50))
for (i in 1:6){
	xFDR10=min(which(extrapolate(FDR_threshold[,grep('FDR_Gene_',cn(FDR_threshold))[i]])<0.05))
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]])[xFDR10],col=acol2[i],cex=1.3,pch=16)	
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]])[xFDR10],cex=1.3,pch=1)	
}
legend(12,2500,lty=1,lwd=2,legend=c("1Mb","500kb","100kb","50kb","10kb","5kb"),col=acol2,bty='n')
legend(12,3500,lty=1:2,lwd=2,legend=c("observed","false positives (estimated)"),col='grey',bty='n')
legend(14,3150,pch=16,legend=c("detected at 5%FDR"),col='grey',bty='n')
dev.off()



###################################
###		Figure 2A	gene		###
###################################
par(mar=c(4,4,1,1))
thresholdsDetail=10^-c(seq(3,13,0.1),15,20,30,50)
extrapolate=function(x){approx(-log10(thresholds),x,-log10(thresholdsDetail))$y}
plot(-log10(thresholds),FDR_threshold[,'NbAssoc_Gene_1000 kb'],col='#00000000',axes=F,xlim=c(0,50),ylab='Number of genes with sQTL',xlab=expression(-log[10](P-value)))
for (i in 1:6){
	lines(-log10(thresholds),FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]],col=acol2[i],lwd=2)
	lines(-log10(thresholds),FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]]*(1-FDR_threshold[,grep('FDR_Gene_',cn(FDR_threshold))[i]]),col=acol2[i],lwd=2,lty=3)
}
for (i in 1:6){
	xFDR10=min(which(extrapolate(FDR_threshold[,grep('FDR_Gene_',cn(FDR_threshold))[i]])<0.05))
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]]*(1-FDR_threshold[,grep('FDR_Gene_',cn(FDR_threshold))[i]]))[xFDR10],col=acol2[i],cex=1.3,pch=16)	
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]]*(1-FDR_threshold[,grep('FDR_Gene_',cn(FDR_threshold))[i]]))[xFDR10],cex=1.3,pch=1)
}
for (i in 1:6){
	xFDR10=min(which(extrapolate(FDR_threshold[,grep('FDR_Gene_',cn(FDR_threshold))[i]])<0.05))
	cat(cn(FDR_threshold)[grep('NbAssoc_Gene',cn(FDR_threshold))[i]],extrapolate(FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]])[xFDR10],'\n')
	cat('truepos: ',cn(FDR_threshold)[grep('NbAssoc_Gene',cn(FDR_threshold))[i]],extrapolate(FDR_threshold[,grep('NbAssoc_Gene',cn(FDR_threshold))[i]]*(1-FDR_threshold[,grep('FDR_Gene_',cn(FDR_threshold))[i]]))[xFDR10],'\n')
	}

axis(2,las=2)
axis(1,at=c(3,5,7,10,20,30,40,50),labels=c(3,5,7,10,20,30,40,50))
legend(12,3000,lty=1,lwd=2,legend=c("1Mb","500kb","100kb"),col=acol2[6:4],bty='n')
legend(30,3000,lty=1,lwd=2,legend=c("50kb","10kb","5kb"),col=acol2[3:1],bty='n')

legend(12,4500,lty=c(1,3),lwd=2,legend=c("observed","true positives (estimated)"),col='grey',bty='n')
legend(15.4,3700,pch=16,legend=c("true positives at 5%FDR"),col='grey',bty='n')
###################################
###			Figure END			###
###################################

###################################
###		Figure 2A	event		###
###################################
par(mar=c(4,4,1,1))
thresholdsDetail=10^-c(seq(3,13,0.1),15,20,30,50)
extrapolate=function(x){approx(-log10(thresholds),x,-log10(thresholdsDetail))$y}
plot(-log10(thresholds),FDR_threshold[,'NbAssoc_Event_1000 kb'],col='#00000000',axes=F,xlim=c(0,50),ylim=c(0,4500),ylab='Number of splice QTL',xlab=expression(-log[10](P-value)))
for (i in 1:6){
	lines(-log10(thresholds),FDR_threshold[,grep('NbAssoc_Event',cn(FDR_threshold))[i]],col=acol2[i],lwd=2)
	lines(-log10(thresholds),FDR_threshold[,grep('NbAssoc_Event',cn(FDR_threshold))[i]]*(1-FDR_threshold[,grep('FDR_Event_',cn(FDR_threshold))[i]]),col=acol2[i],lwd=2,lty=3)
}
for (i in 1:6){
	xFDR10=min(which(extrapolate(FDR_threshold[,grep('FDR_Event_',cn(FDR_threshold))[i]])<0.05))
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(FDR_threshold[,grep('NbAssoc_Event',cn(FDR_threshold))[i]]*(1-FDR_threshold[,grep('FDR_Event_',cn(FDR_threshold))[i]]))[xFDR10],col=acol2[i],cex=1.3,pch=16)	
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(FDR_threshold[,grep('NbAssoc_Event',cn(FDR_threshold))[i]]*(1-FDR_threshold[,grep('FDR_Event_',cn(FDR_threshold))[i]]))[xFDR10],cex=1.3,pch=1)
}
for (i in 1:6){
	xFDR10=min(which(extrapolate(FDR_threshold[,grep('FDR_Event_',cn(FDR_threshold))[i]])<0.05))
	cat(cn(FDR_threshold)[grep('NbAssoc_Event',cn(FDR_threshold))[i]],extrapolate(FDR_threshold[,grep('NbAssoc_Event',cn(FDR_threshold))[i]])[xFDR10],'\n')
	cat('truepos: ',cn(FDR_threshold)[grep('NbAssoc_Event',cn(FDR_threshold))[i]],extrapolate(FDR_threshold[,grep('NbAssoc_Event',cn(FDR_threshold))[i]]*(1-FDR_threshold[,grep('FDR_Event_',cn(FDR_threshold))[i]]))[xFDR10],'\n')
	}

axis(2,las=2)
axis(1,at=c(3,5,7,10,20,30,40,50),labels=c(3,5,7,10,20,30,40,50))
legend(12,3000,lty=1,lwd=2,legend=c("1Mb","500kb","100kb"),col=acol2[6:4],bty='n')
legend(30,3000,lty=1,lwd=2,legend=c("50kb","10kb","5kb"),col=acol2[3:1],bty='n')

legend(12,4500,lty=c(1,3),lwd=2,legend=c("observed","true positives (estimated)"),col='grey',bty='n')
legend(15.4,3700,pch=16,legend=c("true positives at 5%FDR"),col='grey',bty='n')
dev.off()

###################################
###			Figure END			###
###################################

RESobs$FDR_1Mb=FDR_compute(RESobs$pvalue,FDR_threshold[['FDR_Gene_1000 kb']])
RESobs$FDR_5kb=FDR_compute(RESobs$pvalue,FDR_threshold[['FDR_Gene_5 kb']])
RESobs$FDR_1Mb_event=FDR_compute(RESobs$pvalue,FDR_threshold[['FDR_Event_1000 kb']])

luq(RESobs$gene[RESobs$FDR_5kb<0.05 & abs(RESobs$CisDist)<5e3 ]) # 1565
luq(RESobs$gene[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]) # 993
luq(RESobs$gene[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable_all]) # 689

luq(RESobs$event_id[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]) # 1464
luq(RESobs$event_id[RESobs$FDR_1Mb_event<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]) # 1549
luq(RESobs$event_id[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable_all]) # 1029
luq(RESobs$event_id[RESobs$FDR_1Mb_event<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable_all]) # 1089
luq(RESobs$gene[RESobs$FDR_1Mb_event<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]) # 1549
luq(RESobs$gene[RESobs$FDR_1Mb_event<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable_all]) # 1089


# enrichment of immune genes in sQTL
resGO=GOSeq(unique(RESobs$gene[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$toKeep & RESobs$isTestable]),unique(RESobs$gene[RESobs$toKeep & RESobs$isTestable])))
# enrichment of immune genes in sQTL when removing chr 6, remains significant.
resGO=GOSeq(unique(RESobs$gene[RESobs$FDR_1Mb<0.05 & substr(RESobs$locus,1,1) !=6 & abs(RESobs$CisDist)<1e6 & RESobs$toKeep & RESobs$isTestable]),unique(RESobs$gene[RESobs$toKeep & RESobs$isTestable & substr(RESobs$locus,1,1)!=6])))

PSI_Annot$has_sQTL=NA
PSI_Annot$has_sQTL[toKeep]=PSI_Annot$event_id[toKeep]%in%RESobs$event_id[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]
by(PSI_Annot$has_sQTL, PSI_Annot$isImmune,mean,na.rm=T)
#PSI_Annot$isImmune: FALSE
#[1] 0.0845399
#PSI_Annot$isImmune: TRUE
#[1] 0.145
odds.ratio(table(PSI_Annot$has_sQTL, PSI_Annot$isImmune))
#           LowerCI       OR  UpperCI alpha           P
# odds ratio 1.572096 1.836455 2.138528  0.05 6.72101e-14


RESobs_nodup_1Mb=RESobs[which(RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable),]
RESobs_nodup_1Mb=RESobs_nodup_1Mb[order(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$pval),]
RESobs_nodup_1Mb=RESobs_nodup_1Mb[!duplicated(RESobs_nodup_1Mb$event_id),]

RESobs_nodup_1Mb_cond=RESobs[which(RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable),]
RESobs_nodup_1Mb_cond=RESobs_nodup_1Mb_cond[order(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond,RESobs_nodup_1Mb_cond$pval),]
RESobs_nodup_1Mb_cond=RESobs_nodup_1Mb_cond[!duplicated(paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond)),]

###################################
###		Figure 2A	distance	###
###################################

#RESobs_nodup_1Mb_cond_expr=RESobs[which(RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isExpressed),]
#RESobs_nodup_1Mb_cond_expr=RESobs_nodup_1Mb_cond_expr[order(RESobs_nodup_1Mb_cond_expr$event_id,RESobs_nodup_1Mb_cond_expr$cond,RESobs_nodup_1Mb_cond_expr$pval),]
#RESobs_nodup_1Mb_cond_expr=RESobs_nodup_1Mb_cond_expr[!duplicated(paste(RESobs_nodup_1Mb_cond_expr$event_id,RESobs_nodup_1Mb_cond_expr$cond)),]

#plot(RESobs_nodup_1Mb$CisDist/1e6,-log10(RESobs_nodup_1Mb$pval),pch=16,col=gsub('AA','44',colERC5[1]),xlim=c(-1,1),las=1,axes=F,ylab=expression(-log[10](P-value)),xlab='Distance from Event boundaries');axis(2,las=1);axis(1,at=c(-1,-0.5,0,0.5,1),labels=c('1Mb','500kb','0','500kb','1Mb'))
#plot(RESobs_nodup_1Mb$CisDist/1e6,-log10(RESobs_nodup_1Mb$pval),pch=16,col=paste(colPSI[RESobs_nodup_1Mb$event_type],'AA',sep=''),xlim=c(-1,1),las=1,axes=F,ylab=expression(-log[10](P-value)),xlab='Distance from Event boundaries');axis(2,las=1);axis(1,at=c(-1,-0.5,0,0.5,1),labels=c('1Mb','500kb','0','500kb','1Mb'))
#legend('topright',pch=16,col=paste(colPSI,'AA',sep=''),legend=names(colPSI),bty='n')
#
#xx=c(seq(-1,-0.05,0.05),c(-0.005,0.005),seq(0.05,1,0.05))
#lines(0.5*(xx[-1]+xx[-length(xx)]),table(cut(RESobs_nodup_1Mb$CisDist/1e6,xx))/550*70;,col='red')

###################################
###		Figure 2A	END			###
###################################

###########################################################
###		Enrichment of sQTLs in immune genes 			###
###		remains when sQTLs from chr6 have been removed	###
###########################################################

resGO=GOSeq(unique(RESobs_nodup_1Mb$gene),unique(PSI_Annot$gene_id[toKeep]),overOnly=F)
plotGO(resGO[which(resGO$over_represented_pvalue<0.05),])
plotGO(resGO[which(resGO$under_represented_pvalue<0.05),],under=T)
plotGO=function(resGO,mar=18,under=F,...){
	splitname=function(x,sep=' ',nmax=40){y=strsplit(x,sep)
										y=sapply(y,function(z){
														countchar=cumsum(nchar(z));
														group=countchar%/%nmax;
														paste(By(z,group,paste,collapse=' '),collapse='\n')
														})
										x[nchar(x)>(nmax+10)]=y[nchar(x)>(nmax+10)]
										x
			}														
	params=par()
	par(mar=c(4,mar,1,1))
	if(!under){resGO=resGO[nrow(resGO):1,]}
	barplot(-log10(resGO$FDR), main="", horiz=TRUE, names.arg=splitname(resGO$Term),col=ifelse(resGO$ontology=='BP',colPSI[1],ifelse(resGO$ontology=='MF',colPSI[2],colPSI[3])),xlab=expression(-log[10](P[adj])),las=1,...)
	legend("bottomright",fill=colPSI[1:3],legend=c('BP','MF','CC'),bty='n')
	par(mar=params$mar)
	}
# enrichment of immune genes in sQTL when removing chr 6, remains significant.
resGO=GOSeq(unique(RESobs$gene[RESobs$FDR_1Mb<0.05 & substr(RESobs$locus,1,1) !=6 & abs(RESobs$CisDist)<1e6 & RESobs$toKeep & RESobs$isTestable]),unique(RESobs$gene[RESobs$toKeep & RESobs$isTestable & substr(RESobs$locus,1,1)!=6])))

###########################################################
###		Enrichment of sQTLs in immune genes 			###
###						DONE							###
###########################################################


#barplot(table(cut(RESobs_nodup_1Mb$CisDist,c(-1e6,-5e5,-1e5,-5e4,-1e4,-5e3,-1e3,1e3,5e3,1e4,5e4,1e5,5e5,1e6))),las=2)
#barplot(table(cut(RESobs_nodup_1Mb$CisDist[RESobs_nodup_1Mb$event_id%in%setdiff(RESobs_nodup_1Mb$event_id,RESobs_nodup_5kb$event_id)],c(-1e6,-5e5,-1e5,-5e4,-1e4,-5e3,-1e3,1e3,5e3,1e4,5e4,1e5,5e5,1e6))),las=2,col=colERC5[2],add=T)

RESobs_has_sQTL_1Mb=RESobs[which(abs(RESobs$CisDist)<1e6 & paste(RESobs$event_id,RESobs$snps)%in%paste(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$snps)),]
#Detected_sQTL_P5pct=By(RESobs_has_sQTL_1Mb$cond,paste(RESobs_has_sQTL_1Mb$event_id, RESobs_has_sQTL_1Mb$snps),function(x){paste(sort(x),collapse='')})
#RESobs_nodup_1Mb$detected_cond=Detected_sQTL_P5pct[match(paste(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$snps),names(Detected_sQTL_P5pct))]
#RESobs_nodup_1Mb$NbTestableCond=PSI_Annot$NbTestableCond[match(RESobs_nodup_1Mb$event_id,PSI_Annot$event_id)]

#RESobs_has_sQTL_1Mb=RESobs[which(abs(RESobs$CisDist)<1e6 & paste(RESobs$event_id,RESobs$snps)%in%paste(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$snps) & RESobs$isTestable_all),]
Pval_sQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='pvalue');
Pval_sQTL[is.na(Pval_sQTL)]=1
Pval_sQTL=Pval_sQTL[match(RESobs_nodup_1Mb$event_id,Pval_sQTL$event_id),-1]
colnames(Pval_sQTL)=paste('Pvalue',condIndex,sep='_')
Beta_sQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='beta')
Beta_sQTL[is.na(Beta_sQTL)]=0
Beta_sQTL=Beta_sQTL[match(RESobs_nodup_1Mb$event_id,Beta_sQTL$event_id),-1]
colnames(Beta_sQTL)=paste('Beta',condIndex,sep='_')

Testable_sQTL=TESTABLE[match(RESobs_nodup_1Mb$event_id,PSI_Annot[toKeep,'event_id']),]
Expressed_sQTL=(PCTNA<0.05 & SUPPORT>10)[match(RESobs_nodup_1Mb$event_id,PSI_Annot[toKeep,'event_id']),]
colnames(Expressed_sQTL)=paste('Expressed',condIndex,sep='_')
Testable_DIFF_sQTL=TESTABLE_DIFF[match(RESobs_nodup_1Mb$event_id,PSI_Annot[toKeep,'event_id']),-1]

RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,Pval_sQTL,Beta_sQTL,Testable_sQTL,Expressed_sQTL,Testable_DIFF_sQTL)
#RESobs_FDR_nodup[which(RESobs_FDR_nodup$detected_cond==234 & RESobs_FDR_nodup$event_id%in%Testable_all)[1:3],]
#save(RESobs_FDR_nodup,RESobs_hasSQTL,RESobs_FDR,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_V2_detail.Rdata',EVO_IMMUNO_POP))

mm=match(RESobs_nodup_1Mb$snps,Map_imputed$snp.name)
RESobs_nodup_1Mb$FST=Map_imputed$FST_adj[mm]
RESobs_nodup_1Mb$iHS_AFB=Map_imputed$iHS_AFB[mm]
RESobs_nodup_1Mb$iHS_EUB=Map_imputed$iHS_EUB[mm]
RESobs_nodup_1Mb$maf_AFB=Map_imputed$maf_AFB[mm]		
RESobs_nodup_1Mb$maf_EUB=Map_imputed$maf_EUB[mm]

mm=match(RESobs_nodup_1Mb_cond$snps,Map_imputed$snp.name)
RESobs_nodup_1Mb_cond$FST=Map_imputed$FST_adj[mm]
RESobs_nodup_1Mb_cond$iHS_AFB=Map_imputed$iHS_AFB[mm]
RESobs_nodup_1Mb_cond$iHS_EUB=Map_imputed$iHS_EUB[mm]
RESobs_nodup_1Mb_cond$maf_AFB=Map_imputed$maf_AFB[mm]
RESobs_nodup_1Mb_cond$maf_EUB=Map_imputed$maf_EUB[mm]

library(data.table)
MAP_GWAS=as.data.frame(fread(sprintf('gunzip -c %s/Annotation/GWAS/Map_GWAS_LD_allChr.txt.gz',HOME)))
mm=match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)
RESobs_nodup_1Mb$GWAS_Trait_R2_1E5=MAP_GWAS$GWAS_Trait_R2_1E5[mm]
RESobs_nodup_1Mb$GWAS_Trait_R2_1E8=MAP_GWAS$GWAS_Trait_R2_1E8[mm]
RESobs_nodup_1Mb$GWAS_EFO_R2_1E5=MAP_GWAS$GWAS_EFO_R2_1E5[mm]
RESobs_nodup_1Mb$GWAS_EFO_R2_1E8=MAP_GWAS$GWAS_EFO_R2_1E8[mm]

mm=match(RESobs_nodup_1Mb_cond$snps,MAP_GWAS$snp.name)
RESobs_nodup_1Mb_cond$GWAS_Trait_R2_1E5=MAP_GWAS$GWAS_Trait_R2_1E5[mm]
RESobs_nodup_1Mb_cond$GWAS_Trait_R2_1E8=MAP_GWAS$GWAS_Trait_R2_1E8[mm]
RESobs_nodup_1Mb_cond$GWAS_EFO_R2_1E5=MAP_GWAS$GWAS_EFO_R2_1E5[mm]
RESobs_nodup_1Mb_cond$GWAS_EFO_R2_1E8=MAP_GWAS$GWAS_EFO_R2_1E8[mm]

RESobs_nodup_1Mb$RegElt=Map_imputed$RegElt[match(RESobs_nodup_1Mb$snps,Map_imputed$snp.name)]
RESobs_nodup_1Mb_cond$RegElt=Map_imputed$RegElt[match(RESobs_nodup_1Mb_cond$snps,Map_imputed$snp.name)]

Require('allGOterms')
GoToTest=list(immuneResponse='GO:0006955',
	InnateImmuneResponse='GO:0045087',
	AdaptiveImmuneResponse='GO:0002250',
	AntiviralResponse='GO:0051607',
	AntibacterialResponse='GO:0042742',
	TF_DNAbinding='GO:0003700',
	development='GO:0032502',
	olfactory_receptors='GO:0004984')
GoToTest_genes=lapply(GoToTest,function(x){allGOterms$gene[allGOterms$go==x]})
GoToTest_genes$all=unique(allGOterms$gene)

RESobs_nodup_1Mb_cond$isImmune=RESobs_nodup_1Mb_cond$gene%in%GoToTest_genes$immuneResponse
RESobs_nodup_1Mb$isImmune=RESobs_nodup_1Mb$gene%in%GoToTest_genes$immuneResponse

#save(RESobs_nodup_1Mb,RESobs_nodup_1Mb_cond,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup.Rdata',EVO_IMMUNO_POP))
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup.Rdata',EVO_IMMUNO_POP))
write.table(RESobs_nodup_1Mb,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/table/Table_S2A_RESobs_nodup_1Mb.txt',HOME),sep='\t',quote=F)

##########################################################################################
#####		Mechanistic Bases  of splicing	Fig 2B, 2C and Supplementary Figs		  ####
##########################################################################################

library(data.table)
Require('Map_imputed')

RESobs_nodup_1Mb$RegElt=gsub('-validated','',Map_imputed$RegElt[match(RESobs_nodup_1Mb$snps,Map_imputed$snp.name)])
RESobs_nodup_1Mb_cond$RegElt=gsub('-validated','',Map_imputed$RegElt[match(RESobs_nodup_1Mb_cond$snps,Map_imputed$snp.name)])

RESobs_nodup_1Mb$RegElt_Text=RESobs_nodup_1Mb$RegElt
RESobs_nodup_1Mb$RegElt_Text[is.na(RESobs_nodup_1Mb$RegElt_Text)]=''
#RESobs_nodup_1Mb_cond$RegElt[is.na(RESobs_nodup_1Mb_cond$RegElt)]=''

SPIDEX=as.data.frame(fread(sprintf("%s/Annotation/SPIDEX/SNPlist_allSPIDEX.txt",HOME)))

library(GenomicRanges)
GRange_Map=makeGRangesFromDataFrame(Map_imputed[which(Map_imputed$SNPfreq>0.05 ),c('chromosome','position','allele.1','allele.2','minor_allele','ancestral_allele')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_PSI=makeGRangesFromDataFrame(PSI_Annot[toKeep,],seqnames.field='chrom',start.field='start',end.field='end',keep.extra.columns=TRUE)
Cisdist=1e6
GR_cis=flank(flank(GRange_PSI,Cisdist),-2*Cisdist-1)
ooCis=findOverlaps(GR_cis, GRange_Map)
wCis_1Mb=which(Map_imputed$SNPfreq>=0.05)[unique(subjectHits(ooCis))]
wCis=wCis_1Mb

RESobs_nodup_1Mb$SPIDEX_PSI=SPIDEX$DeltaPSI_SPIDEX[match(RESobs_nodup_1Mb$snps,SPIDEX$snp.name)]
RESobs_nodup_1Mb$SPIDEX_Z=SPIDEX$Z_SPIDEX[match(RESobs_nodup_1Mb$snps,SPIDEX$snp.name)]
RESobs_nodup_1Mb$Exonic=SPIDEX$Exon_zone[match(RESobs_nodup_1Mb$snps,SPIDEX$snp.name)]
RESobs_nodup_1Mb$groupSplice=ifelse(is.na(RESobs_nodup_1Mb$SPIDEX_Z),'',ifelse(abs(RESobs_nodup_1Mb$SPIDEX_Z)>2,paste(RESobs_nodup_1Mb$Exonic,'deleterious'),RESobs_nodup_1Mb$Exonic))

RESobs_nodup_1Mb_cond$SPIDEX_PSI=SPIDEX$DeltaPSI_SPIDEX[match(RESobs_nodup_1Mb_cond$snps,SPIDEX$snp.name)]
RESobs_nodup_1Mb_cond$SPIDEX_Z=SPIDEX$Z_SPIDEX[match(RESobs_nodup_1Mb_cond$snps,SPIDEX$snp.name)]
RESobs_nodup_1Mb_cond$Exonic=SPIDEX$Exon_zone[match(RESobs_nodup_1Mb_cond$snps,SPIDEX$snp.name)]
RESobs_nodup_1Mb_cond$groupSplice=ifelse(is.na(RESobs_nodup_1Mb_cond$SPIDEX_Z),'',ifelse(abs(RESobs_nodup_1Mb_cond$SPIDEX_Z)>2,paste(RESobs_nodup_1Mb_cond$Exonic,'deleterious'),RESobs_nodup_1Mb_cond$Exonic))

Geno_sQTL=getGenos(unique(RESobs_nodup_1Mb_cond$snps))
r2_sQTL=cor(t(as.matrix(Geno_sQTL[-(1:5)])),use='p')^2
library(igraph)
cor_graph=graph_from_adjacency_matrix(r2_sQTL>0.8)
wc <- cluster_walktrap(cor_graph)
haplos=membership(wc)
names(haplos)=Geno_sQTL$snp.name
RESobs_nodup_1Mb$haplo=paste('haplo',haplos)[match(RESobs_nodup_1Mb$snps,names(haplos))]
RESobs_nodup_1Mb_cond$haplo=paste('haplo',haplos)[match(RESobs_nodup_1Mb_cond$snps,names(haplos))]

save(RESobs_nodup_1Mb,RESobs_nodup_1Mb_cond,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup.Rdata',EVO_IMMUNO_POP))

###################################
###		add response data		###
###################################

RESobs_expr=RESobs;rm(RESobs)
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_response_V5.Rdata',EVO_IMMUNO_POP))
RESobs_has_sQTL_1Mb=RESobs[which(abs(RESobs$CisDist)<1e6 & paste(RESobs$event_id,RESobs$snps)%in%paste(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$snps)),]
#Detected_sQTL_P5pct=By(RESobs_has_sQTL_1Mb$cond,paste(RESobs_has_sQTL_1Mb$event_id, RESobs_has_sQTL_1Mb$snps),function(x){paste(sort(x),collapse='')})
#RESobs_nodup_1Mb$detected_cond=Detected_sQTL_P5pct[match(paste(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$snps),names(Detected_sQTL_P5pct))]
#RESobs_nodup_1Mb$NbTestableCond=PSI_Annot$NbTestableCond[match(RESobs_nodup_1Mb$event_id,PSI_Annot$event_id)]

Pval_rsQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='pvalue');
Pval_rsQTL=Pval_rsQTL[match(RESobs_nodup_1Mb$event_id,Pval_rsQTL$event_id),-1]
Pval_rsQTL[is.na(Pval_rsQTL)]=1
colnames(Pval_rsQTL)=paste('Pvalue_rsQTL',condIndex[-1],sep='_')
Beta_rsQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='beta')
Beta_rsQTL=Beta_rsQTL[match(RESobs_nodup_1Mb$event_id,Beta_rsQTL$event_id),-1]
Beta_rsQTL[is.na(Beta_rsQTL)]=0
colnames(Beta_rsQTL)=paste('Beta_rsQTL',condIndex[-1],sep='_')

RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,Pval_rsQTL,Beta_rsQTL)
save(RESobs_nodup_1Mb,RESobs_nodup_1Mb_cond,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup.Rdata',EVO_IMMUNO_POP))


###################################
###		add Gerp scores 		###
###################################

library(data.table)
Map_Gerp=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Map_Gerp_snp.txt',HOME))
median(Map_Gerp$V5[Map_Gerp$V5!=0],na.rm=T)
RESobs_nodup_1Mb$GerpRS=Map_Gerp$V5[match(RESobs_nodup_1Mb$snps,Map_Gerp$V1)]

odds.ratio(table(Map_Gerp$V5[Map_Gerp$V5!=0]> 4 , Map_Gerp$V1[Map_Gerp$V5!=0]%in% RESobs_nodup_1Mb$snps))
odds.ratio(table(Map_Gerp$V5[Map_Gerp$V5!=0]> 2 , Map_Gerp$V1[Map_Gerp$V5!=0]%in% RESobs_nodup_1Mb$snps))
odds.ratio(table(Map_Gerp$V5[Map_Gerp$V5!=0]< -2 , Map_Gerp$V1[Map_Gerp$V5!=0]%in% RESobs_nodup_1Mb$snps))
odds.ratio(table(Map_Gerp$V5[Map_Gerp$V5!=0]< -4 , Map_Gerp$V1[Map_Gerp$V5!=0]%in% RESobs_nodup_1Mb$snps))


###################################
###		add eQTL overlap 		###
###################################


load(sprintf('%s/03_Analysis/cis-eQTL/Table_eQTL_Finale_04072016.Rdata',HOME))
luq(PSI_Annot$gene_id[toKeep])
gene_tested=unique(PSI_Annot$gene_id[toKeep])
has_sQTL=gene_tested%in%RESobs_nodup_1Mb$gene
has_eQTL=gene_tested%in%RES_cis_Final$gene
odds.ratio(table(has_sQTL, has_eQTL))
#            LowerCI       OR  UpperCI alpha            P
# odds ratio 2.294884 2.696806 3.167383  0.05 3.600575e-33


#snps_eQTL_EUB=sapply(RES_cis_Final$snps[RES_cis_Final$population=='EUB'],getSNP,pop='EUB')
#snps_eQTL_AFB=sapply(RES_cis_Final$snps[RES_cis_Final$population=='AFB'],getSNP,pop='AFB')

Genos_eQTL=getGenos(unique(RES_cis_Final$snps))
snps_eQTL=t(as.matrix(Genos_eQTL[-(1:5)]))
colnames(snps_eQTL)=Genos_eQTL[[1]]

snps_eQTL_EUB=snps_eQTL[substr(rownames(snps_eQTL),1,3)=='EUB',match(RES_cis_Final$snps[RES_cis_Final$population=='EUB'],colnames(snps_eQTL))]
snps_eQTL_AFB=snps_eQTL[substr(rownames(snps_eQTL),1,3)=='AFB',match(RES_cis_Final$snps[RES_cis_Final$population=='AFB'],colnames(snps_eQTL))]

Genos_sQTL=getGenos(unique(RESobs_nodup_1Mb_cond$snps))
snps_sQTL=t(as.matrix(Genos_sQTL[-(1:5)]))
colnames(snps_sQTL)=Genos_sQTL[[1]]
snps_sQTL=snps_sQTL[,match(RESobs_nodup_1Mb$snps,colnames(snps_sQTL))]

gene_sQTL=mapply(function(gene,cond){FPKM_gene[gene,match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(FPKM_gene))]},RESobs_nodup_1Mb$gene,RESobs_nodup_1Mb$cond)

r2_sQTL_eQTL=mapply(function(gene,snp){
		snps_EUB=RES_cis_Final$snps[RES_cis_Final$population=='EUB' & RES_cis_Final$gene==gene]
		snps_AFB=RES_cis_Final$snps[RES_cis_Final$population=='AFB' & RES_cis_Final$gene==gene]
		LD_EUB=cor(snps_sQTL[match(rownames(snps_eQTL_EUB),rownames(snps_sQTL)),snp],snps_eQTL_EUB[,snps_EUB],use='p')^2
		LD_AFB=cor(snps_sQTL[match(rownames(snps_eQTL_AFB),rownames(snps_sQTL)),snp],snps_eQTL_AFB[,snps_AFB],use='p')^2
		max(c(LD_AFB,LD_EUB),na.rm=T)
 },RESobs_nodup_1Mb$gene,RESobs_nodup_1Mb$snps)


psi_sQTL=mapply(function(event,cond){PSI_prov[match(event,PSI_Annot[toKeep,'event_id']),match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(PSI_prov))]},RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$cond)

P_sQTL_geneExpr=apply(rbind(snps_sQTL,gene_sQTL),2,function(x){
n=length(x)/2
SNP=x[1:n]
GENE=x[n+1:n]
POP=substr(rownames(snps_sQTL),1,3)
summary(lm(GENE~SNP+POP))$coeff[2,4]
})
B_sQTL_geneExpr=apply(rbind(snps_sQTL,gene_sQTL),2,function(x){
n=length(x)/2
SNP=x[1:n]
GENE=x[n+1:n]
POP=substr(rownames(snps_sQTL),1,3)
summary(lm(GENE~SNP+POP))$coeff[2,1]
})

RESobs_nodup_1Mb$r2_sQTL_eQTL=r2_sQTL_eQTL
RESobs_nodup_1Mb$P_sQTL_geneExpr=P_sQTL_geneExpr
RESobs_nodup_1Mb$B_sQTL_geneExpr=B_sQTL_geneExpr


###################################
###		add causality test 		###
###################################

getLikelihoods=function(x,y,z,u=NULL){
	if(!is.null(u)){
		w=which(!is.na(y*x*z*u))
		}else{
		w=which(!is.na(y*x*z))
		}
	x=x[w]
	y=y[w]
	z=z[w]
	if(!is.null(u)){
		u=u[w]
	mod_xCausal=logLik(lm(y~x+u))+logLik(lm(x~z+u))
	mod_yCausal=logLik(lm(x~y+u))+logLik(lm(y~z+u))
	mod_xyInd=logLik(lm(y~z+u))+logLik(lm(x~z+u))
#	mod_xyIndCor=logLik(lm(x~y+z+u))+logLik(lm(y~z+u))
#	mod_yxIndCor=logLik(lm(y~x+z+u))+logLik(lm(x~z+u))
		}else{
	mod_xCausal=logLik(lm(y~x))+logLik(lm(x~z))
	mod_yCausal=logLik(lm(x~y))+logLik(lm(y~z))
	mod_xyInd=logLik(lm(y~z))+logLik(lm(x~z))
#	mod_xyIndCor=logLik(lm(x~y+z))+logLik(lm(y~z))
#	mod_yxIndCor=logLik(lm(y~x+z))+logLik(lm(x~z))
	}
	c(mod_xCausal,mod_yCausal,mod_xyInd)
}

getProbs=function(myPSI,geneFPKM,SNP,u=NULL){
	logLiks=getLikelihoods(myPSI,geneFPKM,SNP,u=u)
	probs=exp(logLiks[1:3])/sum(exp(logLiks[1:3]))
	names(probs)=c('SNP->Splice->Expr','SNP->Expr->Splice','Splice<-SNP->Expr')
	probs
}

getLikelihoods_uncertain=function(x,y,z1,z2,u=NULL){
	if(!is.null(u)){
		w=which(!is.na(y*x*z1*z2*u))
		}else{
		w=which(!is.na(y*x*z1*z2))
		}
	x=x[w]
	y=y[w]
	z=z[w]
	if(!is.null(u)){
		u=u[w]
	mod_xCausal=logLik(lm(y~x+u))+logLik(lm(x~z1+u))
	mod_yCausal=logLik(lm(x~y+u))+logLik(lm(y~z2+u))
	mod_xyInd=logLik(lm(y~z2+u))+logLik(lm(x~z1+u))
#	mod_xyIndCor=logLik(lm(x~y+z+u))+logLik(lm(y~z+u))
#	mod_yxIndCor=logLik(lm(y~x+z+u))+logLik(lm(x~z+u))
		}else{
	mod_xCausal=logLik(lm(y~x))+logLik(lm(x~z1))
	mod_yCausal=logLik(lm(x~y))+logLik(lm(y~z2))
	mod_xyInd=logLik(lm(y~z2))+logLik(lm(x~z1))
#	mod_xyIndCor=logLik(lm(x~y+z))+logLik(lm(y~z))
#	mod_yxIndCor=logLik(lm(y~x+z))+logLik(lm(x~z))
	}
	c(mod_xCausal,mod_yCausal,mod_xyInd)
}

getProbs_uncertain=function(myPSI,geneFPKM,PSI.pred,GENE.pred,u=NULL){
	logLiks=getLikelihoods_uncertain(myPSI,geneFPKM,PSI.pred,GENE.pred,u=u)
	probs=exp(logLiks[1:3])/sum(exp(logLiks[1:3]))
	names(probs)=c('SNP->Splice->Expr','SNP->Expr->Splice','Splice<-SNP->Expr')
	probs
}


ProbsModel=apply(rbind(snps_sQTL,gene_sQTL,psi_sQTL),2,function(x){
	n=length(x)/3
	SNP=x[1:n]
	GENE=x[n+1:n]
	PSI=x[2*n+1:n]
	POP=as.numeric(substr(rownames(snps_sQTL),1,3)=='AFB')
	getProbs(PSI,GENE,SNP,POP)
})

chrom_sQTL=as.numeric(gsub('([0-9]+):[0-9]+Mb','\\1',RESobs_nodup_1Mb$locus))
pos_sQTL=as.numeric(RESobs_nodup_1Mb$pos)
myDist=1e6

ProbsModel_uncertainSNP=apply(rbind(chrom_sQTL,pos_sQTL,myDist,gene_sQTL,psi_sQTL),2,function(x){
	n=(length(x)-3)/2
	chrom=x[1]
	pos=x[2]
	dist=x[3]/2
	Geno=2-as.matrix(getGenos_locus(chrom,pos-dist,pos+dist)[-(1:5)])
	MAF=apply(Geno,1,mean,na.rm=T)
	Geno=Geno[MAF>=0.05,]
	Geno=impute.knn(Geno)$data
	GENE=x[3+1:n]
	PSI=x[3+n+1:n]

	require(elasticnet)
	w=which(!is.na(GENE) & !is.na(PSI))
	cv.GENE=cv.glmnet(t(Geno[,w]),GENE[w],alpha=0.5)
	GENE.mod=glmnet(t(Geno[,w]),GENE[w],alpha=0.5,lambda=cv.GENE$lambda.min)
	GENE.pred=(t(as.matrix(GENE.mod$beta))%*%Geno[,w])[1,]+GENE.mod$a0
	cv.PSI=cv.glmnet(t(Geno[,w]),PSI[w],alpha=0.5)
	PSI.mod=glmnet(t(Geno[,w]),PSI[w],alpha=0.5,lambda=cv.PSI$lambda.min)
	PSI.pred=(t(as.matrix(PSI.mod$beta))%*%Geno[,w])[1,]+PSI.mod$a0
	POP=as.numeric(substr(rownames(snps_sQTL),1,3)=='AFB')

	getProbs_uncertain(PSI,GENE,PSI.pred,GENE.pred,POP)
})

BestModel=apply(ProbsModel,2,function(x){ifelse(any(is.na(x)),NA,which.max(x))})
RESobs_nodup_1Mb$BestModel=rownames(ProbsModel)[BestModel]
RESobs_nodup_1Mb$BestModel_Prob=ProbsModel[cbind(BestModel,1:ncol(ProbsModel))]

save(RESobs_nodup_1Mb,RESobs_nodup_1Mb_cond,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup.Rdata',EVO_IMMUNO_POP))

library(data.table)
#Map_imputed=as.data.frame(fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/selection/resampling_selection/Map_selection.txt',HOME)))
Map_select=as.data.frame(fread(sprintf('%s/Annotation/Select/Map_Select_LD_allChr_V2.txt',HOME)))
#cn(Map_select)
mm=match(RESobs_nodup_1Mb$snps,Map_select$snp.name)

RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,Map_select[mm,-(1:3)])
mm=match(RESobs_nodup_1Mb_cond$snps,Map_select$snp.name)
RESobs_nodup_1Mb_cond=cbind(RESobs_nodup_1Mb_cond,Map_select[mm,-(1:3)])

save(RESobs_nodup_1Mb,RESobs_nodup_1Mb_cond,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup_withSelect.Rdata',EVO_IMMUNO_POP))






Require('Map_imputed')
tab0=table(gsub('-validated','',Map_imputed$RegElt[wCis]),exclude='')
tab1=table(gsub('-validated','',RESobs_nodup_1Mb$RegElt),RESobs_nodup_1Mb$event_type,exclude='')
#tab2=table(gsub('-validated','',RESobs_nodup_5kb$RegElt),RESobs_nodup_5kb$event_type,exclude='')
barplot(cbind(tab0/sum(tab0),t(t(tab1)/apply(tab1,2,sum))))
#barplot(cbind(tab0/sum(tab0),tab2/sum(tab2)))
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/SPIDEX/RegElt_barplot.pdf',HOME))
###################################
###	Figure 2B	RegElt detail 	###
###################################
par(xpd=T)
barplot(cbind(All=tab0[1:4]/sum(tab0),t(t(tab1[1:4,])/apply(tab1,2,sum))),col=c(colPSI[c(1:3,6)]),las=1,ylab='% of splice QTL',ylim=c(0,0.4))
legend('topright',fill=c(colPSI[c(1:3,6)]),c('CTCF','Enhancer','Promoter','Promoter Flanks'),bty='n')
par(xpd=F)
###################################
###	Figure 2B RegElt detail END ###
###################################

###################################
###	Figure 2B RegElt grouped	###
###################################
par(xpd=T)
barplot(cbind(All=tab0[1:4]/sum(tab0),sQTL=apply(tab1[1:4,],1,sum)/sum(tab1),rsQTL=apply(tab2[1:4,],1,sum)/sum(tab2)),col=c(colPSI[c(1:3,6)]),las=1,ylab='% of splice QTL',las=3)
legend(fill=c(colPSI[c(1:3,6)]),legend=c('CTCF','Enhancer','Promoter','Promoter Flanks'),bty='n',ncol=3)
par(xpd=F)
# las =3 for horizontal
###################################
###	Figure 2B RegElt grouped END###
###################################

tabAll=table(gsub('-validated','',Map_imputed$RegElt[wCis]),exclude='');names(tabAll)[is.na(names(tabAll))]='x'
X=cbind(tabAll,sapply(c(1e-5,1e-10,1e-30),function(th){tab=table(gsub('-validated','',RESobs_nodup_1Mb$RegElt[RESobs_nodup_1Mb$pval<th]),exclude='');names(tab)[is.na(names(tab))]='x';tab=tab[names(tabAll)];tab[is.na(tab)]=0;tab}))
Proportions_RegElt=t(t(X[-5,])/apply(X,2,sum))
colnames(Proportions_RegElt)=c('All','sQTL, P<10-5','sQTL, P<10-10','sQTL, P<10-30')
###################################
###		Figure 2B 				###
###################################
barplot(Proportions_RegElt,col=colPSI[c(1:3,6)],las=1,ylim=c(0,0.4),ylab='% of splice QTL')
barplot(Proportions_RegElt,col=colPSI[c(1:3,6)],las=1,ylim=c(0,0.4),ylab='% of splice QTL',las=2)
###################################
###	Figure 2B Significance		###
###################################

OR=NULL
for (j in 1:4){
	for (i in 1:7){
		tab=rbind(c(sum(tab0[-j])-sum(tab1[-j,i]),sum(tab1[-j,i])),
			c(tab0[j]-tab1[j,i],tab1[j,i]))
		myOR=odds.ratio(tab)
		myOR$event_type=colnames(tab1)[i]
		myOR$RegElt_type=rownames(tab1)[j]
		OR=rbind(OR,myOR)
	}
	tab=rbind(c(sum(tab0[-j])-sum(tab1[-j,]),sum(tab1[-j,])),
		c(tab0[j]-sum(tab1[j,]),sum(tab1[j,])))
	myOR=odds.ratio(tab)
	myOR$event_type="all sQTL"
	myOR$RegElt_type=rownames(tab1)[j]
	OR=rbind(OR,myOR)
	tab=rbind(c(sum(tab0[-j])-sum(tab2[-j,]),sum(tab2[-j,])),
		c(tab0[j]-sum(tab2[j,]),sum(tab2[j,])))
	myOR=odds.ratio(tab)
	myOR$event_type="rsQTL"
	myOR$RegElt_type=rownames(tab2)[j]
	OR=rbind(OR,myOR)
}
OR=as.data.frame(OR)

OR$FDR=p.adjust(OR$P,'fdr')

###################################
###		Figure 2B END			###
###################################

Map_imputed$SPIDEX_Z=SPIDEX$Z_SPIDEX[match(Map_imputed$snp.name,SPIDEX$snp.name)]
Map_imputed$Exonic=SPIDEX$Exon_zone[match(Map_imputed$snp.name,SPIDEX$snp.name)]
Map_imputed$groupSplice=ifelse(is.na(Map_imputed$SPIDEX_Z),'',ifelse(abs(Map_imputed$SPIDEX_Z)>2,paste(Map_imputed$Exonic,'deleterious'),Map_imputed$Exonic))

acol= c(rev(brewer.pal(4,'YlOrRd')),brewer.pal(4,'YlGnBu'))
rcol = SetAlpha(acol, 0.5)
#expectedCTCF=mean(grepl('CTCF',Map_imputed$RegElt[which(Map_imputed$SNPfreq>0.05)]),na.rm=T)
#observedCTCF=mean(grepl('CTCF',x$RegElt),na.rm=T)
#
#expectedEnhancer=mean(grepl('Flank|Enhancer',Map_imputed$RegElt[which(Map_imputed$SNPfreq>0.05)]),na.rm=T)
#observedEnhancer=mean(grepl('Flank|Enhancer',x$RegElt),na.rm=T)
tab=table(Map_imputed$groupSplice[wCis])
X=cbind(tab,sapply(c(1e-5,1e-10,1e-30),function(th){tab=table(RESobs_nodup_1Mb$groupSplice[RESobs_nodup_1Mb$pval<th])}))
Proportions_Exon=t(t(X[-1,])/apply(X,2,sum))
colnames(Proportions_Exon)=c('All','sQTL, P<10-5','sQTL, P<10-10','sQTL, P<10-30')

###################################
###		Figure 2C 				###
###################################
par(mar=c(7,5,1,3))
barplot(Proportions_Exon,col=acol[c(3,1,6,7)],las=3,ylim=c(0,0.4),ylab='% of splice QTL')
#barplot(Proportions_Exon,col=acol[c(3,1,6,7)],las=1,ylim=c(0,0.4),ylab='% of splice QTL',las=2)
###################################
###		Figure 2C END			###
###################################

#barplot(Proportions_Exon,col=acol[c(3,1,6,7)],las=1,ylim=c(0,0.6),ylab='% of splice QTL')
par(xpd=T)
barplot(cbind(All=tabE0[1+c(2,1,4,3)]/sum(tabE0),apply(tabE1,1,sum)[1+c(2,1,4,3)]/sum(tabE1)),col=acol[c(1,3,7,6)],las=1,ylim=c(0,0.6),ylab='% of splice QTL')
legend(0,0.67,fill=acol[c(7,6,3,1)], legend=c('intronic, <300 nt from Splice site',"intronic, |SPIDEX's Z|>2",'exonic, <300 nt from Splice site',"exonic, |SPIDEX's Z|>2"),bty='n')
par(xpd=T)
dev.off()
tabE0=table(Map_imputed$groupSplice[wCis],exclude=NA)
tabE1=table(RESobs_nodup_1Mb$groupSplice,RESobs_nodup_1Mb$event_type,exclude=NA)
OR=NULL
for (j in 2:5){
	for (i in 1:7){
		tab=rbind(c(sum(tabE0[-j])-sum(tabE1[-j,i]),sum(tabE1[-j,i])),
			c(tabE0[j]-tabE1[j,i],tabE1[j,i]))
		myOR=odds.ratio(tab)
		myOR$event_type=colnames(tabE1)[i]
		myOR$Splice_type=rownames(tabE1)[j]
		OR=rbind(OR,myOR)
	}
	tab=rbind(c(sum(tabE0[-j])-sum(tabE1[-j,]),sum(tabE1[-j,])),
		c(tabE0[j]-sum(tabE1[j,]),sum(tabE1[j,])))
	myOR$event_type="All"
	myOR$Splice_type=rownames(tabE1)[j]
	OR=rbind(OR,myOR)
}
for (i in 1:7){
	tab=rbind(c(tabE0[1]-tabE1[1,i],tabE1[1,i]),
	c(sum(tabE0[-1])-sum(tabE1[-1,i]),sum(tabE1[-1,i])))
	myOR=odds.ratio(tab)
	myOR$event_type=colnames(tabE1)[i]
	myOR$Splice_type='All'
	OR=rbind(OR,myOR)
}
tab=rbind(c(sum(tabE0[c(1,4:5)])-sum(tabE1[c(1,4:5),]),sum(tabE1[c(1,4:5),])),
		c(sum(tabE0[2:3])-sum(tabE1[2:3,]),sum(tabE1[2:3,])))
	myOR=odds.ratio(tab)
	myOR$event_type='All'
	myOR$Splice_type='All Exonic'
	OR=rbind(OR,myOR)
tab=rbind(c(sum(tabE0[c(1,2:3)])-sum(tabE1[c(1,2:3),]),sum(tabE1[c(1,2:3),])),
				c(sum(tabE0[4:5])-sum(tabE1[4:5,]),sum(tabE1[4:5,])))
	myOR=odds.ratio(tab)
	myOR$event_type='All'
	myOR$Splice_type='All Intronic'
	OR=rbind(OR,myOR)
tab=rbind(c(sum(tabE0[c(2,4)])-sum(tabE1[c(2,4),]),sum(tabE1[c(2,4),])),
				c(sum(tabE0[c(3,5)])-sum(tabE1[c(3,5),]),sum(tabE1[c(3,5),])))
	myOR=odds.ratio(tab)
	myOR$event_type='All'
	myOR$Splice_type='deleterious'

	OR=rbind(OR,myOR)
	for (j in 2:4){
		tab=rbind(c(X[1,1]-X[1,j],sum(X[1,j])),
			c(sum(X[-1,1])-sum(X[-1,j]),sum(X[-1,j])))
		myOR=odds.ratio(tab)
		myOR$event_type="All"
		myOR$Splice_type=colnames(Proportions_Exon)[j]
		OR=rbind(OR,myOR)
		tab=rbind(c(sum(X[c(2,4),1])-sum(X[c(2,4),j]),sum(X[c(2,4),j])),
				c(sum(X[c(3,5),1])-sum(X[c(3,5),j]),sum(X[c(3,5),j])))
		myOR=odds.ratio(tab)
		myOR$event_type='All'
		myOR$Splice_type=paste('deleterious',colnames(Proportions_Exon)[j])
		OR=rbind(OR,myOR)
	}
write.table(OR,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/SPIDEX/SpliceElt_OR.txt',HOME),sep='\t',quote=F,row.names=F)

##########################################################################################
#####					frequency and sharing of sQTLs								  ####
##########################################################################################

###### Figure 2D with 1Mb
tab_type=table(RESobs_nodup_1Mb[,'event_type'])
tab_type_tested=table(PSI_Annot$event_type[toKeep])

barplot(tab_type/tab_type_tested,beside=T,col=colPSI,las=1)
OR=NULL
for (i in 1:length(colPSI)){
	y=c(tab_type_tested[i]-tab_type[i],tab_type[i])
	X=rbind(apply(cbind(tab_type_tested-tab_type,tab_type),2,sum)-y,y)
	OR=rbind(OR,odds.ratio(X))
}
rownames(OR)=names(tab_type_tested)
OR$FDR=p.adjust(OR$P,'fdr')
write.table(OR,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/OR_byType_V3_1Mb.txt',HOME),sep='\t',quote=F,row.names=F)

###### Figure 2D with 5kb
tab_type=table(RESobs_nodup_5kb[,'event_type'])
tab_type_tested=table(PSI_Annot$event_type[toKeep])

barplot(tab_type/tab_type_tested,beside=T,col=colPSI,las=1)
OR=NULL
for (i in 1:length(colPSI)){
	y=c(tab_type_tested[i]-tab_type[i],tab_type[i])
	X=rbind(apply(cbind(tab_type_tested-tab_type,tab_type),2,sum)-y,y)
	OR=rbind(OR,odds.ratio(X))
}
rownames(OR)=names(tab_type_tested)
OR$FDR=p.adjust(OR$P,'fdr')
write.table(OR,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/OR_byType_V3_5kb.txt',HOME),sep='\t',quote=F,row.names=F)


###################################
###		Figure 2D				###
###################################
par(xpd=F)
plot(c(0,2.5),c(1,7.5),col='#00000000',axes=F,ylab='',xlab='Odds ratio')
axis(1)
axis(2,at=0.5+1:7,names(colPSI),las=2)
abline(v=1,col='grey')
points(OR$OR,0.5+1:7,pch=16,cex=1.3,col=colPSI)
for(i in 1:7){segments(OR$LowerCI[i],0.5+i,OR$UpperCI[i],0.5+i,lwd=2,col=colPSI[i])}
###################################
###		Figure 2D	END			###
###################################

###################################
###		Figure 2E	numbers		###
###################################
tab_NbCond=table(nchar(RESobs_nodup_5kb$detected_cond_testable))
tab_NbCond_AllTestable=table(nchar(RESobs_nodup_5kb$detected_cond_testable[RESobs_nodup_5kb$NbTestableCond==5 & RESobs_nodup_5kb$pval<1e-6]))
tab_NbCond_AllTestable=table(nchar(RESobs_nodup_5kb$detected_cond_testable[RESobs_nodup_5kb$NbTestableCond==5]))
tab_NbCond
tab_NbCond_AllTestable
(tab_NbCond)/sum(tab_NbCond)
(tab_NbCond_AllTestable)/sum(tab_NbCond_AllTestable)

tab_NbCond=table(nchar(RESobs_nodup_1Mb$detected_cond_testable))
tab_NbCond_AllTestable=table(nchar(RESobs_nodup_1Mb$detected_cond_testable[RESobs_nodup_1Mb$NbTestableCond==5]))
tab_NbCond
tab_NbCond_AllTestable
(tab_NbCond)/sum(tab_NbCond)
(tab_NbCond_AllTestable)/sum(tab_NbCond_AllTestable)

###################################
###		Figure 2E				###
###################################

par(xpd=T)
X=cbind(all=rev(tab_NbCond),shared=rev(tab_NbCond_AllTestable))
barplot(X,ylab='number of sQTLs',col=acol[c(1,3,4,6,7)],las=1,ylim=c(0,max(apply(X,2,sum))*1.1))
text(0:1*1.2+0.7,apply(X,2,sum)+max(X)*0.09,labels=apply(X,2,sum))
#legend(1.3,1100,bty='n',paste(1:5,'cond'),fill=rev(acol[c(1,3,4,6,7)]))
legend(1.3,1600,bty='n',paste(1:5,'cond'),fill=rev(acol[c(1,3,4,6,7)]))

###################################
###		Figure 2E horizontal	###
###################################

par(xpd=T)
X=cbind(all=rev(tab_NbCond),shared=rev(tab_NbCond_AllTestable))
barplot(X,ylab='number of sQTLs',col=acol[c(1,3,4,6,7)],las=3,ylim=c(0,max(apply(X,2,sum))*1.15))
text(0:1*1.2+0.7,apply(X,2,sum)+max(X)*0.3,labels=apply(X,2,sum),srt=90)
#legend(1.3,1100,bty='n',paste(1:5,'cond'),fill=rev(acol[c(1,3,4,6,7)]))
#legend(1.3,1700,bty='n',paste(1:5,'cond'),fill=rev(acol[c(1,3,4,6,7)]))

#barplot(t(t(X)/apply(X,2,sum)),ylab='%s Events with sQTL',col=acol,las=1,ylim=c(0,1.1))
#text(0:1*1.2+0.7,c(1,1)*1.05,labels=paste('N =',apply(X,2,sum)))

###################################
###		Figure 2E	END			###
###################################

##########################################################################################
#####					GO enrichment of sQTLs										  ####
##########################################################################################
# sQTls are enriched in defense response genes
resGO=GOSeq(unique(RESobs_nodup_1Mb$gene),unique(PSI_Annot$gene_id[toKeep]),overOnly=F,addCI=T)
resGO_cond=lapply(1:5,function(i){GOSeq(unique(RESobs_nodup_1Mb_cond$gene[RESobs_nodup_1Mb_cond$cond==i]),unique(PSI_Annot$gene_id[toKeep][PSI_Annot[toKeep,grep('Testable_',cn(PSI_Annot))[i]]]),overOnly=F,addCI=T)})
resGO_cond_table=do.call(rbind,lapply(1:5,function(i){cbind(resGO_cond[[i]],cond=rep(condIndex[i],nrow(resGO_cond[[i]])))}))
resGO$cond='ALL_POOLED'
TableS2C=rbind(resGO_cond_table,resGO)
write.table(TableS2C,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/Table_GO_sQTL.txt',HOME),quote=F,sep='\t',row.names=F)
pdf('',width=9.5,height=5.6)
plotGO(resGO[1:15,])
TableS2C=read.table(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/Table_GO_sQTL.txt',HOME),sep='\t',header=T)
plotGO(TableS2C[which(TableS2C$cond=='ALL_POOLED')[1:15],])
plotGO(TableS2C[rev(which(TableS2C$cond=='ALL_POOLED'))[1:15],])
gene_list=intersect(GoToTest_genes$all,unique(PSI_Annot$gene_id[toKeep]))
odds.ratio(table(gene_list%in%RESobs_nodup_1Mb$gene,gene_list%in%GoToTest_genes$immuneResponse))


###################################
###		Figure 2F				###
###################################

hist(pmax(0,r2_sQTL_eQTL),main='',xlab='r2 with eQTL',col='darkgrey',br=seq(0,1,0.1),border='darkgrey',las=2)
hist(r2_sQTL_eQTL,add=T,br=seq(0,1,0.1),col='lightgrey',border='darkgrey')
hist(r2_sQTL_eQTL[r2_sQTL_eQTL>0.8],add=T,br=seq(0,1,0.1),col=colERC5[2],border='darkgrey')
barplot(By(r2_sQTL_eQTL>0.8,RESobs_nodup_1Mb$event_type,mean),col=colPSI)

###################################
###		Figure 2F	END			###
###################################


x=load('/Volumes/evo_immuno_pop/Maxime/SNP_annotation/GWAScat_12062015.Rdata')
nrow(GWAScat) # 18875
sum(GWAScat$p.Value<1e-8) # 6844

GWAScat_2015=GWAScat
library(gwascat)
library(rtracklayer)
ch=import.chain("/Volumes/@home/Annotation/liftoverFiles/hg38ToHg19.over.chain")
#http://www.genome.gov/admin/gwascatalog.txt
tab = read.delim(sprintf('%s/Annotation/GWAS/gwas_catalog_v1.0.1-associations_e89_r2017-06-26_downloaded05072017_current.tsv',HOME), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tab = gwascat:::fixNonASCII(tab)
tab$CHR_ID[tab$CHR_ID=='23']='X'
tab$CHR_POS=as.numeric(tab$CHR_POS)
tab=tab[!is.na(tab$CHR_POS),]
tab_GR= makeGRangesFromDataFrame(tab,seqnames='CHR_ID', start.field ='CHR_POS',end.field='CHR_POS',keep.extra.columns=TRUE)
seqlevelsStyle(tab_GR) = "UCSC"  # necessary
cur19 = liftOver(tab_GR, ch)
GWAScat=as.data.frame(cur19)
GWAScat$seqnames=gsub('chr','',as.character(GWAScat$seqnames))
GWAScat$strand=as.character(GWAScat$strand)
nrow(GWAScat) #40796
sum(GWAScat$P.VALUE<1e-8) # 14004
save(GWAScat,file=sprintf('%s/Annotation/GWAS/gwas_catalog_v1_downloaded05072017.RData',HOME))


###################################
###		Figure 2G				###
###################################


###################################
###		Figure 2G	END			###
###################################



###################################
###		Table S2A				###
###################################
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V3_nodup.Rdata',EVO_IMMUNO_POP))
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V3.Rdata',EVO_IMMUNO_POP))
TableS2=RESobs_nodup_1Mb
Pval_sQTL=sapply(1:5,function(i){RESobs$pval[match(paste(RESobs_nodup_1Mb$snps,RESobs_nodup_1Mb$event_id,i),paste(RESobs$snps,RESobs$event_id,RESobs$cond))]})
#FDR_sQTL=sapply(1:5,function(i){RESobs$FDR_1Mb[match(paste(RESobs_nodup_1Mb$snps,RESobs_nodup_1Mb$event_id,i),paste(RESobs$snps,RESobs$event_id,RESobs$cond))]})
R2_sQTL=sapply(1:5,function(i){RESobs$R2[match(paste(RESobs_nodup_1Mb$snps,RESobs_nodup_1Mb$event_id,i),paste(RESobs$snps,RESobs$event_id,RESobs$cond))]})

beta_sQTL=sapply(1:5,function(i){RESobs$beta[match(paste(RESobs_nodup_1Mb$snps,RESobs_nodup_1Mb$event_id,i),paste(RESobs$snps,RESobs$event_id,RESobs$cond))]})
colnames(Pval_sQTL)=paste('Pval',condIndex,sep='_')
#colnames(FDR_sQTL)=paste('FDR',condIndex,sep='_')
colnames(R2_sQTL)=paste('R2',condIndex,sep='_')
colnames(beta_sQTL)=paste('beta',condIndex,sep='_')

Testable=apply(PSI_Annot[match(RESobs_nodup_1Mb$event_id,PSI_Annot$event_id),grep('Testable_',colnames(PSI_Annot))],2,ifelse,'yes','')
colnames(Testable)=paste('Alternatively_spliced',condIndex,sep='_')
TableS2=cbind(TableS2,Pval_sQTL,R2_sQTL,beta_sQTL,Testable)
TableS2$isImmune=ifelse(TableS2$isImmune,'yes','')
TableS2$event_coord=gsub('ENSG[0-9]+;.*:[0-9XY]+:(.*):.*','\\1',TableS2$event_id)
TableS2$FDR_1Mb[TableS2$FDR_1Mb<=1e-6]='<1E-6'
TableS2$cond=condIndex[TableS2$cond]
By(sign(TableS2$B_sQTL_geneExpr)*sign(TableS2$beta)>0, paste(TableS2$event_type, TableS2$r2_sQTL_eQTL>0.8),function(x){binom.test(sum(x),length(x))})
# Exon inclusion increase => expression decreases
# Intron retention increase => expression increase

write.table(TableS2,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/TableS2_V3.2_OneLinePerEvent.txt',HOME),quote=F,sep='\t',row.names=F)

###################################
###		Table S2A END			###
###################################

RES_NCF1=RESobs[which(RESobs$symbol=='NCF1'),]
RES_MVP=RESobs[which(RESobs$symbol=='MVP'),]
RES_ITGAL=RESobs[which(RESobs$symbol=='ITGAL'),]
RES_RABGAP1=RESobs[which(RESobs$symbol=='RABGAP1'),]
RES_CYP3A5=RESobs[which(RESobs$symbol=='CYP3A5'),]
RES_BTN3A3=RESobs[which(RESobs$symbol=='BTN3A3'),]
RES_CD226=RESobs[which(RESobs$symbol=='CD226'),]
save(RES_CD226,RES_BTN3A3,RES_CYP3A5,RES_RABGAP1,RES_ITGAL,RES_MVP,RES_NCF1,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_V3_targetGenes.Rdata',EVO_IMMUNO_POP))

# /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/plotCoverage/sQTLplot_runner.sh 7 74192282 74202975 rs12113435 NCF1 'ENSG00000158517;RI:7:74197282:74197404-74197868:74197975:+'
# /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/plotCoverage/sQTLplot_runner.sh 16 29840718 29852025 rs10204 MVP	'ENSG00000013364;A5:16:29845718-29847025:29845387-29847025:+'
# /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/plotCoverage/sQTLplot_runner.sh 16 30479219 30490517 rs11574938 ITGAL 'ENSG00000005844;A3:16:30484219-30485286:30484219-30485517:+'
# /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/plotCoverage/sQTLplot_runner.sh 9 125827703 125840831 rs686661 RABGAP1 'ENSG00000011454;SE:9:125832703-125833691:125833780-125835831:+'
# /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/plotCoverage/sQTLplot_runner.sh 7 992665302 99277156 rs776746 CYP3A5 'ENSG00000106258;SE:7:99270302-99270407:99270538-99272156:-'
# /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/plotCoverage/sQTLplot_runner.sh 6 26435876 26449185 rs9467757 BTN3A3 'ENSG00000111801;SE:6:26440876-26443610:26443670-26444185:+'
# /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/plotCoverage/sQTLplot_runner.sh 18 67609768 26454185 rs1823778 CD226 'ENSG00000150637;A3:18:67614768-67623991:67614669-67623991:-'


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################						END  V3					 								################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
















































































max(RESobs$pval[RESobs$FDR<0.05]) # 2.476815e-07
#save(RES,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_perm.Rdata',EVO_IMMUNO_POP))

RESobs_FDR=RESobs[which(RESobs$FDR_5kb <0.05),]
max(RESobs$pval[which(RESobs$FDR_5kb<0.05)]) # 3.863972e-05
#resGO=GOSeq(RESobs$gene[RESobs$FDR<0.05],unique(PSI_Annot$gene_id),addGenes=T)


Require('Map_imputed')
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO.Rdata',EVO_IMMUNO_POP))

mm=match(RESobs$snps,Map_imputed$snp.name)
RESobs$pos=Map_imputed$position[mm]
RESobs$chrom=Map_imputed$chrom[mm]
RESobs$locus=Map_imputed$locus[mm]
RESobs$daf_AFB=Map_imputed$daf_char_AFB[mm]
RESobs$daf_EUB=Map_imputed$daf_char_EUB[mm]
RESobs$SNPfreq=Map_imputed$SNPfreq[mm]
mm=match(RESobs$event_id,PSI_Annot$event_id)
RESobs$event_start=as.numeric(PSI_Annot$start[mm])
RESobs$event_end=as.numeric(PSI_Annot$end[mm])
RESobs$CisDist=ifelse(RESobs$pos<RESobs$event_start,RESobs$pos-RESobs$event_start,ifelse(RESobs$pos<RESobs$event_end,0,RESobs$pos-RESobs$event_end))
RESobs=RESobs[which(!is.na(RESobs$pos)),]
max(RESobs$CisDist)
RESobs$gene=PSI_Annot$gene_id[mm]
RESobs$symbol=PSI_Annot$symbol[mm]
RESobs$event_type=PSI_Annot$event_type[mm]
RESobs$toKeep=RESobs$event_id%in% PSI_Annot$event_id[toKeep]

NbAssocEvent_1=sapply(thresholds,function(th){luq(RESobs$event_id[!is.na(RESobs$p.k) & RESobs$pval<th])})
NbAssocGene_1=sapply(thresholds,function(th){luq(RESobs$gene[!is.na(RESobs$p.k) & RESobs$pval<th])})

RES_OAS1=RESobs[which(RESobs$symbol=='OAS1'),]
RES_NADSYN1=RESobs[which(RESobs$symbol=='NADSYN1'),]
RES_CTSH=RESobs[which(RESobs$symbol=='CTSH'),]
RES_CAST=RESobs[which(RESobs$symbol=='CAST'),]
save(RES_OAS1,RES_NADSYN1,RES_CTSH,RES_CAST,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_targetGenes.Rdata',EVO_IMMUNO_POP))

RESobs=RESobs[which(RESobs$toKeep),]
save(RESobs,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_V2.Rdata',EVO_IMMUNO_POP))

load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_targetGenes.Rdata',EVO_IMMUNO_POP))
NbAssocEvent_2=sapply(thresholds,function(th){luq(RESobs$event_id[!is.na(RESobs$p.k) & RESobs$pval<th & RESobs$toKeep])})
NbAssocGene_2=sapply(thresholds,function(th){luq(RESobs$gene[!is.na(RESobs$p.k) & RESobs$pval<th & RESobs$toKeep])})
NbAssocEvent_3=sapply(thresholds,function(th){luq(RESobs$event_id[!is.na(RESobs$p.k) & RESobs$pval<th & RESobs$toKeep & RESobs$SNPfreq>0.05 ])})
NbAssocGene_3=sapply(thresholds,function(th){luq(RESobs$gene[!is.na(RESobs$p.k) & RESobs$pval<th & RESobs$toKeep & RESobs$SNPfreq>0.05 ])})
NbAssocEvent_4=sapply(thresholds,function(th){luq(RESobs$event_id[!is.na(RESobs$p.k) & RESobs$pval<th & RESobs$toKeep & RESobs$SNPfreq>0.05 & RESobs$CisDist<1e5])})
NbAssocGene_4=sapply(thresholds,function(th){luq(RESobs$gene[!is.na(RESobs$p.k) & RESobs$pval<th & RESobs$toKeep & RESobs$SNPfreq>0.05 & RESobs$CisDist<1e5])})
save(NbAssocEvent_1,NbAssocEvent_2,NbAssocEvent_3,NbAssocEvent_4,NbAssocGene_1,NbAssocGene_2,NbAssocGene_3,NbAssocGene_4,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_V2_Counts.Rdata',EVO_IMMUNO_POP))

load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_perm.Rdata',EVO_IMMUNO_POP))

mm=match(RES$snps,Map_imputed$snp.name)
RES$pos=Map_imputed$position[mm]
RES$chrom=Map_imputed$chrom[mm]
RES$locus=Map_imputed$locus[mm]
RES$daf_AFB=Map_imputed$daf_char_AFB[mm]
RES$daf_EUB=Map_imputed$daf_char_EUB[mm]
RES$SNPfreq=Map_imputed$SNPfreq[mm]
mm=match(RES$event_id,PSI_Annot$event_id)
RES$event_start=as.numeric(PSI_Annot$start[mm])
RES$event_end=as.numeric(PSI_Annot$end[mm])
RES$CisDist=ifelse(RES$pos<RES$event_start,RES$pos-RES$event_start,ifelse(RES$pos<RES$event_end,0,RES$pos-RES$event_end))
RES=RES[which(!is.na(RES$pos)),]
RES$toKeep=RES$event_id%in% PSI_Annot$event_id[toKeep]

NbAssocEventPerm_1=sapply(thresholds,function(th){luq(RES$event_id[!is.na(RES$p.k) & RES$pval<th])})
NbAssocGenePerm_1=sapply(thresholds,function(th){luq(RES$gene[!is.na(RES$p.k) & RES$pval<th])})

RES=RES[which(RES$toKeep),]
save(RES,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_perm_V2.Rdata',EVO_IMMUNO_POP))

NbAssocEventPerm_2=sapply(thresholds,function(th){luq(RES$event_id[!is.na(RES$p.k) & RES$pval<th & RES$toKeep])})
NbAssocGenePerm_2=sapply(thresholds,function(th){luq(RES$gene[!is.na(RES$p.k) & RES$pval<th & RES$toKeep])})
NbAssocEventPerm_3=sapply(thresholds,function(th){luq(RES$event_id[!is.na(RES$p.k) & RES$pval<th & RES$toKeep & RES$SNPfreq>0.05 ])})
NbAssocGenePerm_3=sapply(thresholds,function(th){luq(RES$gene[!is.na(RES$p.k) & RES$pval<th & RES$toKeep & RES$SNPfreq>0.05 ])})
NbAssocEventPerm_4=sapply(thresholds,function(th){luq(RES$event_id[!is.na(RES$p.k) & RES$pval<th & RES$toKeep & RES$SNPfreq>0.05 & RES$CisDist<1e5])})
NbAssocGenePerm_4=sapply(thresholds,function(th){luq(RES$gene[!is.na(RES$p.k) & RES$pval<th & RES$toKeep & RES$SNPfreq>0.05 & RES$CisDist<1e5])})

save(NbAssocEventPerm_1,NbAssocEventPerm_2,NbAssocEventPerm_3,NbAssocEventPerm_4,NbAssocGenePerm_1,NbAssocGenePerm_2,NbAssocGenePerm_3,NbAssocGenePerm_4,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_perm_V2_Counts.Rdata',EVO_IMMUNO_POP))

#ndup=which(!duplicated(substr(RES$event_id,20,40)))
#event_start=gsub('([0-9XYMT]+):([0-9]+)-?:?(.*)','\\2',substr(RES$event_id[ndup],20,40))
#RES$event_start=as.numeric(event_start[match(substr(RES$event_id,20,40),substr(RES$event_id[ndup],20,40))])
#Nchar=nchar(RES$event_id)
#ndup=which(!duplicated(substr(RES$event_id,Nchar-20,Nchar)))
#event_end=gsub('.*[-:]([0123456789]+)$','\\1',substr(RES$event_id[ndup],Nchar[ndup]-20,Nchar[ndup]-2))
#RES$event_end=as.numeric(event_end[match(substr(RES$event_id,Nchar-20,Nchar-2),substr(RES$event_id[ndup],Nchar[ndup]-20,Nchar[ndup]-2))])

Require('Map_imputed')
#load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO.Rdata',EVO_IMMUNO_POP))
#load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_perm.Rdata',EVO_IMMUNO_POP))
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_perm_V2.Rdata',EVO_IMMUNO_POP))
thresholds=10^-c(seq(3,13),15,20,30,50)
NbAssocGenePerm=list()
NbAssocGenePerm[['1Mb']]=sapply(thresholds,function(th){luq(RES$gene[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<1e6 & RESobs$SNPfreq>0.05)])})
NbAssocGenePerm[['500kb']]=sapply(thresholds,function(th){luq(RES$gene[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<5e5 & RESobs$SNPfreq>0.05)])})
NbAssocGenePerm[['100kb']]=sapply(thresholds,function(th){luq(RES$gene[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<1e5 & RESobs$SNPfreq>0.05)])})
NbAssocGenePerm[['50kb']]=sapply(thresholds,function(th){luq(RES$gene[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<5e4 & RESobs$SNPfreq>0.05)])})
NbAssocGenePerm[['10kb']]=sapply(thresholds,function(th){luq(RES$gene[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<1e4 & RESobs$SNPfreq>0.05)])})
NbAssocGenePerm[['5kb']]=sapply(thresholds,function(th){luq(RES$gene[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<5e3 & RESobs$SNPfreq>0.05)])})
NbAssocGenePerm[['5kbTo500kb']]=sapply(thresholds,function(th){luq(RES$gene[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<1e6 & abs(RES$CisDist)>5e3 & RESobs$SNPfreq>0.05)])})

NbAssocEventPerm=list()
NbAssocEventPerm[['1Mb']]=sapply(thresholds,function(th){luq(RES$event_id[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<1e6 & RESobs$SNPfreq>0.05)])})
NbAssocEventPerm[['500kb']]=sapply(thresholds,function(th){luq(RES$event_id[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<5e5 & RESobs$SNPfreq>0.05)])})
NbAssocEventPerm[['100kb']]=sapply(thresholds,function(th){luq(RES$event_id[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<1e5 & RESobs$SNPfreq>0.05)])})
NbAssocEventPerm[['50kb']]=sapply(thresholds,function(th){luq(RES$event_id[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<5e4 & RESobs$SNPfreq>0.05)])})
NbAssocEventPerm[['10kb']]=sapply(thresholds,function(th){luq(RES$event_id[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<1e4 & RESobs$SNPfreq>0.05)])})
NbAssocEventPerm[['5kb']]=sapply(thresholds,function(th){luq(RES$event_id[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<5e3 & RESobs$SNPfreq>0.05)])})
NbAssocEventPerm[['5kbTo500kb']]=sapply(thresholds,function(th){luq(RES$event_id[which(!is.na(RES$p.k) & RES$pval<th & abs(RES$CisDist)<1e6 & abs(RES$CisDist)>5e3 & RESobs$SNPfreq>0.05)])})

rm(RES);gc()


load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_V2.Rdata',EVO_IMMUNO_POP))
thresholds=10^-c(seq(3,13),15,20,30,50)
NbAssocGene=list()
NbAssocGene[['1Mb']]=sapply(thresholds,function(th){luq(RESobs$gene[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e6 & RESobs$SNPfreq>0.05)])})
NbAssocGene[['500kb']]=sapply(thresholds,function(th){luq(RESobs$gene[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<5e5 & RESobs$SNPfreq>0.05)])})
NbAssocGene[['100kb']]=sapply(thresholds,function(th){luq(RESobs$gene[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e5 & RESobs$SNPfreq>0.05)])})
NbAssocGene[['50kb']]=sapply(thresholds,function(th){luq(RESobs$gene[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<5e4 & RESobs$SNPfreq>0.05)])})
NbAssocGene[['10kb']]=sapply(thresholds,function(th){luq(RESobs$gene[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e4 & RESobs$SNPfreq>0.05)])})
NbAssocGene[['5kb']]=sapply(thresholds,function(th){luq(RESobs$gene[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<5e3 & RESobs$SNPfreq>0.05)])})
NbAssocGene[['5kbTo1Mb']]=sapply(thresholds,function(th){luq(RESobs$gene[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e6 & abs(RESobs$CisDist)>5e3 & RESobs$SNPfreq>0.05)])})
NbAssocGene[['5kbTo500kb']]=sapply(thresholds,function(th){luq(RESobs$gene[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e6 & abs(RESobs$CisDist)>5e3 & RESobs$SNPfreq>0.05)])})

NbAssocEvent=list()
NbAssocEvent[['1Mb']]=sapply(thresholds,function(th){luq(RESobs$event_id[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e6 & RESobs$SNPfreq>0.05)])})
NbAssocEvent[['500kb']]=sapply(thresholds,function(th){luq(RESobs$event_id[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<5e5 & RESobs$SNPfreq>0.05)])})
NbAssocEvent[['100kb']]=sapply(thresholds,function(th){luq(RESobs$event_id[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e5 & RESobs$SNPfreq>0.05)])})
NbAssocEvent[['50kb']]=sapply(thresholds,function(th){luq(RESobs$event_id[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<5e4 & RESobs$SNPfreq>0.05)])})
NbAssocEvent[['10kb']]=sapply(thresholds,function(th){luq(RESobs$event_id[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e4 & RESobs$SNPfreq>0.05)])})
NbAssocEvent[['5kb']]=sapply(thresholds,function(th){luq(RESobs$event_id[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<5e3 & RESobs$SNPfreq>0.05)])})
NbAssocEvent[['5kbTo1Mb']]=sapply(thresholds,function(th){luq(RESobs$event_id[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e6 & abs(RESobs$CisDist)>5e3 & RESobs$SNPfreq>0.05)])})
NbAssocEvent[['5kbTo500kb']]=sapply(thresholds,function(th){luq(RESobs$event_id[which(!is.na(RESobs$p.k) & RESobs$pval<th & abs(RESobs$CisDist)<1e6 & abs(RESobs$CisDist)>5e3 & RESobs$SNPfreq>0.05)])})

save(NbAssocGenePerm,NbAssocEventPerm,NbAssocGene,NbAssocEvent,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_v2_countsDist.Rdata',EVO_IMMUNO_POP))

NbAssocGene=NbAssocGene[c("1Mb","500kb","100kb","50kb","10kb","5kb","5kbTo500kb","5kbTo1Mb")]
NbAssocGenePerm=NbAssocGenePerm[c("1Mb","500kb","100kb","50kb","10kb","5kb","5kbTo500kb","5kbTo1Mb")]
NbAssocEvent=NbAssocEvent[c("1Mb","500kb","100kb","50kb","10kb","5kb","5kbTo500kb","5kbTo1Mb")]
NbAssocEventPerm=NbAssocEventPerm[c("1Mb","500kb","100kb","50kb","10kb","5kb","5kbTo500kb","5kbTo1Mb")]

library(RColorBrewer)
library(DescTools)
acol= c(rev(brewer.pal(4,'YlOrRd')),brewer.pal(4,'YlGnBu'))
acol=acol[1:3,6:8]

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/Nb_detected_sQTLs_gene_MISO.pdf',HOME))
thresholdsDetail=10^-c(seq(3,13,0.1),15,20,30,50)
extrapolate=function(x){approx(-log10(thresholds),x,-log10(thresholdsDetail))$y}
plot(-log10(thresholds),NbAssocGene[[1]],col='#00000000',axes=F,xlim=c(0,50),ylab='Nb Genes with sQTL')
for (i in 1:6){
	lines(-log10(thresholds),NbAssocGene[[i]],col=acol[i],lwd=2)
	lines(-log10(thresholds),NbAssocGenePerm[[i]],col=acol[i],lwd=2,lty=2)
}
axis(2,las=2)
axis(1,at=c(3,5,7,10,20,30,40,50),labels=c(3,5,7,10,20,30,40,50))
for (i in 1:6){
	xFDR10=min(which(extrapolate(NbAssocGenePerm[[i]]/NbAssocGene[[i]])<0.05))
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(NbAssocGene[[i]])[xFDR10],col=acol[i],cex=1.3,pch=16)	
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(NbAssocGene[[i]])[xFDR10],cex=1.3,pch=1)	
}
legend(12,2500,lty=1,lwd=2,legend=c("1Mb","500kb","100kb","50kb","10kb","5kb"),col=acol,bty='n')
legend(12,3500,lty=1:2,lwd=2,legend=c("observed","false positives (estimated)"),col='grey',bty='n')
legend(14,3150,pch=16,legend=c("detected at 5%FDR"),col='grey',bty='n')


thresholdsDetail=10^-c(seq(3,13,0.1),15,20,30,50)
extrapolate=function(x){approx(-log10(thresholds),x,-log10(thresholdsDetail))$y}
plot(-log10(thresholds),NbAssocGene[[1]],col='#00000000',axes=F,xlim=c(0,50),ylab='Nb Genes with sQTL')
for (i in 1:6){
	lines(-log10(thresholds),NbAssocGene[[i]],col=acol[i],lwd=2)
	lines(-log10(thresholds),NbAssocGene[[i]]-NbAssocGenePerm[[i]],col=acol[i],lwd=2,lty=2)
}
for (i in 1:6){
	xFDR10=min(which(extrapolate(NbAssocGenePerm[[i]]/NbAssocGene[[i]])<0.05))
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(NbAssocGene[[i]]-NbAssocGenePerm[[i]])[xFDR10],col=acol[i],cex=1.3,pch=16)	
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(NbAssocGene[[i]]-NbAssocGenePerm[[i]])[xFDR10],cex=1.3,pch=1)	
}
axis(2,las=2)
axis(1,at=c(3,5,7,10,20,30,40,50),labels=c(3,5,7,10,20,30,40,50))
legend(12,2500,lty=1,lwd=2,legend=c("1Mb","500kb","100kb","50kb","10kb","5kb"),col=acol,bty='n')
legend(12,3500,lty=1:2,lwd=2,legend=c("observed","true positives (estimated)"),col='grey',bty='n')
legend(14,3150,pch=16,legend=c("true positives at 5%FDR"),col='grey',bty='n')
dev.off()

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/Nb_detected_sQTLs_events_MISO.pdf',HOME))
thresholdsDetail=10^-c(seq(3,13,0.1),15,20,30,50)
extrapolate=function(x){approx(-log10(thresholds),x,-log10(thresholdsDetail))$y}
plot(-log10(thresholds),NbAssocEvent[[1]],col='#00000000',axes=F,xlim=c(0,50),ylab='Nb Genes with sQTL')
for (i in 1:6){
	lines(-log10(thresholds),NbAssocEvent[[i]],col=acol[i],lwd=2)
	lines(-log10(thresholds),NbAssocEventPerm[[i]],col=acol[i],lwd=2,lty=2)
}
axis(2,las=2)
axis(1,at=c(3,5,7,10,20,30,40,50),labels=c(3,5,7,10,20,30,40,50))
for (i in 1:6){
	xFDR10=min(which(extrapolate(NbAssocEventPerm[[i]]/NbAssocEvent[[i]])<0.05))
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(NbAssocEvent[[i]])[xFDR10],col=acol[i],cex=1.3,pch=16)	
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(NbAssocEvent[[i]])[xFDR10],cex=1.3,pch=1)	
}
legend(12,4000,lty=1,lwd=2,legend=c("1Mb","500kb","100kb","50kb","10kb","5kb"),col=acol,bty='n')
legend(12,5500,lty=1:2,lwd=2,legend=c("observed","false positives (estimated)"),col='grey',bty='n')
legend(14.5,4900,pch=16,legend=c("detected at 5%FDR"),col='grey',bty='n')


thresholdsDetail=10^-c(seq(3,13,0.1),15,20,30,50)
extrapolate=function(x){approx(-log10(thresholds),x,-log10(thresholdsDetail))$y}
plot(-log10(thresholds),NbAssocEvent[[1]],col='#00000000',axes=F,xlim=c(0,50),ylab='Nb Genes with sQTL')
for (i in 1:6){
	lines(-log10(thresholds),NbAssocEvent[[i]],col=acol[i],lwd=2)
	lines(-log10(thresholds),NbAssocEvent[[i]]-NbAssocEventPerm[[i]],col=acol[i],lwd=2,lty=2)
}
for (i in 1:6){
	xFDR10=min(which(extrapolate(NbAssocEventPerm[[i]]/NbAssocEvent[[i]])<0.05))
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(NbAssocEvent[[i]]-NbAssocEventPerm[[i]])[xFDR10],col=acol[i],cex=1.3,pch=16)	
	points(-log10(thresholdsDetail)[xFDR10],extrapolate(NbAssocEvent[[i]]-NbAssocEventPerm[[i]])[xFDR10],cex=1.3,pch=1)	
}
axis(2,las=2)
axis(1,at=c(3,5,7,10,20,30,40,50),labels=c(3,5,7,10,20,30,40,50))
legend(12,4000,lty=1,lwd=2,legend=c("1Mb","500kb","100kb","50kb","10kb","5kb"),col=acol,bty='n')
legend(12,5500,lty=1:2,lwd=2,legend=c("observed","true positives (estimated)"),col='grey',bty='n')
legend(14.5,4900,pch=16,legend=c("true positives at 5%FDR"),col='grey',bty='n')
dev.off()

FDR_compute=function(pval,FP,TP){y=approxfun(c(0,-log10(thresholds),500),c(pmin(-log10(c(1,FP/TP)),6),500))(-log10(pval)); 10^-y}

RESobs$FDR_5kb=NA
RESobs$FDR_100kb=NA
RESobs$FDR_1Mb=NA
RESobs$FDR_5kb_500kb=NA
RESobs$FDR_5kb[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e3 & RESobs$SNPfreq>0.05 )]=FDR_compute(RESobs$pval[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e3 & RESobs$SNPfreq>0.05)],NbAssocEventPerm[['5kb']],NbAssocEvent[['5kb']])
RESobs$FDR_100kb[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e4 & RESobs$SNPfreq>0.05)]=FDR_compute(RESobs$pval[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e4 & RESobs$SNPfreq>0.05)],NbAssocEventPerm[['50kb']],NbAssocEvent[['50kb']])
RESobs$FDR_1Mb[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e5 & RESobs$SNPfreq>0.05)]=FDR_compute(RESobs$pval[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e5 & RESobs$SNPfreq>0.05)],NbAssocEventPerm[['500kb']],NbAssocEvent[['500kb']])
RESobs$FDR_5kb_500kb[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e5 & RESobs$SNPfreq>0.05)]=FDR_compute(RESobs$pval[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e5 & RESobs$SNPfreq>0.05)],NbAssocEventPerm[['5kbTo500kb']],NbAssocEvent[['5kbTo500kb']])
RESobs$FDR_5kb_500kb[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e3 & RESobs$SNPfreq>0.05)]=FDR_compute(RESobs$pval[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e3 & RESobs$SNPfreq>0.05)],NbAssocEventPerm[['5kb']],NbAssocEvent[['5kb']])
save(RESobs,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_V2.Rdata',EVO_IMMUNO_POP))

Testable_all=PSI_Annot$event_id[PSI_Annot$NbTestableCond==5 & PSI_Annot$JuncCovered &  PSI_Annot$noOverlap]
RESobs_FDR=RESobs[which(RESobs$FDR_5kb <0.05),]
RESobs_FDR=RESobs_FDR[order(RESobs_FDR$event_id,RESobs_FDR$cond,RESobs_FDR$pvalue),]
RESobs_FDR_nodup=RESobs_FDR[!duplicated(paste(RESobs_FDR$event_id,RESobs_FDR$cond)),]
RESobs_FDR_nodup=RESobs_FDR_nodup[order(RESobs_FDR_nodup$pvalue),]
RESobs_hasSQTL=RESobs[which(!is.na(RESobs$p.k) & abs(RESobs$CisDist)<5e3 & RESobs$SNPfreq>0.05 & paste(RESobs$event_id,RESobs$snps)%in%paste(RESobs_FDR_nodup$event_id,RESobs_FDR_nodup$snps)),]
Detected_sQTL_P5pct=By(RESobs_hasSQTL$cond,paste(RESobs_hasSQTL$event_id, RESobs_hasSQTL$snps),function(x){paste(sort(x),collapse='')})
RESobs_FDR_nodup$detected_cond=Detected_sQTL_P5pct[match(paste(RESobs_FDR_nodup$event_id,RESobs_FDR_nodup$snps),names(Detected_sQTL_P5pct))]
RESobs_FDR_nodup[which(RESobs_FDR_nodup$detected_cond==234 & RESobs_FDR_nodup$event_id%in%Testable_all)[1:3],]
save(RESobs_FDR_nodup,RESobs_hasSQTL,RESobs_FDR,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_V2_detail.Rdata',EVO_IMMUNO_POP))

# get sQTL_1Mb
RESobs_FDR_1Mb=RESobs[which(RESobs$FDR_1Mb <0.05),]
RESobs_FDR_1Mb=RESobs_FDR_1Mb[order(RESobs_FDR_1Mb$event_id,RESobs_FDR_1Mb$cond,RESobs_FDR_1Mb$pvalue),]
RESobs_FDR_1Mb_nodup=RESobs_FDR_1Mb[!duplicated(paste(RESobs_FDR_1Mb$event_id,RESobs_FDR_1Mb$cond)),]
RESobs_FDR_1Mb_nodup=RESobs_FDR_1Mb_nodup[order(RESobs_FDR_1Mb_nodup$pvalue),]
tab5kb=table(substr(setdiff(RESobs_FDR_nodup$event_id, RESobs_FDR_1Mb_nodup$event_id),17,18))
# A3  A5  AF  AL  MX  RI  SE 
# 58  43  94  72   6  66 185 
sum(tab5kb)
# Total 524

tab1Mb=table(substr(setdiff(RESobs_FDR_1Mb_nodup$event_id, RESobs_FDR_nodup$event_id),17,18))
#A3 A5 AF AL MX RI SE 
# 3  4  3  6  1  5 10 
sum(tab1Mb)
# Total 32
chisq.test(cbind(tab1Mb,tab5kb)) # p=0.60 pas de différences

tab_NbCond=table(nchar(RESobs_FDR_nodup$detected_cond[!duplicated(RESobs_FDR_nodup$event_id)]))
tab_NbCond_AllTestable=table(nchar(RESobs_FDR_nodup$detected_cond[RESobs_FDR_nodup$event_id%in%Testable_all & !duplicated(RESobs_FDR_nodup$event_id)]))

#table(RESobs_FDR_nodup$detected_cond[!duplicated(RESobs_FDR_nodup$event_id)]))


pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/sQTL_sharing_across_Conditions.pdf',HOME),width=4)
X=cbind(all=rev(tab_NbCond),shared=rev(tab_NbCond_AllTestable))
barplot(X,ylab='NbEvents with sQTL',col=acol,las=1,ylim=c(0,max(apply(X,2,sum))*1.1))
text(0:1*1.2+0.7,apply(X,2,sum)+max(X)*0.05,labels=apply(X,2,sum))
legend("topright",bty='n',paste(1:5,'cond'),fill=rev(acol[1:5]))
barplot(t(t(X)/apply(X,2,sum)),ylab='%s Events with sQTL',col=acol,las=1,ylim=c(0,1.1))
text(0:1*1.2+0.7,c(1,1)*1.05,labels=paste('N =',apply(X,2,sum)))
dev.off()

tab_type=table(RESobs_FDR[!duplicated(RESobs_FDR$event_id),'event_type'])
tab_type_tested=table(PSI_Annot$event_type[PSI_Annot$noOverlap & PSI_Annot$NbTestableCond>0 ])

barplot(tab_type/tab_type_tested,beside=Tcol=colPSI,las=1)
OR=NULL
for (i in 1:length(colPSI)){
	y=c(tab_type_tested[i]-tab_type[i],tab_type[i])
	X=rbind(apply(cbind(tab_type_tested-tab_type,tab_type),2,sum)-y,y)
	OR=rbind(OR,odds.ratio(X))
}
rownames(OR)=names(tab_type_tested)

write.table(OR,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/OR_byType.txt',HOME),sep='\t',quote=F,row.names=F)

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/sQTL_by_event_type.pdf',HOME),height=4)
plot(c(0,2.5),c(0,7),col='#00000000',axes=F,ylab='',xlab='Odds ratio')
axis(1)
axis(2,at=1:7,names(colPSI),las=2)
abline(v=1,col='grey')
points(OR$OR,1:7,pch=16,cex=1.3,col=colPSI)
for(i in 1:7){segments(OR$LowerCI[i],i,OR$UpperCI[i],i,lwd=2,col=colPSI[i])}

barplot(tab_type_tested,beside=T,col=colPSI,las=1,ylim=c(0,max(tab_type_tested))*1.15,ylab='Nb event')
barplot(tab_type,beside=T,col='#00000088',add=T,axes=F)
text(0:6*1.2+0.7,tab_type_tested+0.1*max(tab_type_tested),labels=tab_type_tested)
text(0:6*1.2+0.7,tab_type+0.05*max(tab_type_tested),labels=tab_type)
barplot(tab_type/tab_type_tested,beside=T,col=colPSI,las=1,ylab='% of events with a splicing QTL')
abline(h=0.1444589,lty=2,col='lightgrey')
dev.off()

w=which(PSI_Annot$JuncCovered & PSI_Annot$noOverlap & PSI_Annot$NbTestableCond>0)

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/SignalToNoise_by_event_type.pdf',HOME),height=7)
boxplot(log2(PSI_Annot$Signal2Noise_global[w])~PSI_Annot$event_type[w],col=colPSI,notch=T,pch=16,outcol=colERC5[1],ylab='signal-to-noise',las=1)
abline(h=median(log2(PSI_Annot$Signal2Noise_global[w])),col='grey',lty=2)
dev.off()

write.table(RESobs_FDR_nodup,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/sQTL_full_list.txt',HOME),sep='\t',quote=F,row.names=F)
	# overlap with sQTL obtained by leafcutter ?
	# leafcutter specific events ?
	# what is new ?
w=which(PSI_Annot$JuncCovered & PSI_Annot$noOverlap & PSI_Annot$NbTestableCond>0)
resGO=GOSeq(unique(RESobs_FDR_nodup$gene),unique(PSI_Annot$gene[w]))

resGO=GOSeq(unique(RESobs_FDR_nodup$gene),unique(PSI_Annot$gene[w]))
resGO$cond='all'

resGOByCond=lapply(1:5,function(i){GOSeq(unique(RESobs_FDR_nodup$gene[RESobs_FDR_nodup$cond==i]),unique(PSI_Annot$gene[which(PSI_Annot$JuncCovered & PSI_Annot$noOverlap & PSI_Annot[,grep('Testable',colnames(PSI_Annot))[i]])]))})

lapply(resGOByCond,d)
resGOByCond=do.call(rbind,sapply(1:5,function(i){cbind(resGOByCond[[i]],cond=rep(i,nrow(resGOByCond[[i]])))}))

###################################################
###			overlap with eQTL by type	 V2		###
###################################################

Require('RES_cis_Full')
cis_eQTL=RES_cis_Full[which(RES_cis_Full$pop=='ALL' & RES_cis_Full$pval<1e-7 & RES_cis_Full$type=='eQTL'),]
cis_eQTL=cis_eQTL[order(cis_eQTL$gene,cis_eQTL$cond,cis_eQTL$pval),]
cis_eQTL=cis_eQTL[which(!duplicated(paste(cis_eQTL$gene,cis_eQTL$cond))),]
cis_eQTL=cis_eQTL[order(cis_eQTL$pval),]
rm(RES_cis_Full)

RESobs_FDR_nodup=RESobs_FDR_nodup[order(RESobs_FDR_nodup$event_id,RESobs_FDR_nodup$pval),]
RESobs_FDR_nodupCond=RESobs_FDR_nodup[!duplicated(RESobs_FDR_nodup$event_id),]

geneInCommon=intersect(paste(cis_eQTL$gene,cis_eQTL$cond),paste(RESobs_FDR_nodupCond$gene,RESobs_FDR_nodupCond$cond))
has_eQTL=paste(RESobs_FDR_nodupCond$gene,RESobs_FDR_nodupCond$cond)%in%paste(cis_eQTL$gene,cis_eQTL$cond)

geneInCommon_eQTL=cis_eQTL$snps[match(paste(RESobs_FDR_nodupCond$gene,RESobs_FDR_nodupCond$cond)[has_eQTL],paste(cis_eQTL$gene,cis_eQTL$cond))]
geneInCommon_sQTL=RESobs_FDR_nodupCond$snps[paste(RESobs_FDR_nodupCond$gene,RESobs_FDR_nodupCond$cond)%in%geneInCommon]
geneInCommon_type=RESobs_FDR_nodupCond$event_type[paste(RESobs_FDR_nodupCond$gene,RESobs_FDR_nodupCond$cond)%in%geneInCommon]
luq(geneInCommon) # 286
luq(paste(RESobs_FDR_nodupCond$gene,RESobs_FDR_nodupCond$cond)) # 1038
luq(paste(cis_eQTL$gene,cis_eQTL$cond)) # 8395

#704/2593: 27%  also have an eQTL
#286/1038: 27%  also have an eQTL

By(paste(RESobs_FDR_nodupCond$gene,RESobs_FDR_nodupCond$cond)%in%paste(cis_eQTL$gene,cis_eQTL$cond),RESobs_FDR_nodupCond$event_type,mean)

getGenotypes=function( snp_list){
	SNPs=getGenos(snp_list)
	snp.names=SNPs[,1]
	SNPs=as.matrix(SNPs[-(1:5)])
	rownames(SNPs)=snp.names
	SNPs=SNPs[match(snp_list,snp.names),order(colnames(SNPs))]
}
geneInCommon_eQTL_geno=getGenotypes(geneInCommon_eQTL)
geneInCommon_sQTL_geno=getGenotypes(geneInCommon_sQTL)

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/SplicingQTL_PSI_eQTL_Corr_bestSNP_V2.pdf',HOME),height=5,width=5)
r2_sQTL_eQTL=diag(cor(t(geneInCommon_eQTL_geno),t(geneInCommon_sQTL_geno),use='p')^2)
hist(r2_sQTL_eQTL,xlab='r2',main='r2, best sQTL/eQTL',br=50,col='grey')
Pct_HighLD=c()
for(i in 1:5){
	hist(r2_sQTL_eQTL[grep(paste('',i),geneInCommon)],xlab='r2',main='r2, best sQTL/eQTL',col=colERC5[i],breaks=50)
	Pct_HighLD[i]=mean(r2_sQTL_eQTL[grep(paste('',i),geneInCommon)]>=0.8)
}
names(Pct_HighLD)=condIndex
barplot(Pct_HighLD,col=colERC5,ylab='% of genes with r2>0.8 between sQTL and eQTL')

Pct_highLD_byType=By(r2_sQTL_eQTL>0.8,geneInCommon_type, mean)



has_eQTL_LD=has_eQTL
has_eQTL_LD[has_eQTL]=r2_sQTL_eQTL[match(geneInCommon_sQTL,rn(geneInCommon_sQTL_geno))]>0.8
LD_eQTL=rep(0,length(has_eQTL))
LD_eQTL[has_eQTL]=r2_sQTL_eQTL[match(geneInCommon_sQTL,rn(geneInCommon_sQTL_geno))]

sum(has_eQTL_LD) # 101
sum(r2_sQTL_eQTL>0.8) # 101
sum(LD_eQTL>0.8) # 101

tab_has_eQTL_LD_type=table(has_eQTL_LD, RESobs_FDR_nodupCond$event_type)


OR=NULL
for (i in 1:length(colPSI)){
	y=tab_has_eQTL_LD_type[,i]
	X=rbind(apply(cbind(tab_has_eQTL_LD_type[,-i]),1,sum),y)
	OR=rbind(OR,odds.ratio(X))
}
rownames(OR)=colnames(tab_has_eQTL_LD_type)
write.table(OR,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/OR_eQTL_seQTL_overlap_byType_V2.txt',HOME),sep='\t',quote=F,row.names=F)

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/Overlpa_sQTL_eQTL_by_event_type_V2.pdf',HOME),height=4)
barplot( By(has_eQTL_LD, RESobs_FDR_nodupCond$event_type,mean),col=colPSI,ylim=c(0,0.2),las=1,ylab='% of sQTL overlapping an eQTL')
abline(h=mean(has_eQTL_LD),col='grey',lty=2)
plot(c(0,2.5),c(1,7),col='#00000000',axes=F,ylab='',xlab='Odds ratio')
axis(1)
axis(2,at=1:7,names(colPSI),las=2)
abline(v=1)
points(OR$OR,1:7,pch=16,cex=1.3,col=colPSI)
for(i in 1:7){segments(OR$LowerCI[i],i,OR$UpperCI[i],i,lwd=2,col=colPSI[i])}
barplot(table(RESobs_FDR_nodup$event_type),beside=T,col=colPSI,las=1,ylim=c(0,max(table(RESobs_FDR_nodup$event_type)))*1.15,ylab='Nb sQTL')
barplot(tab_has_eQTL_LD_type[2,],beside=T,col='#00000088',add=T,axes=F)
text(0:6*1.2+0.7,table(RESobs_FDR_nodup$event_type)+0.1*max(table(RESobs_FDR_nodup$event_type)),labels=table(RESobs_FDR_nodup$event_type))
text(0:6*1.2+0.7,tab_has_eQTL_LD_type[2,]+0.05*max(table(RESobs_FDR_nodup$event_type)),labels=tab_has_eQTL_LD_type[2,])
dev.off()


pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/Overlpa_sQTL_eQTL_by_event_type_v2.pdf',HOME),height=4,width=4)
barplot( By(LD_eQTL>0.8, RESobs_FDR_nodupCond$event_type,mean),col=colPSI,ylim=c(0,0.2),las=1,ylab='% of sQTL overlapping an eQTL')
abline(h=mean(LD_eQTL>0.8),col='grey',lty=2)
hist(LD_eQTL,br=30,col=grey(0.3),las=1,xlab='LD with eQTL peak SNP')
hist(LD_eQTL[has_eQTL],br=30,col=grey(0.8),las=1,add=T)
legend('topright',fill=grey(c(0.3,0.8)),legend=c('gene has no eQTL','gene has eQTL'),bty='n')
dev.off()






###########################################################
##		Test enrichment in FST/iHS/ GWAS / aSNP			###
###########################################################

Require('local_eQTL_list')
boxplot(abs(Map_imputed$signedFST_EUB[unique(mm)]),abs(Map_imputed$signedFST_EUB[mm_EUB]))
mm=match(RESobs_FDR_nodup$snps,Map_imputed$snp.name)
mm_EUB=match(keep_EUB,Map_imputed$snp.name)
mm_AFB=match(keep_AFB,Map_imputed$snp.name)

RESobs_FDR_nodup$iHS_EUB=Map_imputed$signedIHS_EUB[mm]
RESobs_FDR_nodup$iHS_AFB=Map_imputed$signedIHS_AFB[mm]
RESobs_FDR_nodup$FST=Map_imputed$FST_adj[mm]
RESobs_FDR_nodup$Max_FST_R2_EUB=Map_imputed$Max_FST_R2_EUB[mm]
RESobs_FDR_nodup$Max_FST_R2_AFB=Map_imputed$Max_FST_R2_AFB[mm]

RESobs_FDR_nodup$P_CLS2_EUB=Map_imputed$P_CLS2_EUB[mm]
RESobs_FDR_nodup$P_CLS2_AFB=Map_imputed$P_CLS2_AFB[mm]
RESobs_FDR_nodup$aSNP=Map_imputed$aSNP[mm]
RESobs_FDR_nodup$aSNP_R2_EUB=Map_imputed$aSNP_R2_EUB[mm]
RESobs_FDR_nodup$GWAS_Trait_R2_EUB=Map_imputed$GWAS_Trait_R2_EUB[mm]
RESobs_FDR_nodup$GWAS_Trait_R2_AFB=Map_imputed$GWAS_Trait_R2_AFB[mm]
RESobs_FDR_nodup$iHS_EUB_R2=Map_imputed$Max_signedIHS_R2_EUB[mm]
RESobs_FDR_nodup$aiHS_EUB_R2=Map_imputed$Max_IHS_R2_EUB[mm]
RESobs_FDR_nodup$iHS_AFB_R2=Map_imputed$Max_signedIHS_R2_AFB[mm]
RESobs_FDR_nodup$aiHS_AFB_R2=Map_imputed$Max_IHS_R2_AFB[mm]
RESobs_FDR_nodup$maf_EUB= Map_imputed$maf_EUB[mm]
RESobs_FDR_nodup$maf_AFB= Map_imputed$maf_AFB[mm]

Q99_FST=quantile(Map_imputed$FST_adj,0.99)
Q99_IHS_EUB=quantile(-Map_imputed$signedIHS_EUB[mm_EUB],0.99,na.rm=T)
Q99_IHS_AFB=quantile(-Map_imputed$signedIHS_AFB[mm_AFB],0.99,na.rm=T)

#quantile(Map_imputed$FST_adj[Map_imputed$daf_EUB>0.05],0.99)


###########################
##		FIRST ROUND		###
###########################

# get target MAF distribution
mmSNPfreq=Map_imputed$SNPfreq[mm] 
mmSNPfreq_bin=round(mmSNPfreq/2.5,2)*2.5 
SNPfreq_bin=round(Map_imputed$SNPfreq[Map_imputed$SNPfreq>0.05]/2.5,2)*2.5
binCounts=table(mmSNPfreq_bin) # Nb of SNP to sample for each MAF bin 

sum(Map_imputed[sample_global,'GWAS_Trait_R2_EUB']!='',na.rm=T)

# defined SNP_support (set of sampled SNPs)
SNP_Support=Map_imputed$SNPfreq>0.05
NSAMP=200
# get statistics for matched set of SNPs
meanFST=meanFST_overQ99=meanIHS_EUB=meanIHS_AFB=meanIHS_EUB_overQ99=meanIHS_AFB_overQ99=NbaSNP_R2_EUB=NbGWAS_SNP_EUB=c()
for (samp in 1:NSAMP){
	tic=Sys.time()
	cat(samp)
	sample_global=By(which(SNP_support),SNPfreq_bin,sample,max(binCounts)) # over sample each bin
	for (i in 1:length(sample_global)){
		sample_global[[i]]=sample_global[[i]][1:binCounts[i]] # sub sample to right amount
		}
	sample_global=unlist(sample_global)
	meanFST[samp]=mean(Map_imputed[sample_global,'FST_adj'],na.rm=T)
	meanFST_overQ99[samp]=mean(Map_imputed[sample_global,'FST_adj']>Q99_FST,na.rm=T)
	meanIHS_EUB[samp]=mean(Map_imputed[sample_global,'signedIHS_EUB'],na.rm=T)
	meanIHS_AFB[samp]=mean(Map_imputed[sample_global,'signedIHS_AFB'],na.rm=T)
	meanIHS_EUB_overQ99[samp]=mean(-Map_imputed[sample_global,'signedIHS_EUB']>Q99_IHS_EUB,na.rm=T)
	meanIHS_AFB_overQ99[samp]=mean(-Map_imputed[sample_global,'signedIHS_AFB']>Q99_IHS_AFB,na.rm=T)
	NbaSNP_R2_EUB[samp]=sum(Map_imputed[sample_global,'aSNP_R2_EUB'],na.rm=T)
	NbGWAS_SNP_EUB[samp]=sum(Map_imputed[sample_global,'GWAS_Trait_R2_EUB']!='',na.rm=T)
	toc=Sys.time()
	print(toc-tic)
	flush.console()
	}

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/Enrichments_GWAS_aSNP_Selection_FirstRound.pdf',HOME))
hist(NbGWAS_SNP_EUB,xlim=c(80,300),col='grey',las=1)
points(sum(RESobs_FDR_nodup$GWAS_Trait_R2_EUB!='',na.rm=T),0,pch=16,col='red')

hist(NbaSNP_R2_EUB,xlim=c(20,100),col='grey',las=1)
points(sum(RESobs_FDR_nodup$aSNP_R2_EUB,na.rm=T),0,pch=16,col='red')

hist(meanIHS_EUB,col='grey',las=1)
points(mean(RESobs_FDR_nodup$iHS_EUB,na.rm=T),0,pch=16,col='red')

hist(meanIHS_AFB,col='grey',las=1)
points(mean(RESobs_FDR_nodup$iHS_AFB,na.rm=T),0,pch=16,col='red')

hist(meanIHS_EUB_overQ99,col='grey',las=1)
points(mean(RESobs_FDR_nodup$iHS_EUB>Q99_IHS_EUB,na.rm=T),0,pch=16,col='red')

hist(meanIHS_AFB_overQ99,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_AFB>Q99_IHS_AFB,na.rm=T),0,pch=16,col='red')

hist(meanFST_overQ99,col='grey',las=1)
points(mean(RESobs_FDR_nodup$FST>Q99_FST,na.rm=T),0,pch=16,col='red')
dev.off()


NSAMP=1000
Resamp=list()
###########################
##		EUB ONLY		###
###########################

for (cond in 1:5){
mm=match(RESobs_FDR_nodup$snps[RESobs_FDR_nodup$cond==cond & RESobs_FDR_nodup$maf_EUB>0.05],Map_imputed$snp.name)
mmSNPfreq=Map_imputed$SNPfreq[mm] 
mmSNPfreq_bin=round(mmSNPfreq/2.5,2)*2.5 

#mmLD=cut(Map_imputed$NbLD_SNP_EUB[mm],c(0,1,2,5,10,20,50,100,1e6)))
#mmSNPfreq_bin=round(mmSNPfreq/2.5,2)*2.5 
binCounts=table(mmSNPfreq_bin) # Nb of SNP to sample for each MAF bin 

# defined SNP_support (set of sampled SNPs)
SNP_Support=which(Map_imputed$SNPfreq>0.05 & Map_imputed$maf_EUB>0.05 & Map_imputed$snp.name%in%keep_EUB)
SNPfreq_bin=round(Map_imputed$SNPfreq[SNP_Support]/2.5,2)*2.5 #cut(Map_imputed$NbLD_SNP_EUB[SNP_Support],c(0,1,2,5,10,20,50,100,1e6)))

Q95_FST=quantile(Map_imputed$FST_adj[SNP_Support],0.95)
Q95_FST_R2_EUB=quantile(Map_imputed$Max_FST_R2_EUB[SNP_Support],0.95)
Q95_IHS_EUB=quantile(-Map_imputed$signedIHS_EUB[SNP_Support],0.95,na.rm=T)
Q95_IHS_R2_EUB=quantile(-Map_imputed$Max_signedIHS_R2_EUB[SNP_Support],0.95,na.rm=T)
Q95_aIHS_EUB=quantile(abs(Map_imputed$signedIHS_EUB[SNP_Support]),0.95,na.rm=T)
Q95_aIHS_R2_EUB=quantile(Map_imputed$Max_IHS_R2_EUB[SNP_Support],0.95,na.rm=T)

Q99_FST=quantile(Map_imputed$FST_adj[SNP_Support],0.99)
Q99_FST_R2_EUB=quantile(Map_imputed$Max_FST_R2_EUB[SNP_Support],0.99)
Q99_IHS_EUB=quantile(-Map_imputed$signedIHS_EUB[SNP_Support],0.99,na.rm=T)
Q99_IHS_R2_EUB=quantile(-Map_imputed$Max_signedIHS_R2_EUB[SNP_Support],0.99,na.rm=T)
Q99_aIHS_EUB=quantile(abs(Map_imputed$signedIHS_EUB[SNP_Support]),0.99,na.rm=T)
Q99_aIHS_R2_EUB=quantile(Map_imputed$Max_IHS_R2_EUB[SNP_Support],0.99,na.rm=T)


# get statistics for matched set of SNPs
if(cond==1){Resamp[['EUB']]=list()}
Resamp[['EUB']][[cond]]=list()
Resamp[['EUB']][[cond]]$meanFST=c()
Resamp[['EUB']][[cond]]$meanFST_overQ95=c()
Resamp[['EUB']][[cond]]$meanFST_overQ99=c()
Resamp[['EUB']][[cond]]$meanIHS_EUB=c()
Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ95=c()
Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ99=c()
Resamp[['EUB']][[cond]]$meanaIHS_EUB=c()
Resamp[['EUB']][[cond]]$meanaIHS_EUB_overQ95=c()
Resamp[['EUB']][[cond]]$meanaIHS_EUB_overQ99=c()

Resamp[['EUB']][[cond]]$meanFST_R2=c()
Resamp[['EUB']][[cond]]$meanFST_R2_overQ95=c()
Resamp[['EUB']][[cond]]$meanFST_R2_overQ99=c()
Resamp[['EUB']][[cond]]$meanIHS_EUB_R2=c()
Resamp[['EUB']][[cond]]$meanIHS_EUB_R2_overQ95=c()
Resamp[['EUB']][[cond]]$meanIHS_EUB_R2_overQ99=c()
Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2=c()
Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2_overQ95=c()
Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2_overQ99=c()

Resamp[['EUB']][[cond]]$NbaSNP_R2_EUB=c()
Resamp[['EUB']][[cond]]$NbGWAS_SNP_EUB=c()

for (samp in 1:NSAMP){
	tic=Sys.time()
	cat(samp)
	sample_global=By(SNP_Support,SNPfreq_bin,sample,max(binCounts)) # over sample each bin
	for (i in 1:length(sample_global)){
		sample_global[[i]]=sample_global[[i]][1:binCounts[i]] # sub sample to right amount
		}

	sample_global=unlist(sample_global)
	Resamp[['EUB']][[cond]]$meanFST[samp]=mean(Map_imputed[sample_global,'FST_adj'],na.rm=T)
	Resamp[['EUB']][[cond]]$meanFST_overQ95[samp]=mean(Map_imputed[sample_global,'FST_adj']>Q95_FST,na.rm=T)
	Resamp[['EUB']][[cond]]$meanFST_overQ99[samp]=mean(Map_imputed[sample_global,'FST_adj']>Q99_FST,na.rm=T)
	Resamp[['EUB']][[cond]]$meanFST_R2[samp]=mean(Map_imputed[sample_global,'Max_FST_R2_EUB'],na.rm=T)
	Resamp[['EUB']][[cond]]$meanFST_R2_overQ95[samp]=mean(Map_imputed[sample_global,'Max_FST_R2_EUB']>Q95_FST,na.rm=T)
	Resamp[['EUB']][[cond]]$meanFST_R2_overQ99[samp]=mean(Map_imputed[sample_global,'Max_FST_R2_EUB']>Q99_FST,na.rm=T)
	Resamp[['EUB']][[cond]]$meanIHS_EUB[samp]=mean(Map_imputed[sample_global,'signedIHS_EUB'],na.rm=T)
	Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ95[samp]=mean(-Map_imputed[sample_global,'signedIHS_EUB']>Q95_IHS_EUB,na.rm=T)
	Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ99[samp]=mean(-Map_imputed[sample_global,'signedIHS_EUB']>Q99_IHS_EUB,na.rm=T)
	Resamp[['EUB']][[cond]]$meanIHS_EUB_R2[samp]=mean(Map_imputed[sample_global,'Max_signedIHS_R2_EUB'],na.rm=T)
	Resamp[['EUB']][[cond]]$meanIHS_EUB_R2_overQ95[samp]=mean(-Map_imputed[sample_global,'Max_signedIHS_R2_EUB']>Q95_IHS_R2_EUB,na.rm=T)
	Resamp[['EUB']][[cond]]$meanIHS_EUB_R2_overQ99[samp]=mean(-Map_imputed[sample_global,'Max_signedIHS_R2_EUB']>Q99_IHS_R2_EUB,na.rm=T)
	Resamp[['EUB']][[cond]]$meanaIHS_EUB[samp]=mean(abs(Map_imputed[sample_global,'signedIHS_EUB']),na.rm=T)
	Resamp[['EUB']][[cond]]$meanaIHS_EUB_overQ95[samp]=mean(abs(Map_imputed[sample_global,'signedIHS_EUB'])>Q95_aIHS_EUB,na.rm=T)
	Resamp[['EUB']][[cond]]$meanaIHS_EUB_overQ99[samp]=mean(abs(Map_imputed[sample_global,'signedIHS_EUB'])>Q99_aIHS_EUB,na.rm=T)
	Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2[samp]=mean(abs(Map_imputed[sample_global,'Max_IHS_R2_EUB']),na.rm=T)
	Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2_overQ95[samp]=mean(abs(Map_imputed[sample_global,'Max_IHS_R2_EUB'])>Q95_aIHS_R2_EUB,na.rm=T)
	Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2_overQ99[samp]=mean(abs(Map_imputed[sample_global,'Max_IHS_R2_EUB'])>Q99_aIHS_R2_EUB,na.rm=T)
	Resamp[['EUB']][[cond]]$NbaSNP_R2_EUB[samp]=sum(Map_imputed[sample_global,'aSNP_R2_EUB'],na.rm=T)
	Resamp[['EUB']][[cond]]$NbGWAS_SNP_EUB[samp]=sum(Map_imputed[sample_global,'GWAS_Trait_R2_EUB']!='',na.rm=T)
	toc=Sys.time()
	print(toc-tic)
	flush.console()
	}
}	

###########################
##		AFB ONLY		###
###########################
for (cond in 1:5){
mm=match(RESobs_FDR_nodup$snps[RESobs_FDR_nodup$cond==cond & RESobs_FDR_nodup$maf_EUB>0.05],Map_imputed$snp.name)
mmSNPfreq=Map_imputed$SNPfreq[mm] 
mmSNPfreq_bin=round(mmSNPfreq/2.5,2)*2.5 

#mmLD=cut(Map_imputed$NbLD_SNP_AFB[mm],c(0,1,2,5,10,20,50,100,1e6)))
#mmSNPfreq_bin=round(mmSNPfreq/2.5,2)*2.5 
binCounts=table(mmSNPfreq_bin) # Nb of SNP to sample for each MAF bin 

# defined SNP_support (set of sampled SNPs)
SNP_Support=which(Map_imputed$SNPfreq>0.05 & Map_imputed$maf_AFB>0.05 & Map_imputed$snp.name%in%keep_AFB)
SNPfreq_bin=round(Map_imputed$SNPfreq[SNP_Support]/2.5,2)*2.5 #cut(Map_imputed$NbLD_SNP_AFB[SNP_Support],c(0,1,2,5,10,20,50,100,1e6)))

Q95_FST=quantile(Map_imputed$FST_adj[SNP_Support],0.95)
Q95_FST_R2_AFB=quantile(Map_imputed$Max_FST_R2_AFB[SNP_Support],0.95)
Q95_IHS_AFB=quantile(-Map_imputed$signedIHS_AFB[SNP_Support],0.95,na.rm=T)
Q95_IHS_R2_AFB=quantile(-Map_imputed$Max_signedIHS_R2_AFB[SNP_Support],0.95,na.rm=T)
Q95_aIHS_AFB=quantile(abs(Map_imputed$signedIHS_AFB[SNP_Support]),0.95,na.rm=T)
Q95_aIHS_R2_AFB=quantile(Map_imputed$Max_IHS_R2_AFB[SNP_Support],0.95,na.rm=T)

Q99_FST=quantile(Map_imputed$FST_adj[SNP_Support],0.99)
Q99_FST_R2_AFB=quantile(Map_imputed$Max_FST_R2_AFB[SNP_Support],0.99)
Q99_IHS_AFB=quantile(-Map_imputed$signedIHS_AFB[SNP_Support],0.99,na.rm=T)
Q99_IHS_R2_AFB=quantile(-Map_imputed$Max_signedIHS_R2_AFB[SNP_Support],0.99,na.rm=T)
Q99_aIHS_AFB=quantile(abs(Map_imputed$signedIHS_AFB[SNP_Support]),0.99,na.rm=T)
Q99_aIHS_R2_AFB=quantile(Map_imputed$Max_IHS_R2_AFB[SNP_Support],0.99,na.rm=T)

# get statistics for matched set of SNPs
if(cond==1){Resamp[['AFB']]=list()}

Resamp[['AFB']][[cond]]=list()
Resamp[['AFB']][[cond]]$meanFST=c()
Resamp[['AFB']][[cond]]$meanFST_overQ95=c()
Resamp[['AFB']][[cond]]$meanFST_overQ99=c()
Resamp[['AFB']][[cond]]$meanIHS_AFB=c()
Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ95=c()
Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ99=c()
Resamp[['AFB']][[cond]]$meanaIHS_AFB=c()
Resamp[['AFB']][[cond]]$meanaIHS_AFB_overQ95=c()
Resamp[['AFB']][[cond]]$meanaIHS_AFB_overQ99=c()

Resamp[['AFB']][[cond]]$meanFST_R2=c()
Resamp[['AFB']][[cond]]$meanFST_R2_overQ95=c()
Resamp[['AFB']][[cond]]$meanFST_R2_overQ99=c()
Resamp[['AFB']][[cond]]$meanIHS_AFB_R2=c()
Resamp[['AFB']][[cond]]$meanIHS_AFB_R2_overQ95=c()
Resamp[['AFB']][[cond]]$meanIHS_AFB_R2_overQ99=c()
Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2=c()
Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2_overQ95=c()
Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2_overQ99=c()

Resamp[['AFB']][[cond]]$NbGWAS_SNP_AFB=c()

for (samp in 1:NSAMP){
	tic=Sys.time()
	cat(samp)
	sample_global=By(SNP_Support,SNPfreq_bin,sample,max(binCounts)) # over sample each bin
	for (i in 1:length(sample_global)){
		sample_global[[i]]=sample_global[[i]][1:binCounts[i]] # sub sample to right amount
		}

	sample_global=unlist(sample_global)
	Resamp[['AFB']][[cond]]$meanFST[samp]=mean(Map_imputed[sample_global,'FST_adj'],na.rm=T)
	Resamp[['AFB']][[cond]]$meanFST_overQ95[samp]=mean(Map_imputed[sample_global,'FST_adj']>Q95_FST,na.rm=T)
	Resamp[['AFB']][[cond]]$meanFST_overQ99[samp]=mean(Map_imputed[sample_global,'FST_adj']>Q99_FST,na.rm=T)
	Resamp[['AFB']][[cond]]$meanFST_R2[samp]=mean(Map_imputed[sample_global,'Max_FST_R2_AFB'],na.rm=T)
	Resamp[['AFB']][[cond]]$meanFST_R2_overQ95[samp]=mean(Map_imputed[sample_global,'Max_FST_R2_AFB']>Q95_FST,na.rm=T)
	Resamp[['AFB']][[cond]]$meanFST_R2_overQ99[samp]=mean(Map_imputed[sample_global,'Max_FST_R2_AFB']>Q99_FST,na.rm=T)
	Resamp[['AFB']][[cond]]$meanIHS_AFB[samp]=mean(Map_imputed[sample_global,'signedIHS_AFB'],na.rm=T)
	Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ95[samp]=mean(-Map_imputed[sample_global,'signedIHS_AFB']>Q95_IHS_AFB,na.rm=T)
	Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ99[samp]=mean(-Map_imputed[sample_global,'signedIHS_AFB']>Q99_IHS_AFB,na.rm=T)
	Resamp[['AFB']][[cond]]$meanIHS_AFB_R2[samp]=mean(Map_imputed[sample_global,'Max_signedIHS_R2_AFB'],na.rm=T)
	Resamp[['AFB']][[cond]]$meanIHS_AFB_R2_overQ95[samp]=mean(-Map_imputed[sample_global,'Max_signedIHS_R2_AFB']>Q95_IHS_R2_AFB,na.rm=T)
	Resamp[['AFB']][[cond]]$meanIHS_AFB_R2_overQ99[samp]=mean(-Map_imputed[sample_global,'Max_signedIHS_R2_AFB']>Q99_IHS_R2_AFB,na.rm=T)
	Resamp[['AFB']][[cond]]$meanaIHS_AFB[samp]=mean(abs(Map_imputed[sample_global,'signedIHS_AFB']),na.rm=T)
	Resamp[['AFB']][[cond]]$meanaIHS_AFB_overQ95[samp]=mean(abs(Map_imputed[sample_global,'signedIHS_AFB'])>Q95_aIHS_AFB,na.rm=T)
	Resamp[['AFB']][[cond]]$meanaIHS_AFB_overQ99[samp]=mean(abs(Map_imputed[sample_global,'signedIHS_AFB'])>Q99_aIHS_AFB,na.rm=T)
	Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2[samp]=mean(abs(Map_imputed[sample_global,'Max_IHS_R2_AFB']),na.rm=T)
	Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2_overQ95[samp]=mean(abs(Map_imputed[sample_global,'Max_IHS_R2_AFB'])>Q95_aIHS_R2_AFB,na.rm=T)
	Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2_overQ99[samp]=mean(abs(Map_imputed[sample_global,'Max_IHS_R2_AFB'])>Q99_aIHS_R2_AFB,na.rm=T)
	Resamp[['AFB']][[cond]]$NbGWAS_SNP_AFB[samp]=sum(Map_imputed[sample_global,'GWAS_Trait_R2_AFB']!='',na.rm=T)
	toc=Sys.time()
	print(toc-tic)
	flush.console()
	}
}

save(Resamp,RESobs_FDR_nodup,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/Enrichments_GWAS_aSNP_Selection_%s.Rdata',HOME,cond))

y=RESobs_FDR_nodup[order(RESobs_FDR_nodup$pval),]
y[which(y$aSNP_R2_EUB & !duplicated(y$event_id)),]
MakePlots_EUB=function(cond){
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/Enrichments_GWAS_aSNP_Selection_EUB_only_%s.pdf',HOME,cond))
wEUB=which(RESobs_FDR_nodup$maf_EUB>0.05 & RESobs_FDR_nodup$cond==cond)
hist(Resamp[['EUB']][[cond]]$NbGWAS_SNP_EUB,xlim=c(0,80),col='grey',las=1)
points(sum(RESobs_FDR_nodup$GWAS_Trait_R2_EUB[wEUB]!='',na.rm=T),0,pch=16,col=colERC5[cond])

#points(sum(RESobs_FDR_nodup$GWAS_Trait_R2_EUB[]!='',na.rm=T),0,pch=16,col='red')

hist(Resamp[['EUB']][[cond]]$NbaSNP_R2_EUB,xlim=c(0,80),col='grey',las=1)
points(sum(RESobs_FDR_nodup$aSNP_R2_EUB[wEUB],na.rm=T),0,pch=16,col=colERC5[cond])
barplot(table(Resamp[['EUB']][[cond]]$NbaSNP_R2_EUB),xlim=c(0,0.7+10*1.2),las=1)
points(sum(RESobs_FDR_nodup$aSNP_R2_EUB[wEUB],na.rm=T)*1.2+0.7,0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB,col='grey',las=1)
points(mean(RESobs_FDR_nodup$iHS_EUB[wEUB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ95,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_EUB[wEUB]>Q95_IHS_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ99,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_EUB[wEUB]>Q99_IHS_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_R2,col='grey',las=1,xlim=c(-1,0))
points(mean(RESobs_FDR_nodup$iHS_EUB_R2[wEUB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_R2_overQ95,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_EUB_R2[wEUB]>Q95_IHS_R2_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_R2_overQ99,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_EUB_R2[wEUB]>Q99_IHS_R2_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanaIHS_EUB,col='grey',las=1)
points(mean(abs(RESobs_FDR_nodup$iHS_EUB[wEUB]),na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ95,col='grey',las=1)
points(mean(abs(RESobs_FDR_nodup$iHS_EUB[wEUB])>Q95_aIHS_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ99,col='grey',las=1)
points(mean(abs(RESobs_FDR_nodup$iHS_EUB[wEUB])>Q99_aIHS_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2,col='grey',las=1,xlim=c(0.95,1.5))
points(mean(RESobs_FDR_nodup$aiHS_EUB_R2[wEUB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2_overQ95,col='grey',las=1,xlim=c(0,0.15))
points(mean(RESobs_FDR_nodup$aiHS_EUB_R2[wEUB]>Q95_aIHS_R2_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanaIHS_EUB_R2_overQ99,col='grey',las=1)
points(mean(RESobs_FDR_nodup$aiHS_EUB_R2[wEUB]>Q99_aIHS_R2_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB,col='grey',las=1)
points(mean(RESobs_FDR_nodup$iHS_EUB[wEUB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ95,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_EUB[wEUB]>Q95_IHS_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanIHS_EUB_overQ99,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_EUB[wEUB]>Q99_IHS_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanFST,col='grey',las=1,xlim=c(0,1))
points(mean(RESobs_FDR_nodup$FST[wEUB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanFST_R2,col='grey',las=1,xlim=c(0,1))
points(mean(RESobs_FDR_nodup$Max_FST_R2_EUB[wEUB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanFST_R2_overQ95,col='grey',las=1,xlim=c(0,0.15))
points(mean(RESobs_FDR_nodup$Max_FST_R2_EUB[wEUB]>Q95_FST_R2_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['EUB']][[cond]]$meanFST_R2_overQ99,col='grey',las=1,xlim=c(0,0.1))
points(mean(RESobs_FDR_nodup$Max_FST_R2_EUB[wEUB]>Q99_FST_R2_EUB,na.rm=T),0,pch=16,col=colERC5[cond])

dev.off()
}

MakePlots_AFB=function(cond){
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/Enrichments_GWAS_aSNP_Selection_AFB_only_%s.pdf',HOME,cond))
wAFB=which(RESobs_FDR_nodup$maf_AFB>0.05 & RESobs_FDR_nodup$cond==cond)
hist(Resamp[['AFB']][[cond]]$NbGWAS_SNP_AFB,xlim=c(0,80),col='grey',las=1)
points(sum(RESobs_FDR_nodup$GWAS_Trait_R2_AFB[wAFB]!='',na.rm=T),0,pch=16,col=colERC5[cond])

#points(sum(RESobs_FDR_nodup$GWAS_Trait_R2_AFB[]!='',na.rm=T),0,pch=16,col='red')

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB,col='grey',las=1)
points(mean(RESobs_FDR_nodup$iHS_AFB[wAFB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ95,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_AFB[wAFB]>Q95_IHS_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ99,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_AFB[wAFB]>Q99_IHS_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_R2,col='grey',las=1,xlim=c(-1,0))
points(mean(RESobs_FDR_nodup$iHS_AFB_R2[wAFB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_R2_overQ95,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_AFB_R2[wAFB]>Q95_IHS_R2_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_R2_overQ99,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_AFB_R2[wAFB]>Q99_IHS_R2_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanaIHS_AFB,col='grey',las=1)
points(mean(abs(RESobs_FDR_nodup$iHS_AFB[wAFB]),na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ95,col='grey',las=1)
points(mean(abs(RESobs_FDR_nodup$iHS_AFB[wAFB])>Q95_aIHS_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ99,col='grey',las=1)
points(mean(abs(RESobs_FDR_nodup$iHS_AFB[wAFB])>Q99_aIHS_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2,col='grey',las=1,xlim=c(0.95,1.5))
points(mean(RESobs_FDR_nodup$aiHS_AFB_R2[wAFB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2_overQ95,col='grey',las=1,xlim=c(0,0.15))
points(mean(RESobs_FDR_nodup$aiHS_AFB_R2[wAFB]>Q95_aIHS_R2_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanaIHS_AFB_R2_overQ99,col='grey',las=1)
points(mean(RESobs_FDR_nodup$aiHS_AFB_R2[wAFB]>Q99_aIHS_R2_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB,col='grey',las=1)
points(mean(RESobs_FDR_nodup$iHS_AFB[wAFB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ95,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_AFB[wAFB]>Q95_IHS_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanIHS_AFB_overQ99,col='grey',las=1)
points(mean(-RESobs_FDR_nodup$iHS_AFB[wAFB]>Q99_IHS_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanFST,col='grey',las=1,xlim=c(0,1))
points(mean(RESobs_FDR_nodup$FST[wAFB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanFST_R2,col='grey',las=1,xlim=c(0,1))
points(mean(RESobs_FDR_nodup$Max_FST_R2_AFB[wAFB],na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanFST_R2_overQ95,col='grey',las=1,xlim=c(0,0.15))
points(mean(RESobs_FDR_nodup$Max_FST_R2_AFB[wAFB]>Q95_FST_R2_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

hist(Resamp[['AFB']][[cond]]$meanFST_R2_overQ99,col='grey',las=1,xlim=c(0,0.1))
points(mean(RESobs_FDR_nodup$Max_FST_R2_AFB[wAFB]>Q99_FST_R2_AFB,na.rm=T),0,pch=16,col=colERC5[cond])

dev.off()
}

for(cond in 1:5){
MakePlots_EUB(cond)
MakePlots_AFB(cond)
}


# get target MAF distribution
mm=match(RESobs_FDR_nodup$snps,Map_imputed$snp.name)
mmSNPfreq=(Map_imputed$maf_AFB[mm]+Map_imputed$maf_EUB[mm])/2
mmSNPfreq_bin=round(mmSNPfreq/2.5,2)*2.5 
SNPfreq_bin=round((Map_imputed$maf_AFB+Map_imputed$maf_EUB)[Map_imputed$SNPfreq>0.05]/2/2.5,2)*2.5
binCounts=table(mmSNPfreq_bin) # Nb of SNP to sample for each MAF bin 

sum(Map_imputed[sample_global,'GWAS_Trait_R2_EUB']!='',na.rm=T)


NSAMP=1000

Resamp[['ALL']]=list()

# defined SNP_support (set of sampled SNPs)
SNP_Support=Map_imputed$SNPfreq>0.05

Q95_FST=quantile(Map_imputed$FST_adj[SNP_Support],0.95)
Q95_IHS_EUB=quantile(-Map_imputed$signedIHS_EUB[SNP_Support],0.95,na.rm=T)
Q95_aIHS_EUB=quantile(abs(Map_imputed$signedIHS_EUB[SNP_Support]),0.95,na.rm=T)

Q99_FST=quantile(Map_imputed$FST_adj[SNP_Support],0.99)
Q99_IHS_EUB=quantile(-Map_imputed$signedIHS_EUB[SNP_Support],0.99,na.rm=T)
Q99_aIHS_EUB=quantile(abs(Map_imputed$signedIHS_EUB[SNP_Support]),0.99,na.rm=T)

library(GenomicRanges)
GR_SNP=Make
GR_Exons=

# get statistics for matched set of SNPs
if(cond==1){Resamp[['ALL']]=list()}
Resamp[['ALL']][[cond]]=list()
Resamp[['ALL']][[cond]]$meanFST=c()
Resamp[['ALL']][[cond]]$meanFST_overQ95=c()
Resamp[['ALL']][[cond]]$meanFST_overQ99=c()
Resamp[['ALL']][[cond]]$meanIHS_EUB=c()
Resamp[['ALL']][[cond]]$meanIHS_EUB_overQ95=c()
Resamp[['ALL']][[cond]]$meanIHS_EUB_overQ99=c()
Resamp[['ALL']][[cond]]$meanaIHS_EUB=c()
Resamp[['ALL']][[cond]]$meanaIHS_EUB_overQ95=c()
Resamp[['ALL']][[cond]]$meanaIHS_EUB_overQ99=c()

Resamp[['ALL']][[cond]]$meanFST_R2=c()
Resamp[['ALL']][[cond]]$meanFST_R2_overQ95=c()
Resamp[['ALL']][[cond]]$meanFST_R2_overQ99=c()
Resamp[['ALL']][[cond]]$meanIHS_EUB_R2=c()
Resamp[['ALL']][[cond]]$meanIHS_EUB_R2_overQ95=c()
Resamp[['ALL']][[cond]]$meanIHS_EUB_R2_overQ99=c()
Resamp[['ALL']][[cond]]$meanaIHS_EUB_R2=c()
Resamp[['ALL']][[cond]]$meanaIHS_EUB_R2_overQ95=c()
Resamp[['ALL']][[cond]]$meanaIHS_EUB_R2_overQ99=c()

Resamp[['ALL']][[cond]]$NbaSNP_R2_EUB=c()
Resamp[['ALL']][[cond]]$NbGWAS_SNP_EUB=c()



# get statistics for matched set of SNPs
for (samp in 1:NSAMP){
	tic=Sys.time()
	cat(samp)
	sample_global=By(which(SNP_support),SNPfreq_bin,sample,max(binCounts)) # over sample each bin
	for (i in 1:length(sample_global)){
		sample_global[[i]]=sample_global[[i]][1:binCounts[i]] # sub sample to right amount
		}
	sample_global=unlist(sample_global)
	meanFST[samp]=mean(Map_imputed[sample_global,'FST_adj'],na.rm=T)
	meanFST_overQ99[samp]=mean(Map_imputed[sample_global,'FST_adj']>Q99_FST,na.rm=T)
	meanIHS_EUB[samp]=mean(Map_imputed[sample_global,'signedIHS_EUB'],na.rm=T)
	meanIHS_AFB[samp]=mean(Map_imputed[sample_global,'signedIHS_AFB'],na.rm=T)
	meanIHS_EUB_overQ99[samp]=mean(-Map_imputed[sample_global,'signedIHS_EUB']>Q99_IHS_EUB,na.rm=T)
	meanIHS_AFB_overQ99[samp]=mean(-Map_imputed[sample_global,'signedIHS_AFB']>Q99_IHS_AFB,na.rm=T)
	NbaSNP_R2_EUB[samp]=sum(Map_imputed[sample_global,'aSNP_R2_EUB'],na.rm=T)
	NbGWAS_SNP_EUB[samp]=sum(Map_imputed[sample_global,'GWAS_Trait_R2_EUB']!='',na.rm=T)
	toc=Sys.time()
	print(toc-tic)
	flush.console()
	}


dev.off()

######################################
##			leafcutter aSNPs		##
######################################

load(sprintf("%s/Maxime/Splicing/sQTL/leafcutter/Cis-sQTL_ALL_allCond_allChr_nodup.Rdata",EVO_IMMUNO_POP))
mm=match(res$SNP,Map_imputed$snp.name)
res$aSNP=Map_imputed$aSNP[mm]
res$aSNP_R2=Map_imputed$aSNP_R2_EUB[mm]
res$FST=Map_imputed$FST_adj[mm]
res$daf_EUB=Map_imputed$daf_EUB[mm]
res$daf_AFB=Map_imputed$daf_AFB[mm]
res$iHS_AFB=Map_imputed$signedIHS_AFB[mm]
res$iHS_EUB=Map_imputed$signedIHS_EUB[mm]

res=res[order(res$pvalue),]
aSNP_res=res[which(res$aSNP_R2 & res$FDR<0.05),]
aSNP_res=aSNP_res[!duplicated(aSNP_res$clustID),]

write.table(aSNP_res,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/aSNP_res_leafcutter.txt',HOME),quote=F,sep='\t',row.names=F)

RESobs_FDR_nodup

aSNP_res=RESobs_FDR_nodup[which(RESobs_FDR_nodup$aSNP_R2),]
aSNP_res=aSNP_res[!duplicated(aSNP_res$event_id),]
write.table(aSNP_res,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/aSNP_res_MISO.txt',HOME),quote=F,sep='\t',row.names=F)


mm=match(res$SNP,Map_imputed$snp.name)

mm_EUB=match(keep_EUB,Map_imputed$snp.name)
boxplot(Map_imputed$signedIHS_EUB[mm],Map_imputed$signedIHS_EUB[mm_EUB])



load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO.Rdata',EVO_IMMUNO_POP))
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_counts.Rdata',EVO_IMMUNO_POP))
# todo: plot P-values long-range cis sQTL / Short range cis-sQTL.
RESobs_nodup=RESobs[order(RESobs$event_id,RESobs$cond,RESobs$pval),]
RESobs_nodup=RESobs_nodup[which(!duplicated(paste(RESobs_nodup$event_id,RESobs_nodup$cond))),]
RESobs_nodup=RESobs_nodup[which(RESobs_nodup$FDR_5kb_500kb<0.05),]
save(RESobs_nodup,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-pisQTL_MISO_nodup.Rdata',EVO_IMMUNO_POP))

By(RESobs_nodup$gene,RESobs_nodup$event_type,luq)
#A3  A5  AF  AL  MX  RI  SE 
#257 239 392 262 138 154 629 
By(RESobs_nodup$gene,RESobs_nodup$cond,luq)
#  1   2   3   4   5 
#842 695 757 647 609 
By(RESobs_nodup$gene,paste(RESobs_nodup$event_type,abs(RESobs_nodup$CisDist)<5e3),luq)

By(RESobs_nodup$gene,paste(RESobs_nodup$event_type,abs(RESobs_nodup$CisDist)<5e4),luq)
NbGene_tested_byType=By(PSI_Annot[PSI_Annot$NbCondTestable>0,'gene_id'],PSI_Annot$event_type[PSI_Annot$NbCondTestable>0],luq)
#  A3   A5   AF   AL   MX   RI   SE 
#2343 2230 2902 1465  925 1250 3564 
type='AF'

x=c(257,239,392, 262,138,154 ,629) 
y=c(2343,2230,2902,1465,925,1250,3564)

resGO=GOSeq(unique(RESobs_nodup[which(!duplicated(RESobs_nodup$gene)),"gene"]),unique(PSI_Annot[PSI_Annot$NbCondTestable>0,'gene_id']))

resGO=GOSeq(RESobs_nodup[which(!duplicated(RESobs_nodup$gene) & RESobs_nodup$gene),"gene"],PSI_Annot[PSI_Annot$NbCondTestable>0,'gene_id'])

luq(RESobs$gene[which(!is.na(RESobs$p.k) & (RESobs$pval<1e-4 & abs(RESobs$CisDist)<5e3) |  (RESobs$pval<1e-10 & abs(RESobs$CisDist)<5e5))]) # 2073
luq(RESobs$gene[which(!is.na(RESobs$p.k) & (RESobs$pval<1e-4 & abs(RESobs$CisDist)<5e3) |  (RESobs$pval<1e-11 & abs(RESobs$CisDist)<5e5))]) # 2064

NbEvent=By(RESobs$event_id[RESobs$pval<1e-7],RESobs$event_type[RESobs$pval<1e-7],luq)
NbGene=By(RESobs$gene[RESobs$pval<1e-7],RESobs$event_type[RESobs$pval<1e-7],luq)
NbGene_tested_byType=By(PSI_Annot[PSI_Annot$NbCondTestable>0,'gene_id'],PSI_Annot$event_type[PSI_Annot$NbCondTestable>0],luq)
NbGeneByCond=By(RESobs$gene[RESobs$pval<1e-7],RESobs$cond[RESobs$pval<1e-7],luq)
NbEvent_tested=table(PSI_Annot$event_type[PSI_Annot$NbCondTestable>0])
NbGene_tested_byCond=apply(PSI_Annot[,grep('Testable',colnames(PSI_Annot))[1:5]],2,function(x){luq(PSI_Annot$gene_id[x])})	

pdf(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/Figures/Cis_psiQTL_ALL_miso_1e7.pdf',sep=''))
barplot(NbEvent,col=colPSI[names(NbEvent)],ylab='Nb Event with psiQTL')
barplot(NbEvent/NbEvent_tested,col=colPSI[names(NbEvent)],ylab='% Event with psiQTL')
barplot(NbGene,col=colPSI[names(NbGene)],ylab='Nb Gene with psiQTL')
barplot(NbGene/NbGene_tested_byType,col=colPSI[names(NbGene)],ylab='% Gene with psiQTL')
barplot(NbGeneByCond,col=colERC5,ylab='Nb Gene with psiQTL By Condition')
dev.off()

x=RESobs[which(!is.na(RESobs$p.k) & RESobs$FDR<0.05),]
# threshold=1e-7 (5% FDR)
x=x[order(x$gene,x$pop,x$cond,x$pvalue),]
x_nodup=x[!duplicated(paste(x$gene,x$cond,x$pop)),]
x_nodup=x_nodup[order(x_nodup$pvalue),]
luq(RESobs$snps[RESobs$gene%in%x_nodup$gene & RESobs$cond%in%x_nodup$cond & RESobs$pop%in%x_nodup$pop])

Cis_psiQTL=RESobs[paste(RESobs$gene,RESobs$cond,RESobs$pop)%in%paste(x_nodup$gene,x_nodup$cond,x_nodup$pop),]
Cis_psiQTL_1e3=RESobs
Cis_psiQTL_nodup=x_nodup
save(Cis_psiQTL,Cis_psiQTL_1e3,Cis_psiQTL_nodup,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/sQTL/MatrixEQTL-cis/Cis_psiQTL_ALL_miso.Rdata',sep=''))
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/sQTL/MatrixEQTL-cis/Cis_psiQTL_ALL_miso.Rdata',sep=''))
#load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/sQTL/MatrixEQTL-cis/Cis_psiQTL_ALL.Rdata',sep=''))
Cis_psiQTL$beta=-Cis_psiQTL$beta
Cis_psiQTL_1e3$beta=-Cis_psiQTL_1e3$beta
Cis_psiQTL_nodup$beta=-Cis_psiQTL_nodup$beta

x=Cis_psiQTL_1e3[order(Cis_psiQTL_1e3$event_id,Cis_psiQTL_1e3$pop,Cis_psiQTL_1e3$cond,Cis_psiQTL_1e3$pvalue),]
x_nodup=x[!duplicated(paste(x$event_id,x$cond,x$pop)),]
Cis_psiQTL_1e3$R2max=x_nodup$R2[match(paste(Cis_psiQTL_1e3$event_id,Cis_psiQTL_1e3$cond,Cis_psiQTL_1e3$pop),paste(x_nodup$event_id,x_nodup$cond,x_nodup$pop))]
Cis_psiQTL$R2max=x_nodup$R2[match(paste(Cis_psiQTL$event_id,Cis_psiQTL$cond,Cis_psiQTL$pop),paste(x_nodup$event_id,x_nodup$cond,x_nodup$pop))]

Cis_psiQTL$PeakScore=Cis_psiQTL$R2/Cis_psiQTL$R2max
Cis_psiQTL_1e3$PeakScore=Cis_psiQTL_1e3$R2/Cis_psiQTL_1e3$R2max
save(Cis_psiQTL,Cis_psiQTL_1e3,Cis_psiQTL_nodup,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/sQTL/MatrixEQTL-cis/Cis_psiQTL_ALL_v2.Rdata',sep=''))



MeanExpr=log2(1+GeneAnnot[,10:14])
lFC=MeanExpr[,-1]-MeanExpr[,1]%o%rep(1,4)
MeanExpr_log=t(apply(FPKM_gene,1,By,SampleAnnot$cond,mean))
lFC_log=MeanExpr_log[,-1]-MeanExpr_log[,1]%o%rep(1,4)
max_alFC=apply(abs(lFC),1,max)
max_lFC=apply(lFC,1,max)
gene_list=intersect(names(max_lFC),unique(PSI_Annot$gene_id[toKeep]))
odds.ratio(table(max_lFC[gene_list]>1,gene_list%in%RESobs_nodup$gene_id)





################ PIECHART sQTL SHARING
sQTL_sharing_1E7=apply(RESobs_nodup_1Mb[,paste('Pvalue_',condIndex,sep='')]<4.74e-7,1,sum)

sQTL_sharing=apply(RESobs_nodup_1Mb[,paste('Pvalue_',condIndex,sep='')]<1e-3,1,sum)

expressed_sharing=apply(RESobs_nodup_1Mb[,paste('Expressed_',condIndex,sep='')],1,sum)
mean(sQTL_sharing>1) # 0.8784153 
sum(sQTL_sharing>1) # 1286
mean(sQTL_sharing==1) # 0.1215847
sum(sQTL_sharing==1) # 178


mean(sQTL_sharing==5) # 0.4610656
sum(sQTL_sharing==5) # 675

sum(sQTL_sharing<5) # 789
mean(expressed_sharing[sQTL_sharing<5]<5) # 0.4486692
sum(expressed_sharing[sQTL_sharing<5]<5) # 354

mean(sQTL_sharing[expressed_sharing==5]==5) # 0.546875
sum(sQTL_sharing[expressed_sharing==5]==5) # 525

mean(sQTL_sharing[expressed_sharing==5]==1) # 0.103125
sum(sQTL_sharing[expressed_sharing==5]==1) # 99

par(xpd=T)
X=cbind(all=rev(table(sQTL_sharing_1E7)),shared=rev(table(sQTL_sharing[expressed_sharing==5])))
counts=apply(X,2,sum)
X=t(t(X)/counts)
x=barplot(X,ylab='Percentage of sQTLs',col=acol[c(1,3,4,6,7)],las=1,ylim=c(0,max(apply(X,2,sum))*1.1),space=0.5)
text(x,apply(X,2,sum)+max(X)*0.09,labels=counts)
#legend(1.3,1100,bty='n',paste(1:5,'cond'),fill=rev(acol[c(1,3,4,6,7)]))
legend(1.3,1600,bty='n',paste(1:5,'cond'),fill=rev(acol[c(1,3,4,6,7)]))


par(mar=c(3,3,3,3))
tab=rev(table(sQTL_sharing[expressed_sharing==5]))
tabpct=paste(round(100*tab/sum(tab), 1),'%')
pie(tab,col=acol[c(1,3,4,6,7)],init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab,col=acol[c(1,3,4,6,7)],init.angle=90,labels=rep(' ',5)) # 4 x 4 inches


Genos_sQTL=read.table(file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/SNP_genos_sQTL_bestSNPperCondition.txt',sep=''),sep='\t',header=T)
snps_sQTL=as.matrix(Genos_sQTL[-1])
rownames(snps_sQTL)=gsub('.',':',Genos_sQTL[[1]],fixed=T)
snps_sQTL=t(snps_sQTL[match(RESobs_nodup_1Mb$snps,rownames(snps_sQTL)),])

psi_sQTL=mapply(function(event,cond){PSI_prov[match(event,PSI_Annot[toKeep,'event_id']),match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(PSI_prov))]},RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$cond)
gene_sQTL=mapply(function(gene,cond){FPKM_gene[gene,match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(FPKM_gene))]},RESobs_nodup_1Mb$gene,RESobs_nodup_1Mb$cond)


myPSIs=t(PSI_prov[match(RESobs_nodup_1Mb$event_id[expressed_sharing==5],PSI_Annot[toKeep,'event_id']),])
mySNPs=snps_sQTL[,match(RESobs_nodup_1Mb$snps[expressed_sharing==5],colnames(snps_sQTL))]

# we reverse the SNP to encode allele.1 as the major allele
mySNPs_cond=mySNPs[substr(rownames(myPSIs),1,6),]

allModels=cbind(NS=rep(rep(0:1,e=16),1),
				LPS=rep(rep(0:1,e=8),2),
				PAM3CSK4=rep(rep(0:1,e=4),4),
				R848=rep(rep(0:1,e=2),8),
				IAV=rep(rep(0:1,e=1),16))

LL_full=apply(rbind(mySNPs_cond,myPSIs),2,function(SNP_and_PSI){
    n=length(SNP_and_PSI)/2
    mySNP_cond=SNP_and_PSI[1:n]
    myPSI_cond=SNP_and_PSI[n+1:n]
    LL_psi=apply(allModels,1,function(x,myPSI,SNP_cond){
		mySNP=SNP_cond*rep(x,table(SampleAnnot$cond))
		LL=logLik(lm(myPSI~mySNP+SampleAnnot$CondPop))
	},myPSI_cond,mySNP_cond)
	LL_psi
    })

LL_probs=apply(LL_full, 2,function(LL_psi){exp(LL_psi-min(LL_psi))/sum(exp(LL_psi-min(LL_psi)))})
WhichCond=apply(LL_full, 2,function(LL_psi){allModels[which.max(LL_psi),]})
table(apply(WhichCond,2,paste,collapse='')
#00001 00010 00011 00100 00101 00110 00111 01000 01001 01010 01011 01100 01101 01110 01111 10000 10001 10010 10011 10100 10101 10110 10111 11000 11001 11010 11011 11100 11101 11110 11111 
#   26     7     8     3     3     2     3    10     2     1     2     2     1    22    35     5     8     1     4     4     4     1     4     1     4     1    14     4    16    78   684 


    
par(mar=c(3,3,3,3))
tab=rev(table(apply(WhichCond,2,sum)))
tabpct=paste(round(100*tab/sum(tab), 1),'%')
pie(tab,col=acol[c(1,3,4,6,7)],init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab,col=acol[c(1,3,4,6,7)],init.angle=90,labels=rep(' ',5)) # 4 x 4 inches


