#################################################################
###		Differential splicing between Populations : Test 	  ###
#################################################################

Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V5_withCounts.Rdata',sep=''))
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME))

library(impute)
library(VennDiagram)

PSI_cond=PSI_prov
PSI_Annot_cond=PSI_Annot[toKeep,]
rownames(PSI_cond)=PSI_Annot_cond$event_id

Pval=matrix(1,length(toKeep),5)
rownames(Pval)=rn(PSI_Annot)[toKeep]
for (cond in 1:5){
	Pval[,cond]=apply(PSI_cond,1,function(x){P=try(wilcox.test(x[grep(paste('AFB[0-9]+-',cond,sep=''),cn(PSI_cond))],x[grep(paste('EUB[0-9]+-',cond,sep=''),cn(PSI_cond))])$p.value);if(class(P)=='try-error'){NA}else{P}})
	}

FDR=Pval
FDR=matrix(p.adjust(as.numeric(Pval),'fdr'),nrow(Pval),ncol(Pval))

DELTA_PSI=matrix(0,length(toKeep),5)
rownames(DELTA_PSI)=rn(PSI_Annot)[toKeep]
for(cond in 1:5){
	wtest=match(rn(DELTA_PSI),rn(PSI_Annot))
	DELTA_PSI[,cond]=apply(PSI_cond,1,function(x){mean(x[grep(paste('AFB[0-9]+-',cond,sep=''),cn(PSI_cond))])-mean(x[grep(paste('EUB[0-9]+-',cond,sep=''),cn(PSI_cond))])})
}

TESTABLE=as.matrix(PSI_Annot[toKeep,paste('Testable',condIndex,sep='_')])
GeneSymbols=PSI_Annot$symbol[toKeep]

Require('IAV_RNAseq')
pop=grepl(paste('AFB[0-9]+-',cond,sep=''),IAV_RNAseq$ID)
mRNA_IAV=IAV_RNAseq$flu_FPKM
PvalAdjIAV=apply(PSI_cond[,match(IAV_RNAseq$ID,cn(PSI_cond))],1,function(x){P=try(summary(lm(x~pop+mRNA_IAV))$coeff[2,4]);if(class(P)=='try-error'){NA}else{P}})

w=which(TESTABLE,arr=T)
PopDiffSpliceFull=cbind(PSI_Annot_cond[w[,1],],Pval=Pval[w],DeltaPSI=DELTA_PSI[w],cond=w[,2],PvalAdjIAV=cbind(NA,NA,NA,NA,PvalAdjIAV)[w])
PopDiffSpliceFull$DELTA_PSI_R=DELTA_PSI[w]-DELTA_PSI[cbind(w[,1],1)]
PopDiffSpliceFull$FDR=FDR[w]
PopDiffSpliceFull=PopDiffSpliceFull[order(PopDiffSpliceFull$Pval),]
rownames(PopDiffSpliceFull)=NULL
#PopDiffSpliceFull=PopDiffSpliceFull[!duplicated(paste(PopDiffSpliceFull$gene_id,PopDiffSpliceFull$cond)),]

lFC=population.diff[,grep('Beta_.*_EUBvsAFB',cn(population.diff))]
rownames(lFC)=rn(FPKM_gene)
PopDiffSpliceFull$lFC_gene=NA
for(i in 1:5){
	PopDiffSpliceFull$lFC_gene[PopDiffSpliceFull$cond==i]=lFC[match(PopDiffSpliceFull$gene[PopDiffSpliceFull$cond==i],rownames(lFC)),i]
}

MeanExpr=log2(1+GeneAnnot[,10:14])
lFC_cond=MeanExpr[,-1]-MeanExpr[,1]%o%rep(1,4)
MeanExpr_log=t(apply(FPKM_gene,1,By,SampleAnnot$cond,mean))
lFC_log=MeanExpr_log[,-1]-MeanExpr_log[,1]%o%rep(1,4)
lFC_event=cbind(0,lFC_cond[match(PSI_Annot_cond[,'gene_id'],rn(lFC_cond)),])

for(i in 2:5){
	PopDiffSpliceFull$lFC_cond[PopDiffSpliceFull$cond==i]=lFC_log[match(PopDiffSpliceFull$gene[PopDiffSpliceFull$cond==i],rownames(lFC_log)),i-1]
}


table(apply(FDR<0.05 & TESTABLE,1,sum))

testable_DSG=apply(TESTABLE,2,function(x){sort(unique(GeneSymbols[x]))})
DSG=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE,2,function(x){sort(unique(GeneSymbols[x]))})

luq(unlist(DSG)) # 515
isDAS=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE,1,any)
luq(setdiff(unlist(DSG),DSG[[1]])) # 309 (testable and differential >0.05 in STIM but NOT in NS)

mean(Pval[isDAS,1]>0.05) # 0.2911392

isTestable_DSG=apply(TESTABLE,2,function(x){DSG%in%GeneSymbols[x]}) # in which condition is a DSG testable
isDSG=apply(TESTABLE & Pval<0.05,2,function(x){DSG%in%GeneSymbols[x]}) # in which condition is a DSG marginally significant
apply((isDSG & isTestable_DSG) == isTestable_DSG,1,all) # is a DSG shared across all conditions where it is observed
mean(apply((isDSG & isTestable_DSG) == isTestable_DSG,1,all)) # 34.4% 

##### Figure 3B NbCond DSG VS Nb Cond Spliced 
#library(DescTools)
#tab=table(apply((isDSG & isTestable_DSG),1,sum), apply(isTestable_DSG,1,sum))
#rcol = SetAlpha(acol[ncol(tab)+nrow(tab):1], 0.5)
#PlotCirc( tab[5:1,],acol = rep(acol[c(1,3,6,7,8)],2),rcol = rcol[rev(c(1,3,6,7,8))])


DSG=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE,2,function(x){sort(unique(GeneSymbols[x]))})
#DSG=apply(FDR<0.05 & abs(DELTA_PSI)>0.05,2,function(x){sort(unique(GeneSymbols[x]))})
#DSG=apply(FDR<0.05,2,function(x){sort(unique(GeneSymbols[x]))})
library(VennDiagram)
names(DSG)=condIndex#rep("",5)
VD=venn.diagram(DSG,col=colERC[2*c(0,1,2,3,4)+1],fill=colERC[2*c(0,1,2,3,4)+1],filename=NULL,margin=0.05,main='DSG deltaPSI=0.05')
grid.newpage()
grid.draw(VD)
barplot(sapply(DSG,length),col=colERC5,ylab='DIUG |deltaPSI|>0.05',las=2)

DAS=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE,2,function(x){rn(DELTA_PSI)[x]))})
library(VennDiagram)
names(DAS)=condIndex#rep("",5)
VD=venn.diagram(DAS,col=colERC[2*c(0,1,2,3,4)+1],fill=colERC[2*c(0,1,2,3,4)+1],filename=NULL,margin=0.05,main='DAS deltaPSI=0.05')
grid.newpage()
grid.draw(VD)
barplot(sapply(DAS,length),col=colERC5,ylab='DIUG |deltaPSI|>0.05',las=2)


Pval_r=matrix(1,length(toKeep),5)
rownames(Pval_r)=rn(PSI_Annot)[toKeep]
for (cond in 2:5){
	Pval_r[,cond]=apply(PSI_cond,1,function(x){P=try(wilcox.test(x[grep(paste('AFB[0-9]+-',cond,sep=''),cn(PSI_cond))]-x[grep(paste('AFB[0-9]+-',1,sep=''),cn(PSI_cond))],x[grep(paste('EUB[0-9]+-',cond,sep=''),cn(PSI_cond))]-x[grep(paste('EUB[0-9]+-',1,sep=''),cn(PSI_cond))])$p.value);if(class(P)=='try-error'){NA}else{P}})
	}

Pval_r=Pval_r[,-1]
FDR_r=Pval_r
FDR_r=matrix(p.adjust(as.numeric(Pval_r),'fdr'),nrow(Pval_r),ncol(Pval_r))

DELTA_PSI_r=matrix(0,length(toKeep),4)
rownames(DELTA_PSI_r)=rn(PSI_Annot)[toKeep]
for(cond in 2:5){
	wtest=match(rn(DELTA_PSI_r),rn(PSI_Annot))
	DELTA_PSI_r[,cond-1]=apply(PSI_cond,1,function(x){mean(x[grep(paste('AFB[0-9]+-',cond,sep=''),cn(PSI_cond))]-x[grep(paste('AFB[0-9]+-',1,sep=''),cn(PSI_cond))])-mean(x[grep(paste('EUB[0-9]+-',cond,sep=''),cn(PSI_cond))]-x[grep(paste('EUB[0-9]+-',1,sep=''),cn(PSI_cond))])})
}

save(DELTA_PSI,DELTA_PSI_r,Pval,Pval_r,PvalAdjIAV,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/PopDiffSpliceData_V5.Rdata',HOME))



load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup.Rdata',EVO_IMMUNO_POP))
PopDiffSpliceFull$has_sQTL=paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond)%in%paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond)
PopDiffSpliceFull$has_stim_sQTL=PopDiffSpliceFull$has_sQTL & !PopDiffSpliceFull$Testable_NS

PopDiffSpliceFull$Gene_has_sQTL=PopDiffSpliceFull$gene_id%in%RESobs_nodup_1Mb$gene

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

#load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V3.Rdata',EVO_IMMUNO_POP))


PopDiffSpliceFull$isImmune=PopDiffSpliceFull$gene %in% GoToTest_genes$immuneResponse
PopDiffSplice=PopDiffSpliceFull[PopDiffSpliceFull$FDR<0.05 &  abs(PopDiffSpliceFull$DeltaPSI)>0.05,]

RESobs_nodup_1Mb_cond$isPopDAS=paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond)%in%paste(PopDiffSplice$event_id,PopDiffSplice$cond)
RESobs_nodup_1Mb_cond$isPopDSG=paste(RESobs_nodup_1Mb_cond$gene,RESobs_nodup_1Mb_cond$cond)%in%paste(PopDiffSplice$gene_id,PopDiffSplice$cond)

has_popDS=By(PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05,PopDiffSpliceFull$gene_id,any)
isImmune=By(PopDiffSpliceFull$isImmune,PopDiffSpliceFull$gene_id,any)
has_sQTL=By(PopDiffSpliceFull$has_sQTL,PopDiffSpliceFull$gene_id,any)
isPopDiff=By(abs(PopDiffSpliceFull$lFC_gene)>0.2,PopDiffSpliceFull$gene_id,any)
odds.ratio(table(has_popDS, isImmune))
#             LowerCI       OR  UpperCI alpha         P
#odds ratio 0.9306314 1.252016 1.663517  0.05 0.1261912
odds.ratio(table(has_popDS, has_sQTL))
#            LowerCI       OR  UpperCI alpha             P
#odds ratio 6.745278 8.231732 10.06622  0.05 4.744636e-101
odds.ratio(table(has_popDS, has_sQTL & isImmune))
#           LowerCI       OR  UpperCI alpha            P
#odds ratio 2.60077 3.828613 5.563901  0.05 3.488089e-11
odds.ratio(table(has_popDS[has_sQTL],isImmune[has_sQTL]))
#             LowerCI        OR UpperCI alpha         P
#odds ratio 0.6215742 0.9304796 1.37647  0.05 0.7741337
odds.ratio(table(has_popDS, isPopDiff))
#            LowerCI       OR  UpperCI alpha           P
#odds ratio 1.890735 2.283787 2.760462  0.05 2.38239e-18
maxFC=apply(lFC_log,1,max)[match(unique(PopDiffSpliceFull$gene),rn(lFC_log))]
odds.ratio(table(has_popDS, maxFC>1))

#has_popDS_ns=By((PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05)[PopDiffSpliceFull$cond==1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==1],any)
#has_sQTL_ns=By(PopDiffSpliceFull$has_sQTL[PopDiffSpliceFull$cond==1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==1],any)
#has_popDS_stim=By((PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05)[PopDiffSpliceFull$cond>1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond>1],any)
#has_sQTL_stim=By(PopDiffSpliceFull$has_sQTL[PopDiffSpliceFull$cond>1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond>1],any)

#has_popDS_flu=By((PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05)[PopDiffSpliceFull$cond==5],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==5],any)

has_sQTL_cond=sapply(1:5,function(i){By(PopDiffSpliceFull$has_sQTL[PopDiffSpliceFull$cond==i],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==i],any)})
has_popDS_cond=sapply(1:5,function(i){By((PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05)[PopDiffSpliceFull$cond==i],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==i],any)})

for(i in 1:5){
	print(c(all=mean(has_popDS_cond[[i]]),sQTL=mean(has_popDS_cond[[i]][has_sQTL_cond[[i]]]),mean(has_sQTL_cond[[i]])))
}

# Fig 7A :

for(i in 1:5){
	PopDiffSpliceFull[,paste('Pval',condIndex[i],sep='_')]=Pval[match(PopDiffSpliceFull$event_id,PSI_Annot_cond$event_id),i]
	}
FDR_th=max(Pval[FDR<0.05])



### V1 	
has_popDS_ns=By((PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05)[PopDiffSpliceFull$cond==1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==1],any)
has_sQTL_ns=By(PopDiffSpliceFull$has_sQTL[PopDiffSpliceFull$cond==1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==1],any)
has_popDS_stim=By((PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05 & !(PopDiffSpliceFull$gene_id%in%names(has_popDS_ns)[has_popDS_ns]))[PopDiffSpliceFull$cond>1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond>1],any)
has_sQTL_stim=By(PopDiffSpliceFull$has_stim_sQTL[PopDiffSpliceFull$cond>1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond>1],any)

odds.ratio(table(has_popDS_ns,has_sQTL_ns))
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 11.82845 16.28278 22.61931  0.05 1.821606e-74
odds.ratio(table(has_popDS_stim,has_sQTL_stim))
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 3.779465 5.304344 7.374372  0.05 8.056178e-20

#####################################
###	Figure 7A	V1	popDSG	5x5inches ###
###################################@##
par(mar=c(5.1,5.5,4.1,0.5))
barplot(c(all_ns=mean(has_popDS_ns),all_stim=mean(has_popDS_stim), sQTL_ns=mean(has_popDS_ns[has_sQTL_ns]),sQTL_stim=mean(has_popDS_stim[has_sQTL_stim]))*100,col=acol[c(6,7,3,1)],space=c(0.3,0.1,0.4,0.1)*1.2,las=2,ylab='Percentage of differentially spliced\n genes between populations')

### V2	
has_popDS_ns=By((PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05)[PopDiffSpliceFull$cond==1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==1],any)
has_sQTL_ns=By(PopDiffSpliceFull$has_sQTL[PopDiffSpliceFull$cond==1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==1],any)
has_popDS_stim=By((PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05)[PopDiffSpliceFull$cond>1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond>1],any)
has_sQTL_stim=By(PopDiffSpliceFull$has_sQTL[PopDiffSpliceFull$cond>1],PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond>1],any)

odds.ratio(table(has_popDS_ns,has_sQTL_ns))
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 11.82845 16.28278 22.61931  0.05 1.821606e-74
odds.ratio(table(has_popDS_stim,has_sQTL_stim))
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 6.571674 8.120051 10.05133  0.05 3.547894e-88


PopDiffSpliceFull$isPopDAS_ns= PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05 & PopDiffSpliceFull$cond==1
PopDiffSpliceFull$isPopDAS_stim= PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05 & PopDiffSpliceFull$cond>1 & !(PopDiffSpliceFull$gene_id%in%names(has_popDS_ns)[has_popDS_ns])

#####################################
###	Figure 7A V2	popDSG	5x5inches ###
###################################@##
par(mar=c(5.1,5.5,4.1,0.5))
barplot(c(all_ns=mean(has_popDS_ns),all_stim=mean(has_popDS_stim), sQTL_ns=mean(has_popDS_ns[has_sQTL_ns]),sQTL_stim=mean(has_popDS_stim[has_sQTL_stim]))*100,col=acol[c(6,7,3,1)],space=c(0.3,0.1,0.4,0.1)*1.2,las=2,ylab='Percentage of differentially spliced\n genes between populations')
#####################################
###	Figure 7A V3	popDSG	5x5inches ###
###################################@##
barplot(c(nosQTL_ns=mean(has_popDS_ns[!has_sQTL_ns]),nosQTL_stim=mean(has_popDS_stim[!has_sQTL_stim]), sQTL_ns=mean(has_popDS_ns[has_sQTL_ns]),sQTL_stim=mean(has_popDS_stim[has_sQTL_stim]))*100,col=acol[c(6,7,3,1)],space=c(0.3,0.1,0.4,0.1)*1.2,las=2,ylab='Percentage of differentially spliced\n genes between populations')

PopDiffSpliceFull$isPopDAS_ns= PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05 & PopDiffSpliceFull$cond==1
PopDiffSpliceFull$isPopDAS_stim= PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05 & PopDiffSpliceFull$cond>1 & !(PopDiffSpliceFull$gene_id%in%names(has_popDS_ns)[has_popDS_ns])


x=PopDiffSpliceFull[which(PopDiffSpliceFull$isPopDAS_ns & !PopDiffSpliceFull$has_sQTL),]
x[order(-abs(x$DeltaPSI))[1:4],c('symbol','event_type','cond','DeltaPSI','Pval')]

x=PopDiffSpliceFull[which(PopDiffSpliceFull$isPopDAS_stim & !PopDiffSpliceFull$has_sQTL),]
x[order(-abs(x$DeltaPSI))[1:4],c('symbol','event_type','cond','DeltaPSI','Pval')]

x=PopDiffSpliceFull[which(PopDiffSpliceFull$isPopDAS_ns & PopDiffSpliceFull$has_sQTL),]
x[order(-abs(x$DeltaPSI))[1:4],c('symbol','event_type','cond','DeltaPSI','Pval')]

x=PopDiffSpliceFull[which(PopDiffSpliceFull$isPopDAS_stim & PopDiffSpliceFull$has_sQTL),]
x[order(-abs(x$DeltaPSI))[1:4],c('symbol','event_type','cond','DeltaPSI','Pval')]


#x=PopDiffSpliceFull[which(PopDiffSpliceFull$isPopDAS_stim & PopDiffSpliceFull$has_sQTL),]
#x[order(abs(x$Pval))[1:20],c('symbol','event_type','cond','DeltaPSI','Pval')]
#
#x=PopDiffSpliceFull[which(PopDiffSpliceFull$isPopDAS_ns & PopDiffSpliceFull$has_sQTL),]
#x[order(abs(x$Pval))[1:20],c('symbol','event_type','cond','DeltaPSI','Pval')]
#
#x=PopDiffSpliceFull[which(PopDiffSpliceFull$isPopDAS_ns & !PopDiffSpliceFull$has_sQTL),]
#x[order(abs(x$Pval))[1:20],c('symbol','event_type','cond','DeltaPSI','Pval')]
#
#x=PopDiffSpliceFull[which(PopDiffSpliceFull$isPopDAS_stim & !PopDiffSpliceFull$has_sQTL),]
#x[order(abs(x$Pval))[1:20],c('symbol','event_type','cond','DeltaPSI','Pval')]


###################################
###		Figure 4A				###
###################################
# V1: without immune sQTL
barplot(c(all=mean(has_popDS),immune=mean(has_popDS[isImmune]),sQTL=mean(has_popDS[has_sQTL]))*100,col=acol[c(6,3,1)],las=2,ylab='Percentage of differentially spliced\n genes between populations')
# V2: with immune sQTL
barplot(c(all=mean(has_popDS),immune=mean(has_popDS[isImmune]),sQTL=mean(has_popDS[has_sQTL]),immune_sQTL=mean(has_popDS[has_sQTL & isImmune]))*100,col=acol[c(6,7,3,1)],las=2,ylab='Percentage of differentially spliced\n genes between populations')

barplot(c(all=mean(has_popDS),sQTL=mean(has_popDS[has_sQTL]))*100,col=acol[c(6,3)],las=2,ylab='Percentage of differentially spliced\n genes between populations')

#barplot(c(all_ns=mean(has_popDS_ns),all_stim=mean(has_popDS_stim),sQTL_ns=mean(has_popDS_ns[has_sQTL_ns]),sQTL_stim=mean(has_popDS_stim[has_sQTL_stim]))*100,col=acol[c(6,7,3,1)],las=2,ylab='Percentage of differentially spliced\n genes between populations')
barplot(c(all_ns=mean(has_popDS_ns),sQTL_ns=mean(has_popDS_ns[has_sQTL_ns]),all_stim=mean(has_popDS_stim),sQTL_stim=mean(has_popDS_stim[has_sQTL_stim]))*100,col=acol[c(6,7,3,1)],space=c(0.3,0.1,0.4,0.1)*1.2,las=2,ylab='Percentage of differentially spliced\n genes between populations')

barplot(c(all_ns=mean(has_popDS_ns),sQTL_ns=mean(has_popDS_ns[has_sQTL_ns]),all_stim=mean(has_popDS_stim),sQTL_stim=mean(has_popDS_stim[has_sQTL_stim]))*100,col=colERC[1:4],space=c(0.3,0.1,0.4,0.1)*1.2,las=2,ylab='Percentage of differentially spliced\n genes between populations')



###################################
###		Figure 4A	END			###
###################################

load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup_withSelect_RBP.Rdata',EVO_IMMUNO_POP))
allSQTL=sapply(unique(RESobs_nodup_1Mb_cond$snps),getSNP)

RESobs_nodup_1Mb_cond$isPopDAS=paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond)%in%paste(PopDiffSplice$event_id,PopDiffSplice$cond)
RESobs_nodup_1Mb_cond$isPopDSG=paste(RESobs_nodup_1Mb_cond$gene,RESobs_nodup_1Mb_cond$cond)%in%paste(PopDiffSplice$gene_id,PopDiffSplice$cond)
luq(paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond)[RESobs_nodup_1Mb_cond$isPopDAS])

library(mediation)
PopDiffSpliceFull$mediation_coeff=0
PopDiffSpliceFull$mediation_coeff_low=0
PopDiffSpliceFull$mediation_coeff_high=0

w=which(PopDiffSpliceFull$has_sQTL & (PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05))
snp_ids=rep('',nrow(PopDiffSpliceFull))
snp_ids[w]=RESobs_nodup_1Mb_cond$snps[match(paste(PopDiffSpliceFull$event_id[w],PopDiffSpliceFull$cond[w]),paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond))]

PopDiffSpliceFull$snp_ids=snp_ids
PopDiffSpliceFull$snp_beta=0
PopDiffSpliceFull$snp_beta[w]=RESobs_nodup_1Mb_cond$beta[match(paste(PopDiffSpliceFull$event_id[w],PopDiffSpliceFull$cond[w]),paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond))]
PopDiffSpliceFull$snp_R2=0
PopDiffSpliceFull$snp_R2[w]=RESobs_nodup_1Mb_cond$R2[match(paste(PopDiffSpliceFull$event_id[w],PopDiffSpliceFull$cond[w]),paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond))]

for (i in w){
	myPSI=PSI_cond[PopDiffSpliceFull$event_id[i],SampleAnnot$cond==PopDiffSpliceFull$cond[i]]
	pop=ifelse(substr(names(myPSI),1,3)=='EUB',0,1)
	SNP=allSQTL[match(substr(names(myPSI),1,6),rn(allSQTL)),snp_ids[i]]
	mod.y=lm(myPSI~pop+SNP)
	mod.snp=lm(SNP~pop)
	mm=mediate(mod.snp,mod.y, treat='pop',mediator='SNP')
	PopDiffSpliceFull$mediation_coeff[i]=mm$n0
	PopDiffSpliceFull$mediation_coeff_low[i]=mm$n0.ci[1]
	PopDiffSpliceFull$mediation_coeff_high[i]=mm$n0.ci[2]
}

By(PopDiffSpliceFull$snp_beta[w]>0, PopDiffSpliceFull$event_type[w],mean)
#       A3        A5        AF        AL        MX        RI        SE 
#0.5753425 0.4605263 0.4415584 0.5333333 0.4000000 0.6190476 0.4250000 
wilcox.test(PopDiffSpliceFull$snp_beta[w][PopDiffSpliceFull$event_type[w]=='AF'])  # p-value = 0.474


w=which(PopDiffSpliceFull$has_sQTL & (PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05))
boxplot(abs(PopDiffSpliceFull$DeltaPSI[w])~cut(abs(PopDiffSpliceFull$snp_beta[w]),c(0,0.05,0.1,0.15,0.2,1)))
boxplot(abs(PopDiffSpliceFull$mediation_coeff[w])~cut(abs(PopDiffSpliceFull$snp_beta[w]),c(0,0.05,0.1,0.15,0.2,1)),add=T)

MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),c(0,0.05,0.1,0.15,0.2,1)),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),c(0,0.05,0.1,0.15,0.2,1)),mean)

cuts=c(0,0.05,0.1,0.15,0.2,1)
cuts=quantile(abs(PopDiffSpliceFull$snp_beta[w]),seq(0,1,l=5));cuts[1]=0;cuts[length(cuts)]=Inf
cuts=c(0,0.1,0.2,1)
MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w],cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),paste(PopDiffSpliceFull$cond[w],cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w],cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)

barplot(MeanDPSI,col=colERC[c(0:4*2+1)[rep(1:5,e=length(cuts)-1)]],las=2)
barplot(MeanMEDPSI,add=T,col=colERC[c(0:4*2+2)[rep(1:5,e=length(cuts)-1)]], axisnames=F,axes=F)

MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w],cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),paste(PopDiffSpliceFull$cond[w],cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w],cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)

MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
barplot(MeanDPSI,col=colERC[c(0:1*2+1)[rep(1:2,e=length(cuts)-1)]],las=2,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
x=barplot(MeanMEDPSI,add=T,col=colERC[c(0:1*2+2)[rep(1:2,e=length(cuts)-1)]], axisnames=F,axes=F,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
par(xpd=T)
text(x,MeanDPSI+max(MeanDPSI)*0.2,labels=paste(round(MeanMED*100),'%'),srt=90)

cuts=quantile(abs(PopDiffSpliceFull$snp_R2[w]),seq(0,1,l=5));cuts[1]=0;cuts[length(cuts)]=Inf
cuts=c(0,0.2,0.4,0.6,1)
MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_R2[w]),cuts)),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_R2[w]),cuts)),mean)
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_R2[w]),cuts)),mean)
barplot(MeanDPSI,col=colERC[c(0:1*2+1)[rep(1:2,e=length(cuts)-1)]],las=2,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
x=barplot(MeanMEDPSI,add=T,col=colERC[c(0:1*2+2)[rep(1:2,e=length(cuts)-1)]], axisnames=F,axes=F,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
par(xpd=T)
text(x,MeanDPSI+max(MeanDPSI)*0.2,labels=paste(round(MeanMED*100),'%'),srt=90)

########### Fig 3B mediation effect by sQTL effect size
w=which(PopDiffSpliceFull$has_sQTL & (PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05))
cuts=c(0,0.05,0.1,0.15,0.2,1)
MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts),mean)
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts),mean)


w=which(PopDiffSpliceFull$has_sQTL & (PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05))
cuts=c(0,0.05,0.1,0.15,0.2,1)
MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts),mean)
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts),mean)


########### Fig 3B mediation effect by sQTL effect size
w=which(PopDiffSpliceFull$has_sQTL & (PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05))
barplot(MeanDPSI,col=colERC[1],las=2)
x=barplot(MeanMEDPSI,add=T,col=colERC[2], axisnames=F,axes=F)
par(xpd=T)
text(x,MeanDPSI+max(MeanDPSI)*0.2,labels=paste(round(MeanMED*100),'%'),srt=90)

By(abs(PopDiffSpliceFull$DeltaPSI[w]*PopDiffSpliceFull$mediation_coeff[w]),cut(abs(PopDiffSpliceFull$snp_beta[w]),c(0,0.05,0.1,0.15,0.2,1)),mean)

By(PopDiffSpliceFull$mediation_coeff[w],cut(abs(PopDiffSpliceFull$snp_beta[w]),c(0,0.05,0.1,0.15,0.2,1)),mean)

barplot(PopDiffSpliceFull$mediation_coeff[w[order(-PopDiffSpliceFull$mediation_coeff[w])]])
barplot((PopDiffSpliceFull$DeltaPSI[w]*PopDiffSpliceFull$mediation_coeff[w])[order(-PopDiffSpliceFull$DeltaPSI[w])[1:20]]])


mean(PopDiffSpliceFull$mediation_coeff])
By(PopDiffSpliceFull$mediation_coeff[w],PopDiffSpliceFull$cond[w],mean)
PopDiffSplice=PopDiffSpliceFull[which(PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05),]
save(PopDiffSplice,PopDiffSpliceFull,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/Table_S3_PopDiffSpliceFull_V5.Rdata',HOME))
write.table(PopDiffSpliceFull,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/Table_S3_PopDiffSpliceFull_V5.txt',HOME),quote=F,row.names=F,sep='\t')

RESobs_nodup_1Mb_cond
PopDiffSplice=PopDiffSplice[order()]
barplot(By(pmax(0,pmin(1,PopDiffSplice[which(PopDiffSplice$has_sQTL),'mediation_coeff'])),cut(PopDiffSplice[which(PopDiffSplice$has_sQTL),'DeltaPSI'],c(-0.6,0.2,-0.1,0.1,0.2,0.6),by=0.2),mean))


barplot(PopDiffSplice[which(PopDiffSplice$has_sQTL),'DeltaPSI']*pmax(0,pmin(1,PopDiffSplice[which(PopDiffSplice$has_sQTL),'mediation_coeff'])),col=colERC5[PopDiffSplice[which(PopDiffSplice$has_sQTL),'cond']],border='#00000000', width=0.5)



########### enrichement of splice QTL among up-regulated genes

RESobs_nodup_1Mb$maxFC=max(lFC_log[match(RESobs_nodup_1Mb$gene,rownames(lFC_log)),])
RESobs_nodup_1Mb$FC_maxR2Cond=lFC_log[cbind(match(RESobs_nodup_1Mb$gene,rownames(lFC_log)),RESobs_nodup_1Mb$condition-1)]
maxFC=apply(lFC_log,1,max)[match(unique(PSI_Annot$gene_id[toKeep]),rn(lFC_log))]
hasQTL=unique(PSI_Annot$gene_id[toKeep])%in%RESobs_nodup_1Mb$gene
odds.ratio(table(maxFC>1, hasQTL))
#            LowerCI       OR  UpperCI alpha            P
# odds ratio 1.698845 2.050053 2.468974  0.05 7.613463e-14

for(i in 2:5){
	maxFC=lFC_log[match(unique(PSI_Annot$gene_id[toKeep]),rn(lFC_log)),i-1]
	hasQTL=unique(PSI_Annot$gene_id[toKeep])%in%RESobs_nodup_1Mb_cond$gene[RESobs_nodup_1Mb_cond$cond==i]
	print(i)
	print(odds.ratio(table(maxFC>1, hasQTL)))
}
#[1] 2
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 1.558671 2.198156 3.055684  0.05 8.989865e-06
#[1] 3
#            LowerCI     OR  UpperCI alpha            P
#odds ratio 1.972795 2.7127 3.692409  0.05 1.844609e-09
#[1] 4
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 2.814402 3.653722 4.718305  0.05 6.166981e-21
#[1] 5
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 3.557666 4.607005 5.940807  0.05 1.847335e-28

for(i in 2:5){
	TESTABLE=as.matrix(PSI_Annot[toKeep,paste('Testable',condIndex,sep='_')])
	maxFC=lFC_log[match(unique(PSI_Annot$gene_id[toKeep][TESTABLE[,i]]),rn(lFC_log)),i-1]
	hasQTL=unique(PSI_Annot$gene_id[toKeep][TESTABLE[,i]])%in%RESobs_nodup_1Mb_cond$gene[RESobs_nodup_1Mb_cond$cond==i]
	print(i)
	print(odds.ratio(table(maxFC>1, hasQTL)))
}
#[1] 2
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 1.429926 2.034348 2.856486  0.05 6.138164e-05
#[1] 3
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 1.605016 2.216995 3.033657  0.05 1.601423e-06
#[1] 4
#            LowerCI       OR  UpperCI alpha          P
#odds ratio 2.091701 2.732691 3.553744  0.05 3.0818e-13
#[1] 5
#            LowerCI      OR  UpperCI alpha            P
#odds ratio 2.199733 2.86093 3.705977  0.05 1.027836e-14

###################################
###		Figure 4B	START		###
###################################

# by R2
cuts=quantile(abs(PopDiffSpliceFull$snp_R2[w]),seq(0,1,l=5));cuts[1]=0;cuts[length(cuts)]=Inf
cuts=c(0,0.2,0.4,0.6,1)
w=which(PopDiffSpliceFull$has_sQTL & (PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05))
MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_R2[w]),cuts)),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_R2[w]),cuts)),mean)
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_R2[w]),cuts)),mean)
barplot(MeanDPSI,col=colERC[c(0:1*2+1)[rep(1:2,e=length(cuts)-1)]],las=2,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
x=barplot(MeanMEDPSI,add=T,col=colERC[c(0:1*2+2)[rep(1:2,e=length(cuts)-1)]], axisnames=F,axes=F,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
par(xpd=T)
text(x,MeanDPSI+max(MeanDPSI)*0.2,labels=paste(round(MeanMED*100),'%'),srt=90)
# by beta value
cuts=c(0,0.05,0.1,0.15,0.2,1)
w=which(PopDiffSpliceFull$has_sQTL & (PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05))
MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
barplot(MeanDPSI,col=colERC[c(0:1*2+1)[rep(1:2,e=length(cuts)-1)]],las=2,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
x=barplot(MeanMEDPSI,add=T,col=colERC[c(0:1*2+2)[rep(1:2,e=length(cuts)-1)]], axisnames=F,axes=F,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
par(xpd=T)
text(x,MeanDPSI+max(MeanDPSI)*0.2,labels=paste(round(MeanMED*100),'%'),srt=90)

# with FancyColors
barplot(MeanDPSI,col=acol[c(6,3)[rep(1:2,e=length(cuts)-1)]],las=2,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
x=barplot(MeanMEDPSI,add=T,col=acol[c(7,1)[rep(1:2,e=length(cuts)-1)]], axisnames=F,axes=F,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
par(xpd=T)
text(x,MeanDPSI+max(MeanDPSI)*0.2,labels=paste(round(MeanMED*100),'%'),srt=90)

###################################
###		Figure 4B	END			###
###################################
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/Table_S3_PopDiffSpliceFull_V5.Rdata',HOME))
#odds.ratio(table(PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DELTA_PSI)>0.05,PopDiffSpliceFull$isImmune))



###### NO GO enrichment at 20% FDR 
resGO=lapply(1:5,function(i){GOSeq(unique(PopDiffSplice$gene_id[PopDiffSplice$cond==i]),unique(PopDiffSpliceFull$gene_id[PopDiffSpliceFull$cond==i]),FDR=1)})
###### NO GO enrichment 

###################################
###		Table S3A				###
###################################

TableS3=PSI_Annot_cond[c(1:8,grep('MeanPSI', cn(PSI_Annot_cond)),grep('Support', cn(PSI_Annot_cond)))]
colnames(TableS3)[grep('Support',colnames(TableS3))]=paste('gene_FPKM',condIndex,sep='_')
DSE=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE,2,ifelse,'yes','')
for (i in 1:5){DSE[!(TESTABLE[,i]),i]='n.d.'}
colnames(DSE)=paste('Differentially_spliced_AFB-EUB_',condIndex,sep='')
TESTABLE_yes=apply(TESTABLE,2,ifelse,'yes','')
colnames(TESTABLE_yes)=paste('Alternatively_spliced_',condIndex,sep='')
colnames(Pval)=paste('P-value_AFB-EUB_',condIndex,sep='')
colnames(FDR)=paste('FDR_AFB-EUB_',condIndex,sep='')
colnames(DELTA_PSI)=paste('Delta_PSI_AFB-EUB_',condIndex,sep='')
TableS3=cbind(TableS3,TESTABLE_yes,Pval,FDR,DELTA_PSI,DSE)
has_sQTL=apply(outer(TableS3$event_id,1:5,paste),2,function(x){ifelse(x%in%paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond),'yes','')})
has_sQTL_snp=apply(outer(TableS3$event_id,1:5,paste),2,function(x){ifelse(x%in%paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond),RESobs_nodup_1Mb_cond$snps[match(x,paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond))],'')})
for (i in 1:5){has_sQTL[!(TESTABLE[,i]),i]='n.d.'}

load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/Table_S3_PopDiffSpliceFull_V5.Rdata',HOME))

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
TableS3$isImmune=ifelse(TableS3$gene_id%in%GoToTest_genes$immuneResponse,'yes','')
TableS3$event_coord=gsub('ENSG[0-9]+;.*:[0-9XY]+:(.*):.*','\\1',TableS3$event_id)
mediation_coeff=matrix(0,nrow(TableS3),5)
mediation_coeff=apply(outer(TableS3$event_id,1:5,paste),2,function(x){ifelse(x%in%paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond),round(PopDiffSpliceFull[match(x,paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond)),'mediation_coeff'],2),'')})
mediation_coeff_CI=apply(outer(TableS3$event_id,1:5,paste),2,function(x){w=match(x,paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond));ifelse(is.na(w) | !PopDiffSpliceFull$has_sQTL[w] | PopDiffSpliceFull$FDR[w]>0.05 | abs(PopDiffSpliceFull$DeltaPSI[w])<0.05,'',paste(round(PopDiffSpliceFull[w,'mediation_coeff'],2),'[',round(PopDiffSpliceFull$mediation_coeff_low[w],2),'-',round(PopDiffSpliceFull$mediation_coeff_high[w],2),']'))})
mediation_coeff_CI_max=apply(outer(TableS3$event_id,1:5,paste),2,function(x){w=match(x,paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond));ifelse(is.na(w) | !PopDiffSpliceFull$has_sQTL[w] | PopDiffSpliceFull$FDR[w]>0.05 | abs(PopDiffSpliceFull$DeltaPSI[w])<0.05,'',paste(round(pmax(0,pmin(1,PopDiffSpliceFull[w,'mediation_coeff'])),2),'[',round(pmax(0,pmin(1,PopDiffSpliceFull[w,'mediation_coeff_low'])),2),'-',round(pmax(0,pmin(1,PopDiffSpliceFull[w,'mediation_coeff_high'])),2),']'))})
mediation_Pct=apply(outer(TableS3$event_id,1:5,paste),2,function(x){w=match(x,paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond));ifelse(is.na(w) | !PopDiffSpliceFull$has_sQTL[w] | PopDiffSpliceFull$FDR[w]>0.05 | abs(PopDiffSpliceFull$DeltaPSI[w])<0.05,'',paste(100*round(pmax(0,pmin(1,PopDiffSpliceFull[w,'mediation_coeff'])),2),'[',100*round(pmax(0,pmin(1,PopDiffSpliceFull[w,'mediation_coeff_low'])),2),'-',100*round(pmax(0,pmin(1,PopDiffSpliceFull[w,'mediation_coeff_high'])),2),']'))})
colnames(mediation_Pct)=paste('mediation_Pct',condIndex)
colnames(mediation_coeff)=paste('mediation_coeff',condIndex)
colnames(mediation_coeff_CI)=paste('mediation_coeff_CI',condIndex)
colnames(has_sQTL)=paste('has_sQTL',condIndex)
colnames(has_sQTL_snp)=paste('has_sQTL_snp',condIndex)

TableS3=cbind(TableS3,has_sQTL,has_sQTL_snp,mediation_coeff,mediation_coeff_CI,mediation_Pct)
write.table(TableS3,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/TableS3_V5_OneLinePerEvent.txt',HOME),quote=F,sep='\t',row.names=F)
TableS3A=read.table(file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/TableS3_V5_OneLinePerEvent.txt',HOME),sep='\t',quote='')
TableS3A$empty=''
TableS3A$DiffSplicedAny=apply(TableS3A[,grep('Differentially_spliced_AFB.EUB',colnames(TableS3A))]=='yes',1,any)
TableS3A$FDR=apply(TableS3A[,grep('FDR_AFB.EUB_',colnames(TableS3A))],1,min)
mycols=c('event_id','symbol','event_type','FDR','DiffSplicedAny',sapply(condIndex,function(cond){c('empty',paste(c('Alternatively_spliced_','P.value_AFB.EUB_','FDR_AFB.EUB_','Delta_PSI_AFB.EUB_','Differentially_spliced_AFB.EUB_','has_sQTL_snp.','mediation_Pct.'),cond,sep=''))}))
write.table(TableS3A[,mycols],file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/TableS3A_V5_OneLinePerEvent.txt',HOME),quote=F,sep='\t',row.names=F)

TableS3A=as.data.frame(fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/TableS3A_V5_OneLinePerEvent.txt',HOME))
###################################
###		Table S3A	END			###
###################################


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################						END  V5					 								################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################




TESTABLE=as.matrix(PSI_Annot[toKeep,paste('Testable',condIndex,sep='_')])
GeneSymbols=PSI_Annot$symbol[toKeep]

Require('IAV_RNAseq')
pop=grepl(paste('AFB[0-9]+-',cond,sep=''),IAV_RNAseq$ID)
mRNA_IAV=IAV_RNAseq$flu_FPKM
PvalAdjIAV=apply(PSI_cond[,match(IAV_RNAseq$ID,cn(PSI_cond))],1,function(x){P=try(summary(lm(x~pop+mRNA_IAV))$coeff[2,4]);if(class(P)=='try-error'){NA}else{P}})
w=which(TESTABLE,arr=T)
PopDiffSpliceFull=cbind(PSI_Annot_cond[w[,1],],Pval=Pval[w],DeltaPSI=DELTA_PSI[w],cond=w[,2],PvalAdjIAV=cbind(NA,NA,NA,NA,PvalAdjIAV)[w])
PopDiffSpliceFull$DELTA_PSI_R=DELTA_PSI[w]-DELTA_PSI[cbind(w[,1],1)]
PopDiffSpliceFull$FDR=FDR[w]
PopDiffSpliceFull=PopDiffSpliceFull[order(PopDiffSpliceFull$Pval),]
rownames(PopDiffSpliceFull)=NULL
#PopDiffSpliceFull=PopDiffSpliceFull[!duplicated(paste(PopDiffSpliceFull$gene_id,PopDiffSpliceFull$cond)),]

lFC=population.diff[,grep('Beta_.*_EUBvsAFB',cn(population.diff))]
rownames(lFC)=rn(FPKM_gene)
PopDiffSpliceFull$lFC_gene=NA
for(i in 1:5){
	PopDiffSpliceFull$lFC_gene[PopDiffSpliceFull$cond==i]=lFC[match(PopDiffSpliceFull$gene[PopDiffSpliceFull$cond==i],rownames(lFC)),i]
}
PopDiffSpliceFull$has_sQTL=paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond)%in%paste(RESobs_FDR_nodup$event_id,RESobs_FDR_nodup$cond)
PopDiffSpliceFull$isImmune=PopDiffSpliceFull$gene %in% GoToTest_genes$immuneResponse
PopDiffSplice=PopDiffSpliceFull[PopDiffSpliceFull$FDR<0.05 &  abs(PopDiffSpliceFull$DeltaPSI)>0.05,]
chisq.test(table(PopDiffSpliceFull$FDR<0.05 &  abs(PopDiffSpliceFull$DeltaPSI)>0.05,PopDiffSpliceFull$has_sQTL))

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/DIU_POP_PSI_MISO.pdf',HOME),height=5.8,width=5.8)
DIUG=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE,2,function(x){sort(unique(GeneSymbols[x]))})
library(VennDiagram)
names(DIUG)=condIndex#rep("",5)
VD=venn.diagram(DIUG,col=colERC[2*0:4+2],fill=colERC[2*0:4+1],filename=NULL,margin=0.05,main='DIUG deltaPSI=0.05')
grid.newpage()
grid.draw(VD)
barplot(sapply(DIUG,length),col=colERC5,ylab='DIUG |deltaPSI|>0.05',las=2)

DIUG=apply(FDR<0.05 &TESTABLE,2,function(x){sort(unique(GeneSymbols[x]))})
library(VennDiagram)
names(DIUG)=condIndex#rep("",5)
VD=venn.diagram(DIUG,col=colERC[2*0:4+2],fill=colERC[2*0:4+1],filename=NULL,margin=0.05,main='DIUG 5% FDR')
grid.newpage()
grid.draw(VD)
barplot(sapply(DIUG,length),col=colERC5,ylab='DIUG 5% FDR',las=2)
dev.off()


X=sapply(c(0,0.05,0.1),function(i){w=which(PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>i);By(PopDiffSpliceFull$has_sQTL[w],PopDiffSpliceFull$cond[w],mean)})
X0=By(PopDiffSpliceFull$has_sQTL,PopDiffSpliceFull$cond,mean)
rownames(X)=condIndex
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/DiffPSI_Populations_sQTL_MISO.pdf',HOME)
barplot(t(X)*100,col=paste(rep(substr(colERC5,1,7),e=3),rep(c('55','99','DD'),5),sep=''),beside=T,las=1,ylim=c(0,100),ylab='% of splicing event under genetic control')
for(i in 1:5){
	rect(1+4*(i-1),X0[i]*100,4+4*(i-1),0,col='#00000088',lwd=2)
	rect(1+4*(i-1),X0[i]*100,4+4*(i-1),0,col='#DDDDDD',lwd=2,density=10)
}
par(xpd=T)
legend(1,110,ncol=4,fill=c(paste(rep(substr(colERC5[1],1,7),e=3),c('55','99','DD'),sep=''),'#00000088'),border=c(rep('#00000000',3),'#DDDDDD'),legend=c('FDR<5%','FDR<5%\n|DeltaPSI|>5%','FDR<5%\n|DeltaPSI|>10%','all'),bty="n")
legend(1,110,ncol=4,col=c(rep('#00000000',3),'#DDDDDD'),fill=c(rep('#00000000',3),'#DDDDDD'),border=c(rep('#00000000',3),'#DDDDDD'),density=20,legend=c('FDR<5%','FDR<5%\n|DeltaPSI|>5%','FDR<5%\n|DeltaPSI|>10%','all'),bty="n",text.col='#00000000')
par(xpd=F)
dev.off()

mygroups=as.factor(ifelse(!PopDiffSpliceFull$has_sQTL,ifelse(!PopDiffSpliceFull$isImmune,'non-immune','immune'),ifelse(!PopDiffSpliceFull$isImmune,'non-immune\nsQTL','immune\nsQTL')))
mygroups=reorder(mygroups, new.order=c(3,1, 4, 2))
tab=table(mygroups,PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05)

barplot(c(mean(has_popDS),mean(has_popDS[isImmune]),mean(has_popDS[has_sQTL]))*100,col=acol[c(4,3,1)],las=1,ylab='Percentage of differentially spliced genes')

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/DiffPSI_Populations_sQTL_immune_MISO.pdf',HOME))
barplot(tab[,2]/apply(tab,1,sum)*100,col=acol[c(4,5,3,1)],las=1,ylab='Percentage of population Differential Splicing')
dev.off()

RESobs_FDR_nodup$DiffPSI_pop=paste(RESobs_FDR_nodup$event_id,RESobs_FDR_nodup$cond)%in%paste(PopDiffSplice$event_id,PopDiffSplice$cond)
RESobs_FDR_nodup$isImmune=RESobs_FDR_nodup$gene %in% GoToTest_genes$immuneResponse

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/DIU_POP_PSI_byType_boxplot_MISO.pdf',HOME),height=4,width=6)
boxplot(PopDiffSpliceFull$DeltaPSI~ paste(PopDiffSpliceFull$event_type,PopDiffSpliceFull$cond),col=rep(colPSI,e=5),pch=16,cex=0.5,notch=T,ylim=c(-0.2,0.2),las=2)
abline(h=0,col='grey',lty=2)
boxplot(PopDiffSpliceFull$DeltaPSI~ paste(PopDiffSpliceFull$event_type,PopDiffSpliceFull$cond),col=rep(colPSI,e=5),pch=16,cex=0.5,notch=T,ylim=c(-0.06,0.06),las=2,outcol='#00000000')
abline(h=0,col='red')
dev.off()

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/DiffPSI_Populations_sQTL_iHS.pdf',HOME))
layout(matrix(1:2,1))
wEUB=which(!duplicated(RESobs_FDR_nodup$event_id) & RESobs_FDR_nodup$maf_EUB>0.05)
P=wilcox.test(RESobs_FDR_nodup[wEUB,'iHS_EUB']~RESobs_FDR_nodup$DiffPSI_pop[wEUB]) # p-value = 0.02157

boxplot(RESobs_FDR_nodup[wEUB,'iHS_EUB']~RESobs_FDR_nodup$DiffPSI_pop[wEUB],notch=TRUE,las=1,pch=16,cex=0.7,axes=F,col=colERC5[1:2],ylim=c(-3,3),ylab='iHS EUB',sub=paste('P=',Format(P$p.value,2)))
axis(2,las=2);axis(1,at=1:2,labels=c('sQTL','popDS-sQTL'))

wAFB=which(!duplicated(RESobs_FDR_nodup$event_id) & RESobs_FDR_nodup$maf_AFB>0.05)
P=wilcox.test(RESobs_FDR_nodup[wAFB,'iHS_AFB']~RESobs_FDR_nodup$DiffPSI_pop[wAFB]) # p-value = 0.2541
boxplot(RESobs_FDR_nodup[wAFB,'iHS_AFB']~RESobs_FDR_nodup$DiffPSI_pop[wAFB],notch=TRUE,las=1,pch=16,cex=0.7,axes=F,col=colERC5[1:2],ylim=c(-3,3),ylab='iHS AFB',sub=paste('P=',Format(P$p.value,2)))
axis(2,las=2);axis(1,at=1:2,labels=c('sQTL','popDS-sQTL'))
dev.off()

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/DiffPSI_Populations_sQTL_Fst.pdf',HOME),width=4)
wEUB=which(!duplicated(RESobs_FDR_nodup$event_id))
P=wilcox.test(RESobs_FDR_nodup[wEUB,'FST']~RESobs_FDR_nodup$DiffPSI_pop[wEUB]) # p-value = 0.02157
boxplot(RESobs_FDR_nodup[wEUB,'FST']~RESobs_FDR_nodup$DiffPSI_pop[wEUB],notch=TRUE,las=1,pch=16,cex=0.7,axes=F,col=colERC5[1:2],ylim=c(0,1),ylab='FST',sub=paste('P=',P$p.value))
axis(2,las=2);axis(1,at=1:2,labels=c('sQTL','popDS-sQTL'))
dev.off()

wEUB=which(!duplicated(RESobs_FDR_nodup$event_id) & RESobs_FDR_nodup$maf_EUB>0.05 & RESobs_FDR_nodup$daf_EUB>RESobs_FDR_nodup$daf_AFB)
wilcox.test(RESobs_FDR_nodup[wEUB,'iHS_EUB']~RESobs_FDR_nodup$DiffPSI_pop[wEUB]) # p-value = 0.02157

wAFB=which(!duplicated(RESobs_FDR_nodup$event_id) & RESobs_FDR_nodup$maf_AFB>0.05 & RESobs_FDR_nodup$daf_AFB>RESobs_FDR_nodup$daf_EUB)
wilcox.test(RESobs_FDR_nodup[wAFB,'iHS_AFB']~RESobs_FDR_nodup$DiffPSI_pop[wAFB]) # p-value = 0.05703

axis(2,las=2)
legend(1,110,ncol=1,density=1:10,legend=rep('',10),bty="n",text.col='#00000000')


#rect(5,X0[2]*100,8,0,col='#00000088')
#rect(9,X0[3]*100,12,0,col='#00000088')
#rect(13,X0[4]*100,16,0,col='#00000088')
#rect(17,X0[5]*100,20,0,col='#00000088')

save(PSI_cond,PSI_Annot_cond,FDR,DELTA_PSI,TESTABLE,Pval,PopDiffSplice,PopDiffSpliceFull,GeneSymbols,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_AllTypes_DEbyPop_MISO_clean.Rdata',sep=''))
write.table(PopDiffSpliceFull,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/PopDiffSpliceFull.txt',HOME),quote=F,row.names=F,sep='\t')