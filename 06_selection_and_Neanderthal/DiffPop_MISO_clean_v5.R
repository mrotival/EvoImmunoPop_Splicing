#################################################################
###		Differential splicing between Populations             ###
#################################################################

allGOTerms=as.data.frame(fread('allGOterms_EnsGRC37_13042017.txt'))

GoToTest=list(immuneResponse='GO:0006955',
	InnateImmuneResponse='GO:0045087',
	AdaptiveImmuneResponse='GO:0002250',
	AntiviralResponse='GO:0051607',
	AntibacterialResponse='GO:0042742',
	TF_DNAbinding='GO:0003700',
	development='GO:0032502',
	olfactory_receptors='GO:0004984')
	
GoToTest_genes=lapply(GoToTest,function(x){allGOTerms$'Gene stable ID'[allGOTerms$'GO term accession'==x]})
GoToTest_genes$all=unique(allGOTerms$'Gene stable ID')

# load data from Quach et al, Cell 2016. of differential expression between populations
population.diff=read.table(paste(HOME,'/03_Analysis/DiffExpression/DiffExp_population_ByCondition_geneLevel.txt',sep=''),sep='\t',quote='',header=T,comment='',dec='.')

load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/Ens70_HISAT/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME))

library(impute)
library(VennDiagram)

PSI_Annot=PSI_Annot[toKeep,]
PSI_cond=PSI_prov
rownames(PSI_cond)=PSI_Annot$event_id

Pval=matrix(1,nrow(PSI_Annot),5)
rownames(Pval)=rownames(PSI_Annot)
for (cond in 1:5){
	Pval[,cond]=apply(PSI_cond,1,function(x){P=try(wilcox.test(x[grep(paste('AFB[0-9]+-',cond,sep=''),colnames(PSI_cond))],x[grep(paste('EUB[0-9]+-',cond,sep=''),colnames(PSI_cond))])$p.value);if(class(P)=='try-error'){NA}else{P}})
	}
colnames(Pval)=paste('P-value_AFB-EUB_',condIndex,sep='')


FDR=matrix(p.adjust(as.numeric(Pval),'fdr'),nrow(Pval),ncol(Pval))
colnames(FDR)=paste('FDR_AFB-EUB_',condIndex,sep='')

DELTA_PSI=matrix(0,nrow(PSI_Annot),5)
rownames(DELTA_PSI)=rownames(PSI_Annot)
colnames(DELTA_PSI)=paste('Delta_PSI_AFB-EUB_',condIndex,sep='')

for(cond in 1:5){
	DELTA_PSI[,cond]=apply(PSI_cond,1,function(x){mean(x[grep(paste('AFB[0-9]+-',cond,sep=''),colnames(PSI_cond))])-mean(x[grep(paste('EUB[0-9]+-',cond,sep=''),colnames(PSI_cond))])})
}

TESTABLE=as.matrix(PSI_Annot[,paste('Testable',condIndex,sep='_')])
GeneSymbols=PSI_Annot$symbol

w=which(TESTABLE,arr=T)
PopDiffSpliceFull=cbind(PSI_Annot[ w[, 1], ], Pval = Pval[w], DeltaPSI = DELTA_PSI[w], cond = w[,2])
PopDiffSpliceFull$FDR=FDR[w]
PopDiffSpliceFull=PopDiffSpliceFull[order(PopDiffSpliceFull$Pval),]
rownames(PopDiffSpliceFull)=NULL
#PopDiffSpliceFull=PopDiffSpliceFull[!duplicated(paste(PopDiffSpliceFull$gene_id,PopDiffSpliceFull$cond)),]

lFC=population.diff[,grep('Beta_.*_EUBvsAFB',colnames(population.diff))]
rownames(lFC)=rownames(FPKM_gene)
PopDiffSpliceFull$lFC_gene=NA
for(cond in 1:5){
    matchGene=match(PopDiffSpliceFull$gene[PopDiffSpliceFull$cond==cond],rownames(lFC))
	PopDiffSpliceFull$lFC_gene[PopDiffSpliceFull$cond==cond]=lFC[ matchGene,cond]
}

MeanExpr=log2(1+GeneAnnot[,c("NS_mean", "LPS_mean", "PAM3_mean", "R848_mean", "Flu_mean")])
lFC_cond=MeanExpr[,-1]-MeanExpr[,1]%o%rep(1,4)
lFC_event=cbind(0,lFC_cond[match(PSI_Annot[,'gene_id'],rownames(lFC_cond)),])

PopDiffSpliceFull$lFC_cond=NA
for(cond in 2:5){
    matchGene = match(PopDiffSpliceFull$gene[PopDiffSpliceFull$cond==cond],rownames(lFC_cond))
    lFC_column=cond-1
	PopDiffSpliceFull$lFC_cond[PopDiffSpliceFull$cond==cond]=lFC_cond[matchGene,lFC_column]
}

for(i in 1:5){
	PopDiffSpliceFull[,paste('Pval',condIndex[i],sep='_')]=Pval[match(PopDiffSpliceFull$event_id,PSI_Annot$event_id),i]
	}
FDR_th=max(Pval[FDR<0.05])

table(apply(FDR<0.05 & TESTABLE,1,sum))
#    0     1     2     3     4     5 
#12245  2593   657   317   183   178 
testable_DSG=apply(TESTABLE,2,function(x){sort(unique(GeneSymbols[x]))})

###################################################################################################
########################## Figure S6A Nb DSG per cond (venn diagram) ##############################
###################################################################################################

DSG=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE,2,function(x){sort(unique(GeneSymbols[x]))})

length(unique(unlist(DSG))) # 515
length(unique(setdiff(unlist(DSG),DSG[[1]]))) # 309 (testable and differential >0.05 in STIM but NOT in NS)

library(VennDiagram)
names(DSG)=condIndex#rep("",5)
VD=venn.diagram(DSG,col=colERC[2*c(0,1,2,3,4)+1],fill=colERC[2*c(0,1,2,3,4)+1],filename=NULL,margin=0.05,main='DSG deltaPSI=0.05')
grid.newpage()
grid.draw(VD)
#barplot(sapply(DSG,length),col=colERC5,ylab='DSG |deltaPSI|>0.05',las=2)

###################################################################################################
##########################              END FIGURE S6A               ##############################
###################################################################################################

#load sQTL data 
RESobs_nodup_1Mb_cond=as.data.frame(fread('data/RESobs_nodup_1Mb_cond.txt'))
RESobs_nodup_1Mb=as.data.frame(fread('data/RESobs_nodup_1Mb.txt'))

PopDiffSpliceFull$has_sQTL=paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond)%in%paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond)
PopDiffSpliceFull$has_stim_sQTL=PopDiffSpliceFull$has_sQTL & !PopDiffSpliceFull$Testable_NS

PopDiffSpliceFull$Gene_has_sQTL=PopDiffSpliceFull$gene_id%in%RESobs_nodup_1Mb$gene
PopDiffSpliceFull$isImmune=PopDiffSpliceFull$gene %in% GoToTest_genes$immuneResponse

#########################################
# popDSG are enriched in sQTLs at basal state and upon stimulation  	
has_popDS=By(PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05,PopDiffSpliceFull$gene_id,any)
has_sQTL=By(PopDiffSpliceFull$has_sQTL,PopDiffSpliceFull$gene_id,any)

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

###################################################################################################
######################      FIGURE 6A   popDSG	5x5inches              ############################
###################################################################################################
colors_DSG=c("#A1DAB4", "#41B6C4", "#FECC5C","#E31A1C")
par(mar=c(5.1,5.5,4.1,0.5))
barplot(c(all_ns=mean(has_popDS_ns), all_stim=mean(has_popDS_stim), 
            sQTL_ns=mean(has_popDS_ns[has_sQTL_ns]),sQTL_stim=mean(has_popDS_stim[has_sQTL_stim]))*100,
            col=colors_DSG,space=c(0.3,0.1,0.4,0.1)*1.2,las=2,ylab='Percentage of differentially spliced\n genes between populations')

###################################################################################################
##########################          END FIGURE 6A                    ##############################
###################################################################################################

allSQTL=getGenos(unique(RESobs_nodup_1Mb_cond$snps))
snps_sQTL=t(as.matrix(allSQTL[-(1:5)]))
colnames(snps_sQTL)=allSQTL[[1]]

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
	SNP=allSQTL[match(substr(names(myPSI),1,6),rownames(allSQTL)),snp_ids[i]]
	mod.y=lm(myPSI~pop+SNP)
	mod.snp=lm(SNP~pop)
	mm=mediate(mod.snp,mod.y, treat='pop',mediator='SNP')
	PopDiffSpliceFull$mediation_coeff[i]=mm$n0
	PopDiffSpliceFull$mediation_coeff_low[i]=mm$n0.ci[1]
	PopDiffSpliceFull$mediation_coeff_high[i]=mm$n0.ci[2]
}

###################################
###		Figure 6B	START		###
###################################

#### by beta value
cuts=c(0,0.05,0.1,0.15,0.2,1)
w=which(PopDiffSpliceFull$has_sQTL & (PopDiffSpliceFull$FDR<0.05 & abs(PopDiffSpliceFull$DeltaPSI)>0.05))
# mean change in PSI value
MeanDPSI=By(abs(PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
# mean Percentage mediated
MeanMED=By(abs(PopDiffSpliceFull$mediation_coeff[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)
# mean change in PSI value mediated by the genotype
MeanMEDPSI=By(abs(PopDiffSpliceFull$mediation_coeff[w]*PopDiffSpliceFull$DeltaPSI[w]),paste(PopDiffSpliceFull$cond[w]>1,cut(abs(PopDiffSpliceFull$snp_beta[w]),cuts)),mean)

# make the plot
barplot(MeanDPSI,col=colERC[c(0:1*2+1)[rep(1:2,e=length(cuts)-1)]],las=2,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
x=barplot(MeanMEDPSI,add=T,col=colERC[c(0:1*2+2)[rep(1:2,e=length(cuts)-1)]], axisnames=F,axes=F,space=c(rep(0.2,length(cuts)-1),0.7,rep(0.2,length(cuts)-2)))
par(xpd=T)
text(x,MeanDPSI+max(MeanDPSI)*0.2,labels=paste(round(MeanMED*100),'%'),srt=90)

###################################
###		Figure 6B	END			###
###################################

###################################
###		Table S4A				###
###################################

TableS4=PSI_Annot[c("event_id", "event_type", "chrom", "start", "end", "strand","gene_id", "symbol",
                     "MeanPSI_NS", "MeanPSI_LPS", "MeanPSI_PAM3CSK4", "MeanPSI_R848", "MeanPSI_IAV", 
                     "Support_NS", "Support_LPS", "Support_PAM3CSK4", "Support_R848", "Support_IAV","coding_type")]

colnames(TableS4)[grep('Support',colnames(TableS4))]=paste('gene_FPKM',condIndex,sep='_')
DSE=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE,2,ifelse,'yes','') # Differentially spliced event
for (i in 1:5){DSE[!(TESTABLE[,i]),i]='n.d.'}
colnames(DSE)=paste('Differentially_spliced_AFB-EUB_',condIndex,sep='')

TESTABLE_yes=apply(TESTABLE,2,ifelse,'yes','')
colnames(TESTABLE_yes)=paste('Alternatively_spliced_',condIndex,sep='')

TableS4=cbind(TableS4,TESTABLE_yes,Pval,FDR,DELTA_PSI,DSE)
has_sQTL=apply(outer(TableS4$event_id,1:5,paste),2,function(x){ifelse(x%in%paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond),'yes','')})
has_sQTL_snp=apply(outer(TableS4$event_id,1:5,paste),2,function(x){ifelse(x%in%paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond),RESobs_nodup_1Mb_cond$snps[match(x,paste(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond))],'')})
for (i in 1:5){has_sQTL[!(TESTABLE[,i]),i]='n.d.'}


TableS4$isImmune=ifelse(TableS4$gene_id%in%GoToTest_genes$immuneResponse,'yes','')
TableS4$event_coord=gsub('ENSG[0-9]+;.*:[0-9XY]+:(.*):.*','\\1',TableS4$event_id)

mediation_Pct= apply(outer(TableS4$event_id,1:5,paste), 2, function(x){ 
                                w=match(x,paste(PopDiffSpliceFull$event_id,PopDiffSpliceFull$cond));
                                PCT_mediated = 100 * round( pmax(0, pmin(1, PopDiffSpliceFull[w, 'mediation_coeff'] ) ), 2 )
                                PCI_mediated_lowerBound = 100 * round( pmax( 0, pmin( 1, PopDiffSpliceFull[ w, 'mediation_coeff_low' ] ) ), 2 )
                                PCI_mediated_higherBound = 100 * round( pmax( 0, pmin( 1, PopDiffSpliceFull[ w, 'mediation_coeff_high' ] ) ), 2 )
                                is_DSG_SQTL= !is.na(w) & PopDiffSpliceFull$has_sQTL[w] & PopDiffSpliceFull$FDR[w]<=0.05 & abs(PopDiffSpliceFull$DeltaPSI[w])>=0.05
                                ifelse(is_DSG_SQTL, paste(PCT_mediated,'[',PCI_mediated_lowerBound,'-',PCI_mediated_higherBound,']'), '' ) 
                            } )

colnames(mediation_Pct)=paste('mediation_Pct',condIndex)
colnames(has_sQTL)=paste('has_sQTL',condIndex)
colnames(has_sQTL_snp)=paste('has_sQTL_snp',condIndex)

TableS4A=cbind(TableS4,has_sQTL,has_sQTL_snp,mediation_Pct)
TableS4A$empty=''
TableS4A$DiffSplicedAny=apply(TableS4A[,grep('Differentially_spliced_AFB.EUB',colnames(TableS4A))]=='yes',1,any)
TableS4A$FDR=apply(TableS4A[,grep('FDR_AFB.EUB_',colnames(TableS4A))],1,min)
mycols=c('event_id','symbol','event_type','FDR','DiffSplicedAny',sapply(condIndex,function(cond){c('empty',paste(c('Alternatively_spliced_','P.value_AFB.EUB_','FDR_AFB.EUB_','Delta_PSI_AFB.EUB_','Differentially_spliced_AFB.EUB_','has_sQTL_snp.','mediation_Pct.'),cond,sep=''))}))
write.table(TableS4A[,mycols],file=sprintf('tables/TableS4A_popDSG_OneLinePerEvent.txt',HOME),quote=F,sep='\t',row.names=F)

###################################
###		Table S4A	END			###
###################################

###################################
###		Figure S6B  			###
###################################
TableS4A=as.data.frame(fread(sprintf('tables/TableS4A_popDSG_OneLinePerEvent.txt',HOME))

mediation=list()
for(i in 1:5){
	w=which(TableS4A[,paste('Differentially_spliced_AFB.EUB_',condIndex[i],sep='')]=='yes')
	mediation[[i]]=pop_Dif[w,paste("mediation_Pct.",condIndex[i],sep='')]
	mediation[[i]]=ifelse(mediation[[i]]=="","0",mediation[[i]])
	mediation[[i]]=as.numeric(sapply(strsplit(mediation[[i]],' '),function(x){x[1]}))
}

layout(1:2)
x=unlist(mediation)[unlist(mediation)!=0]
group=as.factor(rep(1:5,sapply(mediation,function(x){sum(x!=0)})))
Pct0=sapply(mediation,function(x){mean(x==0)})
vioPlot_byGroup(x,group,ylim=c(-20,100),col=colERC5)
for(i in 1:5){
	rect(i-0.9/2,-15,i+0.9/2,-5,col='white')
	rect(i-0.45,-15,i-0.45+Pct0[i]*0.9,-5,col=colERC5[i])}
}

###################################
###	        	END  			###
###################################

