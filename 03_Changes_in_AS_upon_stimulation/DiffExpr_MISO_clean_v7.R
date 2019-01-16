Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME))

GeneAnnot=read.table(paste(HOME,'/Annotation/GeneAnnotation_hg37_ens70.txt',sep=''),sep='\t',comment='',quote='',stringsAsFactors=FALSE,header=T,row.names=1)
colnames(GeneAnnot)[1]='Ensembl.Gene.ID'
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='X']=23
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='Y']=24
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='MT']=26
GeneAnnot$Chromosome.Name=as.numeric(GeneAnnot$Chromosome.Name)
GeneAnnot=GeneAnnot[!is.na(GeneAnnot$Chromosome.Name),] # remove alternate chromosomes and unassigned scaffolds
GeneAnnot$Strand=ifelse(GeneAnnot$Strand>0,'+','-') # change strand from -1,1 to '+','-'
 
 # load allGOterms table (containing go ids and matching gene ids)
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

library(impute)
library(VennDiagram)

#################################################################
###			Differential splicing upon Stimulation : Test	  ###
#################################################################

luq=function(x){length(unique(x))}
By=function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x}

PSI_cond=PSI
PSI_Annot_cond=PSI_Annot[toKeep,]
rownames(PSI_cond)=PSI_Annot$event_id


Pval=matrix(1,nrow(PSI_Annot),5)
rownames(Pval)=rownames(PSI_Annot)
for (cond in 2:5){
	Pval[,cond]=apply(PSI_cond,1,function(x){P=try(wilcox.test(x[grep('-1',colnames(PSI_cond))],x[grep(paste('-',cond,sep=''),colnames(PSI_cond))])$p.value);if(class(P)=='try-error'){NA}else{P}})
	}

FDR=matrix(p.adjust(as.numeric(Pval),'fdr'),nrow(Pval),ncol(Pval))

DELTA_PSI=matrix(0,nrow(PSI_Annot),5)
rownames(DELTA_PSI)=rownames(PSI_Annot)
for(cond in 2:5){
	DELTA_PSI[,cond]=PSI_Annot[,paste('MeanPSI_',condIndex[cond],sep='')]-PSI_Annot[,paste('MeanPSI_',condIndex[1],sep='')]
}


DELTA_DIV=matrix(0,nrow(PSI_Annot),5)
rownames(DELTA_PSI)=rownames(PSI_Annot)
for(cond in 2:5){
	DELTA_DIV[,cond]=abs(PSI_Annot[,paste('MeanPSI_',condIndex[cond],sep='')]-0.5)-abs(PSI_Annot[,paste('MeanPSI_',condIndex[1],sep='')]-0.5)
}


TESTABLE=as.matrix(PSI_Annot[,paste('Testable',condIndex,sep='_')])
SUPPORT=as.matrix(PSI_Annot[,paste('Support',condIndex,sep='_')])
MEANPSI=as.matrix(PSI_Annot[,paste('MeanPSI',condIndex,sep='_')])
PCTNA=as.matrix(PSI_Annot[,paste('PctNA',condIndex,sep='_')])
JUNC=as.matrix(PSI_Annot[,paste('JuncCovered',condIndex,sep='_')])

# list fo genes with an AS event that pass our quality filter in a given condition
testable_list=sapply(1:5,function(cond){unique(PSI_Annot[TESTABLE[,cond],'gene_id'])})
# list of AS event that pass our quality filter in a given condition
testable_event_list=sapply(1:5,function(cond){unique(PSI_Annot[TESTABLE[,cond],'event_id'])})
GeneSymbols=PSI_Annot$symbol
GeneIDs=PSI_Annot$gene_id
luq(unlist(testable_event_list)) # 16173 events
luq(unlist(testable_list)) #[1] 4739 genes
luq(setdiff(unlist(testable_event_list[-1]),testable_event_list[[1]])) # 3367 events
luq(setdiff(unlist(testable_list[-1]),testable_list[[1]])) # 630 genes

# repeat the characteristics of condition NS 5 times (frequent AS event, 

TESTABLE_1=t(sapply(TESTABLE[,1],rep,5))
SUPPORT_1=t(sapply(SUPPORT[,1],rep,5))
MEANPSI_1=t(sapply(MEANPSI[,1],rep,5))
PCTNA_1=t(sapply(PCTNA[,1],rep,5))
JUNC_1=t(sapply(JUNC[,1],rep,5))

TESTABLE_DIFF=(PCTNA_1<0.05 & PCTNA<0.05) & (SUPPORT_1>10 & SUPPORT>10) & (TESTABLE| TESTABLE_1)
colnames(TESTABLE_DIFF)=paste('tested_DiffSplice',condIndex,sep='_')

# For each gene and condition, compute shannon entropy
SHANNON=apply(-MEANPSI*log(MEANPSI)-(1-MEANPSI)*log(1-MEANPSI),2,By,GeneSymbols,sum,na.rm=T)
# For each gene, compute shannon entropy in =ns cond (5 times)
SHANNON_1=apply(-MEANPSI_1*log(MEANPSI_1)-(1-MEANPSI_1)*log(1-MEANPSI_1),2,By,GeneSymbols,sum,na.rm=T)

luq(GeneSymbols[apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE_DIFF,1,any)]) # 1919
# list Differentially spliced genes

DSG=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE_DIFF,2,By,GeneSymbols,any)
luq(GeneSymbols[apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE_DIFF,1,any)]) # 1919

# percentage of DSG that increase their Entropy
apply(SHANNON>SHANNON_1 & DSG,2,sum,na.rm=T)/apply(DSG,2,sum,na.rm=T)
#      MeanPSI_NS      MeanPSI_LPS MeanPSI_PAM3CSK4     MeanPSI_R848      MeanPSI_IAV 
#             NaN        0.7350902        0.7055556        0.7276708        0.7496218 

####################################################
##  generate lists of frequent AS events and DSG  ##
#################@##################################

TableS2=PSI_Annot[c(1:8,grep('MeanPSI', colnames(PSI_Annot)),grep('Support', colnames(PSI_Annot)),grep('coding_type', colnames(PSI_Annot)))]
colnames(TableS2)[grep('Support',colnames(TableS2))]=paste('gene_FPKM',condIndex,sep='_')
DSE=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE_DIFF,2,ifelse,'yes','')
for (i in 1:5){DSE[!TESTABLE_DIFF[,i],i]='n.d.'}
colnames(DSE)=paste('Differentially_spliced_NS-',condIndex,sep='')
TESTABLE_yes=apply(TESTABLE,2,ifelse,'yes','')
TESTABLE_DIFF_yes=apply(TESTABLE_DIFF,2,ifelse,'yes','')
colnames(TESTABLE_yes)=paste('Alternatively_spliced_',condIndex,sep='')
colnames(Pval)=paste('P-value_NS-',condIndex,sep='')
colnames(FDR)=paste('FDR_NS-',condIndex,sep='')
colnames(DELTA_PSI)=paste('Delta_PSI_NS-',condIndex,sep='')

lFC=log2(1+TableS2[,grep('gene_FPKM',colnames(TableS2))])
lFC=lFC[,-1]-lFC[,1]%o%rep(1,4)
colnames(lFC)=gsub('gene_FPKM','log2FC',colnames(lFC))

TableS2=cbind(TableS2,TESTABLE_yes,Pval[,-1],FDR[,-1],DELTA_PSI[,-1],DSE[,-1],lFC,PCTNA,TESTABLE_DIFF_yes)
TableS2$isDiffSpliced=ifelse(apply(DSE=='yes',1,any),'yes','')
TableS2$maxlFC=apply(TableS2[,grepl('log2FC',colnames(TableS2))],1,max)

TableS2$coding_type_REF_to_ALT=TableS2$coding_type
TableS2$coding_type_ALT_to_REF=ifelse(TableS2$coding_type=='loss of function','gain of function',ifelse(TableS2$coding_type=='gain of function','loss of function',TableS2$coding_type))

Gene_changes_percentage_of_coding_transcript = unique(TableS2$gene_id) %in% TableS2$gene_id[TableS2$coding_type%in%c('gain of function','loss of function') & TableS2$"isDiffSpliced"=='yes']
Gene_changes_protein_isoform = unique(TableS2$gene_id)%in%TableS2$gene_id[TableS2$coding_type%in%c('modified protein') & TableS2$"isDiffSpliced"=='yes']

table(Gene_changes_percentage_of_coding_transcript,Gene_changes_protein_isoform)

TableS2$isImmune=ifelse(TableS2$gene_id%in%GoToTest_genes$immuneResponse,'yes','')
TableS2$event_coord=gsub('ENSG[0-9]+;.*:[0-9XY]+:(.*):.*','\\1',TableS2$event_id)
TableS2$empty=''

TableS2A=TableS2[,c('event_id','gene_id','symbol','event_type','chrom','event_coord','strand','start','end',paste(c('gene_FPKM','MeanPSI','Alternatively_spliced'),rep(condIndex,e=3),sep='_'),'empty',paste('PctNA',condIndex,sep='_'),'coding_type')]
write.table(TableS2A,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/TableS2A.txt',HOME),quote=F,sep='\t',row.names=F)

TableS2B=TableS2[,c('event_id','symbol','event_type','chrom','event_coord','strand','start','end','empty','isDiffSpliced','Alternatively_spliced_NS',paste(c('Alternatively_spliced_','P-value_NS-','FDR_NS-','Delta_PSI_NS-','Differentially_spliced_NS-'),rep(condIndex[-1],e=5),sep=''),'empty',colnames(lFC),'isImmune','coding_type')]
write.table(TableS2B,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/TableS2B.txt',HOME),quote=F,sep='\t',row.names=F)

##################################################################
###			Differential splicing upon Stimulation By Type		###
###################################################################

w=which(TESTABLE_DIFF,arr=T)
CondDiffSpliceFull=cbind(PSI_Annot[w[,1],],Pval=Pval[w],DeltaPSI=DELTA_PSI[w],cond=w[,2])
rownames(CondDiffSpliceFull)=NULL

#### add info on immune VS non immune
CondDiffSpliceFull$isImmune=CondDiffSpliceFull$gene%in%GoToTest_genes$immuneResponse

#### add info on logFC (based one Gene Annot)
###### TODO: recode this from Support to improve clarity
MeanExpr=log2(1+GeneAnnot[,c("NS_mean", "LPS_mean", "PAM3_mean", "R848_mean", "Flu_mean")])
lFC=MeanExpr[,-1]-MeanExpr[,1]%o%rep(1,4)
CondDiffSpliceFull$lFC_gene=NA
for(i in 2:5){
	CondDiffSpliceFull$lFC_gene[CondDiffSpliceFull$cond==i]=lFC[match(CondDiffSpliceFull$gene[CondDiffSpliceFull$cond==i],rownames(lFC)),i-1]
}

CondDiffSpliceFull=CondDiffSpliceFull[CondDiffSpliceFull$cond!=1,]
CondDiffSpliceFull$FDR=p.adjust(CondDiffSpliceFull$Pval,'fdr')
CondDiffSpliceFull=CondDiffSpliceFull[order(-abs(CondDiffSpliceFull$DeltaPSI)),]

# Impact of Differential splicing on protein coding potential
CondDiffSpliceFull$DeltaPSI_coding=CondDiffSpliceFull$DeltaPSI*ifelse(CondDiffSpliceFull$coding_type=='loss of function',-1,ifelse(CondDiffSpliceFull$coding_type=='gain of function',1,0))

# signifiant AS event only
CondDiffSplice=CondDiffSpliceFull[which(CondDiffSpliceFull$FDR<0.05 & abs(CondDiffSpliceFull$DeltaPSI)>0.05),]

save(PSI_cond,PSI_Annot,FDR,DELTA_PSI,TESTABLE,TESTABLE_DIFF,Pval,CondDiffSplice,CondDiffSpliceFull,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/PSI_events_AllTypes_DEbyCond_MISO_clean_V5.Rdata',HOME))


####################################################
##          Comparison between DEG and DSG        ##
####################################################

#### Number of DSG without DEG
length(unique(CondDiffSpliceFull$gene[CondDiffSpliceFull$FDR<0.05 &  abs(CondDiffSpliceFull$DeltaPSI)>0.05 & abs(CondDiffSpliceFull$lFC_gene)<0.2]))
# 460

#### Enrichment of DSG in DEG
odds.ratio(table(CondDiffSpliceFull$FDR<0.05 &  abs(CondDiffSpliceFull$DeltaPSI)>0.05 ,abs(CondDiffSpliceFull$lFC_gene)>0.5))
#            LowerCI       OR  UpperCI alpha             P
# odds ratio 1.898261 2.008408 2.125205  0.05 1.951038e-135

is_DSG=sapply(2:5,function(i){ifelse(unique(PSI_Annot$gene)%in%CondDiffSpliceFull$gene[CondDiffSpliceFull$FDR<0.05 & abs(CondDiffSpliceFull$DeltaPSI)>0.05 & CondDiffSpliceFull$cond==i],'X','')})
TableS2B=read.table(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/TableS2B.txt',HOME),sep='\t',header=T)
lFC_event_TableS2=TableS2B[,grep('log2FC_',colnames(TableS2B))]
rownames(lFC_event_TableS2)=TableS2B$event_id

is_DEG=apply(lFC[match(unique(PSI_Annot$gene),rownames(lFC)),],2,function(x){y=cut(x,c(-Inf,-1,-0.5,-0.2,0.2,0.5,1,Inf)); levels(y)=c('---','--','-','','+','++','+++');as.character(y)})
rownames(is_DEG)=unique(PSI_Annot$gene)

TableS2C=data.frame(gene=unique(PSI_Annot$gene),symbol=G2S(unique(PSI_Annot$gene)),is_DSG,is_DEG)
colnames(TableS2C)=c('gene','symbol','DSG_LPS','DSG_PAM3CSK4','DSG_R848','DSG_IAV','DEG_LPS','DEG_PAM3CSK4','DEG_R848','DEG_IAV')
TableS2C$isDSG_nonDEG=ifelse((TableS2C$DEG_LPS=='' & TableS2C$DSG_LPS!='') | (TableS2C$DEG_PAM3CSK4=='' & TableS2C$DSG_PAM3CSK4!='') |(TableS2C$DEG_R848=='' & TableS2C$DSG_R848!='') | (TableS2C$DEG_IAV=='' & TableS2C$DSG_IAV!=''),'yes','')
TableS2C$isDSG_nonDEG=gsub('^/|/$','',gsub(' +','/',paste(ifelse(TableS2C$DEG_LPS=='' & TableS2C$DSG_LPS!='','LPS',''),ifelse(TableS2C$DEG_PAM3CSK4=='' & TableS2C$DSG_PAM3CSK4!='','PAM3CSK4',''),ifelse(TableS2C$DEG_R848=='' & TableS2C$DSG_R848!='','R848',''),ifelse(TableS2C$DEG_IAV=='' & TableS2C$DSG_IAV!='','IAV',''))))
write.table(TableS2C,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/TableS2C_DSG_nonDEG.txt',HOME),row.names=F,quote=F,sep='\t')


##############################################
##              get some numbers            ##
##############################################

length(unique(TableS2B[TableS2B$isDiffSpliced=='yes','symbol']))
# 1919

# consequence of the change in AS upon stimulation
table(c(TableS2$coding_type_REF_to_ALT[TableS2$ALT_increases_upon_stimulation], TableS2$coding_type_ALT_to_REF[TableS2$REF_increases_upon_stimulation])
#gain of function   isoform change loss of function       non coding 
#             433             2019              847              308 
table(c(TableS2$coding_type_REF_to_ALT[TableS2$ALT_increases_upon_stimulation], TableS2$coding_type_ALT_to_REF[TableS2$REF_increases_upon_stimulation]))/sum(TableS2$ALT_increases_upon_stimulation|TableS2$REF_increases_upon_stimulation,na.rm=T)
#gain of function   isoform change loss of function       non coding 
#      0.12343216       0.57554162       0.24144812       0.08779932 

#####  change in coding potential upon stimulation (Fig S2B)
w=which(CondDiffSpliceFull$DeltaPSI_coding!=0)
P_wilcox_coding=By(CondDiffSpliceFull$DeltaPSI_coding[w],CondDiffSpliceFull$cond[w],function(x){wilcox.test(x)$p.value})
#            2             3             4             5 
#5.452425e-127 2.459506e-109 7.735900e-106  1.013553e-29 


##################################
## 			Figure 2C			##
##################################

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/sup_figures/Fig1C.pdf',HOME),height=5.5,width=3)
DIUG_5=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE_DIFF,2,function(x){sort(unique(GeneSymbols[x]))})
names(DIUG_5)=paste('NS',condIndex,sep='-')#rep("",5)
DIUG_10=apply(FDR<0.05 & abs(DELTA_PSI)>0.1 & TESTABLE_DIFF,2,function(x){sort(unique(GeneSymbols[x]))})
names(DIUG_10)=paste('NS',condIndex,sep='-')#rep("",5)

par(mar=c(9,4,1,1))
barplot(sapply(DIUG_5[-1],length)/apply(TESTABLE_DIFF[,-1],2,function(x){luq(GeneSymbols[x])})*100,col=colERC[1:4*2+1],ylab='percentage of differentially\n spliced genes',las=3)
barplot(sapply(DIUG_10[-1],length)/apply(TESTABLE_DIFF[,-1],2,function(x){luq(GeneSymbols[x])})*100,col=colERC5[-1],add=T,axisnames=FALSE,axes=F)
dev.off()


##################################
## 			Figure END			##
##################################


########################################################################################
#### 	Comparisons of the percentage of differentially spliced genes per condition	####
########################################################################################

tab=cbind(apply(FDR<0.05 & abs(DELTA_PSI)>0.1 & TESTABLE_DIFF,2,function(x){luq(GeneSymbols[x])}),apply(TESTABLE_DIFF,2,function(x){luq(GeneSymbols[x])}))
tab[,2]=tab[,2]-tab[,1]
tab=tab[-1,]
OR=NULL
for (i in 1:3){
	for (j in (i+1):4){
	OR=rbind(OR,c(odds.ratio(tab[c(i,j),2:1]),paste(condIndex[j+1],'-',condIndex[i+1],'_DeltaPSI=10%',sep='')))
	}
}
tab=cbind(apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE_DIFF,2,function(x){luq(GeneSymbols[x])}),apply(TESTABLE_DIFF,2,function(x){luq(GeneSymbols[x])}))
tab[,2]=tab[,2]-tab[,1]
tab=tab[-1,]
for (i in 1:3){
	for (j in (i+1):4){
	OR=rbind(OR,c(odds.ratio(tab[c(i,j),2:1]),paste(condIndex[j+1],'-',condIndex[i+1],'_DeltaPSI=5%',sep='')))
	}
}
tab=cbind(apply(FDR<0.05 & abs(DELTA_PSI)>0.1 & abs(lFC_event)<0.2 & TESTABLE_DIFF,2,function(x){luq(GeneSymbols[x])}),apply(TESTABLE_DIFF & abs(lFC_event)<0.2,2,function(x){luq(GeneSymbols[x])}))
tab[,2]=tab[,2]-tab[,1]
tab=tab[-1,]
for (i in 1:3){
	for (j in (i+1):4){
	OR=rbind(OR,c(odds.ratio(tab[c(i,j),2:1]),paste(condIndex[j+1],'-',condIndex[i+1],'_DeltaPSI=10%_lFCunder0.2',sep='')))
	}
}
tab=cbind(apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & abs(lFC_event)<0.2 & TESTABLE_DIFF,2,function(x){luq(GeneSymbols[x])}),apply(TESTABLE_DIFF & abs(lFC_event)<0.2,2,function(x){luq(GeneSymbols[x])}))
tab[,2]=tab[,2]-tab[,1]
tab=tab[-1,]
for (i in 1:3){
	for (j in (i+1):4){
	OR=rbind(OR,c(odds.ratio(tab[c(i,j),2:1]),paste(condIndex[j+1],'-',condIndex[i+1],'_DeltaPSI=5%_lFCunder0.2',sep='')))
	}
}

write.table(OR,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/OR_DIU_STIMcomparisons_PSI_MISO_V7.txt',HOME),sep='\t',quote=F,row.names=F)


##################################
## 			Figure 2D			##
##################################
P_wilcox=By(CondDiffSpliceFull$DeltaPSI,paste(CondDiffSpliceFull$event_type,CondDiffSpliceFull$cond),function(x){wilcox.test(x)$p.value})
DPSI_mean=By(CondDiffSpliceFull$DeltaPSI,paste(CondDiffSpliceFull$event_type,CondDiffSpliceFull$cond),mean)
resDeltaPSI=data.frame(cond_type=names(P_wilcox),DPSI_mean,P_wilcox,P_adj=p.adjust(P_wilcox))
write.table(resDeltaPSI,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/sup_figures/Figure_2D_Pvalues_V7.txt',HOME),sep='\t',quote=F,row.names=F)

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/figures/Fig2D_NS-STIM_DPSI_byType_boxplot_V7.pdf',HOME),width=2.5,height=3.75)
nlevels=luq(paste(CondDiffSpliceFull$event_type,CondDiffSpliceFull$cond))
DF=data.frame(event_type=CondDiffSpliceFull$event_type,cond=CondDiffSpliceFull$cond)
DF=DF[order(DF$event_type,DF$cond),]
DF=DF[!duplicated(paste(DF$event_type,DF$cond)),]
DF$myColors=colPSI[match(DF[,1],unique(DF$event_type))]
DF$x_pos=1:nrow(DF)+as.numeric(as.factor(DF$event_type))*0.5
plot(1:nlevels,rep(0,nlevels),ylim=c(-0.2,0.2),xlim=range(DF$x_pos),axes=F,xlab='',ylab='',main='',col='#00000000',las=3)
for (i in 1:nrow(DF)){
			w=which(CondDiffSpliceFull$event_type==DF[i,'event_type'] & CondDiffSpliceFull$cond==DF[i,'cond'])
			if(length(w)>0){
				boxplot(CondDiffSpliceFull$DeltaPSI[w],at=DF$x_pos[i],add=T,col=DF$myColors[i],axes=F,pch='.',outcol='#00000000',notch=T,wex=2,pars = list(boxwex = 1.5))
				}
			rect(DF$x_pos[i]-0.4,-0.2+.008,DF$x_pos[i]+0.4,-0.2-.008,col=colERC5[DF$cond[i]])
#			rect(DF$x_pos[i]-0.4,0.2+.008,DF$x_pos[i]+0.4,0.2-.008,col=grey(1+log10(resDeltaPSI$P_adj[i])/172)) # withPany
#			rect(DF$x_pos[i]-0.4,0.2+.008,DF$x_pos[i]+0.4,0.2-.008,col=grey(ifelse(resDeltaPSI$P_adj[i]<1e-5,ifelse(resDeltaPSI$P_adj[i]<1e-10,ifelse(resDeltaPSI$P_adj[i]<1e-15,0.2,0.5),0.8),1))) # withPthresholds
			rect(DF$x_pos[i]-0.4,0.2+.008,DF$x_pos[i]+0.4,0.2-.008,col=grey((1-pmin(-log10(resDeltaPSI$P_adj[i]),16)/16))) # withP16
			}
par(xpd=F)
abline(h=0,col='red',lty=1)
axis(2,las=3); 
par(xpd=T)
for(j in names(colPSI)){ 
	axis(1,at=range(DF$x_pos[DF$event_type==j])+c(-0.4,0.4),labels=c('',''))
	text(y=-0.27,x=mean(range(DF$x_pos[DF$event_type==j])),labels=j,srt=90)
}
par(xpd=F)
#boxPlot_byGroup(CondDiffSpliceFull$DeltaPSI,CondDiffSpliceFull$event_type,CondDiffSpliceFull$cond,colorBy=1,colors=colPSI)
boxplot(CondDiffSpliceFull$DeltaPSI~ paste(CondDiffSpliceFull$event_type,CondDiffSpliceFull$cond),col=rep(colPSI,e=4),pch=16,cex=0.5,notch=T,ylim=c(-0.2,0.2),outcol='#00000000',las=3)
abline(h=0,col='red',lty=1)

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/figures/FIG2D_scale_Pany_P16.pdf',HOME),width=2.5,height=3.75)
Image(matrix(0:100%o%2-50),col=grey(seq(0,1,l=100)))
rect(-1,1,1,0)
axis(2,at=c(0,5,10,15)/16,labels=c(0,5,10,15))
dev.off()

##################################
## 			Figure END			##
##################################


##################################
## 			Figure 2E			##
##################################
CondDiffSpliceFull_noDup=CondDiffSpliceFull[order(CondDiffSpliceFull$gene_id,-abs(CondDiffSpliceFull$lFC_gene),-CondDiffSpliceFull$Signal2Noise_global),]
CondDiffSpliceFull_noDup=CondDiffSpliceFull_noDup[!duplicated(CondDiffSpliceFull_noDup$gene_id),]
max_lFC=CondDiffSpliceFull_noDup$lFC_gene
names(max_lFC)=CondDiffSpliceFull_noDup$gene_id
isDSG=By(CondDiffSpliceFull$FDR<0.05 & abs(CondDiffSpliceFull$DeltaPSI)>0.05,CondDiffSpliceFull$gene_id,any)
isImmune=names(isDSG)%in% GoToTest_genes$immuneResponse
max_lFC= max_lFC[names(isDSG)]

#isDSG=By(CondDiffSpliceFull$FDR<0.05 &  abs(CondDiffSpliceFull$DeltaPSI)>0.10,CondDiffSpliceFull$gene_id,any)
par(mar=c(2,2,2,3))
category=ifelse(isDSG,ifelse(isImmune, 'immune DSG','DSG'),'')
tab=table(category, cut(max_lFC,c(-10,-1,-0.5,0.5,1,10)))
mosaicplot(t(tab),col=acol[c(6,3,1)],las=1,main='',ylab='')
axis(4,at=seq(-0.036, 0.986,l=5),labels=seq(0,1,l=5)*100,las=2)
dev.off()
odds.ratio(table(abs(max_lFC)>0.5, isDSG))
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 1.918847 2.239849 2.617389  0.05 4.349366e-26

#for 10% DSG
#           LowerCI      OR  UpperCI alpha            P
#odds ratio 1.427738 1.74965 2.154834  0.05 1.715034e-08

x=splitGene_byCond(S2G('STAT2'))$ALL;wilcox.test(x[,5]-x[,1])$p.value
#[1] 4.102244e-34
mean(x[,5]-x[,1],na.rm=T)
# 1.351882
x=splitGene_byCond(S2G('NFKB1'))$ALL;wilcox.test(x[,4]-x[,1])$p.value
#[1] 4.327318e-33
mean(x[,4]-x[,1],na.rm=T)
# 1.516889

##################################
## 	Figure 2E	boxplots		##
##################################
# what gene to plot

layout(matrix(1:2,1))
cond=5
myCond=relevel(as.factor(condIndex[as.numeric(SampleAnnot$cond[SampleAnnot$cond%in%c(1,cond)])]),'NS')
boxplot(FPKM_gene[S2G('CLEC7A'),SampleAnnot$cond%in%c(1,cond)]~myCond,las=2,col=colERC5[c(1,cond)],main='',las=1,axes=T,notch=T,ylab='',outpch=16,cex=0.7)
boxplot(PSI_cond['ENSG00000172243;RI:12:10279170:10279307-10280346:10280444:-',SampleAnnot$cond%in%c(1,cond)]~myCond,las=2,col=colERC5[c(1,cond)],main='',las=1,axes=T,notch=T,ylim=c(0,1),outpch=16,cex=0.7)

 
layout(matrix(1:2,1))
cond=4
myCond=relevel(as.factor(condIndex[as.numeric(SampleAnnot$cond[SampleAnnot$cond%in%c(1,cond)])]),'NS')
boxplot(FPKM_gene[S2G('NFKB1'),SampleAnnot$cond%in%c(1,cond)]~myCond,las=2,col=colERC5[c(1,cond)],main='',las=1,axes=T,notch=T,ylab='',outpch=16,cex=0.7)
boxplot(PSI_cond['ENSG00000109320;SE:4:103422945-103432037:103432106-103446669:+',SampleAnnot$cond%in%c(1,cond)]~myCond,las=2,col=colERC5[c(1,cond)],main='',las=1,axes=T,notch=T,ylim=c(0,1),outpch=16,cex=0.7)

P_wilcox=By(CondDiffSpliceFull$DeltaPSI,paste(CondDiffSpliceFull$event_type,CondDiffSpliceFull$cond),function(x){wilcox.test(x)$p.value})


##################################
## 			Figure S3A			##
##################################
library(VennDiagram)

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/sup_figures/FigS2A.pdf',HOME),height=5.8,width=5.8)
par(mar=c(3,3,3,3))
TESTABLE_1=t(sapply(TESTABLE[,1],rep,5))
DIUG=apply(FDR<0.05 & abs(DELTA_PSI)>0.05 & TESTABLE_DIFF,2,function(x){sort(unique(GeneSymbols[x]))})

names(DIUG)=condIndex#rep("",5)
VD=venn.diagram(DIUG[-1],col=colERC[2*c(0,2,3,1)+4],fill=colERC[2*c(0,1,2,3)+3],filename=NULL,margin=0.05,main='DIUG deltaPSI=0.05')
grid.newpage()
grid.draw(VD)
barplot(sapply(DIUG[-1],length)/sapply(DIUG[-1],length),col=colERC5[-1],ylab='DIUG |deltaPSI|>0.05',las=2)
dev.off()

########################
####   Figure S3B   ####
########################

w=which(CondDiffSpliceFull$DeltaPSI_coding!=0)
DF=CondDiffSpliceFull[w,]
DF$cond=factor(condIndex[DF$cond],levels = condIndex[-1])
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/figures/FigS2X_DeltaPSI_coding_v2.pdf',HOME),width=4,height=2.5)
ggplot(DF,aes(cond,DeltaPSI_coding))+geom_violin(aes(fill=cond))+theme_minimal()+geom_boxplot(width=0.25,outlier.shape = NA,notch=T)+scale_fill_manual(values=colERC5[-1])+coord_cartesian(ylim=c(-0.1,0.1))+geom_hline(aes(yintercept=0))
dev.off()

########################
####   		END 	####
########################

####################################
############ Figure S3C ############
####################################

coding_consequence=c(TableS2$coding_type_REF_to_ALT[TableS2$ALT_increases_upon_stimulation], TableS2$coding_type_ALT_to_REF[TableS2$REF_increases_upon_stimulation])
lFC=c(TableS2$maxlFC[TableS2$ALT_increases_upon_stimulation],TableS2$maxlFC[TableS2$REF_increases_upon_stimulation])

coding_consequence_By_FC=t(as.matrix(table(coding_consequence,cut(lFC,c(-Inf,-.5,.5,Inf)))))/as.N(table(cut(lFC,c(-Inf,-.5,.5,Inf))))
#             coding_consequence
#              gain of function isoform change loss of function non coding
#  (-Inf,-0.5]       0.09530792     0.54692082       0.29325513 0.06451613
#  (-0.5,0.5]        0.11411549     0.55545371       0.24014665 0.09028414
#  (0.5, Inf]        0.16016151     0.58411844       0.16554509 0.09017497

coding_consequence_By_FC=melt(coding_consequence_By_FC,value.name='Pct')
colnames(coding_consequence_By_FC)[1]='log2FC'
levels(coding_consequence_By_FC$Max_log2FC)=c('down','same','up')
up_down_color=c(same="#FEE586",down="#E31A1C",up="#A1DAB4")

ggplot(coding_consequence_By_FC, aes(fill=Max_log2FC, y=Pct, x=coding_consequence)) + geom_bar(position="dodge", stat="identity") + theme_minimal() + scale_fill_manual(values=up_down_color)

####################################
######## END Figure S3C ############
####################################


##################################
CondDiffSpliceFull_noDup=CondDiffSpliceFull[order(CondDiffSpliceFull$gene_id,-abs(CondDiffSpliceFull$lFC_gene),-CondDiffSpliceFull$Signal2Noise_global),]
CondDiffSpliceFull_noDup=CondDiffSpliceFull_noDup[!duplicated(CondDiffSpliceFull_noDup$gene_id),]
max_lFC=CondDiffSpliceFull_noDup$lFC_gene
names(max_lFC)=CondDiffSpliceFull_noDup$gene_id
isDSG=By(CondDiffSpliceFull$FDR<0.05 & abs(CondDiffSpliceFull$DeltaPSI)>0.05,CondDiffSpliceFull$gene_id,any)
isImmune=names(isDSG)%in% GoToTest_genes$immuneResponse
max_lFC= max_lFC[names(isDSG)]

# gene has at least one AS event associated to an increase of non coding transcripts after stimulation
isDSG_lof=By(CondDiffSpliceFull$FDR<0.05 & abs(CondDiffSpliceFull$DeltaPSI)>0.05 & CondDiffSpliceFull$DeltaPSI_coding<0,CondDiffSpliceFull$gene_id,any)
# gene has at least one AS event associated to an decrease of non coding transcripts after stimulation
isDSG_gof=By(CondDiffSpliceFull$FDR<0.05 & abs(CondDiffSpliceFull$DeltaPSI)>0.05 & CondDiffSpliceFull$DeltaPSI_coding>0,CondDiffSpliceFull$gene_id,any)

# note that a gene can have both loss of function and gain of function AS events we do not consider those AS events :
odds.ratio(table(max_lFC>0.5, isDSG_gof & !isDSG_lof))
#             LowerCI        OR   UpperCI alpha           P
#odds ratio 0.5102116 0.6869623 0.9125331  0.05 0.008099409
odds.ratio(table(max_lFC>0.5, isDSG_gof & !isDSG_lof))
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 1.385232 1.914218 2.615792  0.05 6.583811e-05
# Among up-regulated genes we see more change toward gain of function and less toward loss of function.

odds.ratio(table(max_lFC< -0.5, isDSG_lof & !isDSG_gof))
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 1.492128 1.827877 2.247557  0.05 1.215996e-09
odds.ratio(table(max_lFC< -0.5, isDSG_gof & !isDSG_lof))
#             LowerCI        OR  UpperCI alpha          P
#odds ratio 0.6022013 0.7880346 1.032818  0.05 0.07874161
# We see the opposite trend among down-regulated genes 


####################################
############ Figure S3D ############
####################################

##### Fig S3D
#### conservation as a function of coding consequence: SE 
tab=table(PSI_Annot$coding_type, gsub(' (all)', '',PSI_Annot$conserved_site_GerpRS,fixed=T), PSI_Annot$event_type)

tab=table(PSI_Annot$coding_type, PSI_Annot$conserved_site_GerpRS, PSI_Annot$event_type)

color_conserv=c("#A1DAB4","#FEE586","#225EA8")
par(mar=c(12,4,1,1))
tabSE=tab[,c('gain of conserved splice site (all)','gain of conserved splice site','non conserved Alternative splice sites'),'SE']
colnames(tabSE)=c('both conserved','one conserved','non conserved')
barplot(t(tabSE/apply(tabSE,1,sum)),col=color_conserv,las=3)

####################################
######## END Figure S3D ############
####################################

