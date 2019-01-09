#######################################################
###				Aggregate psiQTL Pop Combined		###
#######################################################


# set the path to files that will be needed

allPSI_events=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/Ens70_HISAT/PSI_events_ALL_V7.2_withCounts.Rdata',sep='')
Adjusted_PSI_values=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME)

library(data.table)

# allGOterms
allGOTerms=as.data.frame(fread('/Volumes/evo_immuno_pop/Maxime/Shared/allGOterms_EnsGRC37_13042017.txt'))

# load data 
load(allPSI_events)
load(Adjusted_PSI_values)

####################################################################################
###                      Annotate genetic variants                               ###
####################################################################################

# load Map of genetic variants
Map_imputed=as.data.frame(fread('/Volumes/evo_immuno_pop/Maxime/SNP_annotation/Map_imputed_essential_informations.txt'))
####### add informations to the SNP annoation table (Map_imputed)

###### generic question: What SNPs are in Cis from a splicing event ######
library(GenomicRanges)
GRange_Map=makeGRangesFromDataFrame(Map_imputed[which(Map_imputed$SNPfreq>0.05 ),c('chromosome','position','allele.1','allele.2','minor_allele','ancestral_allele')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_PSI=makeGRangesFromDataFrame(PSI_Annot,seqnames.field='chrom',start.field='start',end.field='end',keep.extra.columns=TRUE)
Cisdist=1e6
GR_cis=flank(flank(GRange_PSI,Cisdist),-2*Cisdist-1)
ooCis=findOverlaps(GR_cis, GRange_Map)
wCis=which(Map_imputed$SNPfreq>=0.05)[unique(subjectHits(ooCis))]
###### END of generic question ########

Map_imputed$RegElt=gsub('-validated','',Map_imputed$RegElt)
Map_imputed$RegElt[is.na(Map_imputed$RegElt_Text)]=''

SPIDEX=as.data.frame(fread(sprintf("%s/Annotation/SPIDEX/SNPlist_allSPIDEX.txt",HOME)))

Map_imputed$SPIDEX_Z=SPIDEX$Z_SPIDEX[match(Map_imputed$snp.name,SPIDEX$snp.name)]
Map_imputed$Exonic=SPIDEX$Exon_zone[match(Map_imputed$snp.name,SPIDEX$snp.name)]
Map_imputed$groupSplice=ifelse(is.na(Map_imputed$SPIDEX_Z),'',ifelse(abs(Map_imputed$SPIDEX_Z)>2,paste(Map_imputed$Exonic,'deleterious'),Map_imputed$Exonic))
###########################################################################################################################

# count unique elements of a vector
luq=function(x){length(unique(x))}

# read all snp-gene associations in Cis
perm=0
thresholds=10^-c(seq(3,13,0.1),15,20,30,50) # thresholds used to compute FDR
RES=list()
for(pop in c('ALL')){
	for(cond in 1:5){
		for( CHR in 1:22){
		cat(pop, cond,CHR,'\n')
			load(paste(HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_KW_LOGIT_misoPSI_Ens70_HISAT_V7.Rdata',sep=''))
#			load(paste(HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_V8.Rdata',sep=''))
			RES[[paste(cond,pop,CHR)]]=RESCIS
		}
	}
}
RES=do.call(rbind,RES)
colnames(RES)[colnames(RES)=='gene']='event_id'
RESobs=RES
rm(RES)
gc()

########
########
condIndex=c('NS','LPS','PAM3CSK4','R848','IAV')
PSI_Annot=PSI_Annot[toKeep,]
TESTABLE=as.matrix(PSI_Annot[,paste('Testable',condIndex,sep='_')])
SUPPORT=as.matrix(PSI_Annot[,paste('Support',condIndex,sep='_')])
MEANPSI=as.matrix(PSI_Annot[,paste('MeanPSI',condIndex,sep='_')])
PCTNA=as.matrix(PSI_Annot[,paste('PctNA',condIndex,sep='_')])
JUNC=as.matrix(PSI_Annot[,paste('JuncCovered',condIndex,sep='_')])

testable_list=sapply(1:5,function(cond){unique(PSI_Annot[TESTABLE[,cond],'gene_id'])})
GeneSymbols=PSI_Annot$symbol

MeanExpr=log2(1+GeneAnnot[,c("NS_mean", "LPS_mean", "PAM3_mean", "R848_mean", "Flu_mean")])
lFC=MeanExpr[,-1]-MeanExpr[,1]%o%rep(1,4)

#lFC_event=cbind(0,lFC[match(PSI_Annot[,'gene_id'],rownames(lFC)),])

TESTABLE_1=t(sapply(TESTABLE[,1],rep,5))
SUPPORT_1=t(sapply(SUPPORT[,1],rep,5))
MEANPSI_1=t(sapply(MEANPSI[,1],rep,5))
PCTNA_1=t(sapply(PCTNA[,1],rep,5))
JUNC_1=t(sapply(JUNC[,1],rep,5))

TESTABLE_DIFF=(PCTNA_1<0.05 & PCTNA<0.05) & (SUPPORT_1>10 & SUPPORT>10) & (TESTABLE | TESTABLE_1)
colnames(TESTABLE_DIFF)=paste('tested_DiffSplice',condIndex,sep='_')

#TESTABLE_all = (apply(PCTNA<0.05 & SUPPORT>10, 1 ,all) & apply(TESTABLE,1,any)) %o%rep(1,5)==1

###### annotate SNP
### chromosomal position
mm=match(RESobs$snps,Map_imputed$snp.name)
RESobs$chrom=Map_imputed$chrom[mm]
RESobs$pos=Map_imputed$position[mm]
# locus is defined as ( chrom:position in Mb, eg. 6:30Mb )
RESobs$locus=Map_imputed$locus[mm]
### Allele frequency
# derived allele frequency (ancestral allele assigned based on 6EPO alignments)
RESobs$daf_AFB=Map_imputed$daf_char_AFB[mm] 
RESobs$daf_EUB=Map_imputed$daf_char_EUB[mm]
# global allele frequency (in European and African combined)
RESobs$SNPfreq=Map_imputed$SNPfreq[mm]
# remove SNPs for which position is unavailable or MAF is too low 
RESobs=RESobs[which(!is.na(RESobs$pos) & RESobs$SNPfreq>=0.05),]

###### annotate AS events
mm=match(RESobs$event_id,PSI_Annot$event_id)
# gene, symbol, type,  
RESobs$gene=PSI_Annot$gene_id[mm]
RESobs$symbol=PSI_Annot$symbol[mm]
RESobs$event_type=PSI_Annot$event_type[mm]
RESobs$event_start=as.numeric(PSI_Annot$start[mm])
RESobs$event_end=as.numeric(PSI_Annot$end[mm])
# is the event frequent (testable) in the condition under consideration
RESobs$isTestable=TESTABLE[cbind(mm,RESobs$cond)]
# compute distance between SNP and AS event 
# (distance is negative if SNP is upstream of the event, positive if the SNP is downstream of the AS event)
RESobs$CisDist = ifelse(PSI_Annot$strand[mm]=='+',1,-1) * ifelse(RESobs$pos < RESobs$event_start, RESobs$pos-RESobs$event_start, ifelse(RESobs$pos<RESobs$event_end, 0, RESobs$pos-RESobs$event_end))

# Compute the number of associated genes for each level of significance 
NbAssocGene=sapply(thresholds,function(th){luq(RESobs$gene[RESobs$pval<th & abs(RESobs$CisDist)<1e6 & RESobs$isTestable])})

# load the permutation results and number of associated genes for each level of significance 
NbAssocGenePerm=list()
for(perm in 1:100){
    cat(perm,'')
	RES=list()
	for( CHR in 1:22){
    	RESCHR=list()
		for(cond in 1:5){
            cat(perm,pop, cond,CHR,'\n')
            load(paste(HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_KW_LOGIT_misoPSI_',Source,'_V7.Rdata',sep=''))
            RESCHR[[paste(cond,pop,CHR)]]=as.data.table(RESCIS)
            }
        RES[[CHR]]=rbindlist(RESCHR)
		}
    RES=rbindlist(RES)
    RES[,gene:=substr(event_id,1,15)]
    # annotate snp
    mm=match(RES$snps,Map_imputed$snp.name)
    RES[,chrom:=Map_imputed$chrom[mm]]
    RES[,pos:=Map_imputed$position[mm]]
    RES[,SNPfreq:=Map_imputed$SNPfreq[mm]]
    # annotate event
    mm=match(RES$event_id,PSI_Annot$event_id)
    RES[,isTestable:=TESTABLE[cbind(mm,RES$cond)]]
    RES[,event_start:=as.numeric(PSI_Annot$start[mm])]
    RES[,event_end:=as.numeric(PSI_Annot$end[mm])]
    # compute distance
    RES[,CisDist:= ifelse(PSI_Annot$strand[mm]=='+',1,-1) * ifelse(RES$pos < RES$event_start, RES$pos-RES$event_start, ifelse(RES$pos<RES$event_end, 0, RES$pos-RES$event_end))]
    # filter out low frequency variants (security check)
    RES=RES[SNPfreq>=0.05,]
    # count significant genes for various thresholds 
    NbAssocGenePerm[[perm]]=sapply(thresholds,function(th){luq(RES$gene[RES$pval<th & abs(RES$CisDist)<1e6 & RES$isTestable])})
    # output a few numbers
    cat(NbAssocGenePerm[[perm]][c(1,11,21,31,41,51,61)],'\n')
}
# GET THE MEAN ACROSS PERMUTATIONS
NbAssocGenePerm=do.call(cbind,NbAssocGenePerm)
NbFP_gene=apply(NbAssocGenePerm,1,mean)

# compute FDR for each threshold and 
FDR_threshold=cbind(thresholds,NbFP_gene = NbAssocGene, NbFP_gene= NbFP_gene, FDR_gene=NbFP_Gene/NbAssocGene)
FDR_threshold=as.data.frame(FDR_threshold)

# assign an FDR for each p-value through linear extrapolation on the log10 scale (for FDR < 1E-6, FDR is extrapolated proportionally to the p-value. ) 
FDR_compute=function(pval,FDR_th){y=approxfun(c(0,-log10(FDR_threshold$thresholds),500),c(pmin(-log10(c(1,FDR_th)),6),500))(-log10(pval)); 10^-y}
RESobs$FDR_1Mb=FDR_compute(RESobs$pvalue,FDR_threshold[['FDR_gene']])

# store observed associations and FDR threshold matrix
save(RESobs,FDR_threshold,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V8.Rdata',EVO_IMMUNO_POP))

luq(RESobs$gene[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]) # 993
luq(RESobs$event_id[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]) # 1464

RESobs_nodup_1Mb=RESobs[which(RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable),]
RESobs_nodup_1Mb=RESobs_nodup_1Mb[order(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$pval),]
RESobs_nodup_1Mb=RESobs_nodup_1Mb[!duplicated(RESobs_nodup_1Mb$event_id),]

###################################
###		add response data		###
###################################

# read all snp-response associations in Cis
perm=0
thresholds=10^-c(seq(3,13,0.1),15,20,30,50) # thresholds used to compute FDR
RES=list()
for(pop in c('ALL')){
	for(cond in 1:5){
		for( CHR in 1:22){
		cat(pop, cond,CHR,'\n')
			load(paste(HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_KW_LOGIT_misoPSI_Ens70_HISAT_response_V7.Rdata',sep=''))
#			load(paste(HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_response_V8.Rdata',sep=''))
			RES[[paste(cond,pop,CHR)]]=RESCIS
		}
	}
}
RES=do.call(rbind,RES)
colnames(RES)[colnames(RES)=='gene']='event_id'
RESobs_resp=RES
rm(RES)
gc()

save(RESobs_resp,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_response_V8.Rdata',EVO_IMMUNO_POP))

RESobs_has_rsQTL_1Mb=RESobs[which(abs(RESobs$CisDist)<1e6 & paste(RESobs$event_id,RESobs$snps)%in%paste(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$snps)),]

Pval_rsQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='pvalue');
Pval_rsQTL=Pval_rsQTL[match(RESobs_nodup_1Mb$event_id,Pval_rsQTL$event_id),-1]
Pval_rsQTL[is.na(Pval_rsQTL)]=1
colnames(Pval_rsQTL)=paste('Pvalue_rsQTL',condIndex[-1],sep='_')
Beta_rsQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='beta')
Beta_rsQTL=Beta_rsQTL[match(RESobs_nodup_1Mb$event_id,Beta_rsQTL$event_id),-1]
Beta_rsQTL[is.na(Beta_rsQTL)]=0
colnames(Beta_rsQTL)=paste('Beta_rsQTL',condIndex[-1],sep='_')

RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,Pval_rsQTL,Beta_rsQTL)

########### Annotate sQTLs 
# retrieve statistics and cast them by cond 
RESobs_has_sQTL_1Mb=RESobs[which(abs(RESobs$CisDist)<1e6 & paste(RESobs$event_id,RESobs$snps)%in%paste(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$snps)),]

#Pvalues all conditions
Pval_sQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='pvalue');
Pval_sQTL[is.na(Pval_sQTL)]=1
Pval_sQTL=Pval_sQTL[match(RESobs_nodup_1Mb$event_id,Pval_sQTL$event_id),-1]
colnames(Pval_sQTL)=paste('Pvalue',condIndex,sep='_')

#Beta_values all conditions
Beta_sQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='beta')
Beta_sQTL[is.na(Beta_sQTL)]=0
Beta_sQTL=Beta_sQTL[match(RESobs_nodup_1Mb$event_id,Beta_sQTL$event_id),-1]
colnames(Beta_sQTL)=paste('Beta',condIndex,sep='_')

#R2 all conditions
R2_sQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='R2');
R2_sQTL[is.na(R2_sQTL)]=0
R2_sQTL=R2_sQTL[match(RESobs_nodup_1Mb$event_id,R2_sQTL$event_id),-1]
colnames(R2_sQTL)=paste('R2',condIndex,sep='_')

Testable_sQTL=TESTABLE[match(RESobs_nodup_1Mb$event_id,PSI_Annot[toKeep,'event_id']),]

Expressed_sQTL=(PCTNA<0.05 & SUPPORT>10)[match(RESobs_nodup_1Mb$event_id,PSI_Annot[toKeep,'event_id']),]
colnames(Expressed_sQTL)=paste('Expressed',condIndex,sep='_')

Testable_DIFF_sQTL=TESTABLE_DIFF[match(RESobs_nodup_1Mb$event_id,PSI_Annot[toKeep,'event_id']),-1]

RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,Pval_sQTL,Beta_sQTL,Testable_sQTL,Expressed_sQTL,Testable_DIFF_sQTL)

###################################
###		cluster linked sQTLs	###
###################################

getGenos=function(snps){
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### TODO write the function that loads a specific set of SNPs (without using the SQL database) ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
}

Geno_sQTL=getGenos(unique(RESobs_nodup_1Mb$snps))
r2_sQTL=cor(t(as.matrix(Geno_sQTL[-(1:5)])),use='p')^2
library(igraph)
cor_graph=graph_from_adjacency_matrix(r2_sQTL>0.8)
wc <- cluster_walktrap(cor_graph)
haplos=membership(wc)
names(haplos)=Geno_sQTL$snp.name
RESobs_nodup_1Mb$haplo=paste('haplo',haplos)[match(RESobs_nodup_1Mb$snps,names(haplos))]


###################################
###		Figure 4A	distance	###
###################################
HLA_region=c('6:28Mb','6:29Mb','6:30Mb','6:31Mb','6:32Mb','6:33Mb')
mean(abs(RESobs_nodup_1Mb$CisDist)<1e6) # 66%
mean(abs(RESobs_nodup_1Mb$CisDist[RESobs_nodup_1Mb$locus!=HLA_region])<1e6) #66% 

plot(RESobs_nodup_1Mb$CisDist/1e6,-log10(RESobs_nodup_1Mb$pval),pch=16,col=gsub('AA','44',colERC5[1]),xlim=c(-1,1),las=1,axes=F,ylab=expression(-log[10](P-value)),xlab='Distance from Event boundaries')
axis(2,las=1);
axis(1,at=c(-1,-0.5,0,0.5,1),labels=c('1Mb','500kb','0','500kb','1Mb'))


######################################################
#########           Annotation of sQTL        ########
######################################################

mm=match(RESobs_nodup_1Mb$snps,Map_imputed$snp.name)
#### selection
RESobs_nodup_1Mb$FST=Map_imputed$FST_adj[mm]
RESobs_nodup_1Mb$iHS_AFB=Map_imputed$iHS_AFB[mm]
RESobs_nodup_1Mb$iHS_EUB=Map_imputed$iHS_EUB[mm]
#RESobs_nodup_1Mb$maf_AFB=Map_imputed$maf_AFB[mm]		
#RESobs_nodup_1Mb$maf_EUB=Map_imputed$maf_EUB[mm]

#### mechanisms
RESobs_nodup_1Mb$RegElt=Map_imputed$RegElt[mm]
RESobs_nodup_1Mb$SPIDEX_Z=Map_imputed$SPIDEX_Z[mm]
RESobs_nodup_1Mb$Exonic=Map_imputed$Exonic[mm]
RESobs_nodup_1Mb$groupSplice=Map_imputed$groupSplice[mm]

#### is Immune gene
RESobs_nodup_1Mb$isImmune=RESobs_nodup_1Mb$gene%in%GoToTest_genes$immuneResponse

##########################################################################################
#####		Mechanistic Bases  of splicing	Fig 4B, 4C and Supplementary Figs		  ####
##########################################################################################


tab=table(Map_imputed$groupSplice[wCis])
ImpactOnSplicing=cbind(tab,sapply(c(1e-5,1e-10,1e-30),function(th){tab=table(RESobs_nodup_1Mb$groupSplice[RESobs_nodup_1Mb$pval<th])}))
Proportion_ImpactOnSplicing=t(t(ImpactOnSplicing[-1,])/apply(ImpactOnSplicing,2,sum))
colnames(Proportion_ImpactOnSplicing)=c('All','sQTL, P<10-5','sQTL, P<10-10','sQTL, P<10-30')

###################################
###		Figure 4B 				###
###################################

library(RBrewerPalette)
acol= c(rev(brewer.pal(4,'YlOrRd')),brewer.pal(4,'YlGnBu'))
rcol = SetAlpha(acol, 0.5)
par(mar=c(7,5,1,3))
barplot(Proportion_ImpactOnSplicing,col=acol[c(3,1,6,7)],las=3,ylim=c(0,0.4),ylab='% of splice QTL')

############################################################################################################################################
###################################          Enrichment in regions predicted to alter splicing           ###################################
############################################################################################################################################

OR=NULL
for (j in 2:4){
	tab=rbind(c(ImpactOnSplicing[1,1]-ImpactOnSplicing[1,j],sum(ImpactOnSplicing[1,j])),
		c(sum(ImpactOnSplicing[-1,1])-sum(ImpactOnSplicing[-1,j]),sum(ImpactOnSplicing[-1,j])))
	myOR=odds.ratio(tab)
	myOR$event_type="All"
	myOR$Splice_type=colnames(Proportions_Exon)[j]
	OR=rbind(OR,myOR)
	tab=rbind(c(sum(ImpactOnSplicing[c(2,4),1])-sum(ImpactOnSplicing[c(2,4),j]),sum(ImpactOnSplicing[c(2,4),j])),
			c(sum(ImpactOnSplicing[c(3,5),1])-sum(ImpactOnSplicing[c(3,5),j]),sum(ImpactOnSplicing[c(3,5),j])))
	myOR=odds.ratio(tab)
	myOR$event_type='All'
	myOR$Splice_type=paste('deleterious',colnames(Proportions_Exon)[j])
	OR=rbind(OR,myOR)
	}
write.table(OR,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/SPIDEX/SpliceElt_OR.txt',HOME),sep='\t',quote=F,row.names=F)

save(RESobs_nodup_1Mb,RESobs_nodup_1Mb_cond,file=sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup.Rdata',EVO_IMMUNO_POP))

##########################################################################################
#####					GO enrichment of sQTLs										  ####
##########################################################################################

# sQTls are enriched in defense response genes
resGO=GOSeq(unique(RESobs_nodup_1Mb$gene),unique(PSI_Annot$gene_id),overOnly=F,addCI=T)
resGO$cond='ALL_POOLED'
resGO_cond=lapply(1:5,function(i){GOSeq(unique(RESobs_nodup_1Mb_cond$gene[RESobs_nodup_1Mb_cond$cond==i]),unique(PSI_Annot$gene_id[TESTABLE[,i]]),overOnly=F,addCI=T)})
resGO_cond_table=do.call(rbind,lapply(1:5,function(i){cbind(resGO_cond[[i]],cond=rep(condIndex[i],nrow(resGO_cond[[i]])))}))
TableS3C=rbind(resGO_cond_table,resGO)
write.table(TableS3C,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/TableS3C.txt',HOME),quote=F,sep='\t',row.names=F)

HLA_region=c('6:28Mb','6:29Mb','6:30Mb','6:31Mb','6:32Mb','6:33Mb')
resGO=GOSeq(unique(RESobs_nodup_1Mb$gene[!RESobs_nodup_1Mb$locus%in%HLA_region)]),unique(PSI_Annot$gene_id),overOnly=F,addCI=T)

MeanExpr=log2(1+GeneAnnot[,c("NS_mean", "LPS_mean", "PAM3_mean", "R848_mean", "Flu_mean")])
Expr_level=apply(MeanExpr,1,mean,na.rm=T)
Expr_level[is.na(Expr_level)]=0
Expr_bin=cut(Expr_level,c(-0.01,1,5,10,50,100,500,1000,Inf))

has_sQTL=unique(PSI_Annot$gene_id[toKeep])%in%unique(RESobs_nodup_1Mb$gene)
isImmune=unique(PSI_Annot$gene_id[toKeep])%in%GoToTest_genes$immuneResponse
Expr_bin_test=Expr_bin[match(unique(PSI_Annot$gene_id[toKeep]),GeneAnnot$Ensembl.Gene.ID)]

PctWithsQTL=c()
for (nsamp in 1:10000){
	if(nsamp%%10000==1){ cat(nsamp,'k done\n')}
	samp=c()
	for (x in levels(Expr_bin_test)){
		samp=c(samp,sample(which(Expr_bin_test==x),sum(isImmune[Expr_bin_test==x]),replace=F))
		}
	PctWithsQTL[nsamp]=mean(has_sQTL[samp])
	}

mean(PctWithsQTL>=mean(has_sQTL[isImmune]))
#1000 resample: <0.001
#10000 resample: <0.0001



###################################
###		add eQTL overlap 		###
###################################

 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### TODO UPDATE THIS SECTION TO INCLUDE THE LAST VERSION    ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

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
###		Figure 4C				###
###################################

hist(pmax(0,r2_sQTL_eQTL),main='',xlab='r2 with eQTL',col='darkgrey',br=seq(0,1,0.1),border='darkgrey',las=2)
hist(r2_sQTL_eQTL,add=T,br=seq(0,1,0.1),col='lightgrey',border='darkgrey')
hist(r2_sQTL_eQTL[r2_sQTL_eQTL>0.8],add=T,br=seq(0,1,0.1),col=colERC5[2],border='darkgrey')
barplot(By(r2_sQTL_eQTL>0.8,RESobs_nodup_1Mb$event_type,mean),col=colPSI)

###################################
###		Figure 4C	END			###
###################################

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


 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### TODO UPDATE THIS SECTION TO INCLUDE THE LAST VERSION    ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######


#####################################################################################################################################################################################
#####           Selection analyses                                                                                                                                              #####
#####################################################################################################################################################################################


library(data.table)
#Map_imputed=as.data.frame(fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/selection/resampling_selection/Map_selection.txt',HOME)))
Map_select=as.data.frame(fread(sprintf('%s/Annotation/Select/Map_Select_LD_allChr_V2.txt',HOME)))
#colnames(Map_select)
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

 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### TODO UPDATE THIS SECTION TO INCLUDE THE LAST VERSION    ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######



##########################################################################################
#####				frequency and sharing of sQTLs	across conditions 			      ####
##########################################################################################


################ PIECHART sQTL SHARING
sQTL_sharing_1E7=apply(RESobs_nodup_1Mb[,paste('Pvalue_',condIndex,sep='')]<4.74e-7,1,sum)

expressed_sharing=apply(RESobs_nodup_1Mb[,paste('Expressed_',condIndex,sep='')],1,sum)

Genos_sQTL=read.table(file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/SNP_genos_sQTL_bestSNPperCondition.txt',sep=''),sep='\t',header=T)
snps_sQTL=as.matrix(Genos_sQTL[-1])
rownames(snps_sQTL)=gsub('.',':',Genos_sQTL[[1]],fixed=T)
snps_sQTL=t(snps_sQTL[match(RESobs_nodup_1Mb$snps,rownames(snps_sQTL)),])

psi_sQTL=mapply(function(event,cond){PSI_prov[match(event,PSI_Annot[,'event_id']),match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(PSI_prov))]},RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$cond)
gene_sQTL=mapply(function(gene,cond){FPKM_gene[gene,match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(FPKM_gene))]},RESobs_nodup_1Mb$gene,RESobs_nodup_1Mb$cond)


myPSIs=t(PSI_prov[match(RESobs_nodup_1Mb$event_id[expressed_sharing==5],PSI_Annot[,'event_id']),])
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

####################################################################################################################################
############################################ END sQTL SHARING ######################################################################
####################################################################################################################################

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





