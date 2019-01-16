
###############################################################################################################
#####			selection and Neanderthal Tables and enrichments 										  #####
###############################################################################################################
# write.table(Map_imputed[,c(1,2,4:6,29,30,8,15,16,19,20,23:28,37:40,45:52,55:64,82,90)],file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/selection/resampling_selection/Map_selection.txt',HOME),sep='\t',quote=F,row.names=F)
library(data.table)
library(GenomicRanges)

################################################################################################################################################################
########    LD_pruned_SNPs_80_EUB and LD_pruned_SNPs_80_AFB are obtained using plink with the following command                                         ########
########    plink --bfile Genotype_${CHR} --keep {POP}.txt --indep-pairwise 1000 50 0.8 --maf 0.05 --chr ${CHR} --out LD_pruned_SNPs_80_${POP}_${CHR}   ########
########    the output of these command is then concatenated into a single file                                                                         ########
################################################################################################################################################################

LD_pruned_SNPs_80_EUB=fread( 'LD_pruned_SNPs_80_EUB.prune.in',header=FALSE )$V1
LD_pruned_SNPs_80_AFB=fread( 'LD_pruned_SNPs_80_AFB.prune.in',header=FALSE )$V1

load('data/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load('data/PSI_adjust_all_16173_V5.Rdata')
PSI_Annot=PSI_Annot[toKeep,]

GRange_Map=makeGRangesFromDataFrame(Map_select[which(Map_select$SNPfreq>0.05 ),c('chromosome','position','allele.1','allele.2','minor_allele','ancestral_allele')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_PSI=makeGRangesFromDataFrame(PSI_Annot,seqnames.field='chrom',start.field='start',end.field='end',keep.extra.columns=TRUE)
Cisdist=1e6
GR_cis=reduce(union(flank(GRange_PSI,Cisdist,start=T),union(flank(GRange_PSI,Cisdist,start=F),GRange_PSI)))
ooCis=findOverlaps(GR_cis, GRange_Map)
wCis=which(Map_select$SNPfreq>0.05)[unique(subjectHits(ooCis))]

# defined SNP_support (set of sampled SNPs)
wCis_EUB=wCis[which(Map_select$maf_EUB[wCis]>0.05 & Map_select$snp.name[wCis]%in%LD_pruned_SNPs_80_EUB)] # The distribution of |iHS| is defined based on independant SNPs (in the population under consideration)
wCis_AFB=wCis[which(Map_select$maf_AFB[wCis]>0.05 & Map_select$snp.name[wCis]%in%LD_pruned_SNPs_80_AFB)] # The distribution of |iHS| is defined based on independant SNPs (in the population under consideration)
rm(GRange_Map,GRange_PSI,GR_cis)

RESobs_nodup_1Mb_cond=as.data.frame(fread('data/RESobs_nodup_1Mb_cond.txt'))
TableS4=as.data.frame(fread('data/RESobs_nodup_1Mb.txt'))

TableS4$FST_Pemp=sapply(TableS4$FST_adj,function(x){mean(Map_select$FST_adj[wCis]>x)})
TableS4$aIHS_AFB_Pemp=sapply(abs(TableS4$iHS_AFB),function(x){mean(abs(Map_select$iHS_AFB[wCis_AFB])>x,na.rm=T)})
TableS4$aIHS_EUB_Pemp=sapply(abs(TableS4$iHS_EUB),function(x){mean(abs(Map_select$iHS_EUB[wCis_EUB])>x,na.rm=T)})

#TableS4=as.data.frame(fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_diffPSI/TableS2_V3.2_OneLinePerEvent.txt',HOME)))

####################################################################################

select_cols=c('maf_EUB','maf_AFB','FST_adj','iHS_AFB','iHS_EUB','NbSNP_100kb','NbSNP_out95_FST_100kb','NbSNP_out95_aIHS_EUB_100kb','NbSNP_out95_aIHS_AFB_100kb','aSNP','aSNP_R2_EUB','haploLength','haploLength_archaic')
mm=match(TableS4$snps,Map_select$snp.name)
for (i in select_cols){
	TableS4[,x]=Map_select[mm,x]
	}

####################################################################################
####################################################################################
####################################################################################

getEmpiricalPval=function(Nb_informative_SNP_train,NbOutliers_train,Nb_informative_SNP_test,NbOutliers_test){
    # compares the ratio of NbOutliers_test/Nb_informative_SNP_test to the distribution of NbOutliers_train/Nb_informative_SNP_train accounting for uncertainty due to variability in Number of SNPs in the window. 
		library(bbmle)
		library(emdbook)
		library(SuppDists)

	mtmp<- function(prob,theta){-sum(dbetabinom(NbOutliers_train,prob,Nb_informative_SNP_train,theta,log=TRUE))}
	m0<-mle2(mtmp,start=list(prob=0.2,theta=9),data=list(f=30)) # f=30 is just there because we need and argument to the function, it in not actually used.

	prob=m0@coef[1]
	theta=m0@coef[2]
	shape1=prob*theta
	shape2=theta-shape1
	Pval=pghyper(NbOutliers_test,a=-shape1, N=-shape1-shape2,k=Nb_informative_SNP_test,low=F)
	# P-values for Outliers Enrichment were obtained by fitting a beta-binomial (to account for over-dispersion) to the GW distribution of the number of outliers within random 100kb window.
}

# here we subsample 100000SNPs of each category to speed up computation
samp_ALL=sample(wCis,100000)
samp_EUB=sample(wCis_EUB,100000)
samp_AFB=sample(wCis_AFB,100000)

mm_EUB=match(TableS4$snps[TableS4$maf_EUB>0.05],Map_select$snp.name)
mm_AFB=match(TableS4$snps[TableS4$maf_AFB>0.05],Map_select$snp.name)
mm_ALL=match(TableS4$snps,Map_select$snp.name)

# FST enrichment
TableS4$Pct_out95_100kb_FST=TableS4$NbSNP_out95_FST_100kb/TableS4$NbSNP_100kb
TableS4$out100kb_FST_95_Pemp=getEmpiricalPval(Map_select$NbSNP_100kb[samp_ALL],Map_select$NbSNP_out95_FST_100kb[samp_ALL],Map_select$NbSNP_100kb[mm_ALL],Map_select$NbSNP_out95_FST_100kb[mm_ALL])

# iHS enrichment AFB
TableS4$Pct_out95_100kb_aIHS_AFB=TableS4$NbSNP_out95_aIHS_AFB_100kb/TableS4$NbSNP_100kb
TableS4$out100kb_aIHS_95_AFB_Pemp=NA
TableS4$out100kb_aIHS_95_EUB_Pemp[TableS4$maf_EUB>0.05]=getEmpiricalPval(Map_select$NbSNP_daf5_EUB_100kb[samp_EUB],Map_select$NbSNP_out95_aIHS_EUB_100kb[samp_EUB],Map_select$NbSNP_daf5_EUB_100kb[mm_EUB],Map_select$NbSNP_out95_aIHS_EUB_100kb[mm_EUB])

# iHS enrichment EUB
TableS4$Pct_out95_100kb_aIHS_EUB=TableS4$NbSNP_out95_aIHS_EUB_100kb/TableS4$NbSNP_100kb
TableS4$out100kb_aIHS_95_AFB_Pemp=NA
TableS4$out100kb_aIHS_95_AFB_Pemp[TableS4$maf_AFB>0.05]=getEmpiricalPval(Map_select$NbSNP_daf5_AFB_100kb[samp_AFB],Map_select$NbSNP_out95_aIHS_AFB_100kb[samp_AFB],Map_select$NbSNP_daf5_AFB_100kb[mm_AFB],Map_select$NbSNP_out95_aIHS_AFB_100kb[mm_AFB])

write.table(TableS4,file=sprintf('tables/TableS4_Selection.txt',HOME),quote=F,sep='\t',row.names=F)


####################################################################################
#################       Table S4B FST outliers              ########################
####################################################################################

colsFST=c('event_id','symbol','event_type','is_sQTL_of_stimulationSpecificGene','is_rsQTL','snps','cond','beta','FDR','daf_char_AFB','daf_char_EUB','FST_adj','FST_Pemp','NbSNP_100kb','Pct_out95_100kb_FST','out100kb_FST_95_Pemp','coding_type')
TableS4B_FST=TableS4[TableS4$FST_Pemp<=0.05 & TableS4$out100kb_FST_95_Pemp<=0.05,colsFST]
TableS4B_FST=TableS4B_FST[order(TableS4B_FST$FST_Pemp),]
write.table(TableS4B_FST,file='tables/TableS4B_FST.txt',quote=F,sep='\t',row.names=F)

####################################################################################
#################       Table S4C IHS outliers (AFB)        ########################
####################################################################################

colsIHS_AFB=c('event_id','symbol','event_type','is_sQTL_of_stimulationSpecificGene','is_rsQTL','snps','cond','beta','FDR','daf_char_AFB','daf_char_EUB','iHS_AFB','aIHS_AFB_Pemp','NbSNP_100kb','outPct_aIHS_95_AFB','out100kb_aIHS_95_AFB_Pemp','coding_type')
TableS4C_IHS_AFB=TableS4[TableS4$aIHS_AFB_Pemp<=0.05 & TableS4$out100kb_aIHS_95_AFB_Pemp<=0.05,colsIHS_AFB]
TableS4C_IHS_AFB=TableS4C_IHS_AFB[order(TableS4C_IHS_AFB$aIHS_AFB_Pemp),]
write.table(TableS4C_IHS_AFB,file='tables/TableS4C_IHS_AFB.txt',quote=F,sep='\t',row.names=F)

####################################################################################
#################       Table S4D IHS outliers (EUB)        ########################
####################################################################################

colsIHS_EUB=c('event_id','symbol','event_type','is_sQTL_of_stimulationSpecificGene','is_rsQTL','snps','cond','beta','FDR','daf_char_AFB','daf_char_EUB','iHS_EUB','aIHS_EUB_Pemp','NbSNP_100kb','outPct_aIHS_95_EUB','out100kb_aIHS_95_EUB_Pemp','coding_type')
TableS4D_IHS_EUB=TableS4[TableS4$aIHS_EUB_Pemp<=0.05 & TableS4$out100kb_aIHS_95_EUB_Pemp<=0.05,colsIHS_EUB]
TableS4D_IHS_EUB=TableS4D_IHS_EUB[order(TableS4D_IHS_EUB$aIHS_EUB_Pemp),]
write.table(TableS4D_IHS_EUB,file='tables/TableS4D_IHS_EUB.txt',quote=F,sep='\t',row.names=F)

###########################################################################################################################################################
#### Fig S7C was obtained from https://popgen.uchicago.edu/ggv/?data=%221000genomes%22&chr=7&pos=99270539 with avoid overlap option enabled  ##############
###########################################################################################################################################################


###############################################################
#### Fig 6C is obtained directly from BAM files  ##############
###############################################################

################################
########    Fig 6D    ##########
################################
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V8.1.Rdata',EVO_IMMUNO_POP))
RES_CYP3A5=RESobs[RESobs$symbol=='CYP3A5',]

library(data.table)
par(mar=c(3,4,0.5,0.5),xpd=F)
w=which(RES_CYP3A5$cond%in%c(3))
plot(RES_CYP3A5$pos[w],-log10(RES_CYP3A5$pvalue[w]),pch=21,cex=0.2+0.1*-log10(RES_CYP3A5$pvalue[w])/3.5,bg=colERC5[RES_CYP3A5$cond[w]],col='#00000000',axes=F,xlab='',ylab=expression(-log[10](P)),ylim=c(0,25),xlim=c(98200072,100272456));axis(2,las=1);
w=which(Map_select$chromosome==7 & Map_select$position>98270002 & Map_select$position<100272456)
axis(1,las=1,at=c(98.2*1e6,99.2*1e6,100.2*1e6),labels=c(98.2,99.2,100.2));
arrows(99270539,-1,99270539,-log10(0.05),code=2,col=colERC[8],lwd=2,length=0.08)

################################
########    Fig 6E    ##########
################################

plot(Map_select$position[w],Map_select$FST_adj[w],xlim=c(98270002,100272456),pch=21,bg="#A2A2A2AA",cex=0.1+Map_select$FST_adj[w]*0.7,ylab=expression(F[ST]),ylim=c(0,1),axes=F);axis(2,las=1);axis(1,las=1,at=c(98.2*1e6,99.2*1e6,100.2*1e6),labels=c(98.2,99.2,100.2));
lines(Map_select$position[w],Map_select$NbSNP_out95_FST_100kb[w]/Map_select$NbSNP_100kb[w],col=substr(colERC5[2],1,7),lwd=2)
arrows(99270539,0.95,99270539,0.85,code=2,col=colERC[8],lwd=2,length=0.08)


############################################################
############### Annotation of archaic sQTLs ################
############################################################


mm=match(RESobs_nodup_1Mb$snps,Map_select$snp.name)
RESobs_nodup_1Mb$Nb_archaic=Map_select$Nb_archaic[mm]
RESobs_nodup_1Mb$haploLength_archaic=Map_select$haploLength_archaic[mm]
RESobs_nodup_1Mb$haploLength_archaic_cM=Map_select$haploLength_archaic_cM[mm]


##########################################################################
#############       Table S4E: Neanderthal-like SNPs         #############
##########################################################################

cols_aSNP_EUB=c("event_id",'symbol','event_type','haploLength_archaic','maf_EUB','snps','conditon','beta','FDR_1Mb','Nb_archaic','Most_frequent_aSNP','GWAS_Trait_R2_1E5','coding_type')
TableS4E_aSNP=TableS4[which(TableS4$aSNP_R2_EUB),cols_aSNP_EUB]
TableS4E_aSNP$haploLength_archaic=as.numeric(TableS4E_aSNP$haploLength_archaic)/1000
write.table(TableS4E_aSNP,file='tables/TableS4E_aSNP.txt',quote=F,sep='\t',row.names=F)

##################################################################################################################################
#############################################  Enrichment in Neanderthal-like SNPs ###############################################
##################################################################################################################################

GRange_Map=makeGRangesFromDataFrame(Map_select[which(Map_select$SNPfreq>0.05 ),c('chromosome','position','allele.1','allele.2','ancestral_allele')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_PSI=makeGRangesFromDataFrame(PSI_Annot, seqnames.field='chrom', start.field='start', end.field='end', keep.extra.columns=TRUE)
Cisdist=1e6
GR_cis=reduce(union(flank(GRange_PSI,Cisdist,start=T),union(flank(GRange_PSI,Cisdist,start=F),GRange_PSI)))
ooCis=findOverlaps(GR_cis, GRange_Map)
wCis=which(Map_select$SNPfreq>0.05)[unique(subjectHits(ooCis))]

save(Map_select,LD_pruned_SNPs_80_EUB,wCis,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/Map_Neanderthal.Rdata',HOME))

NSAMP=1000
Resamp=list()
pop='EUB'
#FILENAME=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/resampling_selection/selection_resamp_V2_%s_cond%s_0.Rdata',HOME,pop,cond)
SNP_Support=wCis[which(Map_select$maf_EUB[wCis]>0.05 & Map_select$snp.name[wCis]%in%LD_pruned_SNPs_80_EUB)]

# define SNP set 

# cond: selection sQTLs detected in a specific condition .
# set cond = 0, for sQTLs across all conditions. 

if(cond>0){
	w=which(RESobs_nodup_1Mb_cond$cond==cond & RESobs_nodup_1Mb_cond$maf_EUB>0.05)
	w=w[!duplicated(RESobs_nodup_1Mb_cond$haplo[w])]
	mm=match(unique(RESobs_nodup_1Mb_cond$snps[w]),Map_select$snp.name)
	}else{
	w=which(RESobs_nodup_1Mb$maf_EUB>0.05)
	w=w[!duplicated(RESobs_nodup_1Mb$haplo[w])]
	mm=match(unique(RESobs_nodup_1Mb$snps[w]),Map_select$snp.name)
}


MAF_bin=cut(pmin(Map_select$SNPfreq[SNP_Support],0.5),seq(0,0.5,0.02)) 
# create bins for sQTLs
mmMAF_bin=cut(pmin(Map_select$SNPfreq[mm],0.5),seq(0,0.5,0.02))
binCounts=table(mmMAF_bin) # Nb of SNP to sample for each MAF bin 
binCounts=binCounts[levels(MAF_bin)] 
names(binCounts)=levels(MAF_bin)
binCounts[is.na(binCounts)]=0

By = function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x}
#Map_select$locus=paste(Map_select$chromosome,':',round(Map_select$position/1e6),'Mb',sep='')

####### Perform the resamplings 
for (samp in 1:NSAMP){
	tic=Sys.time()
	cat(samp,'')
	sample_global=By(SNP_Support,MAF_bin,sample,max(binCounts)) # over sample each bin
	for (i in 1:length(sample_global)){
		if(binCounts[i]>0){
			sample_global[[names(binCounts)[i]]]=sample_global[[names(binCounts)[i]]][1:binCounts[i]] # sub sample to right amount
		}else{
			sample_global[[names(binCounts)[i]]]=c()
		}
	}
	sample_global=unlist(sample_global)
    ##Note : haplotype length= distance maximal distance between 2 SNPs in LD with the target SNP
    ##       archaic haplotype length = distance maximal distance between 2 aSNPs in LD with the target SNP

    ###### get Resampled values 
    # count number of SNP on an archaic haplotype (length haplotype 10,20,30,50 kb) 
	Resamp[[paste('resamp_EUB',cond,'10kb',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>10000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'20kb',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>20000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'30kb',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>30000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'50kb',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>50000,na.rm=T)
    # count number of loci with 1+ sQTL on an archaic haplotype (length haplotype 10,20,30,50 kb) 
	Resamp[[paste('resamp_EUB',cond,'10kb_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>10000)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>20000)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>30000)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>50000)]))
    # count number of loci with 1+ sQTL on an archaic haplotype with 5+ aSNPs (length haplotype 10,20,30,50 kb) 
	Resamp[[paste('resamp_EUB',cond,'10kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>10000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>20000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>30000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>50000 & Map_select[sample_global,'Nb_archaic']>=5)]))
    # count number of SNP on an archaic haplotype with 2+ aSNPs (length archaic haplotype 10,20,30,50 kb) 
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>10000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>20000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>30000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>50000,na.rm=T)
    # count number of loci with 1+ sQTL on an archaic haplotype with 2+ aSNPs (length archaic haplotype 10,20,30,50 kb) 
    Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>10000)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>20000)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>30000)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>50000)]))
    # count number of loci with 1+ sQTL on an archaic haplotype with 5+ aSNPs (length archaic haplotype 10,20,30,50 kb) 
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>10000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>20000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>30000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>50000 & Map_select[sample_global,'Nb_archaic']>=5)]))
    # count number of loci with 1+ sQTL on an archaic haplotype with 10+ aSNPs (length archaic haplotype 10,20,30,50 kb) 
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>10000 & Map_select[sample_global,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>20000 & Map_select[sample_global,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>30000 & Map_select[sample_global,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>50000 & Map_select[sample_global,'Nb_archaic']>=10)]))
	}

w=which(RESobs_nodup_1Mb$maf_EUB>0.05)
w=w[order(-RESobs_nodup_1Mb$haploLength_archaic[w])]
w=w[!duplicated(RESobs_nodup_1Mb$haplo[w])]
mm=match(unique(RESobs_nodup_1Mb$snps[w]),Map_select$snp.name)

###### get observed values (at sQTLs)
# count number of SNP on an archaic haplotype (length haplotype 10,20,30,50 kb) 
Resamp[[paste('resamp_EUB',cond,'10kb',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>10000,na.rm=T)
Resamp[[paste('resamp_EUB',cond,'20kb',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>20000,na.rm=T)
Resamp[[paste('resamp_EUB',cond,'30kb',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>30000,na.rm=T)
Resamp[[paste('resamp_EUB',cond,'50kb',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>50000,na.rm=T)
# count number of loci with 1+ sQTL on an archaic haplotype (length haplotype 10,20,30,50 kb) 
Resamp[[paste('resamp_EUB',cond,'10kb_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>10000)]))
Resamp[[paste('resamp_EUB',cond,'20kb_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>20000)]))
Resamp[[paste('resamp_EUB',cond,'30kb_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>30000)]))
Resamp[[paste('resamp_EUB',cond,'50kb_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>50000)]))
# count number of loci with 1+ sQTL on an archaic haplotype with 5+ aSNPs (length haplotype 10,20,30,50 kb) 
Resamp[[paste('resamp_EUB',cond,'10kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>10000 & Map_select[mm,'Nb_archaic']>=5)]))
Resamp[[paste('resamp_EUB',cond,'20kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>20000 & Map_select[mm,'Nb_archaic']>=5)]))
Resamp[[paste('resamp_EUB',cond,'30kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>30000 & Map_select[mm,'Nb_archaic']>=5)]))
Resamp[[paste('resamp_EUB',cond,'50kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>50000 & Map_select[mm,'Nb_archaic']>=5)]))
    # count number of SNP on an archaic haplotype with 2+ aSNPs (length archaic haplotype 10,20,30,50 kb) 
Resamp[[paste('resamp_EUB',cond,'10kb_archaic',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>10000,na.rm=T)
Resamp[[paste('resamp_EUB',cond,'20kb_archaic',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>20000,na.rm=T)
Resamp[[paste('resamp_EUB',cond,'30kb_archaic',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>30000,na.rm=T)
Resamp[[paste('resamp_EUB',cond,'50kb_archaic',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>50000,na.rm=T)
# count number of loci with 1+ sQTL on an archaic haplotype with 2+ aSNPs (length archaic haplotype 10,20,30,50 kb) 
Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>10000)]))
Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>20000)]))
Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>30000)]))
Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>50000)]))
# count number of loci with 1+ sQTL on an archaic haplotype with 5+ aSNPs (length archaic haplotype 10,20,30,50 kb) 
Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>10000 & Map_select[mm,'Nb_archaic']>=5)]))
Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>20000 & Map_select[mm,'Nb_archaic']>=5)]))
Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>30000 & Map_select[mm,'Nb_archaic']>=5)]))
Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>50000 & Map_select[mm,'Nb_archaic']>=5)]))
# count number of loci with 1+ sQTL on an archaic haplotype with 10+ aSNPs (length archaic haplotype 10,20,30,50 kb) 
Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>10000 & Map_select[mm,'Nb_archaic']>=10)]))
Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>20000 & Map_select[mm,'Nb_archaic']>=10)]))
Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>30000 & Map_select[mm,'Nb_archaic']>=10)]))
Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>50000 & Map_select[mm,'Nb_archaic']>=10)]))

save(Resamp,file='data/Neanderthal_resamples_MAFover5pct.Rdata')

##################################################
######### print the results (Table S4E)  #########
##################################################

dump=sapply(names(Resamp),function(i){x=Resamp[[i]];cat(i,'N=',x$NbaSNP_R2_EUB_sQTL,'E=',mean(x$NbaSNP_R2_EUB,na.rm=T),'P=',mean(x$NbaSNP_R2_EUB>=x$NbaSNP_R2_EUB_sQTL,na.rm=T),'\n')})

