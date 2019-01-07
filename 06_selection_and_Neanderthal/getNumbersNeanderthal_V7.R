library(data.table)

Map_select=as.data.frame(fread(sprintf('%s/Annotation/Select/Map_Select_LD_allChr_V3.txt',HOME)))
Map_aSNP=as.data.frame(fread(sprintf('%s/Martin/Project_Neanderthal/SNPs_all_informations/SNPs_results/CEU/snps.NeanderthalLikeInfo.tsv',EVO_IMMUNO_POP)))

sum(Map_aSNP$aSNP_micha)
#199926
sum(Map_aSNP$NeanderthalLikeSNPs)
# 237422

sum(Map_aSNP$aSNP_micha)/sum(Map_aSNP$NeanderthalLikeSNPs)
# 84%
table(Map_aSNP$genotype_Altai)
#    ./.     0/0     0/1     1/1     1/2 
#  80357 9644041  196683 2390507     915 
mean(Map_aSNP$first_allele==Map_aSNP$REF_Altai_Allele,na.rm=T)
#1
Map_aSNP[which(!Map_aSNP$ALTAI_Alt_Correct)[1:10],]
#     chromosome position first_allele second_allele first_allele_freq        snp_id YRI_Major_SNP YRI_Major_SNP_freq REF_Altai_Allele Alt_Altai_Allele genotype_Altai ALTAI_REF_Correct ALTAI_Alt_Correct NeanderthalLikeSNPs aSNP_micha
#1311           1   899928            G             C         0.0757576  1_899928_G_C             C           0.527778                G                T            0/1              TRUE             FALSE               FALSE      FALSE
#2403           1  1069442            C             T         0.7424240 1_1069442_C_T             C           0.546296                C                G            1/1              TRUE             FALSE               FALSE      FALSE
#2909           1  1137706            C             G         0.9040400 1_1137706_C_G             C           0.560185                C              G,A            1/2              TRUE             FALSE               FALSE      FALSE
#3546           1  1233085            A             T         0.9646460 1_1233085_A_T             T           0.680556                A                G            0/1              TRUE             FALSE               FALSE      FALSE
#4203           1  1354060            G             A         0.9747470 1_1354060_G_A             G           0.981481                G                T            0/1              TRUE             FALSE               FALSE      FALSE
#5903           1  1650529            G             T         0.9949490 1_1650529_G_T             G           1.000000                G                A            0/1              TRUE             FALSE               FALSE      FALSE
#6681           1  1840673            G             A         0.4797980 1_1840673_G_A             G           0.995370                G                T            0/1              TRUE             FALSE               FALSE      FALSE
#8260           1  2127434            C             A         0.6161620 1_2127434_C_A             A           0.731481                C                G            1/1              TRUE             FALSE               FALSE      FALSE
#8832           1  2246193            G             C         0.9848480 1_2246193_G_C             G           1.000000                G                A            1/1              TRUE             FALSE               FALSE      FALSE
#10806          1  2585208            C             A         0.9444440 1_2585208_C_A             C           0.990741                C                G            0/1              TRUE             FALSE               FALSE      FALSE

sum(!Map_aSNP$ALTAI_Alt_Correct,na.rm=T)
# 5645 variant chez ALTAI et en Europe mais avec des alleles différents

table(Map_aSNP$aSNP_micha,Map_aSNP$NeanderthalLikeSNPs)  
#aSNP/nlSNP        FALSE     TRUE
#  FALSE 		12162306    37846
#  TRUE  		     350   199576

Map_aSNP[which(Map_aSNP$aSNP_micha & !Map_aSNP$NeanderthalLikeSNPs)[1:10],]
table(Map_aSNP[which(Map_aSNP$aSNP_micha & !Map_aSNP$NeanderthalLikeSNPs),'genotype_Altai'],exclude='')
# 0/0  1/1 <NA> 
#   1  100  249 
Map_aSNP[which(Map_aSNP$aSNP_micha & !Map_aSNP$NeanderthalLikeSNPs & Map_aSNP$genotype_Altai=='0/0'),]
#        chromosome position first_allele second_allele first_allele_freq         snp_id YRI_Major_SNP YRI_Major_SNP_freq REF_Altai_Allele Alt_Altai_Allele genotype_Altai ALTAI_REF_Correct ALTAI_Alt_Correct NeanderthalLikeSNPs aSNP_micha
#7151031          9 96000622            C             T         0.0454545 9_96000622_C_T             T          0.9583333                C                .            0/0              TRUE              TRUE               FALSE       TRUE


Map_aSNP$NAF_CEU=pmin(1-Map_aSNP$first_allele_freq,Map_aSNP$first_allele_freq)
hist(Map_aSNP$NAF_CEU[Map_aSNP$NeanderthalLikeSNPs])
Map_aSNP[which(Map_aSNP$NeanderthalLikeSNPs & Map_aSNP$NAF_CEU>0.9)[1:10],]
pos_check=paste(Map_aSNP[which(Map_aSNP$NeanderthalLikeSNPs & Map_aSNP$NAF_CEU>0.5),'chromosome'],Map_aSNP[which(Map_aSNP$NeanderthalLikeSNPs & Map_aSNP$NAF_CEU>0.5),'position'])
Map_select[which(paste(Map_select$chromosome,Map_select$position)%in%pos_check)[1:10],]
1/min(Map_aSNP$NAF_CEU[Map_aSNP$NAF_CEU>0])
#0.00505051
N_hapCEU=198
barplot(table(round(Map_aSNP$NAF_CEU[Map_aSNP$NeanderthalLikeSNPs]*N_hapCEU)))
barplot(table(round(Map_aSNP$NAF_CEU[Map_aSNP$NeanderthalLikeSNPs & Map_aSNP$aSNP_micha]*N_hapCEU)),add=T,col=colERC5[2])

mean(round(Map_aSNP$NAF_CEU[Map_aSNP$NeanderthalLikeSNPs & Map_aSNP$aSNP_micha]* N_hapCEU)==1)
# 0.1393153
mean(round(Map_aSNP$NAF_CEU[Map_aSNP$NeanderthalLikeSNPs & !Map_aSNP$aSNP_micha]* N_hapCEU)==1)
# 0.2707552


Map_select$nlSNP=paste(Map_select$chromosome,Map_select$position)%in%paste(Map_aSNP[which(Map_aSNP$NeanderthalLikeSNPs),'chromosome'],Map_aSNP[which(Map_aSNP$NeanderthalLikeSNPs),'position'])
hh=hist(log10(Map_select$haploLength_cM_CEU[Map_select$nlSNP]),br=100)
hist(log10(Map_select$haploLength_cM_CEU[Map_select$aSNP & Map_select$nlSNP]),add=T,col=colERC5[3],br=hh$br)
hist(log10(Map_select$haploLength_cM_CEU[!Map_select$aSNP & Map_select$nlSNP]),add=T,col=colERC5[2],br=hh$br)

hh=hist(log10(Map_select$haploLength_EUB[Map_select$nlSNP]),br=100)
hist(log10(Map_select$haploLength_EUB[Map_select$aSNP & Map_select$nlSNP]),add=T,col=colERC5[3],br=hh$br)
hist(log10(Map_select$haploLength_EUB[!Map_select$aSNP & Map_select$nlSNP]),add=T,col=colERC5[2],br=hh$br)

par(mar=c(10,4,2,2))
layout(matrix(1:2,1))
boxplot(log10(Map_select$haploLength_EUB[Map_select$nlSNP])~ifelse(Map_select$aSNP[Map_select$nlSNP],'overlaps Neanderthal\nhaplotype','no overlaping\nNeanderthalhaplotype'),col=colERC5[1:2],las=2,notch=T,ylab='Nb SNP per Kb')
boxplot(log10(Map_select$haploLength_cM_CEU[Map_select$nlSNP])~ifelse(Map_select$aSNP[Map_select$nlSNP],'overlaps Neanderthal\nhaplotype','no overlaping \nNeanderthal haplotype'),col=colERC5[1:2],las=2,notch=T,ylab='log10 haplo length (cM)')


############### Enrichment in aSNPs ##############

load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/cis-psiQTL_MISO_counts_V7_nodup_withSelect_RBP_coding_intron.Rdata',HOME))
Require('local_eQTL_list')

Map_all=list()
pop='EUB'
for (CHR in 1:22){
    tim=Sys.time()
    RR_CEU=as.data.frame(fread(sprintf('%s/Annotation/RecombinationRate/CEU/CEU-%s-final.txt',HOME,CHR)))
    colnames(RR_CEU)=make.names( colnames(RR_CEU))
	load(sprintf('%s/Annotation/Select/ByCHR/Map_Select_LD_chr%s_V7.Rdata',HOME,CHR))
 	Map=Map[,c('snp.name','chromosome','position','SNPfreq','maf_EUB','maf_AFB','allele.1','allele.2','ancestral_allele','NbLD_SNP_EUB','aSNP_R2_EUB','aSNP','haploLength_EUB')]
    cat('loading Map EUB',CHR,'\n')
    print(Sys.time()-tim)
    LD_table=as.data.frame(fread(paste(EVO_IMMUNO_POP,'/Maxime/GenotypeCalls/Omni5/plink_manip/LD_compute/Imputed_LD_list_80_',pop,'_',CHR,'.ld',sep='')))
    LD_table$maf_1=Map$maf_EUB[match(LD_table$SNP_A,Map$snp.name)]
    LD_table$maf_2=Map$maf_EUB[match(LD_table$SNP_B,Map$snp.name)]
    LD_table=LD_table[which(!is.na(LD_table$maf_2 & LD_table$maf_1)),]
    print(Sys.time()-tim)
# Only SNPs with maf>=0.05 are considered)
# 	}
# Map=do.call(rbind,Map_all)
    cat('haploLength_cM_CEU\n')
    Map$haploLength_cM_CEU=NA
    approx_cM=approx(RR_CEU$Position.bp., y =RR_CEU$Map.cM., Map$position)$y
    
   haploLength=By(c(approx_cM,approx_cM[match(LD_table$SNP_A,Map$snp.name)],approx_cM[match(LD_table$SNP_B,Map$snp.name)]),c(Map$snp.name,LD_table$SNP_B,LD_table$SNP_A),function(x){diff(range(x,na.rm=T))})
    wSNP=match(Map$snp.name,names(haploLength))
    Map$haploLength_cM_CEU=haploLength[wSNP]
    Map$approx_cM=approx_cM
    print(Sys.time()-tim)

    cat('haploLength_EUB_archaic\n')
    Map=as.data.table(Map)
    Map_LD=rbind(cbind(Map[,mget(c('aSNP','position','approx_cM','snp.name'))],matchSNP=Map$snp.name),cbind(Map[match(LD_table$SNP_B,Map$snp.name),mget(c('aSNP','position','approx_cM','snp.name'))],matchSNP=LD_table$SNP_A),cbind(Map[match(LD_table$SNP_A,Map$snp.name),mget(c('aSNP','position','approx_cM','snp.name'))],matchSNP=LD_table$SNP_B))
    haploLength_archaic=Map_LD[aSNP==1,.(Nb_archaic=length(unique(snp.name)),haploLength_archaic=diff(range(position)),haploLength_archaic_cM=diff(range(approx_cM))),by=matchSNP]
    wSNP=match(haploLength_archaic$matchSNP,Map$snp.name)
    Map$haploLength_archaic=0
    Map$haploLength_archaic[wSNP]=haploLength_archaic[,get('haploLength_archaic')]
    Map$haploLength_archaic_cM=0
    Map$haploLength_archaic_cM[wSNP]=haploLength_archaic[,get('haploLength_archaic_cM')]
    Map$Nb_archaic=0
    Map$Nb_archaic[wSNP]=haploLength_archaic[,get('Nb_archaic')]
    Map$locus=paste(Map$chromosome,':',round(Map$position/1e6),'Mb',sep='')
    Map_all[[CHR]]=Map
    print(Sys.time()-tim)
    cat('\n')
}

Map_select=do.call(rbind,Map_all)
tab=table(cut(Map_select$haploLength_archaic[Map_select$aSNP],c(-1,0,1e3,1e4,2e4,5e4,1e5,Inf)),pmin(20,Map_select$Nb_archaic[Map_select$aSNP]),exclude=F)
tab
#                    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20
#  (-1,0]         2721     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
#  (0,1e+03]         0   579   169    71    15    12     0     8     9     0     0     0     0     0     0     0     0     0     0     0
#  (1e+03,1e+04]     0   841  1178  1015   805   558   362   309   215   141   139    93    37    72    66    36    77    36    19   144
#  (1e+04,2e+04]     0   132   341   565   698   716   723   536   600   432   418   261   256   188   120    84   174   194   114   937
#  (2e+04,5e+04]     0    48   131   252   350   670   668   894  1229   810  1008   824   810   867   891   744   659   498   399  5734
#  (5e+04,1e+05]     0     9    19    15    61    77    96   148   157   207   196   311   247   543   480   744   509   545   395 12190
#  (1e+05,Inf]       0     6     4     8    16     1     5    30    19    61    54    26    46    76    53    97   117   174   139 25287

tag_aSNP=intersect(keep_EUB,Map_select$snp.name[Map_select$aSNP_R2_EUB])
tab=table(cut(Map_select$haploLength_archaic[match(tag_aSNP,Map_select$snp.name)],c(-1,0,1e3,1e4,2e4,5e4,1e5,Inf)),pmin(20,Map_select$Nb_archaic[match(tag_aSNP,Map_select$snp.name)]),exclude=F)
tab

Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME))

library(GenomicRanges)
GRange_Map=makeGRangesFromDataFrame(Map_select[which(Map_select$SNPfreq>0.05 ),c('chromosome','position','allele.1','allele.2','ancestral_allele')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_PSI=makeGRangesFromDataFrame(PSI_Annot[toKeep,],seqnames.field='chrom',start.field='start',end.field='end',keep.extra.columns=TRUE)
Cisdist=1e6
#GR_cis=flank(flank(GRange_PSI,Cisdist),-2*Cisdist-1)
GR_cis=union(flank(GRange_PSI,Cisdist,start=T),union(flank(GRange_PSI,Cisdist,start=F),GRange_PSI))
ooCis=findOverlaps(GR_cis, GRange_Map)
wCis=which(Map_select$SNPfreq>0.05)[unique(subjectHits(ooCis))]

save(Map_select,keep_EUB,wCis,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/Map_Neanderthal.Rdata',HOME))
mm=match(RESobs_nodup_1Mb$snps,Map_select$snp.name)
RESobs_nodup_1Mb$Nb_archaic=Map_select$Nb_archaic[mm]
RESobs_nodup_1Mb$haploLength_archaic=Map_select$haploLength_archaic[mm]
RESobs_nodup_1Mb$haploLength_archaic_cM=Map_select$haploLength_archaic_cM[mm]
RESobs_nodup_1Mb$haploLength_EUB=Map_select$haploLength_EUB[mm]
RESobs_nodup_1Mb$haploLength_cM_CEU=Map_select$haploLength_cM_CEU[mm]


MeanExpr=GeneAnnot[grep('_mean',colnames(GeneAnnot))]
MeanExpr[is.na(MeanExpr)]=0
colnames(MeanExpr)=paste('FPKM',condIndex,sep='_')
lFC=log2(1+MeanExpr[,-1])-log2(1+MeanExpr[,1])%o%c(1,1,1,1)
colnames(lFC)=paste('log2FC',condIndex[-1],sep='_')

RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,MeanExpr[RESobs_nodup_1Mb$gene,])
RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,lFC[RESobs_nodup_1Mb$gene,])
RESobs_nodup_1Mb$max_lFC=pmax(RESobs_nodup_1Mb$log2FC_LPS,RESobs_nodup_1Mb$log2FC_PAM3CSK4,RESobs_nodup_1Mb$log2FC_R848,RESobs_nodup_1Mb$log2FC_IAV)

GeneAnnot=cbind(GeneAnnot,lFC,max_lFC=apply(lFC,1,max))
sum(GeneAnnot[unique(PSI_Annot$gene[toKeep]),'NS_mean']<10 & GeneAnnot[unique(PSI_Annot$gene[toKeep]),'max_lFC']>1) # 274

sum(RESobs_nodup_1Mb$max_lFC>1 & RESobs_nodup_1Mb$FPKM_NS<10)
#[1] 108
length(unique(RESobs_nodup_1Mb$gene[RESobs_nodup_1Mb$max_lFC>1 & RESobs_nodup_1Mb$FPKM_NS<10]))
# 74
sum(RESobs_nodup_1Mb$FPKM_NS>10 & RESobs_nodup_1Mb$Pvalue_NS>0.001)


RESobs_nodup_1Mb$is_sQTL_of_stimInducedGene=ifelse(RESobs_nodup_1Mb$max_lFC>1 & RESobs_nodup_1Mb$FPKM_NS<10,'yes','no')
RESobs_nodup_1Mb$rsQTL_tested=ifelse(RESobs_nodup_1Mb$FPKM_NS>10 & RESobs_nodup_1Mb$Pvalue_NS>0.001,'yes','no')
P_val_rsQTL=RESobs_nodup_1Mb[,c('Pvalue_rsQTL_LPS','Pvalue_rsQTL_PAM3CSK4','Pvalue_rsQTL_R848','Pvalue_rsQTL_IAV')]
sum(P_val_rsQTL[RESobs_nodup_1Mb$rsQTL_tested=='yes',]<0.01) # 172

RESobs_nodup_1Mb$is_rsQTL=ifelse(RESobs_nodup_1Mb$rsQTL_tested=='no','na',ifelse(apply(P_val_rsQTL, 1 , min)<0.01,'yes','no'))
table(RESobs_nodup_1Mb$is_rsQTL)
save(RESobs_nodup_1Mb,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/cis-psiQTL_MISO_counts_V7_nodup_withSelect_RBP_coding_intron_aSNP_rsQTL.Rdata',HOME))

table(RESobs_nodup_1Mb$is_rsQTL, RESobs_nodup_1Mb$GWAS_Trait_R2_1E5!='' )
table(RESobs_nodup_1Mb$is_sQTL_of_stimInducedGene, RESobs_nodup_1Mb$GWAS_Trait_R2_1E5!='' )


NSAMP=1000
Resamp=list()
pop='EUB'
#FILENAME=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/resampling_selection/selection_resamp_V2_%s_cond%s_0.Rdata',HOME,pop,cond)
SNP_Support=wCis[which(Map_select$maf_EUB[wCis]>0.05 & Map_select$snp.name[wCis]%in%keep_EUB)]

# define SNP set 

cond=0

if(cond>0){
	w=which(RESobs_nodup_1Mb_cond$cond==cond & RESobs_nodup_1Mb_cond$maf_EUB>0.05)
	w=w[!duplicated(RESobs_nodup_1Mb_cond$haplo[w])]
	mm=match(unique(RESobs_nodup_1Mb_cond$snps[w]),Map_select$snp.name)
	}else{
	w=which(RESobs_nodup_1Mb$maf_EUB>0.05)
	w=w[!duplicated(RESobs_nodup_1Mb$haplo[w])]
	mm=match(unique(RESobs_nodup_1Mb$snps[w]),Map_select$snp.name)
}

MAF_bin=cut(pmin(Map_select$SNPfreq[SNP_Support],0.5),seq(0,0.5,0.02)) #cut(Map_select$NbLD_SNP_EUB[SNP_Support],c(0,1,2,5,10,20,50,100,1e6)))
# create bins for sQTLs
mmMAF_bin=cut(pmin(Map_select$SNPfreq[mm],0.5),seq(0,0.5,0.02))
binCounts=table(mmMAF_bin) # Nb of SNP to sample for each MAF bin 
binCounts=binCounts[levels(MAF_bin)] 
names(binCounts)=levels(MAF_bin)
binCounts[is.na(binCounts)]=0

#Map_select$locus=paste(Map_select$chromosome,':',round(Map_select$position/1e6),'Mb',sep='')
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
	Resamp[[paste('resamp_EUB',cond,'10kb',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>10000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'20kb',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>20000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'30kb',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>30000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'50kb',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>50000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'10kb_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>10000)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>20000)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>30000)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>50000)]))
	Resamp[[paste('resamp_EUB',cond,'10kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>10000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>20000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>30000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_EUB']>50000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>10000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>20000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>30000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic',sep='_')]]$NbaSNP_R2_EUB[samp]=sum(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>50000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>10000)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>20000)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>30000)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>50000)]))
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>10000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>20000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>30000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>50000 & Map_select[sample_global,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>10000 & Map_select[sample_global,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>20000 & Map_select[sample_global,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>30000 & Map_select[sample_global,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB[samp]=length(unique(Map_select$locus[sample_global][which(Map_select[sample_global,'aSNP_R2_EUB'] & Map_select[sample_global,'haploLength_archaic']>50000 & Map_select[sample_global,'Nb_archaic']>=10)]))
	}

	w=which(RESobs_nodup_1Mb$maf_EUB>0.05)
	w=w[order(-RESobs_nodup_1Mb$haploLength_archaic[w])]
	w=w[!duplicated(RESobs_nodup_1Mb$haplo[w])]
	mm=match(unique(RESobs_nodup_1Mb$snps[w]),Map_select$snp.name)
	Resamp[[paste('resamp_EUB',cond,'10kb',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>10000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'20kb',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>20000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'30kb',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>30000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'50kb',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>50000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'10kb_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>10000)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>20000)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>30000)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>50000)]))
	Resamp[[paste('resamp_EUB',cond,'10kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>10000 & Map_select[mm,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>20000 & Map_select[mm,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>30000 & Map_select[mm,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>50000 & Map_select[mm,'Nb_archaic']>=5)]))
#	Resamp[[paste('resamp_EUB',cond,'10kb_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>10000 & Map_select[mm,'Nb_archaic']>=10)]))
#	Resamp[[paste('resamp_EUB',cond,'20kb_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>20000 & Map_select[mm,'Nb_archaic']>=10)]))
#	Resamp[[paste('resamp_EUB',cond,'30kb_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>30000 & Map_select[mm,'Nb_archaic']>=10)]))
#	Resamp[[paste('resamp_EUB',cond,'50kb_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_EUB']>50000 & Map_select[mm,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>10000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>20000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>30000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic',sep='_')]]$NbaSNP_R2_EUB_sQTL=sum(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>50000,na.rm=T)
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>10000)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>20000)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>30000)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>50000)]))
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>10000 & Map_select[mm,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>20000 & Map_select[mm,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>30000 & Map_select[mm,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus_5aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>50000 & Map_select[mm,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'10kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>10000 & Map_select[mm,'Nb_archaic']>=5)]))
	Resamp[[paste('resamp_EUB',cond,'20kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>20000 & Map_select[mm,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'30kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>30000 & Map_select[mm,'Nb_archaic']>=10)]))
	Resamp[[paste('resamp_EUB',cond,'50kb_archaic_locus_10aSNP',sep='_')]]$NbaSNP_R2_EUB_sQTL=length(unique(Map_select$locus[mm][which(Map_select[mm,'aSNP_R2_EUB'] & Map_select[mm,'haploLength_archaic']>50000 & Map_select[mm,'Nb_archaic']>=10)]))

save(Resamp,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/Neanderthal_resamples_MAFover10pct.Rdata',HOME))

dump=sapply(names(Resamp),function(i){x=Resamp[[i]];cat(i,'N=',x$NbaSNP_R2_EUB_sQTL,'E=',mean(x$NbaSNP_R2_EUB,na.rm=T),'P=',mean(x$NbaSNP_R2_EUB>=x$NbaSNP_R2_EUB_sQTL,na.rm=T),'\n')})
#resamp_EUB_0_10kb N= 12 E= 4.128 P= 0.001 
#resamp_EUB_0_20kb N= 12 E= 3.127 P= 0 
#resamp_EUB_0_30kb N= 11 E= 2.378 P= 0 
#resamp_EUB_0_50kb N= 9 E= 1.504 P= 0 
#resamp_EUB_0_10kb_locus N= 10 E= 4.11 P= 0.009 
#resamp_EUB_0_20kb_locus N= 10 E= 3.117 P= 0.001 
#resamp_EUB_0_30kb_locus N= 9 E= 2.373 P= 0.001 
#resamp_EUB_0_50kb_locus N= 7 E= 1.503 P= 0.001 
#resamp_EUB_0_10kb_locus_5aSNP N= 8 E= 2.622 P= 0.003 
#resamp_EUB_0_20kb_locus_5aSNP N= 8 E= 2.205 P= 0 
#resamp_EUB_0_30kb_locus_5aSNP N= 7 E= 1.794 P= 0.002 
#resamp_EUB_0_50kb_locus_5aSNP N= 5 E= 1.199 P= 0.009 


#resamp_EUB_0_10kb_archaic N= 10 E= 2.847 P= 0.001 
#resamp_EUB_0_20kb_archaic N= 10 E= 2.072 P= 0 
#resamp_EUB_0_30kb_archaic N= 7 E= 1.497 P= 0 
#resamp_EUB_0_50kb_archaic N= 5 E= 0.96 P= 0.001 
#resamp_EUB_0_10kb_archaic_locus N= 8 E= 2.835 P= 0.012 
#resamp_EUB_0_20kb_archaic_locus N= 8 E= 2.067 P= 0.002 
#resamp_EUB_0_30kb_archaic_locus N= 6 E= 1.495 P= 0.006 
#resamp_EUB_0_50kb_archaic_locus N= 4 E= 0.959 P= 0.015 
#resamp_EUB_0_10kb_archaic_locus_5aSNP N= 8 E= 2.48 P= 0.003 
#resamp_EUB_0_20kb_archaic_locus_5aSNP N= 8 E= 1.94 P= 0 
#resamp_EUB_0_30kb_archaic_locus_5aSNP N= 6 E= 1.446 P= 0.005 
#resamp_EUB_0_50kb_archaic_locus_5aSNP N= 4 E= 0.951 P= 0.015 
#resamp_EUB_0_10kb_archaic_locus_10aSNP N= 8 E= 1.614 P= 0 
#resamp_EUB_0_20kb_archaic_locus_10aSNP N= 6 E= 1.464 P= 0.003 
#resamp_EUB_0_30kb_archaic_locus_10aSNP N= 5 E= 1.23 P= 0.007 
#resamp_EUB_0_50kb_archaic_locus_10aSNP N= 3 E= 0.891 P= 0.057 

NeanderthallikeHaplo=length(unique(RESobs_nodup_1Mb$locus[RESobs_nodup_1Mb$haploLength_EUB>1e4 & RESobs_nodup_1Mb$aSNP_R2_EUB])

length(unique(RESobs_nodup_1Mb$locus[RESobs_nodup_1Mb$])
mean(Resamp[[paste('resamp_EUB',cond,'10kb',sep='_')]]$NbaSNP_R2_EUB>=8) # 0.043

mean(Resamp[[paste('resamp_EUB',cond,'20kb',sep='_')]]$NbaSNP_R2_EUB>=8) # 0.006

mean(Resamp[[paste('resamp_EUB',cond,'30kb',sep='_')]]$NbaSNP_R2_EUB>=8) # 0.001

mean(Resamp[[paste('resamp_EUB',cond,'50kb',sep='_')]]$NbaSNP_R2_EUB>=8) # 0



    Map=as.data.table(Map)

    Map_LD=rbind( cbind(Map[,mget(c('aSNP','position','approx_cM','snp.name'))],matchSNP=Map$snp.name),
                  cbind(Map[match(LD_table$SNP_B,Map$snp.name),mget(c('aSNP','position','approx_cM','snp.name'))],matchSNP=LD_table$SNP_A),
                  cbind(Map[match(LD_table$SNP_A,Map$snp.name),mget(c('aSNP','position','approx_cM','snp.name'))],matchSNP=LD_table$SNP_B))
    haploLength_archaic=Map_LD[aSNP==1,.(Nb_archaic=length(unique(snp.name)),haploLength_archaic=diff(range(position)),haploLength_archaic_cM=diff(range(approx_cM))),by=matchSNP]
    wSNP=match(haploLength_archaic$matchSNP,Map$snp.name)
    Map$haploLength_archaic=0
    Map$haploLength_archaic[wSNP]=haploLength_archaic[,get('haploLength_archaic')]
    Map$haploLength_archaic_cM=0
    Map$haploLength_archaic_cM[wSNP]=haploLength_archaic[,get('haploLength_archaic_cM')]
    Map$Nb_archaic=0
    Map$Nb_archaic[wSNP]=haploLength_archaic[,get('Nb_archaic')]
    
    Map$locus=paste(Map$chromosome,':',round(Map$position/1e6),'Mb',sep='')

