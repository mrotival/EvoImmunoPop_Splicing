Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/cis-psiQTL_MISO_counts_V7_nodup_withSelect_RBP_coding_intron_aSNP_rsQTL_spliceElt.Rdata',HOME))
library(data.table)
MAP_GWAS=as.data.frame(fread(sprintf('%s/Annotation/GWAS/Map_GWAS_LD_EFOcorrected_allChr.txt',HOME)))

EFO_14=sapply(strsplit(MAP_GWAS$GWAS_EFO_R2_1E5,'//'), function(x){paste(unique(x[!x%in%c('Other trait','Other measurement','Other disease','Biological process')]),collapse='//')})
sum(EFO_14!='') # 229425

#Require('Map_imputed')
library(GenomicRanges)
GRange_Map=makeGRangesFromDataFrame(MAP_GWAS[which(MAP_GWAS$SNPfreq>0.05 ),c('chromosome','position')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_PSI=makeGRangesFromDataFrame(PSI_Annot[toKeep,],seqnames.field='chrom',start.field='start',end.field='end',keep.extra.columns=TRUE)
Cisdist=1e6
GR_cis=flank(flank(GRange_PSI,Cisdist),-2*Cisdist-1)
ooCis=findOverlaps(GR_cis, GRange_Map)
wCis_1Mb=which(MAP_GWAS$SNPfreq>=0.05)[unique(subjectHits(ooCis))]
samp_all=sample(1:nrow(MAP_GWAS),100000)
samp_cis=sample(wCis_1Mb,100000)
wCis=wCis_1Mb



#RESobs_nodup_1Mb$GWAS_EFO_R2_1E5_corrected=MAP_GWAS$GWAS_EFO_R2_1E5[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_EFO_R2_1E5=MAP_GWAS$GWAS_EFO_R2_1E5[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_EFO_R2_1E8=MAP_GWAS$GWAS_EFO_R2_1E8[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_Trait_R2_1E5=MAP_GWAS$GWAS_Trait_R2_1E5[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_Trait_R2_1E8=MAP_GWAS$GWAS_Trait_R2_1E8[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_EFO_Trait_R2_1E5=MAP_GWAS$GWAS_EFO_Trait_R2_1E5[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_EFO_Trait_R2_1E8=MAP_GWAS$GWAS_EFO_Trait_R2_1E8[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
save(RESobs_nodup_1Mb,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/cis-psiQTL_MISO_counts_V7_nodup_withSelect_RBP_coding_intron_aSNP_rsQTL_spliceElt_GWAS.Rdata',HOME))



GWAS_sQTL=RESobs_nodup_1Mb[grepl('asthma|allergy|celiac|type 1 diabetes|Graves|inflammatory bowel disease|multiple sclerosis|psoriasis|rheumatoid arthritis|systemic lupus erythematosus|ulcerative colitis|crohn',tolower(RESobs_nodup_1Mb$GWAS_Trait_R2_1E5)),c('is_rsQTL','symbol','haplo','pvalue','event_type','GWAS_Trait_R2_1E5','GWAS_Trait_R2_1E8','snps','Pvalue_NS','Pvalue_LPS','Pvalue_PAM3CSK4','Pvalue_R848','Pvalue_IAV','Beta_NS','Beta_LPS','Beta_PAM3CSK4','Beta_R848','Beta_IAV','Pvalue_rsQTL_LPS','Pvalue_rsQTL_PAM3CSK4','Pvalue_rsQTL_R848','Pvalue_rsQTL_IAV','event_id',"FPKM_NS","FPKM_LPS","FPKM_PAM3CSK4","FPKM_R848","FPKM_IAV","log2FC_LPS","log2FC_PAM3CSK4","log2FC_R848","log2FC_IAV","max_lFC","is_sQTL_of_stimInducedGene")]
GWAS_sQTL$GWAS_Trait_R2_1E5[GWAS_sQTL$GWAS_Trait_R2_1E5=="Crohn's disease and psoriasis"]="psoriasis/crohn'sdisease"
AI_traits=unique(unlist(strsplit(GWAS_sQTL$GWAS_Trait_R2_1E5,'/')))
AI_traits=AI_traits[grepl('asthma|allergy|celiac|type 1 diabetes|Graves|inflammatory bowel disease|multiple sclerosis|psoriasis|rheumatoid arthritis|systemic lupus erythematosus|ulcerative colitis|crohn',tolower(AI_traits))]
AI_traits_short=c('UC','IBD','MS','T1D','CD','Ath','PS','PS','CEL','SLE','PS','CD','IBD','All','Ath')
names(AI_traits_short)=AI_traits

GWAS_sQTL_trait=strsplit(GWAS_sQTL$GWAS_Trait_R2_1E5,'/')
GWAS_sQTL=GWAS_sQTL[rep(1:length(GWAS_sQTL_trait),sapply(GWAS_sQTL_trait,length)),]
GWAS_sQTL$GWAS_Trait_R2_1E5=unlist(GWAS_sQTL_trait)
GWAS_sQTL$isSignif_8=sapply(1:nrow(GWAS_sQTL),function(x){GWAS_sQTL$GWAS_Trait_R2_1E5[x]%in%strsplit(GWAS_sQTL$GWAS_Trait_R2_1E8[x],'/')[[1]]})
GWAS_sQTL$GWAS_Trait_R2_1E5=AI_traits_short[GWAS_sQTL$GWAS_Trait_R2_1E5]
GWAS_sQTL=GWAS_sQTL[!is.na(GWAS_sQTL$GWAS_Trait_R2_1E5),]
GWAS_sQTL=GWAS_sQTL[order(GWAS_sQTL$pvalue),]
GWAS_sQTL$dup_event=ifelse(duplicated(paste(GWAS_sQTL$symbol,GWAS_sQTL$GWAS_Trait_R2_1E5)),1,0)
write.table(GWAS_sQTL,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/GWAS_sQTL_forGraph.txt',HOME),quote=F,sep='\t',row.names=F)




# get all EFO parents 
EFOcat=unique(unlist(strsplit(RESobs_nodup_1Mb$GWAS_EFO_R2_1E5_corrected,'//')))
EFOcat=EFOcat[EFOcat!='']
w=which(!duplicated(RESobs_nodup_1Mb$haplo)) # all 1271 independent sQTLs 

# defining required Binning (all sQTLs)
mm=match(RESobs_nodup_1Mb$snps[w],MAP_GWAS$snp.name)
mmSNPfreq=MAP_GWAS$SNPfreq[mm]
binSize=0.02
binSize=binSize*100
maf=seq(0,0.5,0.005)
#cbind(maf,floor(maf*100),floor(maf*100/binSize)*binSize)
mmSNPfreq_bin=floor(mmSNPfreq*100/binSize)*binSize 
SNPfreq_bin=floor(MAP_GWAS$SNPfreq[wCis]*100/binSize)*binSize 
binCounts=table(mmSNPfreq_bin) # Nb of SNP to sample for each MAF bin 

NSAMP=1000
Resamp=list()
Resamp[['ALL']]=list()

# defined SNP_support (set of sampled SNPs)
#SNP_support=wCis
isEFO8=list()
isEFO5=list()
for (EFO in EFOcat[-grep('Other|Biological',EFOcat)]){
	isEFO8[[EFO]]=grepl(EFO,MAP_GWAS$GWAS_EFO_R2_1E8)
	isEFO5[[EFO]]=grepl(EFO,MAP_GWAS$GWAS_EFO_R2_1E5)
	Resamp[['ALL']][[paste(EFO,'1e8',sep='_')]]=c()
	Resamp[['ALL']][[paste(EFO,'1e5',sep='_')]]=c()
}

AI_inHouse='celiac|lupus|asthma|psoriasis|inflammatory bowel|Graves|crohn|ulcerative colitis|multiple sclerosis|type 1 diabetes|allergy|rheumatoid arthritis'
AI_inHouse='immune'
isGWAS_all8=MAP_GWAS$GWAS_Trait_R2_1E8!=''
isGWAS_all5=MAP_GWAS$GWAS_Trait_R2_1E5!=''
isGWAS_AI8=grepl(AI_inHouse,tolower(MAP_GWAS$GWAS_Trait_R2_1E8_corrected))
isGWAS_AI5=grepl(AI_inHouse,tolower(MAP_GWAS$GWAS_Trait_R2_1E5_corrected))
isGWAS_other8=isGWAS_all8 & ! isGWAS_AI8
isGWAS_other5=isGWAS_all5 & ! isGWAS_AI5

Resamp[['ALL']][[paste('automImmune_inHouse','1e5',sep='_')]]=c()
Resamp[['ALL']][[paste('automImmune_inHouse','1e8',sep='_')]]=c()
Resamp[['ALL']][[paste('other_inHouse','1e5',sep='_')]]=c()
Resamp[['ALL']][[paste('other_inHouse','1e8',sep='_')]]=c()
Resamp[['ALL']][[paste('GWAS','1e5',sep='_')]]=c()
Resamp[['ALL']][[paste('GWAS','1e8',sep='_')]]=c()


################# resampling EFOs, matcheing on global MAF 
for (samp in 1:NSAMP){
	tic=Sys.time()
	cat(samp)
	sample_global=By(wCis,SNPfreq_bin,sample,max(binCounts)) # over sample each bin
	for (i in 1:length(sample_global)){
		sample_global[[i]]=sample_global[[i]][1:binCounts[i]] # sub sample to right amount
	}
	sample_global=unlist(sample_global)
	for (EFO in EFOcat[-grep('Other|Biological',EFOcat)]){
		Resamp[['ALL']][[paste(EFO,'1e8',sep='_')]][samp]=sum(isEFO8[[EFO]][sample_global],na.rm=T)
		Resamp[['ALL']][[paste(EFO,'1e5',sep='_')]][samp]=sum(isEFO5[[EFO]][sample_global],na.rm=T)
		}
	Resamp[['ALL']][[paste('automImmune_inHouse','1e5',sep='_')]][samp]=sum(isEFO5[[EFO]][sample_global],na.rm=T)
	Resamp[['ALL']][[paste('automImmune_inHouse','1e8',sep='_')]][samp]=sum(isEFO5[[EFO]][sample_global],na.rm=T)
	Resamp[['ALL']][[paste('other_inHouse','1e5',sep='_')]][samp]=sum(isEFO5[[EFO]][sample_global],na.rm=T)
	Resamp[['ALL']][[paste('other_inHouse','1e8',sep='_')]][samp]=sum(isEFO5[[EFO]][sample_global],na.rm=T)
	Resamp[['ALL']][[paste('GWAS','1e5',sep='_')]][samp]=sum(isEFO5[[EFO]][sample_global],na.rm=T)
	Resamp[['ALL']][[paste('GWAS','1e8',sep='_')]][samp]=sum(isEFO5[[EFO]][sample_global],na.rm=T)
		
	toc=Sys.time()
	print(toc-tic)
	flush.console()
}
################# done resampling by EFO

save(Resamp,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/GWAS_resample_EFO_and_type_allsQTL_corrected.Rdata',HOME))




