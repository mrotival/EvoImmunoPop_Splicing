# load( PSI_Annot )
# load( sQTLs )

allPSI_events=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/Ens70_HISAT/PSI_events_ALL_V7.2_withCounts.Rdata',sep='')
load(allPSI_events)

library(data.table)
MAP_GWAS=as.data.frame(fread(sprintf('%s/Annotation/GWAS/Map_GWAS_LD_EFOcorrected_allChr.txt',HOME)))

EFO_14=sapply(strsplit(MAP_GWAS$GWAS_EFO_R2_1E5,'//'), function(x){paste(unique(x[!x%in%c('Other trait','Other measurement','Other disease','Biological process')]),collapse='//')})
sum(EFO_14!='') # 229425

library(GenomicRanges)
GRange_Map=makeGRangesFromDataFrame(MAP_GWAS[which(MAP_GWAS$SNPfreq>0.05 ),c('chromosome','position')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_PSI=makeGRangesFromDataFrame(PSI_Annot[toKeep,],seqnames.field='chrom',start.field='start',end.field='end',keep.extra.columns=TRUE)
Cisdist=1e6
GR_cis=reduce(union(flank(GRange_PSI,Cisdist,start=TRUE),union(flank(GRange_PSI,Cisdist,start=FALSE),GRange_PSI))) # combine, 1Mb upstream,downstream  and region of AS event.
ooCis=findOverlaps(GR_cis, GRange_Map)
wCis_1Mb=which(MAP_GWAS$SNPfreq>=0.05)[unique(subjectHits(ooCis))]
samp_all=sample(1:nrow(MAP_GWAS),100000)
samp_cis=sample(wCis_1Mb,100000)

RESobs_nodup_1Mb$GWAS_EFO_R2_1E5=MAP_GWAS$GWAS_EFO_R2_1E5[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_EFO_R2_1E8=MAP_GWAS$GWAS_EFO_R2_1E8[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_Trait_R2_1E5=MAP_GWAS$GWAS_Trait_R2_1E5[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_Trait_R2_1E8=MAP_GWAS$GWAS_Trait_R2_1E8[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]

#####################################################
##                      Table S3E                  ##
#####################################################

TableS3E = RESobs_nodup_1Mb[,c('event_id','symbol','event_type','snps','coding_type','is_sQTL_of_stimulationSpecificGene','is_rsQTL','GWAS_Trait_R2_1E5')]
write.table(TableS3E,file='table/TableS3E_sQTL_GWAS_overlap.txt',quote=F,sep='\t',row.names=F)

#####################################################
## GWAS Enrichment analysis (Figure 6A)            ##
#####################################################

# get all EFO parents 
EFOcat=unique(unlist(strsplit(RESobs_nodup_1Mb$GWAS_EFO_R2_1E5,'//')))
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
isGWAS_AI8=grepl(AI_inHouse,tolower(MAP_GWAS$GWAS_Trait_R2_1E8))
isGWAS_AI5=grepl(AI_inHouse,tolower(MAP_GWAS$GWAS_Trait_R2_1E5))
isGWAS_other8=isGWAS_all8 & ! isGWAS_AI8
isGWAS_other5=isGWAS_all5 & ! isGWAS_AI5

Resamp[['ALL']][[paste('automImmune_inHouse','1e5',sep='_')]]=c()
Resamp[['ALL']][[paste('automImmune_inHouse','1e8',sep='_')]]=c()
Resamp[['ALL']][[paste('other_inHouse','1e5',sep='_')]]=c()
Resamp[['ALL']][[paste('other_inHouse','1e8',sep='_')]]=c()
Resamp[['ALL']][[paste('GWAS','1e5',sep='_')]]=c()
Resamp[['ALL']][[paste('GWAS','1e8',sep='_')]]=c()

################# resampling EFOs, matching on global MAF 
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

save(Resamp,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/GWAS_resample_EFO_and_type_allsQTL.Rdata',HOME))

# get all EFO parents 
EFOcat=unique(unlist(strsplit(RESobs_nodup_1Mb$GWAS_EFO_R2_1E5,'//')))
EFOcat=EFOcat[EFOcat!='']

w=which(!duplicated(RESobs_nodup_1Mb$haplo)) # all 1271 independent sQTLs 
NbEFO5=NbEFO8=c()
P_EFO5=P_EFO8=c()
isEFO8=list()
isEFO5=list()
for (EFO in EFOcat[-grep('Other|Biological',EFOcat)]){
	isEFO8[[EFO]]=grepl(EFO,MAP_GWAS$GWAS_EFO_R2_1E8)
	isEFO5[[EFO]]=grepl(EFO,MAP_GWAS$GWAS_EFO_R2_1E5)
	NbEFO5[EFO]=sum(grepl(EFO,RESobs_nodup_1Mb$GWAS_EFO_R2_1E5[w]))
	NbEFO8[EFO]=sum(grepl(EFO,RESobs_nodup_1Mb$GWAS_EFO_R2_1E8[w]))
	P_EFO5[EFO]=mean(c(1,Resamp$ALL[[paste(EFO,'1e5',sep='_')]]>NbEFO5[EFO]))
	P_EFO8[EFO]=mean(c(1,Resamp$ALL[[paste(EFO,'1e8',sep='_')]]>NbEFO8[EFO]))
}

##############################################
#####	Fig 3A	Enrichment by EFO   	  ####
##############################################
# 11.25 x 5.5
# 1E5
d=0.3
perm=Resamp$ALL[grep('1e5',names(Resamp$ALL))]
obs=NbEFO5
perm=perm[order(-NbEFO5)]
obs=obs[order(-NbEFO5)]
par(mar=c(4,5,4,1))
nmax=max(max(obs),max(unlist(perm)))
permtab=lapply(perm,function(x){tab=table(x)[as.character(0:nmax)];tab[is.na(tab)]=0;names(tab)=as.character(0:nmax);tab})
library(RColorBrewer)
col13=c(rev(brewer.pal(7,'YlOrRd')),brewer.pal(7,'YlGnBu')[-1])
layout(rbind(rep(1,6)%o%1:14,15,15))
par(mar=0.5*c(6,6,3,1))
plot(0:nmax,permtab[[1]],ylim=c(-1,nmax),xlim=c(-50,450),bty='n',axes=F,col='#00000000',ylab='Number of sQTL overlapping GWAS loci')
axis(2,at=0:nmax,cex.axis=0.8,las=2)
par(mar=0.5*c(6,1,3,1))
for(i in 1:13){
	plot(as.numeric(permtab[[i]]),0:nmax,ylim=c(-1,nmax),xlim=c(-50,450),bty='n',axes=F,col='#00000000')
	axis(1,cex.axis=0.8,las=3)
#	sapply(0:nmax,function(j){segments(table(perm[[i]])[j+1],j,0,j,lwd=2)})
	sapply(0:nmax,function(j){if(permtab[[i]][j+1]>0){polygon(c(permtab[[i]][j+1],permtab[[i]][j+1],0,0),c(j+d,j-d,j-d,j+d),col='grey',border='#000000AA')}})
	points(0,obs[i],bty='n',bg=col13[i],col='black',pch=21,cex=1.5)
#	if(i<5){abline(v=550,lty=2,col='lightgrey')}	
	}
plot.new()
par(xpd=T)
legend(0,1,ncol=4,cex=1,pt.cex=1.2,pt.bg=col13,pch=21,bty='n',legend=names(obs))


P_EFO5
#                Body measurement Lipid or lipoprotein measurement        Digestive system disorder           Immune system disorder 
#                     0.000999001                      0.066933067                      0.000999001                      0.000999001 
#           Neurological disorder                           Cancer        Hematological measurement               Metabolic disorder 
#                     0.001998002                      0.000999001                      0.000999001                      0.001998002 
#      Cardiovascular measurement         Inflammatory measurement         Liver enzyme measurement           Cardiovascular disease 
#                     0.000999001                      0.004995005                      0.092907093                      0.033966034 
#                Response to drug 
#                     0.529470529 

P_EFO5[order(-NbEFO5)]


######################################################################
#####	Fig 3B	overlap between TMBIM1 sQTL and IBD GWAS hit   	  ####
######################################################################
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V8.1.Rdata',EVO_IMMUNO_POP))

SUMSTATS=fread(sprintf("%s/Annotation/GWAS/summaryStats_LDhub/iibdgc-trans-ancestry-filtered-summary-stats/EUR.IBD.gwas_info03_filtered.assoc",HOME))
SUMSTATS=SUMSTATS[match(intersect(res$all$eqtls$snps,SUMSTATS$SNP),SUMSTATS$SNP)]
SUMSTATS$pos=Genos$position[match(SUMSTATS$SNP,Genos$snp.name)]

GWAShits='rs2382817'
sQTL_peak='rs2271543'
mySymbol='TMBIM1'
chrom=2
pos=219142491
event_id="ENSG00000135926;SE:2:219146904-219150486:219150661-219157189:-"
cond=4

Genos=getGenos_locus(chrom,pos-1e6,pos+1e6)
Genos_mat=as.matrix(Genos[,grep('AFB|EUB',colnames(Genos))])
MAF=apply(Genos_mat,1,mean,na.rm=T)
rownames(Genos_mat)=Genos$snp.name

Genos_forLD=Genos_mat[match(SUMSTATS$SNP,rownames(Genos_mat)),match(indiv,colnames(Genos_mat))]
r2_EUB_with_GWAShits=cor(t(Genos_forLD[GWAShits,indiv[grep('EUB',indiv)],drop=F]),t(Genos_forLD[,indiv[grep('EUB',indiv)]]),use='p')^2
r2_EUB_with_sQTL_peak=cor(t(Genos_forLD[sQTL_peak,indiv[grep('EUB',indiv)],drop=F]),t(Genos_forLD[,indiv[grep('EUB',indiv)]]),use='p')^2

RES_TMBIM1=RESobs[RESobs$event_id==event_id & RESobs$cond==cond,]
SUMSTATS$sQTL_Pval=RES_TMBIM1$pvalue[match(SUMSTATS$SNP,RES_TMBIM1$snps)]

color_r2=c("#E31A1C","#FEE586","#A1DAB4","#225EA8")
mergeCols=function(col1,col2,pct2=0.5){colorRampPalette(c(col1,col2))(100)[round(pct2*100)]}
color_r2_dark=sapply(color_r2,mergeCols,'black')

# 5.6x 6.6 inches
layout(1:2)
par(mar=c(3,4,1,1))
getColor_R2=function(r2,color_R2=color_r2){ifelse(r2>0.8,color_R2[1],ifelse(r2>0.5,color_R2[2],ifelse(r2>0.2,color_R2[3],color_R2[4])))}

plot(SUMSTATS$pos/1e6, -log10(SUMSTATS$sQTL_Pval), pch=ifelse(SUMSTATS$SNP%in%sQTL_peak,23,21),
                                                    bg=getColor_R2(r2_EUB_with_sQTL_peak),
                                                    col=getColor_R2(r2_EUB_with_sQTL_peak,color_R2=color_r2_dark),
                                                    cex=0.5+r2_EUB_with_sQTL_peak^2,
                                                    main='', xlab='', ylab='sQTL p-value', las=1, axes=F)
axis(1);axis(2,las=1)

points(SUMSTATS$pos[SUMSTATS$SNP%in%sQTL_peak]/1e6, -log10(SUMSTATS$sQTL_Pval[SUMSTATS$SNP%in%sQTL_peak]) ,pch=23, bg='#FF7A7C', cex=1.5, lwd=2)
legend('topright',pch=c(23,21,21,21,21),pt.bg=c('#FF7A7C',color_r2),col=c('black',color_r2_dark),pt.lwd=c(2,1,1,1,1),pt.cex=c(1.5,1.3,0.85,0.6,0.5),legend=c('               ','','','',''),bty='n')

plot(SUMSTATS$pos/1e6, -log10(SUMSTATS$P), pch=21, bg=getColor_R2(r2_EUB_with_sQTL_peak),
                                                    col=getColor_R2(r2_EUB_with_sQTL_peak,
                                                    color_R2=color_r2_dark),
                                                    cex=0.5+r2_EUB_with_sQTL_peak^2, main='', xlab='', ylab='IBD p-value', axes=F)
axis(1);axis(2,las=1)
points(SUMSTATS$pos[SUMSTATS$SNP%in%sQTL_peak]/1e6,-log10(SUMSTATS$P[SUMSTATS$SNP%in%sQTL_peak]),pch=23,bg='#FF7A7C',cex=1.5,lwd=2,col='black')

######################################################################
#####	                        Fig 3C	   	                      ####
######################################################################

#####################################################
## sQTL overlapping autoimmune GWAS                ##
#####################################################
# select the right subset of columns to use (mainly symbol)

AutoImmune_sQTL=grepl('asthma|allergy|celiac|type 1 diabetes|Graves|inflammatory bowel disease|multiple sclerosis|psoriasis|rheumatoid arthritis|systemic lupus erythematosus|ulcerative colitis|crohn',tolower(RESobs_nodup_1Mb$GWAS_Trait_R2_1E5))

GWAS_sQTL=RESobs_nodup_1Mb[AutoImmune_sQTL,c('event_id','symbol','event_type','haplo','snps','pvalue','is_rsQTL','max_logFC','is_sQTL_of_stimulationSpecificGene','GWAS_Trait_R2_1E5','GWAS_Trait_R2_1E8')]
GWAS_sQTL$GWAS_Trait_R2_1E5[GWAS_sQTL$GWAS_Trait_R2_1E5=="Crohn's disease and psoriasis"]="psoriasis/crohn'sdisease"

# correspondance with short GWAS Names
AI_traits=unique(unlist(strsplit(GWAS_sQTL$GWAS_Trait_R2_1E5,'/')))
AI_traits=AI_traits[grepl('asthma|allergy|celiac|type 1 diabetes|Graves|inflammatory bowel disease|multiple sclerosis|psoriasis|rheumatoid arthritis|systemic lupus erythematosus|ulcerative colitis|crohn',tolower(AI_traits))]
AI_traits_short=c('UC','IBD','MS','T1D','CD','Ath','PS','PS','CEL','SLE','PS','CD','IBD','All','Ath')
names(AI_traits_short)=AI_traits

# split traits to have one trait/sQTL pair per row
GWAS_sQTL_trait=strsplit(GWAS_sQTL$GWAS_Trait_R2_1E5,'/')
GWAS_sQTL=GWAS_sQTL[rep(1:length(GWAS_sQTL_trait),sapply(GWAS_sQTL_trait,length)),]
GWAS_sQTL$GWAS_Trait_R2_1E5=unlist(GWAS_sQTL_trait)

GWAS_sQTL$isSignif_8=sapply(1:nrow(GWAS_sQTL),function(x){GWAS_sQTL$GWAS_Trait_R2_1E5[x]%in%strsplit(GWAS_sQTL$GWAS_Trait_R2_1E8[x],'/')[[1]]}) # is the GWAS hit significant at 1E-8 ?
GWAS_sQTL$GWAS_Trait_R2_1E5=AI_traits_short[GWAS_sQTL$GWAS_Trait_R2_1E5] # get short Names

# remove NA and order by sQTL Pvalue
GWAS_sQTL=GWAS_sQTL[!is.na(GWAS_sQTL$GWAS_Trait_R2_1E5),] 
GWAS_sQTL=GWAS_sQTL[order(GWAS_sQTL$pvalue),]
write.table(GWAS_sQTL,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/GWAS_sQTL_forGraph.txt',HOME),quote=F,sep='\t',row.names=F)

####### plot in cytoscape of related software.
