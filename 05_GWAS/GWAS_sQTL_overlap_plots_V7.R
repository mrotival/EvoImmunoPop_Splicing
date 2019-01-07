load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/GWAS_resample_EFO_and_type_allsQTL_corrected.Rdata',HOME)

Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME))
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/cis-psiQTL_MISO_counts_V7_nodup_withSelect_RBP_coding_introns.Rdata',HOME))
library(data.table)

#EFO_14=sapply(strsplit(MAP_GWAS$GWAS_EFO_R2_1E5,'//'), function(x){paste(unique(x[!x%in%c('Other trait','Other measurement','Other disease','Biological process')]),collapse='//')})
#sum(EFO_14!='') # 229425

MAP_GWAS=as.data.frame(fread(sprintf('%s/Annotation/GWAS/Map_GWAS_LD_EFOcorrected_allChr.txt',HOME)))
RESobs_nodup_1Mb$GWAS_EFO_R2_1E5_corrected=MAP_GWAS$GWAS_EFO_R2_1E5[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]
RESobs_nodup_1Mb$GWAS_EFO_R2_1E8_corrected=MAP_GWAS$GWAS_EFO_R2_1E8[match(RESobs_nodup_1Mb$snps,MAP_GWAS$snp.name)]

# get all EFO parents 
EFOcat=unique(unlist(strsplit(RESobs_nodup_1Mb$GWAS_EFO_R2_1E5_corrected,'//')))
EFOcat=EFOcat[EFOcat!='']

w=which(!duplicated(RESobs_nodup_1Mb$haplo)) # all 1271 independent sQTLs 
NbEFO5=NbEFO8=c()
P_EFO5=P_EFO8=c()
P_EFO5_param=P_EFO8_param=c()
isEFO8=list()
isEFO5=list()
for (EFO in EFOcat[-grep('Other|Biological',EFOcat)]){
	isEFO8[[EFO]]=grepl(EFO,MAP_GWAS$GWAS_EFO_R2_1E8)
	isEFO5[[EFO]]=grepl(EFO,MAP_GWAS$GWAS_EFO_R2_1E5)
	NbEFO5[EFO]=sum(grepl(EFO,RESobs_nodup_1Mb$GWAS_EFO_R2_1E5_corrected[w]))
	NbEFO8[EFO]=sum(grepl(EFO,RESobs_nodup_1Mb$GWAS_EFO_R2_1E8_corrected[w]))
	P_EFO5[EFO]=mean(c(1,Resamp$ALL[[paste(EFO,'1e5',sep='_')]]>NbEFO5[EFO]))
	P_EFO8[EFO]=mean(c(1,Resamp$ALL[[paste(EFO,'1e8',sep='_')]]>NbEFO8[EFO]))
#	P_EFO5_param[EFO]=pnorm(NbEFO5[EFO],mean(Resamp$ALL[[paste(EFO,'1e5',sep='_')]]),sd(Resamp$ALL[[paste(EFO,'1e5',sep='_')]]),low=F)
#	P_EFO8_param[EFO]=pnorm(NbEFO8[EFO],mean(Resamp$ALL[[paste(EFO,'1e8',sep='_')]]),sd(Resamp$ALL[[paste(EFO,'1e8',sep='_')]]),low=F)
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


P_EFO5_param
#                Body measurement Lipid or lipoprotein measurement        Digestive system disorder           Immune system disorder            Neurological disorder                           Cancer        Hematological measurement 
#                    6.120015e-08                     1.564375e-01                     3.479475e-08                     6.740916e-10                     6.117673e-04                     2.544470e-04                     4.738413e-03 
#              Metabolic disorder       Cardiovascular measurement         Inflammatory measurement         Liver enzyme measurement           Cardiovascular disease                 Response to drug 
#                    1.526148e-04                     6.838756e-04                     6.523972e-04                     2.464003e-01                     5.264418e-02                     6.833769e-01 
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
P_EFO5_param[order(-NbEFO5)]

##############################################
#####	Table: Odds ratio EFO 		  	  ####
##############################################

ORs_EFO5=do.call(rbind,lapply(1:length(isEFO5),function(i){
	tab=cbind(table(isEFO5[[i]][wCis]),table(grepl(names(isEFO5)[i],RESobs_nodup_1Mb$GWAS_EFO_R2_1E5_corrected[!duplicated(RESobs_nodup_1Mb$haplo)])))
	myOR=odds.ratio(tab)
	myOR$EFO=names(isEFO5)[i]
	myOR$GWAS_threshold='1e-5'
	myOR
	}))
	
ORs_EFO8=do.call(rbind,lapply(1:length(isEFO8),function(i){
	tab=cbind(table(isEFO8[[i]][wCis]),table(grepl(names(isEFO8)[i],RESobs_nodup_1Mb$GWAS_EFO_R2_1E8_corrected[!duplicated(RESobs_nodup_1Mb$haplo)])))
	myOR=odds.ratio(tab)
	myOR$EFO=names(isEFO8)[i]
	myOR$GWAS_threshold='1e-8'
	myOR
	}))

write.table(rbind(ORs_EFO5,ORs_EFO8),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/OR_GWAS_EFO_corrected.txt',HOME),quote=F,sep='\t',row.names=F)



##################################################
##################################################
##      show colocalization sQTL & GWAS         ##
##################################################
##################################################


Genos_sQTL=read.table(file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/SNP_genos_sQTL_bestSNPperCondition.txt',sep=''),sep='\t',header=T)
snps_sQTL=as.matrix(Genos_sQTL[-1])
rownames(snps_sQTL)=gsub('.',':',Genos_sQTL[[1]],fixed=T)
snps_sQTL=t(snps_sQTL[match(RESobs_nodup_1Mb$snps,rownames(snps_sQTL)),])

psi_sQTL=mapply(function(event,cond){PSI_prov[match(event,PSI_Annot[toKeep,'event_id']),match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(PSI_prov))]},RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$cond)
rownames(psi_sQTL)=rownames(snps_sQTL)
gene_sQTL=mapply(function(gene,cond){FPKM_gene[gene,match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(FPKM_gene))]},RESobs_nodup_1Mb$gene,RESobs_nodup_1Mb$cond)

#Celiac
#Crohns_Disease
#Fasting_Glucose
#Height1
#IBD
#Lupus
#Multiple_sclerosis
#Primary_biliary_cirrhosis
#Rheumatoid_Arthritis
#Type_1_Diabetes
#Ulcerative_Colitis
mySymbol='MVP'
mySymbol='TMBIM1'

j=which(RESobs_nodup_1Mb$symbol==mySymbol)

RESobs_nodup_1Mb[j,c('event_id','GWAS_Trait_R2_1E5','event_type','symbol','pos','locus')]

myTrait='IBD'
chrom=2
pos=RESobs_nodup_1Mb[j,'pos']
event_id=RESobs_nodup_1Mb[j,'event_id']
cond=RESobs_nodup_1Mb[j,'condition']

psi=psi_sQTL[,j]
Genos=getGenos_locus(chrom,pos-1e6,pos+1e6)
Genos_mat=2-as.matrix(Genos[,grep('AFB|EUB',colnames(Genos))])
MAF=apply(Genos_mat,1,mean,na.rm=T)
rownames(Genos_mat)=Genos$snp.name
PCs=read.table(sprintf('%s/02_ERC_Data/SampleCaseReports/EvoImmunoPop_ALL_ancestryPCA.txt',HOME),header=T)
rownames(PCs)=PCs$IID

indiv=SampleAnnot$individu[SampleAnnot$cond==cond & SampleAnnot$pop=='EUB']
indiv=SampleAnnot$individu[SampleAnnot$cond==cond]
Geno_Mat = SlicedData$new()
wFreq=which(MAF>0.05)
Geno_Mat$CreateFromMatrix(Genos_mat[wFreq,match(indiv,colnames(Genos_mat))])

PSI_Mat = SlicedData$new()
psi=psi_sQTL[match(indiv,substr(rownames(psi_sQTL),1,6)),j]
PSI=matrix(qnorm(rank(psi,ties='random')/(length(psi)+1),mean(psi,na.rm=T),sd(psi,na.rm=T)),nrow=1)
PSI_Mat$CreateFromMatrix(PSI)

Cvrt=SlicedData$new()
Cvrt$CreateFromMatrix(t(PCs[match(indiv,PCs$IID),3:4]))
res=Matrix_eQTL_main(Geno_Mat,PSI_Mat,Cvrt,output_file_name=NULL, pvOutputThreshold=1)

SUMSTATS=fread(sprintf("%s/Annotation/GWAS/summaryStats_LDhub/iibdgc-trans-ancestry-filtered-summary-stats/EUR.IBD.gwas_info03_filtered.assoc",HOME))
SUMSTATS=SUMSTATS[match(intersect(res$all$eqtls$snps,SUMSTATS$SNP),SUMSTATS$SNP)]
SUMSTATS$pos=Genos$position[match(SUMSTATS$SNP,Genos$snp.name)]
SUMSTATS$sQTL_Pval=res$all$eqtls$pvalue[match(SUMSTATS$SNP,res$all$eqtls$snps)]
SUMSTATS$sQTL_beta=res$all$eqtls$beta[match(SUMSTATS$SNP,res$all$eqtls$snps)]
SUMSTATS$sQTL_t=res$all$eqtls$statistic[match(SUMSTATS$SNP,res$all$eqtls$snps)]


GWAShits=''
GWAShits='rs2382817'
sQTL_peak='rs2271543'

Genos_forLD=Genos_mat[match(SUMSTATS$SNP,rownames(Genos_mat)),match(indiv,colnames(Genos_mat))]
r2_EUB_with_GWAShits=cor(t(Genos_forLD[GWAShits,indiv[grep('EUB',indiv)],drop=F]),t(Genos_forLD[,indiv[grep('EUB',indiv)]]),use='p')^2
r2_EUB_with_sQTL_peak=cor(t(Genos_forLD[sQTL_peak,indiv[grep('EUB',indiv)],drop=F]),t(Genos_forLD[,indiv[grep('EUB',indiv)]]),use='p')^2

color_r2=c("#E31A1C","#FEE586","#A1DAB4","#225EA8")
mergeCols=function(col1,col2,pct2=0.5){colorRampPalette(c(col1,col2))(100)[round(pct2*100)]}
color_r2_dark=sapply(color_r2,mergeCols,'black')

# 5.6x 6.6 inches
layout(1:2)
par(mar=c(3,4,1,1))
getColor_R2=function(r2,color_R2=color_r2){ifelse(r2>0.8,color_R2[1],ifelse(r2>0.5,color_R2[2],ifelse(r2>0.2,color_R2[3],color_R2[4])))}

plot(SUMSTATS$pos/1e6,-log10(SUMSTATS$sQTL_Pval),pch=ifelse(SUMSTATS$SNP%in%sQTL_peak,23,21),bg=getColor_R2(r2_EUB_with_sQTL_peak),col=getColor_R2(r2_EUB_with_sQTL_peak,color_R2=color_r2_dark),cex=0.5+r2_EUB_with_sQTL_peak^2,main='',xlab='',ylab='sQTL p-value',las=1,axes=F)
axis(1);axis(2,las=1)
points(SUMSTATS$pos[SUMSTATS$SNP%in%sQTL_peak]/1e6,-log10(SUMSTATS$sQTL_Pval[SUMSTATS$SNP%in%sQTL_peak]),pch=23,bg='#FF7A7C',cex=1.5,lwd=2)
#legend('topright',pch=c(23,21,21,21,21),pt.bg=c('#FF7A7C',color_r2),col=c('black',color_r2_dark),pt.lwd=c(2,1,1,1,1),pt.cex=c(1.5,1.3,0.85,0.6,0.5),legend=c('               ','','','',''),bty='n')
plot(SUMSTATS$pos/1e6,-log10(SUMSTATS$P),pch=21,bg=getColor_R2(r2_EUB_with_sQTL_peak),col=getColor_R2(r2_EUB_with_sQTL_peak,color_R2=color_r2_dark),cex=0.5+r2_EUB_with_sQTL_peak^2,main='',xlab='',ylab='IBD p-value',axes=F)
axis(1);axis(2,las=1)
points(SUMSTATS$pos[SUMSTATS$SNP%in%sQTL_peak]/1e6,-log10(SUMSTATS$P[SUMSTATS$SNP%in%sQTL_peak]),pch=23,bg='#FF7A7C',cex=1.5,lwd=2,col='black')



SUMSTATS=fread(sprintf('%s/Annotation/GWAS/summaryStats_LDhub/PASS_%s.sumstats',HOME,myTrait))

SUMSTATS=SUMSTATS[match(intersect(res$all$eqtls$snps,SUMSTATS$SNP),SUMSTATS$SNP)]
SUMSTATS$pos=Genos$position[match(SUMSTATS$SNP,Genos$snp.name)]
SUMSTATS$sQTL_Pval=res$all$eqtls$pvalue[match(SUMSTATS$SNP,res$all$eqtls$snps)]
SUMSTATS$sQTL_beta=res$all$eqtls$beta[match(SUMSTATS$SNP,res$all$eqtls$snps)]
SUMSTATS$sQTL_t=res$all$eqtls$statistic[match(SUMSTATS$SNP,res$all$eqtls$snps)]


#ggplot(SUMSTATS,aes(x=pos,y=
plot(SUMSTATS$pos,-log10(SUMSTATS$sQTL_Pval),pch=16,col=colERC5[1])
points(SUMSTATS$pos,-log10(2*pnorm(abs(SUMSTATS$Z),low=F)),pch=16,col=colERC5[2])

plot(-log10(2*pnorm(abs(SUMSTATS$Z),low=F)),-log10(SUMSTATS$sQTL_Pval),pch=16,col=colERC5[1])

SUMSTATS=fread(sprintf("%s/Annotation/GWAS/summaryStats_LDhub/ucmeta-sumstats.txt",HOME))
SUMSTATS=fread(sprintf("%s/Annotation/GWAS/summaryStats_LDhub/ucmeta-sumstats.txt",HOME))
SUMSTATS=SUMSTATS[match(intersect(res$all$eqtls$snps,SUMSTATS$SNP),SUMSTATS$SNP)]
SUMSTATS$pos=Genos$position[match(SUMSTATS$SNP,Genos$snp.name)]
SUMSTATS$sQTL_Pval=res$all$eqtls$pvalue[match(SUMSTATS$SNP,res$all$eqtls$snps)]
SUMSTATS$sQTL_beta=res$all$eqtls$beta[match(SUMSTATS$SNP,res$all$eqtls$snps)]
SUMSTATS$sQTL_t=res$all$eqtls$statistic[match(SUMSTATS$SNP,res$all$eqtls$snps)]

plot(SUMSTATS$pos,-log10(SUMSTATS$sQTL_Pval),pch=16,col=colERC5[1],main='uc')
points(SUMSTATS$pos,-log10(SUMSTATS$'SCAN-P'),pch=16,col=colERC5[2])

SUMSTATS=fread(sprintf("%s/Annotation/GWAS/summaryStats_LDhub/cd-meta.txt",HOME))
SUMSTATS=SUMSTATS[match(intersect(res$all$eqtls$snps,SUMSTATS$SNP),SUMSTATS$SNP)]
SUMSTATS$pos=Genos$position[match(SUMSTATS$SNP,Genos$snp.name)]
SUMSTATS$sQTL_Pval=res$all$eqtls$pvalue[match(SUMSTATS$SNP,res$all$eqtls$snps)]
SUMSTATS$sQTL_beta=res$all$eqtls$beta[match(SUMSTATS$SNP,res$all$eqtls$snps)]
SUMSTATS$sQTL_t=res$all$eqtls$statistic[match(SUMSTATS$SNP,res$all$eqtls$snps)]





#GWAShits=c('rs2382817','rs11676348','rs13392177','rs11677953')

plot(SUMSTATS$pos,-log10(SUMSTATS$sQTL_Pval),pch=21,bg=colERC5[1],col=ifelse(SUMSTATS$SNP%in%GWAShits,'black','#00000000'),cex=0.5+r2_EUB_with_GWAShits,main='cd')
points(SUMSTATS$pos,-log10(SUMSTATS$P),pch=21,bg=colERC5[2],col=ifelse(SUMSTATS$SNP%in%GWAShits,'black','#00000000'),cex=0.5+r2_EUB_with_sQTL_peak)


GWAShits=''
GWAShits='rs2382817'
sQTL_peak='rs2271543'

Genos_forLD=Genos_mat[match(SUMSTATS$SNP,rownames(Genos_mat)),match(indiv,colnames(Genos_mat))]
r2_EUB_with_GWAShits=cor(t(Genos_forLD[GWAShits,,drop=F]),t(Genos_forLD),use='p')^2
r2_EUB_with_sQTL_peak=cor(t(Genos_forLD[sQTL_peak,,drop=F]),t(Genos_forLD),use='p')^2


color_r2=c("#E31A1C","#FEE586","#A1DAB4","#225EA8")
mergeCols=function(col1,col2,pct2=0.5){colorRampPalette(c(col1,col2))(100)[round(pct2*100)]}
color_r2_dark=sapply(color_r2,mergeCols,'black')

layout(1:2)
getColor_R2=function(r2){ifelse(r2>0.8,color_r2[1],ifelse(r2>0.5,color_r2[2],ifelse(r2>0.2,color_r2[3],color_r2[4])))}

plot(SUMSTATS$pos,-log10(SUMSTATS$sQTL_Pval),pch=ifelse(SUMSTATS$SNP%in%GWAShits,23,21),bg=getColor_R2(r2_EUB_with_GWAShits),cex=0.5+r2_EUB_with_GWAShits^2,main='TMBIM1 sQTL')
points(SUMSTATS$pos[SUMSTATS$SNP%in%GWAShits],-log10(SUMSTATS$sQTL_Pval[SUMSTATS$SNP%in%GWAShits]),pch=23,bg='deeppink',cex=1.8)

plot(SUMSTATS$pos,-log10(SUMSTATS$P),pch=21,bg=getColor_R2(r2_EUB_with_sQTL_peak),cex=0.5+r2_EUB_with_sQTL_peak^2,main='TMBIM1 sQTL')
points(SUMSTATS$pos[SUMSTATS$SNP%in%sQTL_peak],-log10(SUMSTATS$P[SUMSTATS$SNP%in%sQTL_peak]),pch=23,bg='deeppink',cex=1.8)

# 5.6x 6.6 inches
layout(1:2)
par(mar=c(3,4,1,1))
getColor_R2=function(r2,color_R2=color_r2){ifelse(r2>0.8,color_R2[1],ifelse(r2>0.5,color_R2[2],ifelse(r2>0.2,color_R2[3],color_R2[4])))}

plot(SUMSTATS$pos/1e6,-log10(SUMSTATS$sQTL_Pval),pch=ifelse(SUMSTATS$SNP%in%sQTL_peak,23,21),bg=getColor_R2(r2_EUB_with_sQTL_peak),col=getColor_R2(r2_EUB_with_sQTL_peak,color_R2=color_r2_dark),cex=0.5+r2_EUB_with_sQTL_peak^2,main='',xlab='',ylab='sQTL p-value (EUB)',las=1,axes=F)
axis(1);axis(2,las=1)
points(SUMSTATS$pos[SUMSTATS$SNP%in%sQTL_peak]/1e6,-log10(SUMSTATS$sQTL_Pval[SUMSTATS$SNP%in%sQTL_peak]),pch=23,bg='#FF7A7C',cex=1.5,lwd=2)
#legend('topright',pch=c(23,21,21,21,21),pt.bg=c('#FF7A7C',color_r2),col=c('black',color_r2_dark),pt.lwd=c(2,1,1,1,1),pt.cex=c(1.5,1.3,0.85,0.6,0.5),legend=c('               ','','','',''),bty='n')
plot(SUMSTATS$pos/1e6,-log10(SUMSTATS$P),pch=21,bg=getColor_R2(r2_EUB_with_sQTL_peak),col=getColor_R2(r2_EUB_with_sQTL_peak,color_R2=color_r2_dark),cex=0.5+r2_EUB_with_sQTL_peak^2,main='',xlab='',ylab='IBD p-value',axes=F)
axis(1);axis(2,las=1)
points(SUMSTATS$pos[SUMSTATS$SNP%in%sQTL_peak]/1e6,-log10(SUMSTATS$P[SUMSTATS$SNP%in%sQTL_peak]),pch=23,bg='#FF7A7C',cex=1.5,lwd=2,col='black')


plot(SUMSTATS$pos,-log10(SUMSTATS$P),pch=21,bg=colERC5[2],col=ifelse(SUMSTATS$SNP%in%GWAShits,'black','#00000000'),cex=0.5+r2_EUB_with_sQTL_peak^2,main='IBD GWAS')

SUMSTATS_V2=fread(sprintf("%s/Annotation/GWAS/summaryStats_LDhub/iibdgc-trans-ancestry-filtered-summary-stats/IBD_trans_ethnic_association_summ_stats_b37.txt",HOME))
SUMSTATS_V2=SUMSTATS_V2[match(intersect(res$all$eqtls$snps,SUMSTATS_V2$SNP),SUMSTATS_V2$SNP)]
SUMSTATS_V2$pos=Genos$position[match(SUMSTATS_V2$SNP,Genos$snp.name)]
SUMSTATS_V2$sQTL_Pval=res$all$eqtls$pvalue[match(SUMSTATS_V2$SNP,res$all$eqtls$snps)]
SUMSTATS_V2$sQTL_beta=res$all$eqtls$beta[match(SUMSTATS_V2$SNP,res$all$eqtls$snps)]
SUMSTATS_V2$sQTL_t=res$all$eqtls$statistic[match(SUMSTATS_V2$SNP,res$all$eqtls$snps)]

GWAShits=''
GWAShits='rs2382817'
GWAShits=c('rs2382817','rs11676348','rs13392177','rs11677953')
plot(SUMSTATS_V2$pos,-log10(SUMSTATS_V2$sQTL_Pval),pch=21,bg=colERC5[1],col=ifelse(SUMSTATS_V2$SNP%in%GWAShits,'black','#00000000'),cex=ifelse(SUMSTATS_V2$SNP%in%GWAShits,2,1),main='cd')
points(SUMSTATS_V2$pos,SUMSTATS_V2$Mantra_log10BF,pch=21,bg=colERC5[2],col=ifelse(SUMSTATS_V2$SNP%in%GWAShits,'black','#00000000'),cex=ifelse(SUMSTATS_V2$SNP%in%GWAShits,2,1))


load(file=sprintf('%s/Annotation/GWAS/gwas_catalog_v1_downloaded05072017.RData',HOME)) # this WAS lifted over from GRCh38 to GRCh37
library(data.table)

GWAScat[GWAScat$seqnames==2 & GWAScat$start>219e6 & GWAScat$start<220e6 & grepl("(ulcerative colitis|Crohn's disease|inflammatory bowel disease)", GWAScat$MAPPED_TRAIT),]






##################################################
##################################################
##            OR GWAS by traits                 ##
##################################################
##################################################

TraitsToTest=table(unlist(strsplit(MAP_GWAS$GWAS_Trait_1E5[MAP_GWAS$SNPfreq>0.05],'//')))
TraitsToTest=TraitsToTest[TraitsToTest>=20]

ORs_Trait5=do.call(rbind,lapply(1:length(TraitsToTest),function(i){
	if(any(grepl(names(TraitsToTest)[i],RESobs_nodup_1Mb$GWAS_Trait_R2_1E5[!duplicated(RESobs_nodup_1Mb$haplo)]))){
		tab=cbind(table(grepl(names(TraitsToTest)[i],MAP_GWAS$GWAS_Trait_R2_1E5[wCis])),table(grepl(names(TraitsToTest)[i],RESobs_nodup_1Mb$GWAS_Trait_R2_1E5[!duplicated(RESobs_nodup_1Mb$haplo)])))
		myOR=odds.ratio(tab)
		myOR$EFO=names(TraitsToTest)[i]
		myOR$GWAS_threshold='1e-5'	
		myOR
	}else{
		NULL
	}
	}))
	
ORs_Trait8=do.call(rbind,lapply(1:length(TraitsToTest),function(i){
	if(any(grepl(names(TraitsToTest)[i],RESobs_nodup_1Mb$GWAS_Trait_R2_1E8[!duplicated(RESobs_nodup_1Mb$haplo)]))){
	tab=cbind(table(grepl(names(TraitsToTest)[i],MAP_GWAS$GWAS_Trait_R2_1E8[wCis])),table(grepl(names(TraitsToTest)[i],RESobs_nodup_1Mb$GWAS_Trait_R2_1E8[!duplicated(RESobs_nodup_1Mb$haplo)])))
	myOR=odds.ratio(tab)
	myOR$EFO=names(TraitsToTest)[i]
	myOR$GWAS_threshold='1e-8'
	myOR
	}else{
		NULL
	}
	}))

write.table(rbind(ORs_Trait5,ORs_Trait8),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/OR_GWAS_Traits.txt',HOME),quote=F,sep='\t',row.names=F)
ORs_Trait8=ORs_Trait8[order(-ORs_Trait8$LowerCI),]
ORs_Trait8[ORs_Trait8$P<0.05,]

ORs_Trait5=ORs_Trait5[order(-ORs_Trait5$LowerCI),]
ORs_Trait5[ORs_Trait5$P<0.05,]
#              LowerCI        OR    UpperCI alpha            P                                  EFO GWAS_threshold
#odds ratio21 3.9055262 10.674130  23.313305  0.05 2.772945e-05         Chronic lymphocytic leukemia           1e-5
#odds ratio79 3.6055464 11.135845  26.097630  0.05 1.065575e-04              Systolic blood pressure           1e-5
#odds ratio58 2.2134920  6.046027  13.203028  0.05 5.801764e-04                   Multiple sclerosis           1e-5
#odds ratio43 2.1436202  4.700250   8.956148  0.05 1.776196e-04           Inflammatory bowel disease           1e-5
#odds ratio67 2.0792823 10.110691  29.671632  0.05 3.520540e-03            Primary biliary cirrhosis           1e-5
#odds ratio27 1.9828748  4.348089   8.283056  0.05 3.121251e-04                      Crohn's disease           1e-5
#odds ratio10 1.7734832  3.888613   7.408241  0.05 6.863202e-04                      Body mass index           1e-5
#odds ratio17 1.7641778  8.575634  25.165162  0.05 5.541470e-03 Cerebrospinal fluid biomarker levels           1e-5
#odds ratio84 1.5739962  3.924213   8.113993  0.05 2.483708e-03                   Ulcerative colitis           1e-5
#odds ratio37 1.5720298  3.157786   5.671424  0.05 1.010457e-03                               Height           1e-5
#odds ratio56 1.5595895  5.736302  14.734914  0.05 5.724030e-03                             Migraine           1e-5
#odds ratio32 1.3296003  6.460984  18.945619  0.05 1.189759e-02                    Fibrinogen levels           1e-5
#odds ratio31 1.2806025  6.222775  18.245285  0.05 1.314269e-02                           Fibrinogen           1e-5
#odds ratio65 1.2446172  6.047771  17.731144  0.05 1.416994e-02                       Platelet count           1e-5
#odds ratio38 1.2120205  3.740664   8.757207  0.05 1.200310e-02                    Hip circumference           1e-5
#odds ratio77 1.1624959  2.548693   4.854635  0.05 1.058300e-02                        Schizophrenia           1e-5
#odds ratio29 1.0363635  2.829944   6.179221  0.05 2.149883e-02               Educational attainment           1e-5
#odds ratio36 0.9844225  8.149484  29.563054  0.05 2.573283e-02                           Heart rate           1e-5
#odds ratio64 0.9826231  8.131024  29.508865  0.05 2.581955e-02        Pediatric autoimmune diseases           1e-5
#odds ratio88 0.9215114  3.388567   8.699627  0.05 3.228285e-02                  Waist circumference           1e-5
#odds ratio9  0.8634478  4.194821  12.297460  0.05 3.617556e-02                       Blood pressure           1e-5
#odds ratio13 0.8449784 33.650693 191.206590  0.05 2.951767e-02           Bronchopulmonary dysplasia           1e-5
#odds ratio40 0.8076706  6.682856  24.244224  0.05 3.690800e-02                         Hypertension           1e-5
#odds ratio51 0.6926466  5.730682  20.778314  0.05 4.858522e-02              Mean corpuscular volume           1e-5



#GWAS_sQTL=read.table('/Users/mrotival/WORK/Splicing/GWAS_sQTL_cyto_2.txt',header=T,sep='\t',stringsAsFactors=F)
Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))

library(igraph)
GWAS_sQTL=read.table(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/GWAS_sQTL_forGraph.txt',HOME),header=T,sep='\t',stringsAsFactors=F,quote='')
# GWAS_sQTL=read.table('/Users/mrotival/WORK/06_papers/PapierSplicing_20180126/V6.5/GWAS_sQTL_forGraph.txt',header=T,sep='\t',stringsAsFactors=F,quote='')
GWAS_sQTL$LPS_lFC=log2(GWAS_sQTL$LPS_Expr)-log2(GWAS_sQTL$NS_Expr)
GWAS_sQTL$PAM3CSK4_lFC=log2(GWAS_sQTL$PAM3CSK4_Expr)-log2(GWAS_sQTL$NS_Expr)
GWAS_sQTL$R848_lFC=log2(GWAS_sQTL$R848_Expr)-log2(GWAS_sQTL$NS_Expr)
GWAS_sQTL$IAV_lFC=log2(GWAS_sQTL$IAV_Expr)-log2(GWAS_sQTL$NS_Expr)
GWAS_sQTL$Max_lFC=pmax(GWAS_sQTL$LPS_lFC,GWAS_sQTL$PAM3CSK4_lFC,GWAS_sQTL$R848_lFC,GWAS_sQTL$IAV_lFC)

#GWAS_sQTL=GWAS_sQTL[GWAS_sQTL$dup_event==0,]
luq=function(x){length(unique(x))}

#for(i in grep('Testable_',colnames(PSI_Annot))){
#	GWAS_sQTL[[colnames(PSI_Annot)[i]]]=PSI_Annot[match(GWAS_sQTL$event_id,PSI_Annot$event_id),i]
#	}
#GWAS_sQTL$nsSpecific_event=GWAS_sQTL$Testable_NS & !apply(GWAS_sQTL[,grep('Testable_',cn(GWAS_sQTL))[-1]],1,any)
#GWAS_sQTL$stimSpecific_event= !GWAS_sQTL$Testable_NS & apply(GWAS_sQTL[,grep('Testable_',cn(GWAS_sQTL))[-1]],1,any)

colPSI=structure(c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854","#FFD92F", "#E5C494"), .Names = c("A3", "A5", "AF", "AL", "MX","RI", "SE"))

library(igraph)

nodes=data.frame(id=c(unique(GWAS_sQTL$GWAS_Trait_R2_1E5),unique(GWAS_sQTL$event_id),unique(GWAS_sQTL$haplo)),type=rep(c('GWAS','Splice','SNP'),c(luq(GWAS_sQTL$GWAS_Trait_R2_1E5),luq(GWAS_sQTL$event_id),luq(GWAS_sQTL$haplo))))
nodes$color=colPSI[GWAS_sQTL$event_type[match(nodes$id,GWAS_sQTL$event_id)]]
nodes$color[is.na(nodes$color)]='grey'
nodes$color[nodes$type=='SNP']='darkgrey'
nodes$size=ifelse(nodes$type=='GWAS',14,18)
nodes$size[nodes$type=='SNP']=5
nodes$vertex.label.color=ifelse(nodes$type=='GWAS',"black","black")
nodes$vertex.label.color[nodes$type=='SNP']='#00000000'
nodes$vertex.label.family='arial'
nodes$vertex.label.cex=ifelse(nodes$type=='GWAS',0.8,0.6)
nodes$vertex.frame.width= ifelse(GWAS_sQTL$Max_lFC[match(nodes$id,GWAS_sQTL$event_id)]>1, 2, 1)
nodes$vertex.frame.width[is.na(nodes$vertex.frame.width)]=0.01
nodes$vertex.frame.color=substr(colERC5[1],1,7)
#nodes$vertex.frame.color[which(GWAS_sQTL$is_stim_sQTL[match(nodes$id,GWAS_sQTL$symbol)]==1)]=substr(colERC5[4],1,7)
nodes$vertex.frame.color[which(GWAS_sQTL$Max_lFC[match(nodes$id,GWAS_sQTL$event_id)]>1)]='black'
nodes$vertex.frame.color[which(GWAS_sQTL$Max_lFC[match(nodes$id,GWAS_sQTL$event_id)]>1 & GWAS_sQTL$NS_Expr[match(nodes$id,GWAS_sQTL$event_id)]<10)]=colERC[4]
nodes$vertex.label=nodes$id
nodes$vertex.label[nodes$type=='SNP']=''
nodes$vertex.label[nodes$type=='Splice']=GWAS_sQTL$symbol[match(nodes$id[nodes$type=='Splice'],GWAS_sQTL$event_id)]

mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

add.vertex.shape("fcircle", clip=igraph.shape.noclip,
		plot=mycircle, parameters=list(vertex.frame.color=1,
                                  vertex.frame.width=1))

X=GWAS_sQTL[,c('event_id','haplo','event_type','isSignif_8','is_stim_sQTL','Pvalue_NS')]
X$isSignif_8=FALSE
X=X[!duplicated(X),-3]
colnames(X)[1]='id'
X$width=ifelse(X$Pvalue_NS>0.001,3,1)
X$color=colERC[8]
#X$color[which(X$is_stim_sQTL)]=substr(colERC[3],1,7)
X$color[which(X$Pvalue_NS>0.001 )]=substr(colERC[4],1,7)

Y=GWAS_sQTL[,c('GWAS_Trait_R2_1E5','haplo',"isSignif_8",'is_stim_sQTL','Pvalue_NS')]
Y[,4:5]=FALSE
Y=Y[!duplicated(Y),]
colnames(Y)[1]='id'
Y$width=ifelse(Y$isSignif_8,1,1)
Y$color=ifelse(Y$isSignif_8,grey(0.2),grey(0.6))

edges=rbind(X,Y)
#edges$width=ifelse(edges$isSignif_8,1,1)
#edges$color=ifelse(edges$isSignif_8,grey(0.2),grey(0.6))

net <- graph_from_data_frame(d=edges,vertices=nodes,directed=F)
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/GWAS_sQTL_plot_coords_withSNP_v2.Rdata',HOME))
#l3[,2]=(l3[,2]-mean(range(l3[,2])))/diff(range(l3[,2]))*2
#l3[,1]=(l3[,1]-mean(range(l3[,1])))/diff(range(l3[,1]))*2
l3[,1]=(l3[,1]-mean(range(l3[,1])))
l3[,2]=(l3[,2]-mean(range(l3[,2])))
par()$mar
# 13x9 inches
plot(net,layout=l3,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color=grey(0.2),
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label,rescale=F)

# black circle : genes is up regulated after stimulated
# red circle :gene is expressed only upon stimulation
# red: line sQTL is detected only upon stimulation (P<1e-3)


plot(net,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color='black',
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label)


load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/GWAS_sQTL_plot_coords_withSNP_v2.Rdata',HOME))
#l3[,2]=(l3[,2]-mean(range(l3[,2])))/diff(range(l3[,2]))*2
#l3[,1]=(l3[,1]-mean(range(l3[,1])))/diff(range(l3[,1]))*2
l3[,1]=(l3[,1]-mean(range(l3[,1])))
par()$mar
plot(net,layout=l3,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color=grey(0.2),
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label,rescale=F)


tkid <- tkplot(net) #tkid is the id of the tkplot that will open
l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
l2=l
l2[,2]=(l2[,2]-500)/600
l2[,1]=(l2[,1]-750)/600
#save(l2,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/GWAS_sQTL_plot_coords_withSNP.Rdata',HOME))

nodes$vertex.label.color=ifelse(nodes$type=='GWAS',"black","black")
nodes$vertex.label.family='arial'
nodes$vertex.label.cex=ifelse(nodes$type=='GWAS',0.8,0.6)

nodes$vertex.frame.width= ifelse(GWAS_sQTL$is_rsQTL[match(nodes$id,GWAS_sQTL$event_id)]==1, 3, 1)
nodes$vertex.frame.width[is.na(nodes$vertex.frame.width)]=0.01
nodes$vertex.frame.color=colERC[2]
nodes$vertex.frame.color[which(GWAS_sQTL$is_nsSpecific[match(nodes$id,GWAS_sQTL$event_id)]==1)]=substr(colERC[6],1,7)
nodes$vertex.frame.color[which(GWAS_sQTL$is_stimSpecific[match(nodes$id,GWAS_sQTL$event_id)]==1)]=substr(colERC[4],1,7)

nodes$vertex.frame.color="#00000000"
nodes$vertex.frame.width= 1
nodes$vertex.frame.color[which(GWAS_sQTL$nsSpecific_event[match(nodes$id,GWAS_sQTL$event_id)])]=substr(colERC[6],1,7)
nodes$vertex.frame.color[which(GWAS_sQTL$stimSpecific_event[match(nodes$id,GWAS_sQTL$event_id)])]=substr(colERC[4],1,7)

net <- graph_from_data_frame(d=edges,vertices=nodes,directed=F)

plot(net,layout=l2,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color=grey(0.2),
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label,rescale=F)
		
tk_set_coords(tkid,l2*250+500)
l <- tkplot.getcoords(tkid)
l3=l
l3[,2]=(l3[,2]-250)/250
l3[,1]=(l3[,1]-750)/250


net <- graph_from_data_frame(d=edges,vertices=nodes,directed=F)
plot(net,layout=l3,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color=grey(0.2),
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label,rescale=F)
save(l3,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/GWAS_sQTL_plot_coords_withSNP_v2.Rdata',HOME))


plot(net, layout=l3,mark.groups=list(c(1,4,5,8), c(15:17)), mark.col=c("#C5E5E7","#ECD89A"), mark.border=NA)

		