#!/bin/sh

# for CHR in `seq 1 22` ;
#  do sbatch --array=0 -p geh --qos=geh --mem-per-cpu=40000 -J seQTL -o "/pasteur/homes/mrotival/JobOutput/JobOutput_seQTL%a.log" /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/eQTL_ALL/cis_eQTL_tars_runner.sh $CHR
# done;

# sbatch --array=0 -p geh --qos=geh --mem-per-cpu=40000 -J seQTL -o "/pasteur/homes/mrotival/JobOutput/JobOutput_seQTL%a.log" /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/eQTL_ALL/cis_eQTL_tars_runner.sh $CHR

# for CHR in `seq 1 22` ;
#  do sbatch --array=2-100 --mem-per-cpu=40000 -J seQTL -o "/pasteur/homes/mrotival/JobOutput/JobOutput_seQTL%a.log" /pasteur/homes/mrotival/01_scripts/Splicing/PapierSplicing/eQTL_ALL/cis_eQTL_tars_runner.sh $CHR
# done;


source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module load R/3.1.2

# PERMDIR="/pasteur/scratch/mrotival/sQTL"
# mkdir $PERMDIR
LOGDIR="/pasteur/homes/mrotival/JobOutput/"

PERM=${SLURM_ARRAY_TASK_ID}
CHR=$1

FILENAME=${LOGDIR}/scripts/run_eQTLs_chr${CHR}_perm${PERM}_ALL_allCond.R
echo "CHR=${CHR}" > ${FILENAME}
echo "perm=${PERM}" >> ${FILENAME}
cat /pasteur/homes/mrotival/01_scripts/eQTLs/imputation_eQTLs/eQTL_pop-cond-chr_perm_tars_ALL.R >> ${FILENAME}

Rscript ${FILENAME} ||exit 1
echo "R finished Running on Node "
srun hostname || exit 2
exit 0

# nperm=1
#  df=data.frame(chr=rep(rep(1:22,(nperm+1)),e=5),cond=rep(1:5,(nperm+1)*22),perm=rep(0:nperm,each=22*5))
#  write.table(df,file='/pasteur/homes/mrotival/JobOutput/JobList_sQTL.txt',sep='\t',quote=F,row.names=F,col.names=F)

# df=data.frame(cond=rep(1:5,11),perm=rep(0:10,each=5))
# write.table(df,file='/pasteur/homes/mrotival/JobOutput/JobList_ds.txt',sep='\t',quote=F,row.names=F,col.names=F)



#### extract list of 22264 imputed variants that are multi-allelic (multiple annotated variants at the same position)
# these variants have been removed from Map_imputed bringing it down from 19619457 SNPs to 19597193 


############# create eQTL list

NbeQTL=list()
Res_eQTL=list()
for (perm in 0:100){
	cat('\n',perm,':')
	for (CHR in 1:22){
		cat(CHR,'')
		x=load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/permutations/perm',perm,'/Cis-eQTL_ALL_allCond_chr',CHR,'_perm',perm,'_gene.Rdata',sep=''))
		if(perm>0){
			if(CHR==1){
				NbeQTL[[perm]]=NbEQTL
			}else{
				NbeQTL[[perm]]=NbeQTL[[perm]]+NbEQTL
				}
			}
		if(perm==0){
			Res_eQTL[[CHR]]=RESCIS
			}
		}
	}
Res_eQTL=do.call(rbind,Res_eQTL)
NbeQTL=do.call(rbind,NbeQTL)

thresholds=10^-c(seq(3,13),15,20,30,50)
NbEQTL_obs=sapply(thresholds,function(th){luq(Res_eQTL$gene[Res_eQTL$pval < th & pmax(Res_eQTL$CisDist,0) < 1e6])})
FDR_CisGeneNballCondPop=apply(NbeQTL,2,mean)/NbEQTL_obs
FDR_compute=function(pval){ y=approxfun(c(0,-log10(thresholds),500),c(pmin(-log10(c(1,FDR_CisGeneNballCondPop)),6),500))(-log10(pval)); 10^-y}
Res_eQTL$FDR=FDR_compute(Res_eQTL$pvalue)
Res_eQTL_nodup=Res_eQTL[Res_eQTL$FDR<0.05,]
Res_eQTL_nodup=Res_eQTL_nodup[order(Res_eQTL_nodup$gene,Res_eQTL_nodup$cond,Res_eQTL_nodup$pval),]
Res_eQTL_nodup=Res_eQTL_nodup[!duplicated(paste(Res_eQTL_nodup$gene,Res_eQTL_nodup$cond)),]
Res_eQTL_nodup$symbol=G2S(Res_eQTL_nodup$gene)

#w=match(Res_eQTL_nodup$snps,Map_imputed$snp.name)
#Res_eQTL_nodup[['MAF']]=Map_imputed$SNPfreq[w]
#Res_eQTL_nodup[['MAF_AFB']]=Map_imputed$SNPfreq_AF[w]
#Res_eQTL_nodup[['MAF_EUB']]=Map_imputed$SNPfreq_EU[w]
# Res_eQTL_nodup[['FST']]=Map_imputed$FST_adj[w]
# Res_eQTL_nodup[['iHS_AFB']]=Map_imputed$iHS_AFB[w]
# Res_eQTL_nodup[['iHS_EUB']]=Map_imputed$iHS_EUB[w]

Res_eQTL_nodup$beta=-Res_eQTL_nodup$beta
for(i in 1:5){
	w=match(paste(Res_eQTL_nodup$snps,Res_eQTL_nodup$gene,i),paste(Res_eQTL$snps,Res_eQTL$gene,Res_eQTL$cond))
	Res_eQTL_nodup[[paste('Pvalue',condIndex[i],sep='_')]]=Res_eQTL$pvalue[w]
	Res_eQTL_nodup[[paste('Beta',condIndex[i],sep='_')]]=Res_eQTL$beta[w]
}
#Res_eQTL_nodup= Res_eQTL_nodup[,c(1:2,11,3:5,13:25)]
colnames(Res_eQTL_nodup)[5]='MinPvalue'
write.table(Res_eQTL_nodup,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/eQTL_rerun_PopCombined_allCond.txt',sep=''),sep='\t',quote=F,row.names=F)
Res_eQTL_nodup=Res_eQTL_nodup[order(Res_eQTL_nodup$gene,Res_eQTL_nodup$pval),]
Res_eQTL_nodup=Res_eQTL_nodup[!duplicated(Res_eQTL_nodup$gene),]

write.table(Res_eQTL_nodup,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/eQTL_rerun_PopCombined.txt',sep=''),sep='\t',quote=F,row.names=F)

# Res_eQTL_nodup=read.table(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/eQTL_rerun_PopCombined_allCond.txt',sep=''),sep='\t',header=T)

library(snpStats)
SNP_genos=list()
for (CHR in 1:22){
	cat(CHR,'')
	Geno=read.plink(paste(EVO_IMMUNO_POP,'/DATA_FREEZE/ERC_Main_Genotyping_24022015/Imputation/EvoImmunoPop_imputation_200x19619457_chr',CHR,'.bed',sep=''))
	GenoNum=as(Geno$genotype,'numeric') 
	SNP_genos[[CHR]]=GenoNum[,na.omit(match(unique(Res_eQTL_nodup$snps),colnames(GenoNum)))]
}

SNP_genos=do.call(cbind,SNP_genos)
# change to minor allele number
SNP_genos=t(2-SNP_genos)
write.table(data.frame(snp_id=rn(SNP_genos),SNP_genos),file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/SNP_genos_eQTL_rerun_bestSNPperCondition.txt',sep=''),sep='\t',quote=F,row.names=F)

# SNP_genos=read.table(file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/SNP_genos_eQTL_rerun_bestSNPperCondition.txt',sep=''),sep='\t',header=T)
