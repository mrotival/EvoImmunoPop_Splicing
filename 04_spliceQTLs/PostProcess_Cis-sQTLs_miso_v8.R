#######################################################
###				Aggregate psiQTL Pop Combined		###
#######################################################


# set the path to files that will be needed

allPSI_events=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/Ens70_HISAT/PSI_events_ALL_V7.2_withCounts.Rdata',sep='')
Adjusted_PSI_values=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME)

library(data.table)
library(GenomicRanges)

# allGOterms
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

# load data 
load(allPSI_events)
load(Adjusted_PSI_values)

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
			load('data/sQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_V8.Rdata',sep=''))
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

MeanExpr=GeneAnnot[,c("NS_mean", "LPS_mean", "PAM3_mean", "R848_mean", "Flu_mean")]
lFC=log2(1+MeanExpr[,-1])-log2(1+MeanExpr[,1])%o%rep(1,4)

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
RESobs$coding_type=PSI_Annot$coding_type[mm]
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
			load('data/sQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_V8.Rdata',sep='')) 
            RESCHR[[paste(cond,pop,CHR)]]=as.data.table(RESCIS)
            }
        RES[[CHR]]=rbindlist(RESCHR)
		}
    RES=rbindlist(RES)
    colnames(RES)[colnames(RES)=='gene']='event_id'
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

# compute FDR for each threshold 
FDR_threshold=cbind(thresholds,NbFP_gene = NbAssocGene, NbFP_gene= NbFP_gene, FDR_gene=NbFP_gene/NbAssocGene)
FDR_threshold=as.data.frame(FDR_threshold)

# assign an FDR for each p-value through linear extrapolation on the log10 scale (for FDR < 1E-6, FDR is extrapolated proportionally to the p-value. ) 
FDR_compute=function(pval,FDR_th){y=approxfun(c(0,-log10(FDR_threshold$thresholds),500),c(pmin(-log10(c(1,FDR_th)),6),500))(-log10(pval)); 10^-y}
RESobs$FDR_1Mb=FDR_compute(RESobs$pvalue,FDR_threshold[['FDR_gene']])

# store observed associations and FDR threshold matrix
save(RESobs,FDR_threshold,file='data/sQTL/cis-psiQTL_MISO_counts_V8.1.Rdata')

luq(RESobs$gene[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]) # 993
luq(RESobs$event_id[RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable]) # 1464

# keep only the best SNP for each event
RESobs_nodup_1Mb=RESobs[which(RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable),]
RESobs_nodup_1Mb=RESobs_nodup_1Mb[order(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$pval),]
RESobs_nodup_1Mb=RESobs_nodup_1Mb[!duplicated(RESobs_nodup_1Mb$event_id),]

# keep the best SNP for each event/cond (so that we capture multiple sQTL if the best sQTL change upon stimulation)  
RESobs_nodup_1Mb_cond=RESobs[which(RESobs$FDR_1Mb<0.05 & abs(RESobs$CisDist)<1e6 & RESobs$isTestable),]
RESobs_nodup_1Mb_cond=RESobs_nodup_1Mb_cond[order(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond,RESobs_nodup_1Mb_cond$pval),]
RESobs_nodup_1Mb_cond=RESobs_nodup_1Mb_cond[!duplicated(RESobs_nodup_1Mb_cond$event_id,RESobs_nodup_1Mb_cond$cond),]


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

Testable_sQTL=TESTABLE[match(RESobs_nodup_1Mb$event_id,PSI_Annot[,'event_id']),]
colnames(Testable_sQTL)=paste('Alternatively_spliced',condIndex,sep='_')

Expressed_sQTL=(PCTNA<0.05 & SUPPORT>10)[match(RESobs_nodup_1Mb$event_id,PSI_Annot[,'event_id']),]
colnames(Expressed_sQTL)=paste('Expressed',condIndex,sep='_')

Testable_DIFF_sQTL=TESTABLE_DIFF[match(RESobs_nodup_1Mb$event_id,PSI_Annot[,'event_id']),-1]

RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,Pval_sQTL,Beta_sQTL,Expressed_sQTL,Testable_sQTL)

###################################
###		cluster linked sQTLs	###
###################################

Geno_sQTL=getGenos(unique(RESobs_nodup_1Mb$snps))
r2_sQTL=cor(t(as.matrix(Geno_sQTL[-(1:5)])),use='p')^2
library(igraph)
cor_graph=graph_from_adjacency_matrix(r2_sQTL>0.8)
wc <- cluster_walktrap(cor_graph)
haplos=membership(wc)
names(haplos)=Geno_sQTL$snp.name
RESobs_nodup_1Mb$haplo=paste('haplo',haplos)[match(RESobs_nodup_1Mb$snps,names(haplos))]


write.table(RESobs_nodup_1Mb,file='data/RESobs_nodup_1Mb.txt',quote=F,sep='\t',row.names=F)
write.table(RESobs_nodup_1Mb_cond,file='data/RESobs_nodup_1Mb_cond.txt',quote=F,sep='\t',row.names=F)




