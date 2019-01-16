
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




RESobs_nodup_1Mb=fread('data/RESobs_nodup_1Mb_SpliceElt.txt')




###############################################
###		test overlap of sQTLs with eQTLs 	###
###############################################

# load gene expression data
FPKM_gene=as.matrix(fread('data/FPKM_gene.txt'))
SampleAnnot=as.data.frame(fread('data/SampleAnnot.txt'))

# load intronic expression
ByGene_Introns=fread('data/AllSamples_Intronic_Gene_count.txt'))
ByGene_Introns=ByGene_Introns[-1,]
RPKM_cols=colnames(ByGene_Introns)[grep('RPKM',colnames(ByGene_Introns))]
Introns_RPKM=as.matrix(ByGene_Introns[,mget(RPKM_cols)])
rownames(Introns_RPKM)=substr(ByGene_Introns[,intron_id],1,15)
colnames(Introns_RPKM)=substr(colnames(Introns_RPKM),1,8)

# load list of eQTLs
Res_eQTL_nodup_cond=read.table('data/eQTL_rerun_PopCombined_allCond.txt',sep='\t',header=T)

# load genotypes eQTLs
Genos_eQTL=read.table(file='data/SNP_genos_eQTL_rerun_bestSNPperCondition.txt',sep='\t',header=T)
snps_eQTL=as.matrix(Genos_eQTL[-1])
rownames(snps_eQTL)=gsub('.',':',Genos_eQTL[[1]],fixed=T)
snps_eQTL=t(snps_eQTL)


# load genotypes sQTLs
Genos_sQTL=getGenos(RESobs_nodup_1Mb$snps)
snps_sQTL=t(as.matrix(Genos_sQTL[-(1:5)]))
colnames(snps_sQTL)=Genos_sQTL[[1]]
snps_sQTL=snps_sQTL[,match(RESobs_nodup_1Mb$snps,colnames(snps_sQTL))]



# load PSI values sQTL (for the condition with the strongest sQTL)
psi_sQTL=mapply(function(event,cond){
			mmEvent=match(event,PSI_Annot[,'event_id']
			mmSample=match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(PSI_prov))
			PSI_prov[mmEvent,mmSample]},
		RESobs_nodup_1Mb$event_id,
		RESobs_nodup_1Mb$cond)

# load FPKM values of gene with an sQTL (for the condition with the strongest sQTL)
gene_sQTL=mapply(function(gene,cond){
			mmSample=match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(FPKM_gene))
			FPKM_gene[gene,mmSample]},
		RESobs_nodup_1Mb$gene,
		RESobs_nodup_1Mb$cond)

# load FPKM values of gene with an sQTL (for the condition with the strongest sQTL)
intron_sQTL=mapply(function(gene,cond){
			mmSample=match(paste(rownames(snps_sQTL),'-',cond,sep=''),colnames(Introns_RPKM))
			Introns_RPKM[gene,mmSample]},
		RESobs_nodup_1Mb$gene,
		RESobs_nodup_1Mb$cond)

# compute r2 between sQTL and eQTL
r2_sQTL_eQTL=mapply(function(gene,snp,cond){
		mysnps=Res_eQTL_nodup_cond$snps[Res_eQTL_nodup_cond$gene==gene & Res_eQTL_nodup_cond$cond==cond]
		if(length(mysnps)==0){NA)}else{
			LD_both=cor(snps_sQTL[match(rownames(snps_eQTL),rownames(snps_sQTL)),snp],snps_eQTL[,mysnps],use='p')^2
		c(LD_both)}
 },RESobs_nodup_1Mb$gene,RESobs_nodup_1Mb$snps,RESobs_nodup_1Mb$cond)

r2_sQTL_eQTL[is.na(RESobs_nodup_1Mb$r2_sQTL)]=0

## Association of sQTL to gene expression

Expr_effect=apply(rbind(snps_sQTL,gene_sQTL),2,function(x){
	n=length(x)/2
	SNP=x[1:n]
	GENE=x[n+1:n]
	POP=substr(rownames(snps_sQTL),1,3)
	summary(lm(GENE~SNP+POP))$coeff[2,c(1,4)]
})

B_sQTL_geneExpr=Expr_effect[,1]
P_sQTL_geneExpr=Expr_effect[,2]

## Association of sQTL to intronic reads
Intron_effect=apply(rbind(snps_sQTL,intron_sQTL),2,function(x){
	n=length(x)/2
	SNP=x[1:n]
	INTRON=x[n+1:n]
	POP=substr(rownames(snps_sQTL),1,3)
    cat('.')
    if(all(is.na(INTRON))){
		c(NA,NA)
	}else{
		ct=summary(lm(INTRON~SNP+POP))$coeff[2,c(1,4)]
	    ct
    }
})
B_sQTL_intronExpr=Intron_effect[,1]
P_sQTL_intronExpr=Intron_effect[,2]

# number of sQTL that overlap eQTLs or associate to gene expression 
table(RESobs_nodup_1Mb$r2_sQTL_eQTL_sameCond_0>0.8,RESobs_nodup_1Mb$P_sQTL_geneExpr<0.05)

RESobs_nodup_1Mb$r2_sQTL_eQTL=r2_sQTL_eQTL
RESobs_nodup_1Mb$P_sQTL_geneExpr=P_sQTL_geneExpr
RESobs_nodup_1Mb$B_sQTL_geneExpr=B_sQTL_geneExpr
RESobs_nodup_1Mb$P_sQTL_intronExpr=P_sQTL_intronExpr
RESobs_nodup_1Mb$B_sQTL_intronExpr=B_sQTL_intronExpr

###################################
###		Figure 4C				###
###################################

hist(r2_sQTL_eQTL,main='',xlab='r2 with eQTL',col='darkgrey',br=seq(0,1,0.1),border='darkgrey',las=2)
hist(r2_sQTL_eQTL[r2_sQTL_eQTL>0],add=T,br=seq(0,1,0.1),col='lightgrey',border='darkgrey')
hist(r2_sQTL_eQTL[r2_sQTL_eQTL>0.8],add=T,br=seq(0,1,0.1),col=colERC5[2],border='darkgrey')
#barplot(By(RESobs_nodup_1Mb$r2_sQTL_eQTL_sameCond_0>0.8,RESobs_nodup_1Mb$event_type,mean),col=colPSI)

###################################
###		add causality test 		###
###################################

getLikelihoods=function(x,y,z,u=NULL){
	# compute the likelihoods of reactive, causal and independant models
	if(!is.null(u)){
	    if(is.data.frame(u)| is.matrix(u)){
    	    u=as.data.frame(u)
	        w=which(!is.na(y*x*z) & apply(!is.na(u),1,all))
	    }else{
	    	w=which(!is.na(y*x*z*u))
	    	}
		}else{
		w=which(!is.na(y*x*z))
		}
	x=x[w]
	y=y[w]
	z=z[w]
	if(!is.null(u)){
	     if(is.data.frame(u)| is.matrix(u)){
        u=u[w,]
    	mod_xCausal=logLik(lm(y~x+.,data=u))+logLik(lm(x~z+.,data=u))
	    mod_yCausal=logLik(lm(x~y+.,data=u))+logLik(lm(y~z+.,data=u))
     	mod_xyInd=logLik(lm(y~z+.,data=u))+logLik(lm(x~z+.,data=u))
#      	mod_xyIndCor=logLik(lm(y~z+x+.,data=u))+logLik(lm(x~z+.,data=u))
#      	mod_yxIndCor=logLik(lm(x~y+z+.,data=u))+logLik(lm(y~z+.,data=u))
	     }else{
	    u=u[w]
        mod_xCausal=logLik(lm(y~x+u))+logLik(lm(x~z+u))
        mod_yCausal=logLik(lm(x~y+u))+logLik(lm(y~z+u))
        mod_xyInd=logLik(lm(y~z+u))+logLik(lm(x~z+u))
#		mod_xyIndCor=logLik(lm(x~y+z+u))+logLik(lm(y~z+u))
#		mod_yxIndCor=logLik(lm(y~x+z+u))+logLik(lm(x~z+u))
        }
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
	probs=exp(logLiks-min(logLiks))/sum(exp(logLiks-min(logLiks)))
	names(probs)=c('SNP->Splice->Expr','SNP->Expr->Splice','Splice<-SNP->Expr')
	probs
}

ProbsModel=apply(rbind(snps_sQTL,gene_sQTL,psi_sQTL),2,function(x){
	n=(length(x))/3
	SNP=x[1:n]
	GENE=x[n+1:n]
	PSI=x[2*n+1:n]
	POP=as.numeric(substr(rownames(snps_sQTL),1,3)=='AFB')
	P_independent=summary(lm(GENE~PSI+SNP+POP))$coeff[2,4]
	c(getProbs(PSI,GENE,SNP,POP),P_independent)
})


BestModel=apply(ProbsModel,2,function(x){ifelse(any(is.na(x)),NA,ifelse(x[4]<0.05,which.max(x[1:2]),3))})
BestModel_prob=pmax(ProbsModel[1,],ProbsModel[2,])/(ProbsModel[1,]+ProbsModel[2,])


# add these informations to sQTL table
RESobs_nodup_1Mb$BestModel=rownames(ProbsModel)[BestModel]
RESobs_nodup_1Mb$BestModel_Prob=BestModel_prob

############################################
####			Figure 4D				####
############################################

pdf('figures/Fig4D_bestModel_sQTL_P_expr_horiz.pdf',height=5.6,width=3)
w=which(P_sQTL_geneExpr<0.05)
tab=table(RESobs_nodup_1Mb[w,'BestModel'],ifelse(RESobs_nodup_1Mb$r2_sQTL_eQTL[w]>0.8,'yes','no'))
tab=cbind(yes=tab[,2],no=tab[,1])
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)*100
x=barplot(tab,col=getColors(11)[9:11],las=3)
par(xpd=T);
text(x,1.1*max(apply(tab,2,sum)),labels=n_sQTL);par(xpd=F,srt=90)
dev.off()

#                       yes        no
#SNP->Expr->Splice 56.77419  4.678363
#SNP->Splice->Expr 13.54839 34.697856
#Splice<-SNP->Expr 29.67742 60.623782

odds.ratio(table(RESobs_nodup_1Mb$P_sQTL_intronExpr<0.05,ifelse(RESobs_nodup_1Mb$P_sQTL_geneExpr>0.05,'noeQTL', RESobs_nodup_1Mb$BestModel))[,c(1:2)])
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 13.56186 27.70109 63.91699  0.05 1.836052e-35
odds.ratio(table(RESobs_nodup_1Mb$P_sQTL_intronExpr<0.05,ifelse(RESobs_nodup_1Mb$P_sQTL_geneExpr>0.05,'noeQTL', RESobs_nodup_1Mb$BestModel))[,c(3,2)])
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 3.461445 7.483113 18.05425  0.05 1.402946e-09

By(RESobs_nodup_1Mb$P_sQTL_intronExpr<0.05,ifelse(RESobs_nodup_1Mb$P_sQTL_geneExpr>0.05,'noeQTL', RESobs_nodup_1Mb$BestModel),mean,na.rm=T)
#          noeQTL SNP->Expr->Splice SNP->Splice->Expr Splice<-SNP->Expr 
#        0.2534626         0.9042553         0.5562130         0.5690789 


#################################################
######      Table S3C sQTL_eQTL_overlap    ######
#################################################
TableS3C=RESobs_nodup_1Mb
TableS3C$condition_char=condIndex[TableS3C$cond]
TableS3C$has_eQTL=ifelse(TableS3C$r2_sQTL_eQTL!=0,'yes','no')
write.csv(TableS3C[,c('event_id','symbol','event_type','snps','condition_char','has_eQTL','r2_sQTL_eQTL','P_sQTL_geneExpr','BestModel','P_sQTL_intronExpr')],file=sprintf('tables/TableS3C_eQTLoverlap.csv',HOME))



write.table(RESobs_nodup_1Mb,file='data/RESobs_nodup_1Mb_Overlap_eQTL.txt',quote=F,sep='\t',row.names=F)
