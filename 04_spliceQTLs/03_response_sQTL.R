

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
PSI_Annot=PSI_Annot[toKeep,]
load(Adjusted_PSI_values)

# count unique elements of a vector
luq=function(x){length(unique(x))}

RESobs_nodup_1Mb=fread('data/RESobs_nodup_1Mb_Overlap_eQTL.txt')


##########################################################################################
#####				frequency and sharing of sQTLs	across conditions 			      ####
##########################################################################################

################ PIECHART sQTL SHARING
expressed_sharing=apply(RESobs_nodup_1Mb[,paste('Expressed_',condIndex,sep='')],1,sum)

myPSIs=t(PSI_prov[match(RESobs_nodup_1Mb$event_id[expressed_sharing==5],PSI_Annot[,'event_id']),])
mySNPs=snps_sQTL[,match(RESobs_nodup_1Mb$snps[expressed_sharing==5],colnames(snps_sQTL))]
# we assign the genotype of the individual to each sample (duplicating genotype across the 5 conditions)
indiv=substr(rownames(myPSIs),1,6)
mySNPs_cond=mySNPs[indiv,]

# define all 2^5 possible models
allModels=cbind(NS=rep(rep(0:1,e=16),1),
				LPS=rep(rep(0:1,e=8),2),
				PAM3CSK4=rep(rep(0:1,e=4),4),
				R848=rep(rep(0:1,e=2),8),
				IAV=rep(rep(0:1,e=1),16))

# compute the vector of likelihoods
LL_full=apply(rbind(mySNPs_cond,myPSIs),2,function(SNP_and_PSI){
    n=length(SNP_and_PSI)/2
    mySNP_cond=SNP_and_PSI[1:n]
    myPSI_cond=SNP_and_PSI[n+1:n]

    LL_psi=apply(allModels,1,function(x,myPSI,SNP_cond){
		mySNP=SNP_cond*rep(x,table(SampleAnnot$cond))
		# for each condition, assign the 
		LL=logLik(lm(myPSI~mySNP+SampleAnnot$CondPop))
	},myPSI_cond,mySNP_cond)
	LL_psi
    })
# assign a probability to each model
LL_probs=apply(LL_full, 2,function(LL_psi){exp(LL_psi-min(LL_psi))/sum(exp(LL_psi-min(LL_psi)))})

# get the most likely model
WhichCond=apply(LL_full, 2,function(LL_psi){allModels[which.max(LL_psi),]})

#############################################################
#################        Figure S5A         #################
#############################################################
table(apply(WhichCond,2,paste,collapse='')

###### use these numbers in Venn diagrams of Figure S5A
#00001 00010 00011 00100 00101 00110 00111 01000 01001 01010 01011 01100 01101 01110 01111 10000 10001 10010 10011 10100 10101 10110 10111 11000 11001 11010 11011 11100 11101 11110 11111 
#   26     7     8     3     3     2     3    10     2     1     2     2     1    22    35     5     8     1     4     4     4     1     4     1     4     1    14     4    16    78   684 

#############################################################
#################        Figure 4E         #################
#############################################################

par(mar=c(3,3,3,3))
tab=rev(table(apply(WhichCond,2,sum)))
tabpct=paste(round(100*tab/sum(tab), 1),'%')
pie(tab,col=acol[c(1,3,4,6,7)],init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab,col=acol[c(1,3,4,6,7)],init.angle=90,labels=rep(' ',5)) # 4 x 4 inches

####################################################################################################################################
############################################ END sQTL SHARING ######################################################################
####################################################################################################################################

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
			load('data/sQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,'_response_V8.Rdata',sep=''))
			RES[[paste(cond,pop,CHR)]]=RESCIS
		}
	}
}
RES=do.call(rbind,RES)
colnames(RES)[colnames(RES)=='gene']='event_id'
RESobs_resp=RES
rm(RES)
gc()

save(RESobs_resp,file=sprintf('data/cis-psiQTL_MISO_counts_response_V8.Rdata',EVO_IMMUNO_POP))

RESobs_has_rsQTL_1Mb=RESobs_resp[which(abs(RESobs_resp$CisDist)<1e6 & paste(RESobs_resp$event_id,RESobs_resp$snps)%in%paste(RESobs_nodup_1Mb$event_id,RESobs_nodup_1Mb$snps)),]

Pval_rsQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='pvalue');
Pval_rsQTL=Pval_rsQTL[match(RESobs_nodup_1Mb$event_id,Pval_rsQTL$event_id),-1]
Pval_rsQTL[is.na(Pval_rsQTL)]=1
colnames(Pval_rsQTL)=paste('Pvalue_rsQTL',condIndex[-1],sep='_')
Beta_rsQTL=dcast(RESobs_has_sQTL_1Mb, event_id~condition,value.var='beta')
Beta_rsQTL=Beta_rsQTL[match(RESobs_nodup_1Mb$event_id,Beta_rsQTL$event_id),-1]
Beta_rsQTL[is.na(Beta_rsQTL)]=0
colnames(Beta_rsQTL)=paste('Beta_rsQTL',condIndex[-1],sep='_')

RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,Pval_rsQTL,Beta_rsQTL)


#### Annotate Gene Expression and logFC 
MeanExpr=GeneAnnot[,c("NS_mean", "LPS_mean", "PAM3_mean", "R848_mean", "Flu_mean")]
MeanExpr[is.na(MeanExpr)]=0
colnames(MeanExpr)=paste('FPKM',condIndex,sep='_')
lFC=log2(1 + MeanExpr[,-1])-log2(1+MeanExpr[,1])%o%c(1,1,1,1)
colnames(lFC)=paste('log2FC',condIndex[-1],sep='_')

GeneAnnot=cbind(GeneAnnot,lFC,max_logFC=apply(lFC,1,max))

# What fraction of the genes under study are expressed only upon stimulation stimulation specific
sum(GeneAnnot[unique(PSI_Annot$gene),'NS_mean']<10 & GeneAnnot[unique(PSI_Annot$gene),'max_logFC']>1) # 274

#### Annotate Gene Expression and logFC in sQTL table
RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,MeanExpr[RESobs_nodup_1Mb$gene,])
RESobs_nodup_1Mb=cbind(RESobs_nodup_1Mb,lFC[RESobs_nodup_1Mb$gene,])

RESobs_nodup_1Mb$max_logFC=pmax(RESobs_nodup_1Mb$log2FC_LPS,RESobs_nodup_1Mb$log2FC_PAM3CSK4,RESobs_nodup_1Mb$log2FC_R848,RESobs_nodup_1Mb$log2FC_IAV)

#### Annotate sQTL_of_stimulationSpecificGene and rsQTL
RESobs_nodup_1Mb$is_sQTL_of_stimulationSpecificGene=ifelse(RESobs_nodup_1Mb$FPKM_NS<10 & RESobs_nodup_1Mb$max_logFC>1,'yes','no')
RESobs_nodup_1Mb$rsQTL_tested=ifelse(RESobs_nodup_1Mb$FPKM_NS>10 & RESobs_nodup_1Mb$Pvalue_NS>0.001,'yes','no')
RESobs_nodup_1Mb$is_rsQTL= ifelse(RESobs_nodup_1Mb$rsQTL_tested=='yes' & apply(Pval_rsQTL<0.01,1,any), 'yes', 'no')

### How many sQTL of stimulation specific genes (affecting how many genes ?)
sum(RESobs_nodup_1Mb$is_sQTL_of_stimulationSpecificGene=='yes')
#[1] 108
length(unique(RESobs_nodup_1Mb$gene[RESobs_nodup_1Mb$is_sQTL_of_stimulationSpecificGene=='yes']))
# 74

## how many sQTL where tested for rsQTL ?
sum(RESobs_nodup_1Mb$rsQTL_tested=='yes') # 228

sum(Pval_rsQTL[RESobs_nodup_1Mb$rsQTL_tested=='yes',]<0.01) # 172
sum(RESobs_nodup_1Mb$is_rsQTL=='yes') # 129

###############################################
###		Table S3D	response sQTLs			###
###############################################

paste_xy=function(x,y){paste(rep(x,length(y)),rep(y,e=length(x)),sep='')}
TableS3D=RESobs_nodup_1Mb[,c('event_id','symbol','event_type','snps','is_sQTL_of_stimulationSpecificGene','is_rsQTL', 'FPKM_NS',paste_xy(c('FPKM_','log2FC_'),condIndex[-1]),'Pvalue_sQTL_NS','rsQTL_tested',paste_xy(c('Pvalue_rsQTL','Beta_rsQTL'),condIndex[-1]))]

write.table(TableS3D,file='tables/TableS3D_rsQTL.txt',quote=F,sep='\t',row.names=F)

###################################
###		Table S3D END			###
###################################

###################################
###		Figures 4E and 4F		###
###################################

i0=0
library(RcolorBrewer)
for (j in which(RESobs_nodup_1Mb$event_id%in%c("ENSG00000157601;RI:21:42793765:42793798-42793902:42794022:+","ENSG00000185499;A3:1:155162101-155162577:155162074-155162577:-"))){
i0=i0+1
cat(i0,'')
symbol=RESobs_nodup_1Mb$symbol[j]
type=RESobs_nodup_1Mb$event_type[j]
snp=RESobs_nodup_1Mb$snps[j]
map=getMapInfo(snp)
event_id=RESobs_nodup_1Mb$event_id[j]
gene=RESobs_nodup_1Mb$gene[j]
myExpr=2^FPKM_gene[gene,]-1
myPSI=PSI_prov[PSI_Annot[,"event_id"]==event_id,]

ymax=ceiling(max(log2(1+myExpr)))
ymin=min(log2(1+myExpr))

cat(ymax,ymin,'\n')

mySNP=2-getSNP(snp)

mySNP_cond=mySNP[substr(names(myPSI),1,6)]
allele.1=map$allele.2 # minor allele 
allele.2=map$allele.1 # major allele
genotypes_alleles=paste(c(allele.1,allele.1,allele.2),c(allele.1,allele.2,allele.2),sep='')
MeanExprGenoCond=By(log2(1+myExpr),paste(SampleAnnot$condition,mySNP_cond),mean)
MeanExprGenoCond=MeanExprGenoCond[!grepl('NA',names(MeanExprGenoCond))]
myMain=paste(symbol,'-',type,'- NS: FPKM =',round(GeneAnnot[gene,'NS_mean']),', Pval=',Format(RESobs_nodup_1Mb$Pvalue_NS[j],2))
pdf(sprintf('figures/sQTL_%s_PSI_%s_%s_%s.pdf',i0,symbol,type,snp),width=7,height=5)
plot(0:1,0:1,ylim=c(-0.15,1),xlim=c(0.5,20),col='#00000000',axes=F,xlab='',ylab='PSI',main=myMain)
axis(2,las=1)
axis(1,at=c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=4+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=8+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=12+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=16+c(1:3),labels=genotypes_alleles,las=2)
for(i in 1:5){
	w=which(SampleAnnot$condition==unique(SampleAnnot$condition)[i])
	boxplot(myPSI[w][mySNP_cond[w]==0],at=(i-1)*4+1,col=colERC5[i],notch=T,pch=16,cex=0.7,outcol='#00000088',add=T,axes=F,pars=list(boxwex=1.5))
	boxplot(myPSI[w][mySNP_cond[w]==1],at=(i-1)*4+2,col=colERC5[i],notch=T,pch=16,cex=0.7,outcol='#00000088',add=T,axes=F,pars=list(boxwex=1.5))
	boxplot(myPSI[w][mySNP_cond[w]==2],at=(i-1)*4+3,col=colERC5[i],notch=T,pch=16,cex=0.7,outcol='#00000088',add=T,axes=F,pars=list(boxwex=1.5))
}
   x=c(1:3,4+1:3,8+1:3,12+1:3,16+1:3);y=-0.1;dx=0.4;dy=0.04
   colors=grey(1-(MeanExprGenoCond-ymin)/(ymax-ymin))
   rect(x+dx,y-dy,x-dx,y+dy,col=colors)
axis(1,at=4+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=8+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=12+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=16+c(1:3),labels=genotypes_alleles,las=2)	

dev.off()
}


###################################
###		END Figure 4    		###
###################################

write.table(RESobs_nodup_1Mb,file=sprintf('data/RESobs_nodup_1Mb_full.txt',HOME),quote=F,sep='\t',row.names=F)
