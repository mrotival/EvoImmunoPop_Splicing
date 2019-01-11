# select CHR, cond, population 
#CHR=22
#cond=1
#pop='ALL'

# select analysis 
#GENE_STANDARD=TRUE  (0,1,2 linear model)
#GENE_ROBUST=FALSE	 (2 df model on ranks Kruskal wallis -> robust)
#TRANSCRIPT_STANDARD=FALSE (0,1,2 linear model on transcripts)
#TRANSCRIPT_ROBUST=FALSE (2 df model on ranks Kruskal wallis on transcripts)

TEST=FALSE
pop='ALL'
CisDist=1e6
pvCis=1e-3
minFreq=0.05
options(stringsAsFactors=FALSE,max.print=9999)
.libPaths('/pasteur/entites/Geh/Shared_Programs/Rlib')
library(snpStats)
library(MatrixEQTL)
library(fdrtool)
#Geno=read.plink("/pasteur/projets/evo_immuno_pop/Maxime/GenotypeCalls/Omni5/EvoImmunoPop_Omni5_200x3489108.bed")
if(!is.null(CHR)){
	Geno=read.plink(paste(EVO_IMMUNO_POP,'/DATA_FREEZE/ERC_Main_Genotyping_24022015/Imputation/EvoImmunoPop_imputation_200x19619457_chr',CHR,'.bed',sep=''))
	GenoNum=as(Geno$genotype,'numeric') 
	#Map=Geno$map
	#SNPfreq_AF=1-apply(GenoNum[grep('AF',rownames(GenoNum)),],2,mean,na.rm=T)/2
	#SNPfreq_EU=1-apply(GenoNum[grep('EU',rownames(GenoNum)),],2,mean,na.rm=T)/2
	#nAF=1-apply(!is.na(GenoNum[grep('AF',rownames(GenoNum)),]),2,sum,na.rm=T)
	#nEU=1-apply(!is.na(GenoNum[grep('EU',rownames(GenoNum)),]),2,sum,na.rm=T)
	#SNPfreq=(nAF*SNPfreq_AF+nEU*SNPfreq_EU)/(nAF+nEU)
	#Map$chromPos=Map$position+cumsum(c(0,By(Map$position,Map$chromosome,max)[1:23]+chr.dist))[Map$chromosome]
	#Map$SNPfreq=SNPfreq
	#Map$SNPfreq_AF=SNPfreq_AF
	#Map$SNPfreq_EU=SNPfreq_EU
	#Map$FST=1-(nEU*SNPfreq_EU*(1-SNPfreq_EU) + nAF*SNPfreq_AF*(1-SNPfreq_AF))/((nEU+nAF)*(1-SNPfreq)*SNPfreq)
	#Map$locus=paste(Map$chromosome,':',round(Map$position/1e6),'Mb',sep='')
	#Map$posID=paste(Map$chromosome,Map$position,sep=':')
	#save(Map,file=paste('/pasteur/projets/evo_immuno_pop/Maxime/SNP_annotation/imputed_Genotypes/MapFiles/Map_imputation_200x19619457_chr',CHR,'.Rdata',sep=''))
	load(file=paste(EVO_IMMUNO_POP,'/Maxime/SNP_annotation/imputed_Genotypes/MapFiles/Map_imputation_200x19619457_chr',CHR,'.Rdata',sep=''))
	rm(Geno);gc()
	chr.dist=5e7
	By=function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x} # by own version of by (makes a vector with names)
	Chrsep=cumsum(chr.dist+ c(-chr.dist,By(Map$position,Map$chromosome,max)))-chr.dist/2
	toRemove=read.table(paste(EVO_IMMUNO_POP,'/Maxime/SNP_annotation/imputed_Genotypes/multiallelic_Sites.txt',sep=''),header=T)
	GenoNum=GenoNum[,!Map$snp.name%in%toRemove$snp.name]
	Map=Map[which(!Map$snp.name%in%toRemove$snp.name),]
}

load('/pasteur/homes/mrotival/02_ERC_Data/FPKM_Data/allFPKM_data_GC_ComBat_lib-exp-date_filtered.Rdata')
SampleAnnot[nchar(SampleAnnot[,'population'])==2,'population']=paste(SampleAnnot[nchar(SampleAnnot[,'population'])==2,'population'],'B',sep='')

# read GeneAnnot
GeneAnnot=read.table('/pasteur/homes/mrotival/Annotation/GeneAnnotation_hg37_ens70.txt',sep='\t',comment='',quote='',stringsAsFactors=FALSE,header=T,row.names=1)
GeneAnnot=GeneAnnot[which(GeneAnnot$select),]
GeneAnnot=GeneAnnot[match(rownames(FPKM_gene),GeneAnnot$Ensembl.Gene.ID),]
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='X']=23
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='Y']=24
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='MT']=26
GeneAnnot$Chromosome.Name=as.numeric(GeneAnnot$Chromosome.Name)

# read TranscriptAnnot
#TranscriptAnnot=read.table('/pasteur/homes/mrotival/Annotation/TranscriptAnnotation_hg37_ens70.txt',sep='\t',comment='',quote='',stringsAsFactors=FALSE)
#TranscriptAnnot=TranscriptAnnot[match(rownames(FPKM),TranscriptAnnot$Ensembl.Transcript.ID),]
#TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='X']=23
#TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='Y']=24
#TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='MT']=26
#TranscriptAnnot$Chromosome.Name=as.numeric(TranscriptAnnot$Chromosome.Name)

# read FPKM
#TranscriptAnnot=TranscriptAnnot[which(isoAnnot$mF1),]
#FPKM=FPKM[which(isoAnnot$mF1),]
#isoAnnot=isoAnnot[which(isoAnnot$mF1),]
#FPKM=log2(1+pmax(2^FPKM-1,0))
#
# transcript_Annot
#IsoformPos=cbind(transcript=isoAnnot$tracking_id,TranscriptAnnot[match(isoAnnot$tracking_id,TranscriptAnnot$Ensembl.Transcript.ID),c(2,3,4)])
#rownames(IsoformPos)=isoAnnot$tracking_id
#FPKM_ratio=FPKM-FPKM_gene[match(isoAnnot$gene_id,rownames(FPKM_gene)),]
#FPKM_ratio=FPKM_ratio[GeneAnnot[match(isoAnnot$gene_id,rownames(FPKM_gene)),'Nb_transcript']>1,]
#isoAnnot_ratio=isoAnnot[GeneAnnot[match(isoAnnot$gene_id,rownames(FPKM_gene)),'Nb_transcript']>1,]
#IsoformPos_ratio=IsoformPos[match(rownames(FPKM_ratio),rownames(IsoformPos)),]

########################################################################
##				Genes by Pop - Standard - Cis 						  ##
########################################################################

dir.create(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/permutations/perm',perm,sep=''))
RESCIS=list()
#qqbin=seq(0,308.70,0.05)
tim=Sys.time()
indiv=sort(unique(SampleAnnot[,'individu']))
if(perm>0){
	set.seed(perm)
	indiv[grep('AFB',indiv)]=sample(indiv[grep('AFB',indiv)])
	indiv[grep('EUB',indiv)]=sample(indiv[grep('EUB',indiv)])
	}

for(cond in 1:5){
	if(pop=='ALL'){
		wSNP=which(Map$SNPfreq > minFreq)#[1:10000]
		if(TEST){wSNP=wSNP[1:10000]}
		PCs=read.table(paste('/pasteur/homes/mrotival/02_ERC_Data/SampleCaseReports/EvoImmunoPop_',pop,'_ancestryPCA.txt',sep=''),header=T)
		wCondPop=which(SampleAnnot$condition==cond)
		indiv=indiv[indiv%in%SampleAnnot[wCondPop,'individu']]
		wCondPop=match(paste(sort(indiv),cond,sep='-'),colnames(FPKM_gene))
	}
	GenoNum_Mat = SlicedData$new()
	GenoNum_Mat$CreateFromMatrix(t(GenoNum)[wSNP,match(indiv,rownames(GenoNum))])
	GenoNum_Mat$ResliceCombined(sliceSize = 5000)
	FPKM_Mat = SlicedData$new()
	FPKM_Mat$CreateFromMatrix(t(apply(FPKM_gene[,wCondPop],1,function(x){qnorm(rank(x,ties.method='average')/(length(x)+1))})))
	show(FPKM_Mat)
	Cvrt=SlicedData$new()
	Cvrt$CreateFromMatrix(t(PCs[match(indiv,PCs$IID),3:4]))
	res=Matrix_eQTL_main(GenoNum_Mat,FPKM_Mat,Cvrt,output_file_name=NULL, pvOutputThreshold=0,output_file_name.cis=NULL,pvOutputThreshold.cis=pvCis,cisDist=CisDist,snpspos=Map[wSNP,c(2,1,4)],genepos=GeneAnnot[,c(1,2,4,5)],pvalue.hist='qqplot')
	if(nrow(res$cis$eqtls)>0){
		res$cis$eqtls$condition=cond
		res$cis$eqtls$population=pop
		res$cis$eqtls$isCis=TRUE
		# save qqplot data
		RESCIS[[cond]]=res$cis$eqtls
	}
	print(Sys.time()-tim)
}
RESCIS=do.call(rbind,RESCIS)
RESCIS$snps=as.character(RESCIS$snps)
RESCIS$gene=as.character(RESCIS$gene)
wSNP=match(RESCIS$snp,rownames(Map))
wGene=match(RESCIS$gene,rownames(GeneAnnot))
RESCIS$CisDist=ifelse(Map[wSNP,'chromosome']==GeneAnnot[wGene,'Chromosome.Name'],pmax(GeneAnnot[wGene,'Gene.Start..bp.']-Map[wSNP,'position'],Map[wSNP,'position']-GeneAnnot[wGene,'Gene.End..bp.']),Inf)
if(perm==0){
	if(CisDist==1e5 & pvCis==1){
		save(RESCIS,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/permutations/perm',perm,'/Cis-eQTL_ALL_allCond_chr',CHR,'_perm',perm,'_gene_Full_100kb.Rdata',sep=''))
	}else{
		save(RESCIS,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/permutations/perm',perm,'/Cis-eQTL_ALL_allCond_chr',CHR,'_perm',perm,'_gene.Rdata',sep=''))
	}
}else{
	thresholds=10^-c(seq(3,13,0.1),15,20,30,50)
	NbEQTL=sapply(thresholds,function(th){luq(RESCIS$gene[RESCIS$pval < th & abs(RESCIS$CisDist) < 1e6])})
	if(CisDist==1e5 & pvCis==1){
#		save(NbEQTL,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/permutations/perm',perm,'/Cis-eQTL_ALL_allCond_chr',CHR,'_perm',perm,'_gene.Rdata',sep=''))
	}else{
		save(NbEQTL,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/permutations/perm',perm,'/Cis-eQTL_ALL_allCond_chr',CHR,'_perm',perm,'_gene_100kb.Rdata',sep=''))
	}
}

#load(file='/pasteur/projets/evo_immuno_pop/Maxime/transEQTL_Clust/CisTransObserved.Rdata')

q('no')