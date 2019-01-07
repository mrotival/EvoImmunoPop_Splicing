# select CHR, cond, population 
#CHR=22
#cond=1
#pop='ALL'

##################################################################################################
## INPUT:                                                                                       ##
## - CHR: chromosomes for which sQTL are being mapped                                           ##
## -                                                                                            ##
## Calls R packages:                                                                            ##
##   GenomicRanges                                                                              ##
##   data.table                                                                                 ##
##################################################################################################


LOGIT=TRUE
ROBUST=TRUE
# select analysis 
PSI_STANDARD=TRUE  	#(0,1,2 linear model on ratio)
PSI_RESPONSE=TRUE		#(0,1,2 linear model on PSI response)

Source='Ens70_HISAT'
TEST=FALSE # set to TRUE to run on an subset of data (for debugging puposes)
dir.create(paste(HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,sep=''))

minFreq=0.05
options(stringsAsFactors=FALSE,max.print=9999)
.libPaths('/pasteur/entites/Geh/Shared_Programs/Rlib/3.1.2')
library(snpStats)
library(MatrixEQTL)
#Geno=read.plink("/pasteur/projets/evo_immuno_pop/Maxime/GenotypeCalls/Omni5/EvoImmunoPop_Omni5_200x3489108.bed")
cat('ok')
if(!is.null(CHR)){
	Geno=read.plink(paste(EVO_IMMUNO_POP,'/DATA_FREEZE/ERC_Main_Genotyping_24022015/Imputation/EvoImmunoPop_imputation_200x19619457_chr',CHR,'.bed',sep=''))
	GenoNum=2-as(Geno$genotype,'numeric') 
	load(file=paste(EVO_IMMUNO_POP,'/Maxime/SNP_annotation/imputed_Genotypes/MapFiles/Map_imputation_200x19619457_chr',CHR,'.Rdata',sep=''))
	rm(Geno);gc()
	chr.dist=5e7
	By=function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x} # by own version of by (makes a vector with names)
	Chrsep=cumsum(chr.dist+ c(-chr.dist,By(Map$position,Map$chromosome,max)))-chr.dist/2
}
	
cat('ok')
Map$maf_AFB=pmin(Map$SNPfreq_AF,1-Map$SNPfreq_AF)
Map$maf_EUB=pmin(Map$SNPfreq_EU,1-Map$SNPfreq_EU)

Require('Map_imputed')
Require('GeneAnnot')

# read GeneAnnot
#GeneAnnot=read.table('/pasteur/homes/mrotival/Annotation/GeneAnnotation_hg37_ens70.txt',sep='\t',comment='',quote='',stringsAsFactors=FALSE)
#GeneAnnot=GeneAnnot[which(GeneAnnot$select),]
#GeneAnnot=GeneAnnot[match(rownames(FPKM_gene),GeneAnnot$Ensembl.Gene.ID),]
#GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='X']=23
#GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='Y']=24
#GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='MT']=26
#GeneAnnot$Chromosome.Name=as.numeric(GeneAnnot$Chromosome.Name)

# read TranscriptAnnot
#TranscriptAnnot=read.table('/pasteur/homes/mrotival/Annotation/TranscriptAnnotation_hg37_ens70.txt',sep='\t',comment='',quote='',stringsAsFactors=FALSE)
#TranscriptAnnot=TranscriptAnnot[match(rownames(FPKM),TranscriptAnnot$Ensembl.Transcript.ID),]
#TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='X']=23
#TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='Y']=24
#TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='MT']=26
#TranscriptAnnot$Chromosome.Name=as.numeric(TranscriptAnnot$Chromosome.Name)
cat('ok')

# read FPKM
#TranscriptAnnot=TranscriptAnnot[which(isoAnnot$mF1),]
thTrans=0
thCis=0.05
#	for (pop in c('AFB','EUB','ALL')){
pop = 'ALL'

########################################################################
##					PSI	 by Pop - Robust - Cis  		 		 	  ##
########################################################################

if(PSI_STANDARD){
	load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V5_withCounts.Rdata',sep=''))   
	# transcript_Annot
	load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME))	 
	PSI=PSI_prov
	PSI_Annot=PSI_Annot[toKeep,]
	rownames(PSI)=PSI_Annot$event_id
	
	PSI_Pos=data.frame(event=PSI_Annot$event_id,chr=PSI_Annot$chrom,start=PSI_Annot$start,end=PSI_Annot$end)
	rownames(PSI_Pos)=PSI_Annot$event_id
	
	RESCIS=NULL
	RESTRANS=NULL
	#qqbin=seq(0,308.70,0.05)
	tim=Sys.time()
		if(pop=='AFB'){
			wFreq=which(Map$maf_AFB>=minFreq)
			if(TEST){wFreq=wFreq[1:1000]}
			PCs=read.table(paste('/pasteur/homes/mrotival/02_ERC_Data/SampleCaseReports/EvoImmunoPop_',pop,'_ancestryPCA_OmniExome.txt',sep=''),header=T)
			rownames(PCs)=PCs$IID
			wCondPop=which(SampleAnnot$condition==cond & SampleAnnot$population==pop)
			}
		if(pop=='EUB'){
			wFreq=which(Map$maf_EUB>=minFreq)#[1:10000]
			if(TEST){wFreq=wFreq[1:1000]}
			PCs=read.table(paste('/pasteur/homes/mrotival/02_ERC_Data/SampleCaseReports/EvoImmunoPop_',pop,'_ancestryPCA_OmniExome.txt',sep=''),header=T)
			rownames(PCs)=PCs$IID
			wCondPop=which(SampleAnnot$condition==cond & SampleAnnot$population==pop)
		} 
		if(pop=='ALL'){
			wFreq=which(Map$maf_EUB >=minFreq | Map$maf_AFB>=minFreq)#[1:10000]
#			wFreq=which(Map$SNPfreq >=minFreq)
			if(TEST){wFreq=wFreq[1:1000]}
			PCs=read.table(paste('/pasteur/homes/mrotival/02_ERC_Data/SampleCaseReports/EvoImmunoPop_',pop,'_ancestryPCA.txt',sep=''),header=T)
			rownames(PCs)=PCs$IID
			wCondPop=which(SampleAnnot$condition==cond)
		}

		if(perm>0){
			set.seed(perm)
			indiv_0=SampleAnnot[wCondPop,'individu']
			mypop=matrix(ifelse(substr(indiv_0,1,3)=='AFB',1,0),length(indiv_0))
			indiv=indiv_0
			indiv[mypop==0]=sample(indiv[mypop==0])
			indiv[mypop==1]=sample(indiv[mypop==1])
		}else{
			indiv=SampleAnnot[wCondPop,'individu']
			indiv_0=SampleAnnot[wCondPop,'individu']
		}
		GenoNum_Mat = SlicedData$new()
		GenoNum_Mat$CreateFromMatrix( t(GenoNum)[wFreq,match(indiv,rownames(GenoNum))])
		GenoNum_Mat$ResliceCombined(sliceSize = 5000)

		FPKM_Mat = SlicedData$new()
		if(ROBUST){
			FPKM_Mat$CreateFromMatrix(t(apply(PSI[,wCondPop],1,function(x){qnorm(rank(x,ties='random')/(length(x)+1),mean(x),sd(x))})))
		}else{
			FPKM_Mat$CreateFromMatrix(PSI[,wCondPop])
		}
		show(FPKM_Mat)
		Cvrt=SlicedData$new()
		Cvrt$CreateFromMatrix(t(PCs[match(indiv,PCs$IID),3:4]))
		res=Matrix_eQTL_main(GenoNum_Mat,FPKM_Mat,Cvrt,output_file_name=NULL, pvOutputThreshold=thTrans,output_file_name.cis=NULL,pvOutputThreshold.cis=thCis,snpspos=Map[wFreq,c(2,1,4)],genepos=PSI_Pos,pvalue.hist='qqplot')
		res$cis$eqtls$condition=rep(cond,nrow(res$cis$eqtls))
		res$cis$eqtls$population=rep(pop,nrow(res$cis$eqtls))
				
		# save qqplot data
		RESCIS=res$cis$eqtls
	print(Sys.time()-tim)
	RESCIS$snps=as.character(RESCIS$snps)
	RESCIS$gene=as.character(RESCIS$gene)
	RESCIS$type=rep('psiQTL_miso',nrow(RESCIS))
	save(RESCIS,file=paste(HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,ifelse(ROBUST,'_KW','_STD'),ifelse(LOGIT,'_LOGIT',''),'_misoPSI_',Source,'_V7.Rdata',sep=''))
	cat('printed ',HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,ifelse(ROBUST,'_KW','_STD'),ifelse(LOGIT,'_LOGIT',''),'_misoPSI_',Source,'_V7.Rdata',sep='')
	
	#load(file='/pasteur/projets/evo_immuno_pop/Maxime/transEQTL_Clust/CisTransObserved.Rdata')
}

########################################################################
##			Fold Change PSI by Pop - Standard - Cis		 		 	  ##
########################################################################

if(PSI_RESPONSE & cond>1){
	load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V5_withCounts.Rdata',sep=''))  
#	toKeep=which(PSI_Annot[[paste('PctNA',condIndex[1],sep='_')]]<0.8 & PSI_Annot[[paste('PctNA',condIndex[cond],sep='_')]]<0.8 & PSI_Annot$JuncCovered & PSI_Annot$noOverlap & PSI_Annot$NbTestableCond>0)
#	PSI=PSI[toKeep,]
	load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME)) 
	if(LOGIT){
		PSI=logit(PSI_prov)
	}else{
		PSI=PSI_prov
	}
	PSI_Annot=PSI_Annot[toKeep,]
	rownames(PSI)=PSI_Annot$event_id
	
	PSI_Pos=data.frame(event=PSI_Annot$event_id,chr=PSI_Annot$chrom,start=PSI_Annot$start,end=PSI_Annot$end)
	rownames(PSI_Pos)=PSI_Annot$event_id
	
	RESCIS=NULL
	RESTRANS=NULL
	#qqbin=seq(0,308.70,0.05)
	tim=Sys.time()
		wFreq=which(Map$maf_EUB >=minFreq | Map$maf_AFB>=minFreq)#[1:10000]
		if(TEST){wFreq=wFreq[1:2000]}
		PCs=read.table(paste('/pasteur/homes/mrotival/02_ERC_Data/SampleCaseReports/EvoImmunoPop_',pop,'_ancestryPCA.txt',sep=''),header=T)
		wCondPop=which(SampleAnnot$condition==cond)
		w1Pop=which(SampleAnnot$condition==1)
		indiv=intersect(SampleAnnot[wCondPop,'individu'],SampleAnnot[w1Pop,'individu'])
		wCondPop=wCondPop[match(indiv,SampleAnnot[wCondPop,'individu'])]
		w1Pop=w1Pop[match(indiv,SampleAnnot[w1Pop,'individu'])]
		if(perm>0){
			set.seed(perm)
			indiv_0=SampleAnnot[wCondPop,'individu']
			mypop=matrix(ifelse(substr(indiv_0,1,3)=='AFB',1,0),length(indiv_0))
			indiv=indiv_0
			indiv[mypop==0]=sample(indiv[mypop==0])
			indiv[mypop==1]=sample(indiv[mypop==1])
		}else{
			indiv_0=SampleAnnot[wCondPop,'individu']
		}		
		GenoNum_Mat = SlicedData$new()
		GenoNum_Mat$CreateFromMatrix(t(GenoNum)[wFreq,match(indiv,rownames(GenoNum))])
		GenoNum_Mat$ResliceCombined(sliceSize = 5000)
		FPKM_Mat = SlicedData$new()
		if(ROBUST){
			FPKM_Mat$CreateFromMatrix(t(apply(PSI[,wCondPop]-PSI[,w1Pop],1,function(x){qnorm(rank(x,ties='random')/(length(x)+1),mean(x),sd(x))})))
		}else{
			FPKM_Mat$CreateFromMatrix(PSI[,wCondPop]-PSI[,w1Pop])
			}
		Cvrt=SlicedData$new()
		Cvrt$CreateFromMatrix(t(PCs[match(indiv,PCs$IID),3:4]))
		res=Matrix_eQTL_main(GenoNum_Mat,FPKM_Mat,Cvrt,output_file_name=NULL, pvOutputThreshold=thTrans,output_file_name.cis=NULL,pvOutputThreshold.cis=thCis,snpspos=Map[wFreq,c(2,1,4)],genepos=PSI_Pos,pvalue.hist='qqplot')
		res$cis$eqtls$condition=rep(cond,nrow(res$cis$eqtls))
		res$cis$eqtls$population=rep(pop,nrow(res$cis$eqtls))
		# save qqplot data
		RESCIS=res$cis$eqtls
	print(Sys.time()-tim)

	RESCIS$snps=as.character(RESCIS$snps)
	RESCIS$gene=as.character(RESCIS$gene)
	RESCIS$type=rep('resp-psiQTL_miso',nrow(RESCIS))

	save(RESCIS,file=paste(HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,ifelse(ROBUST,'_KW','_STD'),ifelse(LOGIT,'_LOGIT',''),'_misoPSI_response_',Source,'_V7.Rdata',sep=''))
	cat('printed ',HOME,'/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/MatrixEQTL/Perm',perm,'/Cis-sQTL_',pop,'_',cond,'_chr',CHR,ifelse(ROBUST,'_KW','_STD'),ifelse(LOGIT,'_LOGIT',''),'_misoPSI_response_',Source,'_V7.Rdata',sep='')
	}
#}
q('no')