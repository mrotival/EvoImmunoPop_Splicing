####################################
###	  extract genotypes of sQTL  ###
####################################

getGenos=function(snps){
    # extract the list of target SNPs 
    require(snpStats)
    Map_allChr=list()
    Geno_allChr=list()
    for( CHR in 1:22){
    	Geno=read.plink(paste(EVO_IMMUNO_POP,'/DATA_FREEZE/ERC_Main_Genotyping_24022015/Imputation/EvoImmunoPop_imputation_200x19619457_chr',CHR,'.bed',sep=''))
	    GenoNum=as(Geno$genotype,'numeric')
	    Map=Geno$map
	    mm=Map[match(snps,Map$snp.name),]
	    Map_allChr[[CHR]]=Map[mm,,drop=F]
	    Geno_allChr[[CHR]]=2-GenoNum[mm,,drop=F]
	    }
	Map_allChr=rbindlist(Map_allChr)
	Map_allChr$locus=paste(Map_allChr$chromosome,':',round(Map_allChr$position/1e6,0),'Mb',sep='')
	Map_allChr$posID=paste(Map_allChr$chromosome,Map_allChr$position,sep=':')
	Geno_allChr=rbindlist(Geno_allChr)
	Geno=cbind(Map_allChr[c('snp.name','chromosome','position','locus','posID')],Geno_allChr)
	Geno
}

getGenos_locus=function(chrom,start,end){
    # extract all SNPs in a window
    require(snpStats)
    Map_allChr=list()
    Geno_allChr=list()
    CHR=chrom
    Geno=read.plink(paste(EVO_IMMUNO_POP,'/DATA_FREEZE/ERC_Main_Genotyping_24022015/Imputation/EvoImmunoPop_imputation_200x19619457_chr',CHR,'.bed',sep=''))
    GenoNum=as(Geno$genotype,'numeric')
    Map=Geno$map
    mm=which(Map$position>=start & Map$position<=end) 
    Map=Map[mm,,drop=F]
    Geno=2-GenoNum[mm,,drop=F]
	Map$locus=paste(Map$chromosome,':',round(Map$position/1e6,0),'Mb',sep='')
	Map$posID=paste(Map$chromosome,Map$position,sep=':')
	Geno=cbind(Map[c('snp.name','chromosome','position','locus','posID')],Geno)
	Geno
}
