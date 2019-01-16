##################################################################################################################################################################################################
########    LD files Imputed_LD_list_80_${POP}_${CHR}.ld' are obtained using plink with the following command                                                                             ########
########    plink --bfile Genotype_${CHR} --keep {POP}.txt --maf 0.05 --chr ${CHR} --R2 --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0.8 --out Imputed_LD_list_80_${POP}_${CHR}  ########
##################################################################################################################################################################################################


#############################################################################
############### computing length of archaic haplotype in aSNPs ##############
#############################################################################

compute_archaic_haplotype_length=function(Map,CHR){
    #returns Annotated Map
    pop='EUB'
    tim=Sys.time()
    RR_CEU=as.data.frame(fread(sprintf('%s/Annotation/RecombinationRate/CEU/CEU-%s-final.txt',HOME,CHR)))
    colnames(RR_CEU)=make.names( colnames(RR_CEU))

    cat('loading Map EUB',CHR,'\n')
 	Map=Map[,c('snp.name','chromosome','position','SNPfreq','maf_EUB','maf_AFB','allele.1','allele.2','ancestral_allele','aSNP')]

    print(Sys.time()-tim)
    cat('haploLength_cM_CEU\n')
    Map$haploLength_cM_CEU=NA
    approx_cM=approx(RR_CEU$Position.bp., y =RR_CEU$Map.cM., Map$position)$y
    
    haploLength=By(c(approx_cM,approx_cM[match(LD_table$SNP_A,Map$snp.name)],approx_cM[match(LD_table$SNP_B,Map$snp.name)]),c(Map$snp.name,LD_table$SNP_B,LD_table$SNP_A),function(x){diff(range(x,na.rm=T))})
    wSNP=match(Map$snp.name,names(haploLength))
    Map$haploLength_cM_CEU=haploLength[wSNP]
    Map$approx_cM=approx_cM
    print(Sys.time()-tim)

    cat('haploLength_EUB_archaic\n')
    Map=as.data.table(Map)
#   aggregates data across all SNPs in LD
    Map_LD=rbind(cbind(Map[,mget(c('aSNP','position','approx_cM','snp.name','maf_EUB'))],matchSNP=Map$snp.name),
                 cbind(Map[match(LD_table$SNP_B,Map$snp.name),mget(c('aSNP','position','approx_cM','snp.name'))], matchSNP=LD_table$SNP_A),
                 cbind(Map[match(LD_table$SNP_A,Map$snp.name),mget(c('aSNP','position','approx_cM','snp.name'))],matchSNP=LD_table$SNP_B))
    haploLength_archaic=Map_LD[aSNP==1,.(Nb_archaic=length(unique(snp.name)),
                                         haploLength_archaic=diff(range(position)),
                                         haploLength_archaic_cM=diff(range(approx_cM))
                                         best_aSNP=snp.name[which.max(maf_EUB)]),by=matchSNP]

    wSNP=match(haploLength_archaic$matchSNP,Map$snp.name)
#   Add Nb of archaic SNP on haplotype
    Map$Nb_archaic=0
    Map$Nb_archaic[wSNP]=haploLength_archaic[,get('Nb_archaic')]
#   Annotate SNPs linked to 1+ aSNP
    Map$aSNP_R2_EUB=ifelse(Map$Nb_archaic>0,TRUE,FALSE)
#   Add archaic length in Bp
    Map$haploLength_archaic=0
    Map$haploLength_archaic[wSNP]=haploLength_archaic[,get('haploLength_archaic')]
#   Add archaic length in CM
    Map$haploLength_archaic_cM=0
    Map$haploLength_archaic_cM[wSNP]=haploLength_archaic[,get('haploLength_archaic_cM')]
    Map$best_aSNP=''
    Map$best_aSNP[wSNP]=haploLength_archaic[,get('best_aSNP')]
#   Add locus informations
    Map$locus=paste(Map$chromosome,':',round(Map$position/1e6),'Mb',sep='')
    print(Sys.time()-tim)
    cat('\n')
    return(Map)
}


#############################################################################
############### adding informations relative to GWAS           ##############
#############################################################################

library(gwascat)
library(rtracklayer)
library(data.table)

############## loading  GWAS catalog  ######################### 

load_gwas_catalog=function(){
    ## load GWAS catalog 
    GWAScat = read.delim(sprintf("%s/Annotation/GWAS/gwas_catalog_v1.0.1-associations_e89_r2017-06-26_downloaded05072017_current.tsv",HOME), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    GWAScat = gwascat:::fixNonASCII(GWAScat)
    GWAScat$CHR_POS=as.numeric(GWAScat$CHR_POS)
    GWAScat=GWAScat[!is.na(GWAScat$CHR_POS),]
    GWAScat=makeGRangesFromDataFrame(GWAScat,keep.extra.columns=TRUE,seqnames.field='CHR_ID',start.field='CHR_POS',end.field='CHR_POS')
    ## liftover GWAS catalog to hg19
    ch=import.chain("/Volumes/@home/Annotation/liftoverFiles/hg38ToHg19.over.chain")
    seqlevelsStyle(GWAScat) = "UCSC"  # necessary
    cur19 = liftOver(GWAScat, ch)
    GWAScat=as.data.frame(cur19)
    GWAScat$seqnames=gsub('chr','',as.character(GWAScat$seqnames))
    GWAScat$strand=as.character(GWAScat$strand)

    ## annotate GWAS traits
    EFO_map=as.data.frame(fread('data/gwas_catalog_trait-mappings_r2017-07-04_manuallyCurated.tsv') # EFO corrected version
    # for trait that have multiple parent EFO, the original file provided by GWAS catalog reports only one parent EFO category. 
    # The EFO map file was manually curated to report all parent EFOs of traits with multiple parent category.
    EFO_map_0=EFO_map
    colnames(EFO_map_0)=gsub(' ','_',colnames(EFO_map_0))
    GWAScat_0=GWAScat[,c('MAPPED_TRAIT_URI','MAPPED_TRAIT')]
    GWAScat_0=GWAScat_0[!duplicated(GWAScat_0),]
    EFO_map_1=merge(EFO_map_0,as.data.frame(GWAScat_0),by.x='EFO_URI',by.y='MAPPED_TRAIT_URI')
    EFO_map_1[EFO_map_1[,'Parent_term']=='Immune system','Parent_term']='Immune system disorder'
    # concatenate parent EFO
    By = function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x}
    TRAIT2EFO_1=By(EFO_map_1[,'Parent_term'],EFO_map_1$MAPPED_TRAIT,function(x){paste(sort(unique(x)),collapse='//')})
    GWAScat$PARENT.EFO=TRAIT2EFO_1[match(GWAScat$MAPPED_TRAIT,names(TRAIT2EFO_1))]
    GWAScat
}

############## annotating SNPs from GWAS catalog and linked variants ######################### 

annotate_GWAS_hits_with_EFO=function(Map,CHR){

    GWAScat = load_gwas_catalog()
    
    GWAScat=GWAScat[which(GWAScat$seqnames==CHR),]
    GWAScat$posID=paste(GWAScat$seqnames,GWAScat$start,sep=':')
    # retrieve rsids based on position
    GWAScat$snp.name=Map[match(GWAScat$posID,Map$posID),'snp.name']
    # if it fails, check if reported SNP is present in our data
    wNA=which(is.na(GWAScat$snp.name) & GWAScat$SNPS %in% Map$snp.name)
    GWAScat$snp.name[wNA]=GWAScat$SNPS[wNA]
    # if it fails, check if reported SNP_ID match a known rs
    wNA=which(is.na(GWAScat$snp.name) & paste('rs',GWAScat$SNP_ID_CURRENT,sep='') %in% Map$snp.name)
    GWAScat$snp.name[wNA]=paste('rs',GWAScat$SNP_ID_CURRENT[wNA],sep='')	


    ## annotate GWAS SNPs (reported)
    By = function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x}
    # Traits at 10^-5
    GWASSNP_1E5=By(GWAScat$DISEASE.TRAIT,GWAScat$snp.name,function(x){paste(unique(x),collapse='//')})
    Map$GWAS_Trait_1E5=''
    Map$GWAS_Trait_1E5[match(names(GWASSNP_1E5),Map$snp.name)]=GWASSNP_1E5
    # Traits at 10^-8
    GWASSNP_1E8=By(GWAScat$DISEASE.TRAIT[GWAScat$P.VALUE<1e-8],GWAScat$snp.name[GWAScat$P.VALUE<1e-8],function(x){paste(unique(x),collapse='//')})
    Map$GWAS_Trait_1E8=''
    Map$GWAS_Trait_1E8[match(names(GWASSNP_1E8),Map$snp.name)]=GWASSNP_1E8
    # parent EFO at 10^-5
    GWASSNP_EFO_1E5=By(GWAScat$PARENT.EFO,GWAScat$snp.name,function(x){paste(unique(unlist(strsplit(x,'//'))),collapse='//')})
    Map$GWAS_EFO_1E5=''
    Map$GWAS_EFO_1E5[match(names(GWASSNP_EFO_1E5),Map$snp.name)]=GWASSNP_EFO_1E5
    # parent EFO at 10^-8
    GWASSNP_EFO_1E8=By(GWAScat$PARENT.EFO[GWAScat$P.VALUE<1e-8],GWAScat$snp.name[GWAScat$P.VALUE<1e-8],function(x){paste(unique(unlist(strsplit(x,'//'))),collapse='//')})
    Map$GWAS_EFO_1E8=''
    Map$GWAS_EFO_1E8[match(names(GWASSNP_EFO_1E8),Map$snp.name)]=GWASSNP_EFO_1E8

    #### extend these annotation to SNPs in LD with a GWAS SNP in EUB
    combine_traits=function(x){
       # split traits 
       split_traits=unlist(strsplit(x,'//'))
       # remove duplicates and empty fields
       split_traits_clean = setdiff(unique(split_traits),'')
       # collapse traits
       combined_traits = paste(split_traits_clean,collapse='//')
       combined_traits
       }

    Map=as.data.table(Map)
    #annotate pairs of linked SNPs
    Map_LD=rbind(cbind(Map[,mget(c('GWAS_Trait_1E8','GWAS_Trait_1E5','position','snp.name'))],matchSNP=Map$snp.name),
                 cbind(Map[match(LD_table$SNP_B,Map$snp.name),mget(c('aSNP','position','snp.name'))], matchSNP=LD_table$SNP_A),
                 cbind(Map[match(LD_table$SNP_A,Map$snp.name),mget(c('aSNP','position','snp.name'))],matchSNP=LD_table$SNP_B))
    #collapse annotations by SNP
    GWAS_R2=Map_LD[,.(GWAS_Trait_R2_1E5 = combine_traits(GWAS_Trait_1E5),
                     GWAS_Trait_R2_1E8 = combine_traits(GWAS_Trait_1E8),
                     GWAS_EFO_R2_1E5 = combine_traits(GWAS_EFO_1E5),
                     GWAS_EFO_R2_1E8 = combine_traits(GWAS_EFO_1E8)),by=matchSNP]

    ## annotate GWAS SNPs (reported SNPs and linked variants)
    wSNP=match(GWAS_R2$matchSNP,Map$snp.name)
    # Traits at 10^-5
    Map$GWAS_Trait_R2_1E5=''
    Map$GWAS_Trait_R2_1E5[wSNP]=GWAS_R2[,get('GWAS_Trait_R2_1E5')]
    # Traits at 10^-8
    Map$GWAS_Trait_R2_1E8=''
    Map$GWAS_Trait_R2_1E8[wSNP]=GWAS_R2[,get('GWAS_Trait_R2_1E8')]
    # parent EFO at 10^-5
    Map$GWAS_EFO_R2_1E5=''
    Map$GWAS_EFO_R2_1E5[wSNP]=GWAS_R2[,get('GWAS_EFO_R2_1E5')]
    # parent EFO at 10^-8
    Map$GWAS_EFO_R2_1E8=''
    Map$GWAS_EFO_R2_1E8[wSNP]=GWAS_R2[,get('GWAS_EFO_R2_1E8')]

    Map
    }

#############################################################################
###############         computing selection statistics         ##############
#############################################################################

Map_imputed$iHS_EUB[Map_imputed$maf_EUB<0.05]=NA
Map_imputed$iHS_AFB[Map_imputed$maf_AFB<0.05]=NA

wFreq=which(Map_imputed$SNPfreq > 0.05)
q95_aIHS_EUB=quantile(abs(Map_imputed$iHS_EUB[wFreq]),0.95,na.rm=T);q95_aIHS_EUB
q95_aIHS_AFB=quantile(abs(Map_imputed$iHS_AFB[wFreq]),0.95,na.rm=T);q95_aIHS_AFB
q95_FST=quantile(Map_imputed$FST_adj[wFreq],0.95,na.rm=T);q95_FST

annotate_selection=function(Map,CHR){
    GR_Map=makeGRangesFromDataFrame(Map,seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
    myDist=5e4 # 100Kb widow around the SNP
    GR_Window=flank(flank(GR_Map,myDist),-2*myDist-1)
    oo=findOverlaps(GR_Window,GR_Map)
    Map$NbSNP_100kb=By(!is.na(Map$FST_adj)[subjectHits(oo)],queryHits(oo),sum,na.rm=T)
    Map$NbSNP_out95_FST_100kb=By(Map$FST_adj[subjectHits(oo)]>q95_FST,queryHits(oo),sum,na.rm=T)
    Map$NbSNP_out95_aIHS_EUB_100kb=By(abs(Map$iHS_EUB[subjectHits(oo)])>q95_aIHS_EUB,queryHits(oo),sum,na.rm=T)
    Map$NbSNP_out95_aIHS_AFB_100kb=By(abs(Map$iHS_AFB[subjectHits(oo)])>q95_aIHS_AFB,queryHits(oo),sum,na.rm=T)
    Map
}

#############################################################################
############### run annotation script across all chromosomes   ##############
#############################################################################

Map_imputed=as.data.frame(fread('data/Map_imputed_essential_informations.txt'))
Map_imputed=Map_imputed[which(Map_imputed$SNPfreq>=0.05),]

Map_archaic=list()
Map_GWAS=list()
Map_select=list()

for (CHR in 1:22){
    Map=Map_imputed[which(Map_imputed$chrom==CHR),]
    LD_table=as.data.frame(fread(paste('data/LD/Imputed_LD_list_80_',pop,'_',CHR,'.ld',sep='')))
    LD_table$maf_1=Map$maf_EUB[match(LD_table$SNP_A,Map$snp.name)]
    LD_table$maf_2=Map$maf_EUB[match(LD_table$SNP_B,Map$snp.name)]
    LD_table=LD_table[which(!is.na(LD_table$maf_2 & LD_table$maf_1)),]

    Map_archaic[[CHR]]=compute_archaic_haplotype_length(Map,CHR)

    Map_GWAS[[CHR]]=annotate_GWAS_hits_with_EFO(Map,CHR)

    Map_select[[CHR]]=annotate_selection(Map, CHR)
}

Map_GWAS = rbindlist(Map_GWAS)

Map_select = rbindlist(Map_select)

Map_archaic = rbindlist(Map_archaic)

write.table(Map_archaic, file ='data/Map_archaic.txt',quote=F,sep='\t',row.names=F)
write.table(Map_select, file ='data/Map_select.txt',quote=F,sep='\t',row.names=F)
write.table(Map_GWAS, file ='data/Map_GWAS.txt',quote=F,sep='\t',row.names=F)


