#############################################################################
#########      Annotate overlap with regulatory elements      ##############
#############################################################################

####### add informations to the SNP annoation table (Map_imputed)
# load Map of genetic variants
Map_imputed=as.data.frame(fread('/Volumes/evo_immuno_pop/Maxime/SNP_annotation/Map_imputed_essential_informations.txt'))
Map_imputed=Map_imputed[which(Map_imputed$SNPfreq>=0.05),]

###### Quantify effect of SNPs on splicing (Data from Xiong et al.)
SPIDEX=as.data.frame(fread(sprintf("%s/Annotation/SPIDEX/SNPlist_allSPIDEX.txt",HOME)))

###### annotate SNPs that overlap a splice site or splice site flank
SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.4.txt',HOME))
wIntroExpressed=!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand) & SpliceSites$maxFPKM_overlappingGene_withStrand>1
wGerp=!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align
SpliceSites=SpliceSites[which(wIntroExpressed & wGerp),]
## create GRanges containing splice regulatory elements 
Splice_GR=makeGRangesFromDataFrame(SpliceSites,seqnames='chrom_intron',start.field='start_site',end.field='end_site',strand='strand_correct')
SpliceFlank_GR=c(flank(Splice_GR[SpliceSites$isStart,],start=FALSE,8,ignore.strand=TRUE),flank(Splice_GR[SpliceSites$isEnd,],start=TRUE,8,ignore.strand=TRUE))
acceptor_GR=Splice_GR[SpliceSites$type=='acceptor',]
acceptor_GR$WeakAltConst=SpliceSites$WeakAltConst[which(SpliceSites$type=="acceptor")]
branchpoint_expected=flank(GenomicRanges::shift(acceptor_GR,ifelse(strand(acceptor_GR)=='-',16,-16)),21,start=TRUE)

######## Add information from SPIDEX database
Map_imputed$SPIDEX_Z=SPIDEX$Z_SPIDEX[match(Map_imputed$snp.name,SPIDEX$snp.name)]
Map_imputed$Exonic=SPIDEX$Exon_zone[match(Map_imputed$snp.name,SPIDEX$snp.name)]
Map_imputed$groupSplice=ifelse(is.na(Map_imputed$SPIDEX_Z),'',ifelse(abs(Map_imputed$SPIDEX_Z)>2,paste(Map_imputed$Exonic,'deleterious'),Map_imputed$Exonic))

######## Add information splice regulatory elements 
# splice sites
Map_imputed$SpliceSite=''
oo=findOverlaps(Map_GR,Splice_GR)
Map_imputed$SpliceSite[queryHits(oo)]=paste(SpliceSites$type[subjectHits(oo)],'site')

# splice site flanks 
Map_imputed$SpliceSiteFlank=''
oo=findOverlaps(Map_GR,SpliceFlank_GR)
Map_imputed$SpliceSiteFlank[queryHits(oo)]=paste(SpliceSites$type[subjectHits(oo)],'flank')

# branchpoints
oo=findOverlaps(Map_GR,branchpoint_expected)
Map_imputed$Branchpoint=''
Map_imputed$Branchpoint[queryHits(oo)]='Branchpoint'

Map_imputed$SpliceElt=ifelse(grepl('acceptor',Map_imputed$SpliceSite),'acceptor site',
                            ifelse(grepl('donor',Map_imputed$SpliceSite),'donor site',
                            ifelse(grepl('acceptor',Map_imputed$SpliceSiteFlank),'acceptor flank',
                            ifelse(grepl('donor',Map_imputed$SpliceSiteFlank),'donor flank',
                            ifelse(Map_imputed$Branchpoint!='','branchpoint','')))))

rm(SPIDEX, SpliceSites,Splice_GR,SpliceFlank_GR,acceptor_GR);gc()

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ####
##### ##### ##### ##### ##### ##### RUN HOMER ON SNPS TO IDENTIFY RBP  ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ####

size=10
snps_pos=Map_imputed[,c('snp.name','chromosome','position')]
snps_pos$start=snps_pos$position-size
snps_pos$end=snps_pos$position+size
snps_neg=snps_pos
snps_neg$strand='-'
snps_dat=rbind(snps_pos,snps_neg)
snps_dat$chromosome=paste('chr',snps_dat$chromosome,sep='')
snps_dat$snp.name=paste(snps_dat$snp.name,snps_dat$strand,sep='_')
snps_dat$position=NULL
write.table(snps_dat,file='data/snps_data_10bp_1e6snps_both_strands.txt',quote=F,row.names=F,col.names=F,sep='\t')

######	run Homer on all SNPs to discover motifs
# assumes that HOMER is already in the PATH
SNPfile='data/snps_data_10bp_1e6snps_both_strands.txt'
MotifOutputDirectory='data/'

MOTIFS='data/human.rna_o5.motifs' # this file is a modified version of the file provided by HOMER with thresholds set to 5 for all motifs

system('changeNewLine.pl data/snps_data_10bp_1e6snps_both_strands.txt')
cmd = paste('findMotifsGenome.pl %s hg19 %s -find %s -mask -rna -size 20 > %s/Motifs_allSNPs_ERC.txt', SNPfile, MotifOutputDirectory, MOTIFS, MotifOutputDirectory)
system(cmd)

#################################################
######## reading MOTIF information in R	#########
#################################################

library(data.table)
library(TFBSTools)
library(MotIV)

RBP_list=fread('grep -e ">" /Users/maxime/bin/HOMER/data/knownTFs/human.rna.motifs')
SpliceFactor=sapply(strsplit(as.data.frame(RBP_list)$V2,'/'),function(x){x[1]})
TargetMotifName=gsub('>' ,'',as.data.frame(RBP_list)$V1)

motifs <- read.table(file = '/Users/maxime/bin/HOMER/data/knownTFs/human.rna.motifs', fill = T,row.names=NULL, header=F)
table.vector <- cbind(as.vector(motifs[[1]]), as.vector(motifs[[2]]), as.vector(motifs[[3]]), as.vector(motifs[[4]]))
DE <- c(which(grepl('>',table.vector[,1])),nrow(table.vector)+1)
names <- gsub('>' ,'',table.vector[DE[-length(DE)], 1])
listPWM <- list()
for (i in seq(length(DE)-1)) {
   listPWM[[i]] <- matrix(as.numeric(table.vector[(DE[i] + 1):(DE[i+1] - 1), 1:4]), nrow = 4, byrow = T, dimnames = list(c("A","C", "G", "T")))
   colnames(listPWM[[i]]) <- 1:(length(listPWM[[i]])/4)
   names(listPWM)[i] <- table.vector[DE[i], 1]
}
targetSeq=sapply(listPWM,function(x){paste(apply(x,2,function(y){c('A','C','G','U')[which.max(y)]}),collapse='')})
Entropy=sapply(listPWM,function(x){-sum(x*log(x))})
RBP_list=data.frame(TargetMotifName,SpliceFactor,targetSeq,Entropy)
rownames(RBP_list)=1:nrow(RBP_list)
# plot.new()
# MotIV:::seqLogo2(listPWM[[1]]) # lowest entropy
# plot.new()
# MotIV:::seqLogo2(listPWM[[3]]) # highest entropy

####### find binding sites in tested sequence
RBP_motifs=fread('/Users/maxime/bin/HOMER/Z_results/test/Motifs_allSNPs_ERC.txt')
colnames(RBP_motifs)=c('PositionID','Offset','Sequence','Motif Name','Strand','MotifScore')

tab=table(gsub('.*/Homo_sapiens-(.*)-PBM/HughesRNA','\\1',RBP_motifs$V4))
barplot(tab,las=2)
RBP_motifs=RBP_motifs[order(RBP_motifs$PositionID, RBP_motifs$'Motif Name', -RBP_motifs$'MotifScore'),]
RBP_motifs$snp_id=sapply(strsplit(RBP_motifs$PositionID,'_'),function(x){x[1]})
RBP_motifs$Motif_id=gsub('.*/Homo_sapiens-(.*)-PBM/HughesRNA','\\1',RBP_motifs$'Motif Name')
RBP_motifs$SpliceFactor=RBP_list[match(RBP_motifs$Motif_id,RBP_list$TargetMotifName),'SpliceFactor']
NbHitsByMotif=sapply(unique(RBP_motifs$Motif_id),function(id){x=RBP_motifs$MotifScore[RBP_motifs$Motif_id==id];c(NbMotif=length(x),NbMotif_over6=sum(x>6),NbMotif_over7=sum(x>7),NbMotif_over8=sum(x>8))})

#annotate the number of binding sites per motif
RBP_list=cbind(RBP_list,t(NbHitsByMotif)[match(RBP_list$TargetMotifName,colnames(NbHitsByMotif)),])
colnames(RBP_list)=c('TargetMotifName','SpliceFactor','targetSeq','Entropy','NbMotif','NbMotif_over6','NbMotif_over7','NbMotif_over8')
RBP_list$SpliceFactor_short=sapply(strsplit(RBP_list$SpliceFactor,'(',fixed=T),function(x){x[1]})

#########################################################
######## reading SNP overlapping MOTIFs in R	#########
#########################################################

Map_RBP_motifs=Map_imputed[,c('snp.name','chromosome','position','allele.1','allele.2','ancestral_allele','SNPfreq','daf_char_AFB','daf_char_EUB')]

dup=which(duplicated(paste(RBP_motifs$snp_id,RBP_motifs$Motif_id)))
Map_RBP_motifs$RBP_motif=rep('',nrow(Map_RBP_motifs))

j=0
for(i in RBP_list$TargetMotifName){
	j=j+1
	cat(j,':',i,'\n')
	wMot=which(RBP_motifs$Motif_id==i)
	wMot0=setdiff(wMot,dup)
	wMatch=match(RBP_motifs$snp_id[wMot0],Map_RBP_motifs$snp.name)
	mysep=ifelse(j==1,'','/')
	Map_RBP_motifs$RBP_motif[wMatch]=paste(Map_RBP_motifs$RBP_motif[wMatch],i,sep=mysep)
}
#write.table(Map_RBP_motifs,file='/Users/maxime/bin/HOMER/Z_results/Map_RBP_motifs.txt',sep='\t',quote=F,row.names=F,col.names=F)

######################################################################################
###### annotate SNPs that overlap a motif bound by an RBP or splice factor ###########
######################################################################################

Map_RBP_motifs=as.data.frame(fread('/Users/maxime/bin/HOMER/Z_results/Map_RBP_motifs.txt'))
colnames(Map_RBP_motifs)=c('snp.name','chromosome','allele.1','allele.2','ancestral_allele','SNPfreq','daf_char_AFB','daf_char_EUB','RBP_motif','position')
Map_RBP_motifs=Map_RBP_motifs[Map_RBP_motifs$SNPfreq>=0.05,]
Map_imputed$RBP_motif=Map_RBP_motifs$RBP_motif[match(Map_imputed$snp.name,Map_RBP_motifs$snp.name)]


RBP_list=RBP_list[order(RBP_list$SpliceFactor_short,-RBP_list$NbMotif),]
RBP_list=RBP_list[!duplicated(RBP_list$SpliceFactor_short),]

Map_imputed$RBP_motif=sapply(strsplit(Map_imputed$RBP_motif,'/'),function(x){paste(x[x%in%RBP_list$TargetMotifName],collapse='/')})
Map_imputed$SpliceFactor=sapply(strsplit(Map_imputed$RBP_motif,'/'),function(x){paste(sort(RBP_list$SpliceFactor_short[match(x,RBP_list$TargetMotifName)]),collapse='/')})

######################################################################################
#####     annotate SNPs that overlap transcriptional regulatory element       ########
######################################################################################

library(rtracklayer)
library(GenomicRanges)
ch=import.chain("/Volumes/@home/Annotation/liftoverFiles/hg38ToHg19.over.chain")
EnsRegulatory=read.table(paste(HOME,'/Annotation/Ens80_RegFeatures_hg38.txt',sep=''),sep='\t',header=T)
table(EnsRegulatory$Feature.Type)
#   CTCF Binding Site                 Enhancer                 Promoter 	Promoter Flanking Region 
#        117707               127786               	16480             85526 
EnsRegulatory$Chromosome.Name=paste('chr',EnsRegulatory$Chromosome.Name,sep='')
cur=makeGRangesFromDataFrame(EnsRegulatory,seqnames.field='Chromosome.Name',start.field='Start..bp.',end.field='End..bp.',keep.extra.columns=TRUE)
seqlevelsStyle(cur) = "UCSC"  # necessary
EnsRegulatory = as.data.frame(liftOver(cur, ch))
EnsRegulatory$seqnames=gsub('chr','',as.character(EnsRegulatory$seqnames))
EnsRegulatory$strand=as.character(EnsRegulatory$strand)
GRange_Map=makeGRangesFromDataFrame(Map_imputed[,c('chromosome','position','allele.1','allele.2','minor_allele','ancestral_allele')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_RegElt=makeGRangesFromDataFrame(EnsRegulatory[-(5:6)],seqnames.field='seqnames',start.field='start',end.field='end',keep.extra.columns=T)
Overlap=findOverlaps(GRange_Map,GRange_RegElt)

Map_imputed$RegElt=NA
Map_imputed$RegElt[queryHits(Overlap)]=EnsRegulatory$Feature.Type[subjectHits(Overlap)]

###### annotate SNPs that overlap a transcriptional regulatory element
Map_imputed$RegElt=gsub('-validated','',Map_imputed$RegElt)
Map_imputed$RegElt[is.na(Map_imputed$RegElt)]=''


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ###
##### ##### ##### ##### ##### ##### annotate SNPs located at conserved position   ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ###

library(data.table)
library(GenomicRanges)
Gerp_removed=list()
Gerp_chr=list()
for (chr in 1:22){
	cat(chr,'\n')
#	Conservation_chr=fread(sprintf('%s/Annotation/Conservation/chr%s_mammalian_conservation_perBp.txt',HOME,chr))
	hg19_align=fread(sprintf('%s/Annotation/Conservation/SEqAlign/Alignable_Human/chr%s.align.txt',HOME,chr))
	hg19_align$chrom=gsub('hg19.chr','',hg19_align$chrom,fixed=T)
	hg19_align_GR=reduce(makeGRangesFromDataFrame(as.data.frame(hg19_align),seqnames='chrom',start.field='start',end.field='end'))
	
	Gerp_chr[[chr]]=fread(sprintf('%s/Annotation/Conservation/hg19.GERP_scores/chr%s.maf.rates',HOME,chr))
	
	hg19_nalign_GR=setdiff(GRanges(seqnames=chr,range=IRanges(start=1,end=nrow(Gerp_chr[[chr]]))),hg19_align_GR)
	colnames(Gerp_chr[[chr]])=c('GerpExp','GerpRS')
	Gerp_GR=GRanges(seqnames=chr,range=IRanges(start=1:nrow(Gerp_chr[[chr]]),end=1:nrow(Gerp_chr[[chr]])))
	oo=findOverlaps(Gerp_GR,hg19_nalign_GR)
	rm(Gerp_GR)
	Gerp_removed[[chr]]=Gerp_chr[[chr]]$GerpRS[queryHits(oo)]
	Gerp_chr[[chr]]$GerpRS[queryHits(oo)]=NA
	rm(oo)
}
Gerp_all=unlist(lapply(Gerp_chr,function(x){x$GerpRS}))
Gerp_removed=unlist(Gerp_removed)
samp=sample(1:length(Gerp_all),1000000)
Gerp_mean=(Gerp_all[samp]+Gerp_all[samp+1])/2
Gerp_mean=Gerp_mean[!is.na(Gerp_mean)]
Gerp_snp=unlist(sapply(1:22,function(chr){Gerp_chr[[chr]]$GerpRS[as.integer(Map_imputed$position[which(Map_imputed$chrom==chr)])]}))

write.table(Gerp_snp,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_snp.txt',HOME),quote=F,row.names=F,col.names=F,sep='\t')
write.table(Gerp_all[samp],file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_1bp_random.txt',HOME),quote=F,row.names=F,col.names=F,sep='\t')
write.table(Gerp_mean,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_2bp_random.txt',HOME),quote=F,row.names=F,col.names=F,sep='\t')

# here we assume Map_imputed is sorted by chromsome and position
Map_imputed$GerpRS=Gerp_snp
SpliceAnnot=c('snp.name','chromosome','position','SNPfreq','Exonic','SPIDEX_Z','groupSplice','SpliceSite','SpliceSiteFlank','Branchpoint','SpliceElt','RBP_motif','SpliceFactor','RegElt','GerpRS')
write.table(Map_imputed[,SpliceAnnot],file='data/Map_RegElt.txt',quote=F,row.names=F,col.names=F,sep='\t'),




