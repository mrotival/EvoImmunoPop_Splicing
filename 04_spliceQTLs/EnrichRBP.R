
Source='Ens70_HISAT'

load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup_withSelect_RBP.Rdata',EVO_IMMUNO_POP))

branchPoints=fread("/Volumes/@Home/Annotation/intropolis/BranchPointPsoitions_MercerGenRes_2015_DataS2.bed")
colnames(branchPoints)=c('chr','start','end','id','score','strand')
branchPoints$chr=gsub('chr','',branchPoints$chr)
br_GR=flank(makeGRangesFromDataFrame(branchPoints),5,both=T)

Map_RBP_motifs=as.data.frame(fread('/Users/maxime/bin/HOMER/Z_results/Map_RBP_motifs.txt'))
colnames(Map_RBP_motifs)=c('snp.name','chromosome','allele.1','allele.2','ancestral_allele','SNPfreq','daf_char_AFB','daf_char_EUB','RBP_motif','position')

Map_GR=makeGRangesFromDataFrame(Map_RBP_motifs,start.field='position',end.field='position',keep.extra=TRUE)

Map_Gerp=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Map_Gerp_snp.txt',HOME))
colnames(Map_Gerp)=c('snp.name','chromosome','position','SNPfreq','GerpRS')
Map_RBP_motifs$GerpRS=Map_Gerp$GerpRS[match(Map_RBP_motifs$snp.name,Map_Gerp$snp.name)]

Source='Ens70_HISAT'

load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup_withSelect_RBP.Rdata',EVO_IMMUNO_POP))

branchPoints=fread("/Volumes/@Home/Annotation/intropolis/BranchPointPsoitions_MercerGenRes_2015_DataS2.bed")
colnames(branchPoints)=c('chr','start','end','id','score','strand')
branchPoints$chr=gsub('chr','',branchPoints$chr)
br_GR=flank(makeGRangesFromDataFrame(branchPoints),5,both=T)

Map_RBP_motifs=as.data.frame(fread('/Users/maxime/bin/HOMER/Z_results/Map_RBP_motifs.txt'))
colnames(Map_RBP_motifs)=c('snp.name','chromosome','allele.1','allele.2','ancestral_allele','SNPfreq','daf_char_AFB','daf_char_EUB','RBP_motif','position')

Map_GR=makeGRangesFromDataFrame(Map_RBP_motifs,start.field='position',end.field='position',keep.extra=TRUE)

Map_Gerp=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Map_Gerp_snp.txt',HOME))
colnames(Map_Gerp)=c('snp.name','chromosome','position','SNPfreq','GerpRS')
Map_RBP_motifs$GerpRS=Map_Gerp$GerpRS[match(Map_RBP_motifs$snp.name,Map_Gerp$snp.name)]

SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.4.txt',HOME))
wIntroExpressed=!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand) & SpliceSites$maxFPKM_overlappingGene_withStrand>1
wGerp=!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align
Splice_GR=makeGRangesFromDataFrame(SpliceSites[which(wIntroExpressed & wGerp)],seqnames='chrom_intron',start.field='start_site',end.field='end_site',strand='strand_correct')

Map_RBP_motifs$SpliceSite=''
oo=findOverlaps(Map_GR,Splice_GR)
Map_RBP_motifs$SpliceSite[queryHits(oo)]=paste(SpliceSites$type[subjectHits(oo)],SpliceSites$WeakAltConst[subjectHits(oo)],ifelse(SpliceSites$Gerp_site[subjectHits(oo)]>2,'conserved','neutral'))
oo=findOverlaps(Map_GR,c(flank(Splice_GR[SpliceSites$isStart[which(wIntroExpressed & wGerp)],],start=F,8,ignore.strand=TRUE),flank(Splice_GR[SpliceSites$isEnd[which(wIntroExpressed & wGerp)],],start=T,8,ignore.strand=TRUE)))
#SpliceSites[is.na(SpliceSites$WeakAltConst),]
Map_RBP_motifs$SpliceSiteFlank=''
Map_RBP_motifs$SpliceSiteFlank[queryHits(oo)]=paste(SpliceSites$type[subjectHits(oo)],SpliceSites$WeakAltConst[subjectHits(oo)],ifelse(SpliceSites$Gerp_site[subjectHits(oo)]>2,'conserved','neutral'))

acceptor_GR=Splice_GR[SpliceSites$type[which(wIntroExpressed & wGerp)]=='acceptor',]
acceptor_GR$WeakAltConst=SpliceSites$WeakAltConst[which(wIntroExpressed & wGerp & SpliceSites$type=="acceptor")]
acceptor_GR$conserved=ifelse(SpliceSites$Gerp_site[which(wIntroExpressed & wGerp & SpliceSites$type=='acceptor')]>2,'conserved','neutral')
branchpoint_expected=flank(GenomicRanges::shift(acceptor_GR,ifelse(strand(acceptor_GR)=='-',16,-16)),21,start=T)
oo=findOverlaps(Map_GR,branchpoint_expected)
Map_RBP_motifs$SpliceElement=''
Map_RBP_motifs$SpliceElement[queryHits(oo)]=paste(acceptor_GR$WeakAltConst,acceptor_GR$conserved)[subjectHits(oo)]
#oo=findOverlaps(Map_GR,br_GR)
#Map_RBP_motifs$SpliceElement[queryHits(oo)]='observed'

#load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup_withSelect_RBP.Rdata',EVO_IMMUNO_POP)
RESobs_nodup_1Mb$SpliceSiteFlank=Map_RBP_motifs$SpliceSiteFlank[match(RESobs_nodup_1Mb$snps,Map_RBP_motifs$snp.name)]
RESobs_nodup_1Mb$SpliceSite=Map_RBP_motifs$SpliceSite[match(RESobs_nodup_1Mb$snps,Map_RBP_motifs$snp.name)]
RESobs_nodup_1Mb$SpliceElement=Map_RBP_motifs$SpliceElement[match(RESobs_nodup_1Mb$snps,Map_RBP_motifs$snp.name)]

RBP_list=as.data.frame(fread('/Users/maxime/bin/HOMER/Z_results/RBP_list.txt'))
colnames(RBP_list)[c(1:4,12)]=c('TargetMotifName','SpliceFactor','targetSeq','Entropy','SpliceFactor_short')
RESobs_nodup_1Mb$RBP_motif=sapply(strsplit(RESobs_nodup_1Mb$RBP_motif,'/'),function(x){paste(x[x%in%RBP_list$TargetMotifName],collapse='/')})
RBP_list$SpliceFactor_short=sapply(strsplit(RBP_list$SpliceFactor,'(',fixed=T),function(x){x[1]})
RESobs_nodup_1Mb$RBP_name=sapply(strsplit(RESobs_nodup_1Mb$RBP_motif,'/'),function(x){paste(sort(RBP_list$SpliceFactor_short[match(x,RBP_list$TargetMotifName)]),collapse='/')})


RESobs_nodup_1Mb$SpliceElt=ifelse(grepl('acceptor',RESobs_nodup_1Mb$SpliceSite),'acceptor site',
                            ifelse(grepl('donor',RESobs_nodup_1Mb$SpliceSite),'donor site',
                            ifelse(grepl('acceptor',RESobs_nodup_1Mb$SpliceSiteFlank),'acceptor flank',
                            ifelse(grepl('donor',RESobs_nodup_1Mb$SpliceSiteFlank),'donor flank',
                            ifelse(RESobs_nodup_1Mb$SpliceElement!='','branchpoint','')))))

save(RESobs_nodup_1Mb,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/cis-psiQTL_MISO_counts_V7_nodup_withSelect_RBP_coding_intron_aSNP_rsQTL_spliceElt.Rdata',HOME))
load(file='/Volumes/@home/03_Analysis/Splicing/Papier_Splicing/V5/data/Cis_1Mb_PSI_v5_16173_5pctFreq_19597193_snps.Rdata')

j=0
RBP_list$NbMotif=NA
for(i in RBP_list$TargetMotifName){
	cat(j,':',i,'\n')
	j=j+1
	has_motif=grepl(i,Map_RBP_motifs$RBP_motif[wCis])
	RBP_list$NbMotif[j]=sum(has_motif)
	}
RBP_list=RBP_list[order(RBP_list$SpliceFactor_short,-RBP_list$V10),]
RBP_list=RBP_list[!duplicated(RBP_list$SpliceFactor_short),]

OR=NULL
j=0
is_sQTL=Map_RBP_motifs$snp.name[wCis]%in%RESobs_nodup_1Mb$snps[!duplicated(RESobs_nodup_1Mb$haplo)]
for(i in RBP_list$TargetMotifName){
	cat(j,':',i,'\n')
	j=j+1
	has_motif=grepl(i,Map_RBP_motifs$RBP_motif[wCis])
	conserved=Map_RBP_motifs$GerpRS[wCis]>2 	
	if(sum(has_motif)>100){
		myOR=odds.ratio(table(is_sQTL,has_motif))
		myOR$motif=i
		myOR$Elt='RBP'
		myOR$conserved=''
    	OR=rbind(OR,myOR)
		myOR=odds.ratio(table(is_sQTL,has_motif & conserved))
		myOR$motif=i
		myOR$Elt='RBP'
		myOR$conserved='conserved'
		OR=rbind(OR,myOR)
	}
}
	
splice_type=c('weak|alternative|constitutive','weak','alternative','constitutive','weak|alternative')

Require('Map_imputed')
Map_imputed$RegElt[is.na(Map_imputed$RegElt)]=''

RegElt=c("CTCF Binding Site", "Enhancer", "Promoter", "Promoter Flanking Region")
is_sQTL=Map_RBP_motifs$snp.name[wCis]%in%RESobs_nodup_1Mb$snps
for (i in splice_type){
    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceSite[wCis])
    conserved=Map_RBP_motifs$GerpRS[wCis]>2 
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=i
    myOR$conserved=''
    myOR$Elt='splice_site'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & conserved))
    myOR$motif=i
    myOR$conserved='conserved'
    myOR$Elt='splice_site'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & !conserved))
    myOR$motif=i
    myOR$conserved='non conserved'
    myOR$Elt='splice_site'
    OR=rbind(OR,myOR)
    }

for (i in splice_type){
    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceSiteFlank[wCis])
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved=''
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & conserved))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & !conserved))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved='non conserved'
    OR=rbind(OR,myOR)
    }

for (i in splice_type){
    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceElement[wCis])
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=i
    myOR$Elt='branchpoint'
    myOR$conserved=''
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & conserved))
    myOR$motif=i
    myOR$Elt='branchpoint'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & !conserved))
    myOR$motif=i
    myOR$conserved='non conserved'
    myOR$Elt='branchpoint'
    OR=rbind(OR,myOR)
    }

    Splice_site_type=Map_RBP_motifs$GerpRS[wCis]>2 & Map_RBP_motifs$GerpRS[wCis]< 4
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=' 2 < GerpRS < 4'
    myOR$Elt='conservation'
    myOR$conserved='conserved'    
    OR=rbind(OR,myOR)

    Splice_site_type=Map_RBP_motifs$GerpRS[wCis] >= 4
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=' GerpRS >= 4'
    myOR$Elt='conservation'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)


for (i in splice_type){
    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceSiteFlank[wCis])
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved=''
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & conserved))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & !conserved))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved='non conserved'
    OR=rbind(OR,myOR)
    }


for (i in RegElt){
    RegElt_type=gsub('-validated','',Map_imputed$RegElt[wCis])==i
    myOR=odds.ratio(table(is_sQTL,RegElt_type))
    myOR$motif=i
    myOR$Elt='RegElt'
    myOR$conserved=''
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,RegElt_type & conserved))
    myOR$motif=i
    myOR$Elt='RegElt'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,RegElt_type & !conserved))
    myOR$motif=i
    myOR$Elt='RegElt'
    myOR$conserved='non conserved'
    OR=rbind(OR,myOR)
    }
#for (i in c('observed|predicted','observed','predicted')){
#    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceElement[wCis])
#    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
#    myOR$motif=i
#    myOR$Elt='branchpoint'
#    OR=rbind(OR,myOR)
#    }

OR2=merge(OR,RBP_list,by.x='motif',by.y='TargetMotifName',all.x=T)
OR2$name=ifelse(OR2$Elt=='RBP',OR2$SpliceFactor_short,ifelse(OR2$Elt%in%c('conservation','RegElt'),OR2$motif,OR2$Elt))
OR2$Elt[OR2$Elt%in%c('splice_site_Flank','splice_site','branchpoint')]='splicing'

write.table(OR2,file='/Volumes/@home/03_Analysis/Splicing/Papier_Splicing/V7/data/OR_sQTL_functional_elements_V2.txt',quote=F,sep='\t',row.names=F)

#p <- ggplot(OR,aes(name,log2(OR),col=Elt)) + geom_pointrange(ymin= log2(OR$LowerCI),ymax= log2(OR$UpperCI)) + ylim(c(-2,10))
#p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

OR2=OR2[OR2$Elt%in%'conservation'| ((OR2$Elt=='splicing' & grepl('weak\\|alternative\\|constitutive',OR2$motif) | OR2$Elt%in%c('RBP','RegElt')) & !grepl('conserved',OR2$conserved)),]

OR2$Elt=factor(OR2$Elt, levels=c('splicing','conservation','RBP','RegElt'))
OR2=OR2[order(OR2$Elt,OR2$P),]
OR2$name = factor(OR2$name, levels=OR2[order(OR2$Elt,-OR2$LowerCI), "name"])

OR2$fdr=p.adjust(OR2$P,'fdr')
w=which(OR2$fdr<0.05 | OR2$Elt!='RBP')
pdf('/Volumes/@Home/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/FigSX_EnrichmentRegElt.pdf',width=5,heigh=4)
p <- ggplot(OR2[w,],aes(name,log2(OR),col=Elt)) + geom_pointrange( ymin = log2(OR2$LowerCI[w]),ymax= log2(OR2$UpperCI[w])) + ylim(c(-1,6))
p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(aes(yintercept=0))
dev.off()


###################################
###	    OLD :    IGNORE THIS    ###
###################################


tab1=table(gsub('-validated','',RESobs_nodup_1Mb$RegElt),RESobs_nodup_1Mb$event_type,exclude='')
#tab2=table(gsub('-validated','',RESobs_nodup_5kb$RegElt),RESobs_nodup_5kb$event_type,exclude='')
barplot(cbind(tab0/sum(tab0),t(t(tab1)/apply(tab1,2,sum))))
#barplot(cbind(tab0/sum(tab0),tab2/sum(tab2)))

#pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/SPIDEX/RegElt_barplot.pdf',HOME))

###################################
###	Figure 2B	RegElt detail 	###
###################################
par(xpd=T)
barplot(cbind(All=tab0[1:4]/sum(tab0),t(t(tab1[1:4,])/apply(tab1,2,sum))),col=c(colPSI[c(1:3,6)]),las=1,ylab='% of splice QTL',ylim=c(0,0.4))
legend('topright',fill=c(colPSI[c(1:3,6)]),c('CTCF','Enhancer','Promoter','Promoter Flanks'),bty='n')
par(xpd=F)
###################################
###	Figure 2B RegElt detail END ###
###################################

###################################
###	Figure 2B RegElt grouped	###
###################################
par(xpd=T)
barplot(cbind(All=tab0[1:4]/sum(tab0),sQTL=apply(tab1[1:4,],1,sum)/sum(tab1),rsQTL=apply(tab2[1:4,],1,sum)/sum(tab2)),col=c(colPSI[c(1:3,6)]),las=1,ylab='% of splice QTL',las=3)
legend(fill=c(colPSI[c(1:3,6)]),legend=c('CTCF','Enhancer','Promoter','Promoter Flanks'),bty='n',ncol=3)
par(xpd=F)
# las =3 for horizontal
###################################
###	Figure 2B RegElt grouped END###
###################################

tabAll=table(gsub('-validated','',Map_imputed$RegElt[wCis]),exclude='');names(tabAll)[is.na(names(tabAll))]='x'
X=cbind(tabAll,sapply(c(1e-5,1e-10,1e-30),function(th){tab=table(gsub('-validated','',RESobs_nodup_1Mb$RegElt[RESobs_nodup_1Mb$pval<th]),exclude='');names(tab)[is.na(names(tab))]='x';tab=tab[names(tabAll)];tab[is.na(tab)]=0;tab}))
Proportions_RegElt=t(t(X[-5,])/apply(X,2,sum))
colnames(Proportions_RegElt)=c('All','sQTL, P<10-5','sQTL, P<10-10','sQTL, P<10-30')
###################################
###		Figure 2B 				###
###################################
barplot(Proportions_RegElt,col=colPSI[c(1:3,6)],las=1,ylim=c(0,0.4),ylab='% of splice QTL')
barplot(Proportions_RegElt,col=colPSI[c(1:3,6)],las=1,ylim=c(0,0.4),ylab='% of splice QTL',las=2)
###################################
###	Figure 2B Significance		###
###################################


SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.4.txt',HOME))
wIntroExpressed=!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand) & SpliceSites$maxFPKM_overlappingGene_withStrand>1
wGerp=!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align
Splice_GR=makeGRangesFromDataFrame(SpliceSites[which(wIntroExpressed & wGerp)],seqnames='chrom_intron',start.field='start_site',end.field='end_site',strand='strand_correct')

Map_RBP_motifs$SpliceSite=''
oo=findOverlaps(Map_GR,Splice_GR)
Map_RBP_motifs$SpliceSite[queryHits(oo)]=paste(SpliceSites$type[subjectHits(oo)],SpliceSites$WeakAltConst[subjectHits(oo)],ifelse(SpliceSites$Gerp_site[subjectHits(oo)]>2,'conserved','neutral'))
oo=findOverlaps(Map_GR,c(flank(Splice_GR[SpliceSites$isStart[which(wIntroExpressed & wGerp)],],start=F,8,ignore.strand=TRUE),flank(Splice_GR[SpliceSites$isEnd[which(wIntroExpressed & wGerp)],],start=T,8,ignore.strand=TRUE)))
#SpliceSites[is.na(SpliceSites$WeakAltConst),]
Map_RBP_motifs$SpliceSiteFlank=''
Map_RBP_motifs$SpliceSiteFlank[queryHits(oo)]=paste(SpliceSites$type[subjectHits(oo)],SpliceSites$WeakAltConst[subjectHits(oo)],ifelse(SpliceSites$Gerp_site[subjectHits(oo)]>2,'conserved','neutral'))

acceptor_GR=Splice_GR[SpliceSites$type[which(wIntroExpressed & wGerp)]=='acceptor',]
acceptor_GR$WeakAltConst=SpliceSites$WeakAltConst[which(wIntroExpressed & wGerp & SpliceSites$type=="acceptor")]
acceptor_GR$conserved=ifelse(SpliceSites$Gerp_site[which(wIntroExpressed & wGerp & SpliceSites$type=='acceptor')]>2,'conserved','neutral')
branchpoint_expected=flank(GenomicRanges::shift(acceptor_GR,ifelse(strand(acceptor_GR)=='-',16,-16)),21,start=T)
oo=findOverlaps(Map_GR,branchpoint_expected)
Map_RBP_motifs$SpliceElement=''
Map_RBP_motifs$SpliceElement[queryHits(oo)]=paste(acceptor_GR$WeakAltConst,acceptor_GR$conserved)[subjectHits(oo)]
#oo=findOverlaps(Map_GR,br_GR)
#Map_RBP_motifs$SpliceElement[queryHits(oo)]='observed'

#load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup_withSelect_RBP.Rdata',EVO_IMMUNO_POP)
#RESobs_nodup_1Mb$SpliceSiteFlank=Map_RBP_motifs$SpliceSiteFlank[match(RESobs_nodup_1Mb$snps,Map_RBP_motifs$snp.name)]
#RESobs_nodup_1Mb$SpliceSite=Map_RBP_motifs$SpliceSite[match(RESobs_nodup_1Mb$snps,Map_RBP_motifs$snp.name)]
#RESobs_nodup_1Mb$SpliceElement=Map_RBP_motifs$SpliceElement[match(RESobs_nodup_1Mb$snps,Map_RBP_motifs$snp.name)]

load(file='/Volumes/@home/03_Analysis/Splicing/Papier_Splicing/V5/data/Cis_1Mb_PSI_v5_16173_5pctFreq_19597193_snps.Rdata')
RBP_list=as.data.frame(fread('/Users/maxime/bin/HOMER/Z_results/RBP_list.txt'))
colnames(RBP_list)[c(1:4,12)]=c('TargetMotifName','SpliceFactor','targetSeq','Entropy','SpliceFactor_short')

j=0
RBP_list$NbMotif=NA
for(i in RBP_list$TargetMotifName){
	cat(j,':',i,'\n')
	j=j+1
	has_motif=grepl(i,Map_RBP_motifs$RBP_motif[wCis])
	RBP_list$NbMotif[j]=sum(has_motif)
	}
RBP_list=RBP_list[order(RBP_list$SpliceFactor_short,-RBP_list$NbMotif),]
RBP_list=RBP_list[!duplicated(RBP_list$SpliceFactor_short),]

OR=NULL
j=0
is_sQTL=Map_RBP_motifs$snp.name[wCis]%in%RESobs_nodup_1Mb$snps[!duplicated(RESobs_nodup_1Mb$haplo)]
for(i in RBP_list$TargetMotifName){
	cat(j,':',i,'\n')
	j=j+1
	has_motif=grepl(i,Map_RBP_motifs$RBP_motif[wCis])
	conserved=Map_RBP_motifs$GerpRS[wCis]>2 	
	if(sum(has_motif)>100){
		myOR=odds.ratio(table(is_sQTL,has_motif))
		myOR$motif=i
		myOR$Elt='RBP'
		myOR$conserved=''
    	OR=rbind(OR,myOR)
		myOR=odds.ratio(table(is_sQTL,has_motif & conserved))
		myOR$motif=i
		myOR$Elt='RBP'
		myOR$conserved='conserved'
		OR=rbind(OR,myOR)
	}
}
	
splice_type=c('weak|alternative|constitutive','weak','alternative','constitutive','weak|alternative')

Require('Map_imputed')
Map_imputed$RegElt[is.na(Map_imputed$RegElt)]=''

RegElt=c("CTCF Binding Site", "Enhancer", "Promoter", "Promoter Flanking Region")
is_sQTL=Map_RBP_motifs$snp.name[wCis]%in%RESobs_nodup_1Mb$snps
for (i in splice_type){
    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceSite[wCis])
    conserved=Map_RBP_motifs$GerpRS[wCis]>2 
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=i
    myOR$conserved=''
    myOR$Elt='splice_site'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & conserved))
    myOR$motif=i
    myOR$conserved='conserved'
    myOR$Elt='splice_site'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & !conserved))
    myOR$motif=i
    myOR$conserved='non conserved'
    myOR$Elt='splice_site'
    OR=rbind(OR,myOR)
    }

for (i in splice_type){
    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceSiteFlank[wCis])
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved=''
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & conserved))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & !conserved))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved='non conserved'
    OR=rbind(OR,myOR)
    }

for (i in splice_type){
    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceElement[wCis])
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=i
    myOR$Elt='branchpoint'
    myOR$conserved=''
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & conserved))
    myOR$motif=i
    myOR$Elt='branchpoint'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & !conserved))
    myOR$motif=i
    myOR$conserved='non conserved'
    myOR$Elt='branchpoint'
    OR=rbind(OR,myOR)
    }

    Splice_site_type=Map_RBP_motifs$GerpRS[wCis]>2 & Map_RBP_motifs$GerpRS[wCis]< 4
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=' 2 < GerpRS < 4'
    myOR$Elt='conservation'
    myOR$conserved='conserved'    
    OR=rbind(OR,myOR)

    Splice_site_type=Map_RBP_motifs$GerpRS[wCis] >= 4
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=' GerpRS >= 4'
    myOR$Elt='conservation'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)


for (i in splice_type){
    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceSiteFlank[wCis])
    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved=''
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & conserved))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,Splice_site_type & !conserved))
    myOR$motif=i
    myOR$Elt='splice_site_Flank'
    myOR$conserved='non conserved'
    OR=rbind(OR,myOR)
    }


for (i in RegElt){
    RegElt_type=gsub('-validated','',Map_imputed$RegElt[wCis])==i
    myOR=odds.ratio(table(is_sQTL,RegElt_type))
    myOR$motif=i
    myOR$Elt='RegElt'
    myOR$conserved=''
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,RegElt_type & conserved))
    myOR$motif=i
    myOR$Elt='RegElt'
    myOR$conserved='conserved'
    OR=rbind(OR,myOR)
    myOR=odds.ratio(table(is_sQTL,RegElt_type & !conserved))
    myOR$motif=i
    myOR$Elt='RegElt'
    myOR$conserved='non conserved'
    OR=rbind(OR,myOR)
    }
#for (i in c('observed|predicted','observed','predicted')){
#    Splice_site_type=grepl(i,Map_RBP_motifs$SpliceElement[wCis])
#    myOR=odds.ratio(table(is_sQTL,Splice_site_type))
#    myOR$motif=i
#    myOR$Elt='branchpoint'
#    OR=rbind(OR,myOR)
#    }

OR2=merge(OR,RBP_list,by.x='motif',by.y='TargetMotifName',all.x=T)
OR2$name=ifelse(OR2$Elt=='RBP',OR2$SpliceFactor_short,ifelse(OR2$Elt%in%c('conservation','RegElt'),OR2$motif,OR2$Elt))
OR2$Elt[OR2$Elt%in%c('splice_site_Flank','splice_site','branchpoint')]='splicing'

write.table(OR2,file='/Volumes/@home/03_Analysis/Splicing/Papier_Splicing/V7/data/OR_sQTL_functional_elements_V2.txt',quote=F,sep='\t',row.names=F)

#p <- ggplot(OR,aes(name,log2(OR),col=Elt)) + geom_pointrange(ymin= log2(OR$LowerCI),ymax= log2(OR$UpperCI)) + ylim(c(-2,10))
#p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

OR2=OR2[OR2$Elt%in%'conservation'| ((OR2$Elt=='splicing' & grepl('weak\\|alternative\\|constitutive',OR2$motif) | OR2$Elt%in%c('RBP','RegElt')) & !grepl('conserved',OR2$conserved)),]

OR2$Elt=factor(OR2$Elt, levels=c('splicing','conservation','RBP','RegElt'))
OR2=OR2[order(OR2$Elt,OR2$P),]
OR2$name = factor(OR2$name, levels=OR2[order(OR2$Elt,-OR2$LowerCI), "name"])

OR2$fdr=p.adjust(OR2$P,'fdr')
w=which(OR2$fdr<0.05 | OR2$Elt!='RBP')
pdf('/Volumes/@Home/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/FigSX_EnrichmentRegElt.pdf',width=5,heigh=4)
p <- ggplot(OR2[w,],aes(name,log2(OR),col=Elt)) + geom_pointrange( ymin = log2(OR2$LowerCI[w]),ymax= log2(OR2$UpperCI[w])) + ylim(c(-1,6))
p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(aes(yintercept=0))
dev.off()


###################################
###	    OLD :    IGNORE THIS    ###
###################################


tab1=table(gsub('-validated','',RESobs_nodup_1Mb$RegElt),RESobs_nodup_1Mb$event_type,exclude='')
#tab2=table(gsub('-validated','',RESobs_nodup_5kb$RegElt),RESobs_nodup_5kb$event_type,exclude='')
barplot(cbind(tab0/sum(tab0),t(t(tab1)/apply(tab1,2,sum))))
#barplot(cbind(tab0/sum(tab0),tab2/sum(tab2)))

#pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/SPIDEX/RegElt_barplot.pdf',HOME))

###################################
###	Figure 2B	RegElt detail 	###
###################################
par(xpd=T)
barplot(cbind(All=tab0[1:4]/sum(tab0),t(t(tab1[1:4,])/apply(tab1,2,sum))),col=c(colPSI[c(1:3,6)]),las=1,ylab='% of splice QTL',ylim=c(0,0.4))
legend('topright',fill=c(colPSI[c(1:3,6)]),c('CTCF','Enhancer','Promoter','Promoter Flanks'),bty='n')
par(xpd=F)
###################################
###	Figure 2B RegElt detail END ###
###################################

###################################
###	Figure 2B RegElt grouped	###
###################################
par(xpd=T)
barplot(cbind(All=tab0[1:4]/sum(tab0),sQTL=apply(tab1[1:4,],1,sum)/sum(tab1),rsQTL=apply(tab2[1:4,],1,sum)/sum(tab2)),col=c(colPSI[c(1:3,6)]),las=1,ylab='% of splice QTL',las=3)
legend(fill=c(colPSI[c(1:3,6)]),legend=c('CTCF','Enhancer','Promoter','Promoter Flanks'),bty='n',ncol=3)
par(xpd=F)
# las =3 for horizontal
###################################
###	Figure 2B RegElt grouped END###
###################################

tabAll=table(gsub('-validated','',Map_imputed$RegElt[wCis]),exclude='');names(tabAll)[is.na(names(tabAll))]='x'
X=cbind(tabAll,sapply(c(1e-5,1e-10,1e-30),function(th){tab=table(gsub('-validated','',RESobs_nodup_1Mb$RegElt[RESobs_nodup_1Mb$pval<th]),exclude='');names(tab)[is.na(names(tab))]='x';tab=tab[names(tabAll)];tab[is.na(tab)]=0;tab}))
Proportions_RegElt=t(t(X[-5,])/apply(X,2,sum))
colnames(Proportions_RegElt)=c('All','sQTL, P<10-5','sQTL, P<10-10','sQTL, P<10-30')
###################################
###		Figure 2B 				###
###################################
barplot(Proportions_RegElt,col=colPSI[c(1:3,6)],las=1,ylim=c(0,0.4),ylab='% of splice QTL')
barplot(Proportions_RegElt,col=colPSI[c(1:3,6)],las=1,ylim=c(0,0.4),ylab='% of splice QTL',las=2)
###################################
###	Figure 2B Significance		###
###################################

