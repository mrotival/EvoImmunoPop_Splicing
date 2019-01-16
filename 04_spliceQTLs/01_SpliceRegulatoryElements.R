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



##########################################################################
#####		Mechanistic Bases  of splicing	Fig 4B, and Fig S5A		  ####
##########################################################################
RESobs_nodup_1Mb=fread('data/RESobs_nodup_1Mb.txt')
RESobs_nodup_1Mb_cond=fread('data/RESobs_nodup_1Mb_cond.txt')


###################################
###		Figure 4A	distance	###
###################################
HLA_region=c('6:28Mb','6:29Mb','6:30Mb','6:31Mb','6:32Mb','6:33Mb')
mean(abs(RESobs_nodup_1Mb$CisDist)<1e6) # 66%
mean(abs(RESobs_nodup_1Mb$CisDist[RESobs_nodup_1Mb$locus!=HLA_region])<1e6) #66% 

plot(RESobs_nodup_1Mb$CisDist/1e6,-log10(RESobs_nodup_1Mb$pval),pch=16,col=gsub('AA','44',colERC5[1]),xlim=c(-1,1),las=1,axes=F,ylab=expression(-log[10](P-value)),xlab='Distance from Event boundaries')
axis(2,las=1);
axis(1,at=c(-1,-0.5,0,0.5,1),labels=c('1Mb','500kb','0','500kb','1Mb'))


# annotate sQTLs 
RESobs_nodup_1Mb$SpliceElt=Map_imputed$SpliceElt[mm]
RESobs_nodup_1Mb$RBP_motif=Map_imputed$RBP_motif[mm]
RESobs_nodup_1Mb$SpliceFactor=Map_imputed$SpliceFactor[mm]
RESobs_nodup_1Mb$RegElt=Map_imputed$RegElt[mm]
RESobs_nodup_1Mb$GerpRS=Map_imputed$GerpRS[mm]

####### define the set of common variants (MF>5%) in Cis from a splicing event ######
library(GenomicRanges)
Map_GR=makeGRangesFromDataFrame(Map_imputed[,c('chromosome','position','allele.1','allele.2','ancestral_allele')],seqnames.field='chromosome',start.field='position',end.field='position',keep.extra.columns=TRUE)
GRange_PSI=makeGRangesFromDataFrame(PSI_Annot,seqnames.field='chrom',start.field='start',end.field='end',keep.extra.columns=TRUE)
Cisdist=1e6
GR_cis=reduce(union(flank(GRange_PSI,Cisdist,start=T),union(flank(GRange_PSI,Cisdist,start=F),GRange_PSI)))
ooCis=findOverlaps(GR_cis, Map_GR)
wCis=unique(subjectHits(ooCis))
###### END definition of Cis variants 


###################################
###		Figure 4B 				###
###################################
tab=table(Map_imputed$groupSplice[wCis])
ImpactOnSplicing=cbind(tab,sapply(c(1e-5,1e-10,1e-30),function(th){tab=table(RESobs_nodup_1Mb$groupSplice[RESobs_nodup_1Mb$pval<th])}))
Proportion_ImpactOnSplicing=t(t(ImpactOnSplicing[-1,])/apply(ImpactOnSplicing,2,sum))
colnames(Proportion_ImpactOnSplicing)=c('All','sQTL, P<10-5','sQTL, P<10-10','sQTL, P<10-30')

library(RBrewerPalette)
acol= c(rev(brewer.pal(4,'YlOrRd')),brewer.pal(4,'YlGnBu'))
rcol = SetAlpha(acol, 0.5)
par(mar=c(7,5,1,3))
barplot(Proportion_ImpactOnSplicing,col=acol[c(3,1,6,7)],las=3,ylim=c(0,0.4),ylab='% of splice QTL')

############################################################################################################################################
###########################          Enrichment in regions predicted to alter splicing (Pvalues Fig 4B)   ##################################
############################################################################################################################################

##### SPIDEX ENRICHMENTS
OR=NULL
for (j in 2:4){
	tab=rbind(c(ImpactOnSplicing[1,1]-ImpactOnSplicing[1,j],sum(ImpactOnSplicing[1,j])),
		c(sum(ImpactOnSplicing[-1,1])-sum(ImpactOnSplicing[-1,j]),sum(ImpactOnSplicing[-1,j])))
	myOR=odds.ratio(tab)
	myOR$event_type="All"
	myOR$Splice_type=colnames(Proportion_ImpactOnSplicing)[j]
	OR=rbind(OR,myOR)
	tab=rbind(c(sum(ImpactOnSplicing[c(2,4),1])-sum(ImpactOnSplicing[c(2,4),j]),sum(ImpactOnSplicing[c(2,4),j])),
			c(sum(ImpactOnSplicing[c(3,5),1])-sum(ImpactOnSplicing[c(3,5),j]),sum(ImpactOnSplicing[c(3,5),j])))
	myOR=odds.ratio(tab)
	myOR$event_type='All'
	myOR$Splice_type=paste('deleterious',colnames(Proportion_ImpactOnSplicing)[j])
	OR=rbind(OR,myOR)
	}
write.table(OR,file='data/SpliceElt_OR.txt',sep='\t',quote=F,row.names=F)

############################################################################################################################################
####################################### REGULATORY REGIONS ENRICHMENT               (Fig S4A)             ##################################
############################################################################################################################################

is_sQTL=Map_imputed$snp.name[wCis]%in%RESobs_nodup_1Mb$snps
#is_sQTL=Map_imputed$snp.name[wCis]%in%RESobs_nodup_1Mb$snps[!duplicated(RESobs_nodup_1Mb$haplo)]

Splice_site_type=Map_imputed$SpliceSite[wCis]!=''
myOR=odds.ratio(table(is_sQTL,Splice_site_type))
myOR$motif='Donor/acceptor_site'
myOR$Elt='Splicing'
OR=rbind(OR,myOR)

Splice_site_type=Map_RBP_motifs$SpliceSiteFlank[wCis]!=''
myOR=odds.ratio(table(is_sQTL,Splice_site_type))
myOR$motif='Splice site flank'
myOR$Elt='Splicing'
OR=rbind(OR,myOR)

Splice_site_type=Map_imputed$Branchpoint[wCis]!=''
myOR=odds.ratio(table(is_sQTL,Splice_site_type))
myOR$motif='Branchpoint'
myOR$Elt='Splicing'

conserved=Map_imputed$GerpRS[wCis]>2 & Map_imputed$GerpRS[wCis]< 4
myOR=odds.ratio(table(is_sQTL,conserved))
myOR$motif=' 2 < GerpRS < 4'
myOR$Elt='Conservation'
OR=rbind(OR,myOR)

conserved=Map_imputed$GerpRS[wCis] >= 4
myOR=odds.ratio(table(is_sQTL,conserved))
myOR$motif=' GerpRS >= 4'
myOR$Elt='Conservation'
OR=rbind(OR,myOR)

for(i in RBP_list$TargetMotifName){
	cat(j,':',i,'\n')
	j=j+1
	has_motif=grepl(i,Map_imputed$RBP_motif[wCis])
	conserved=Map_imputed$GerpRS[wCis]>2
	if(sum(has_motif)>100){
		myOR=odds.ratio(table(is_sQTL,has_motif))
		myOR$motif=i
		myOR$Elt='RNA binding protein'
    	OR=rbind(OR,myOR)
	}
}

RegElt=c("CTCF Binding Site", "Enhancer", "Promoter", "Promoter Flanking Region")
for (i in RegElt){
    RegElt_type=Map_imputed$RegElt[wCis]==i
    myOR=odds.ratio(table(is_sQTL,RegElt_type))
    myOR$motif=i
    myOR$Elt='RegElt'
    OR=rbind(OR,myOR)
    }

OR_table=merge(OR,RBP_list,by.x='motif',by.y='TargetMotifName',all.x=T)
OR_table$name=ifelse(OR_table$Elt=='RNA binding protein',OR_table$SpliceFactor_short,OR_table$motif)

write.table(OR_table,file='/Volumes/@home/03_Analysis/Splicing/Papier_Splicing/V7/data/OR_sQTL_functional_elements_V2.txt',quote=F,sep='\t',row.names=F)

OR_table$Elt=factor(OR_table$Elt, levels=c('Splicing','Conservation','RNA binding protein','RegElt'))
OR_table=OR_table[order(OR_table$Elt,OR_table$P),]
OR_table$name = factor(OR2$name, levels=OR2[order(OR2$Elt,-OR2$LowerCI), "name"])

OR_table$fdr=p.adjust(OR2$P,'fdr')

######################################################
#########           Figure S4A                ########
######################################################

pdf('figures/FigS4A_EnrichmentRegElt.pdf',width=5,heigh=4)
w=which(OR_table$fdr<0.05 | OR_table$Elt!='RNA binding protein')
p <- ggplot(OR_table[w,],aes(name,log2(OR),col=Elt)) + geom_pointrange( ymin = log2(OR_table$LowerCI[w]),ymax= log2(OR_table$UpperCI[w])) + ylim(c(-1,6))
p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(aes(yintercept=0))
dev.off()

######################################################
#########           Annotation of sQTL        ########
######################################################

mm=match(RESobs_nodup_1Mb$snps,Map_imputed$snp.name)
#### selection
RESobs_nodup_1Mb$FST=Map_imputed$FST_adj[mm]
RESobs_nodup_1Mb$iHS_AFB=Map_imputed$iHS_AFB[mm]
RESobs_nodup_1Mb$iHS_EUB=Map_imputed$iHS_EUB[mm]
#RESobs_nodup_1Mb$maf_AFB=Map_imputed$maf_AFB[mm]		
#RESobs_nodup_1Mb$maf_EUB=Map_imputed$maf_EUB[mm]

#### mechanisms
RESobs_nodup_1Mb$RegElt=Map_imputed$RegElt[mm]
RESobs_nodup_1Mb$SPIDEX_Z=Map_imputed$SPIDEX_Z[mm]
RESobs_nodup_1Mb$Exonic=Map_imputed$Exonic[mm]
RESobs_nodup_1Mb$groupSplice=Map_imputed$groupSplice[mm]

RESobs_nodup_1Mb$SpliceElt=Map_imputed$SpliceElt[mm]
RESobs_nodup_1Mb$SpliceFactor=Map_imputed$SpliceFactor[mm]

#### is Immune gene
RESobs_nodup_1Mb$isImmune=RESobs_nodup_1Mb$gene%in%GoToTest_genes$immuneResponse

###################################
###		Table S3A				###
###################################
paste_xy=function(x,y){paste(rep(x,length(y)),rep(y,e=length(x)),sep='')}
TableS3A=RESobs_nodup_1Mb[,c('event_id','symbol','event_type',paste('Alternatively_spliced_',condIndex),sep=''),'isImmune','cond','snps','pval','pos','CisDist','daf_char_AFB','daf_char_EUB',paste_xy(c('Pvalue_','Beta_','R2_'),condIndex),'groupSplice','spliceFactor','spliceElt','RegElt','coding_type'))]

TableS3A$isImmune=ifelse(TableS3A$isImmune,'yes','')
TableS3A$cond=condIndex[TableS3A$cond]

write.table(TableS3A,file=sprintf('tables/TableS3A_allsQTL.txt',HOME),quote=F,sep='\t',row.names=F)

###################################
###		Table S3A END			###
###################################





##########################################################################################
#####					GO enrichment of sQTLs	(Table S3B)							  ####
##########################################################################################

# sQTls are enriched in defense response genes
resGO=GOSeq(unique(RESobs_nodup_1Mb$gene),unique(PSI_Annot$gene_id),overOnly=F,addCI=T)
resGO$cond='ALL_POOLED'
resGO_cond=lapply(1:5,function(i){GOSeq(unique(RESobs_nodup_1Mb_cond$gene[RESobs_nodup_1Mb_cond$cond==i]),unique(PSI_Annot$gene_id[TESTABLE[,i]]),overOnly=F,addCI=T)})
resGO_cond_table=do.call(rbind,lapply(1:5,function(i){cbind(resGO_cond[[i]],cond=rep(condIndex[i],nrow(resGO_cond[[i]])))}))
TableS3B=rbind(resGO_cond_table,resGO)
write.table(TableS3B,file='tables/TableS3B_GOenrich.txt',quote=F,sep='\t',row.names=F)

HLA_region=c('6:28Mb','6:29Mb','6:30Mb','6:31Mb','6:32Mb','6:33Mb')
resGO=GOSeq(unique(RESobs_nodup_1Mb$gene[!RESobs_nodup_1Mb$locus%in%HLA_region)]),unique(PSI_Annot$gene_id),overOnly=F,addCI=T)

MeanExpr=GeneAnnot[,c("NS_mean", "LPS_mean", "PAM3_mean", "R848_mean", "Flu_mean")]
Expr_level=apply(MeanExpr,1,mean,na.rm=T)
Expr_level[is.na(Expr_level)]=0
Expr_bin=cut(Expr_level,c(-0.01,1,5,10,50,100,500,1000,Inf))

has_sQTL=unique(PSI_Annot$gene_id)%in%unique(RESobs_nodup_1Mb$gene)
isImmune=unique(PSI_Annot$gene_id)%in%GoToTest_genes$immuneResponse

### adjusting on gene expression
Expr_bin_test=Expr_bin[match(unique(PSI_Annot$gene_id),GeneAnnot$Ensembl.Gene.ID)]

PctWithsQTL=c()
for (nsamp in 1:10000){
	if(nsamp%%1000==1){ cat(nsamp-1,'k done\n')}
	samp=c()
	for (x in levels(Expr_bin_test)){
		samp=c(samp,sample(which(Expr_bin_test==x),sum(isImmune[Expr_bin_test==x]),replace=F))
		}
	PctWithsQTL[nsamp]=mean(has_sQTL[samp])
	}

mean(PctWithsQTL>=mean(has_sQTL[isImmune]))
#1000 resample: <0.001
#10000 resample: <0.0001

### adjusting on gene expression at the splicing event level
has_sQTL_event=unique(PSI_Annot$event_id)%in%unique(RESobs_nodup_1Mb$event_id)
isImmune_event=PSI_Annot$gene_id%in%GoToTest_genes$immuneResponse
Expr_bin_test_event=Expr_bin[match(PSI_Annot$gene_id,GeneAnnot$Ensembl.Gene.ID)]

PctEventsWithsQTL=c()
for (nsamp in 1:10000){
	if(nsamp%%1000==0){ cat(nsamp%/%1000,'k done\n')}
	samp=c()
	for (x in levels(Expr_bin_test_event)){
		samp=c(samp,sample(which(Expr_bin_test_event==x),sum(isImmune_event[Expr_bin_test_event==x]),replace=F))
		}
	PctEventsWithsQTL[nsamp]=mean(has_sQTL_event[samp])
	}
mean(PctEventsWithsQTL>=mean(has_sQTL_event[isImmune_event])) 
#1000 resample: <0.001
#10000 resample: <0.0001




write.table(RESobs_nodup_1Mb,file='data/RESobs_nodup_1Mb_SpliceElt.txt',quote=F,sep='\t',row.names=F)



