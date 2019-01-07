library(data.table)
library(ggplot2)

### load gene annotation
GeneAnnot=read.table(paste(HOME,'/Annotation/GeneAnnotation_hg37_ens70.txt',sep=''),sep='\t',comment='',quote='',stringsAsFactors=FALSE,header=T,row.names=1)
colnames(GeneAnnot)[1]='Ensembl.Gene.ID'
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='X']=23
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='Y']=24
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='MT']=26
GeneAnnot$Chromosome.Name=as.numeric(GeneAnnot$Chromosome.Name)
GeneAnnot=GeneAnnot[!is.na(GeneAnnot$Chromosome.Name),] # remove alternate chromosomes and unassigned scaffolds
GeneAnnot$Strand=ifelse(GeneAnnot$Strand>0,'+','-') # change strand from -1,1 to '+','-'


MeanExpr=GeneAnnot[grep('_mean',colnames(GeneAnnot))]
MeanExpr[is.na(MeanExpr)]=0
colnames(MeanExpr)=paste('FPKM',condIndex,sep='_')
lFC=log2(1+MeanExpr[,-1])-log2(1+MeanExpr[,1])%o%c(1,1,1,1)
colnames(lFC)=paste('log2FC',condIndex[-1],sep='_')

### load splice site and junction informations
SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.4.txt',HOME))
Junc_annot=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Junction_annot_V7.2.txt',HOME))

### extract GO informations
allGOTerms=as.data.frame(fread('allGOterms_EnsGRC37_13042017.txt'))
GoToTest=list(immuneResponse='GO:0006955',
	defenseResponse='GO:0006952',
	InnateImmuneResponse='GO:0045087',
	AdaptiveImmuneResponse='GO:0002250',
	AntiviralResponse='GO:0051607',
	AntibacterialResponse='GO:0042742',
	TF_DNAbinding='GO:0003700',
	development='GO:0032502',
	olfactory_receptors='GO:0004984',
	nonsense_mediated_decay='GO:0000184')
GoToTest_genes=lapply(GoToTest,function(x){allGOterms$gene[allGOterms$go==x]})
GoToTest_genes$all=unique(allGOterms$gene)


FoldChange_col=c('< -1'="#41B6C480", '[-1,1]'="#FECC5CBB", '>1'="#E31A1CBB")
wIntro_Expressed=!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand) & !grepl('//',SpliceSites$overlappingGene_withStrand) & SpliceSites$maxFPKM_overlappingGene_withStrand>1
wGerp=!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align


#############################################################
#### Analysis of Number of cryptic splice site per gene  ####
#############################################################
  
Nb_cryptic_SpliceSites=SpliceSites[wGerp & wIntro_Expressed,.(NbCryptic_ANY=sum(!inEnsembl & Gerp_site<2 & WeakAltConst=='weak'),
                         NbCryptic_NS=sum(!inEnsembl & Gerp_site<2 & WeakAltConst=='weak' & Total_count_NS>0),
                         NbCryptic_LPS=sum(!inEnsembl & Gerp_site<2 & WeakAltConst=='weak' & Total_count_LPS>0),
                         NbCryptic_PAM3CSK4=sum(!inEnsembl & Gerp_site<2 & WeakAltConst=='weak' & Total_count_PAM3CSK4>0),
                         NbCryptic_R848=sum(!inEnsembl & Gerp_site<2 & WeakAltConst=='weak' & Total_count_R848>0),
                         NbCryptic_IAV=sum(!inEnsembl & Gerp_site<2 & WeakAltConst=='weak' & Total_count_IAV>0))
,by=.(overlappingGene_withStrand,overlappingSymbol_withStrand)]


# add information on immune VS non-immune genes (and logFC)
Nb_cryptic_SpliceSites$isImmune=ifelse(Nb_cryptic_SpliceSites$overlappingGene_withStrand%in% GoToTest_genes$immuneResponse,'yes','')
Nb_cryptic_SpliceSites=cbind(Nb_cryptic_SpliceSites,lFC[match(Nb_cryptic_SpliceSites$overlappingGene_withStrand,rownames(lFC)),])
Nb_cryptic_SpliceSites$max_logFC=apply(lFC[match(Nb_cryptic_SpliceSites$overlappingGene_withStrand,rownames(lFC)),],1,max)

# reshape to plot Change in Nb_cryptic_SpliceSites as a function of FC by condition

# 1. melt table  (one line per metric/condition/gene)
Nb_cryptic_SpliceSites_melted=melt(Nb_cryptic_SpliceSites,id.vars=c('overlappingGene_withStrand','overlappingSymbol_withStrand','isImmune','NbCryptic_ANY','NbCryptic_NS'),
 variable.name='variable',measure=c("log2FC_LPS","log2FC_PAM3CSK4","log2FC_R848","log2FC_IAV", "NbCryptic_LPS","NbCryptic_PAM3CSK4","NbCryptic_R848","NbCryptic_IAV"))
Nb_cryptic_SpliceSites_melted[,metric:=gsub('(NbCryptic|log2FC)_(LPS|PAM3CSK4|R848|IAV)','\\1',Nb_cryptic_SpliceSites_melted$variable)]
Nb_cryptic_SpliceSites_melted[,condition:=gsub('(NbCryptic|log2FC)_(LPS|PAM3CSK4|R848|IAV)','\\2',Nb_cryptic_SpliceSites_melted$variable)]
Nb_cryptic_SpliceSites_melted[,variable:=NULL]

# 2. cast by condition  (one line per gene/condition)
Nb_cryptic_SpliceSites_perCond= dcast(Nb_cryptic_SpliceSites_melted,overlappingGene_withStrand + overlappingSymbol_withStrand + isImmune + NbCryptic_ANY + NbCryptic_NS + condition ~ metric,value.var='value')
Nb_cryptic_SpliceSites_perCond[,DeltaCryptic:=NbCryptic-NbCryptic_NS]
Nb_cryptic_SpliceSites_perCond[,FoldChange:=cut(log2FC,c(-Inf,-1,1,Inf),labels=c('< -1','[-1,1]','>1'))]
Nb_cryptic_SpliceSites_perCond[,condition:=factor(condition,levels=c('LPS','PAM3CSK4','R848','IAV'))]

write.table(Nb_cryptic_SpliceSites_perCond,file='Nb_cryptic_SpliceSites_perCond',sep='\t',row.names=F,quote=F)

#########################################################
#####       Compute levels of splicing noise        #####
#########################################################


#all constitutive sites to consider
Const_coding_site=SpliceSites$site_id[wIntro_Expressed & SpliceSites$inEnsembl & SpliceSites$coding & SpliceSites$WeakAltConst=='constitutive']
# all cryptic sites to consider
cryptic_site=gsub(' (acceptor|donor)','',SpliceSites$site_id[wIntro_Expressed & wGerp & !SpliceSites$inEnsembl & SpliceSites$WeakAltConst=='weak' & SpliceSites$Gerp_site<2])

### Noisy junctions 
#from start site
codingStart_to_crypic=Junc_annot[Junc_annot$start_id%in%gsub(' (acceptor|donor)','',Const_coding_site) & Junc_annot$end_id%in%cryptic_site, 
                                     .(Nbread_cryptic_NS=sum(Total_count_NS),
                                      Nbread_cryptic_LPS=sum(Total_count_LPS),
                                      Nbread_cryptic_PAM3CSK4=sum(Total_count_PAM3CSK4),
                                      Nbread_cryptic_R848=sum(Total_count_R848),
                                      Nbread_cryptic_IAV=sum(Total_count_IAV)),by=start_id]
                                      

# from end site
codingEnd_to_crypic=Junc_annot[Junc_annot$end_id%in%gsub(' (acceptor|donor)','',Const_coding_site) & Junc_annot$start_id%in%cryptic_site, 
                                    .(Nbread_cryptic_NS=sum(Total_count_NS),
                                      Nbread_cryptic_LPS=sum(Total_count_LPS),
                                      Nbread_cryptic_PAM3CSK4=sum(Total_count_PAM3CSK4),
                                      Nbread_cryptic_R848=sum(Total_count_R848),
                                      Nbread_cryptic_IAV=sum(Total_count_IAV)),by=end_id]
### count all junctions 
#from start site
codingStart_total=Junc_annot[Junc_annot$start_id%in%gsub(' (acceptor|donor)','',Const_coding_site), 
                                    .(Nbread_NS=sum(Total_count_NS),
                                      Nbread_LPS=sum(Total_count_LPS),
                                      Nbread_PAM3CSK4=sum(Total_count_PAM3CSK4),
                                      Nbread_R848=sum(Total_count_R848),
                                      Nbread_IAV=sum(Total_count_IAV)),by=start_id]
                                      
# from end site
codingEnd_total=Junc_annot[Junc_annot$end_id%in%gsub(' (acceptor|donor)','',Const_coding_site),
                                    .(Nbread_NS=sum(Total_count_NS),
                                      Nbread_LPS=sum(Total_count_LPS),
                                      Nbread_PAM3CSK4=sum(Total_count_PAM3CSK4),
                                      Nbread_R848=sum(Total_count_R848),
                                      Nbread_IAV=sum(Total_count_IAV)),by=end_id]

###### merge and compute percentage of noise
codingEnd=merge(codingEnd_to_crypic, codingEnd_total,all=T)
codingEnd[is.na(codingEnd)]=0
codingEnd[,site_id:=paste(end_id,ifelse(grepl('\\+',end_id),'acceptor','donor'))]
codingEnd[,end_id:=NULL]

codingStart=merge(codingStart_to_crypic, codingStart_total,all=T)
codingStart[is.na(codingStart)]=0
codingStart[,site_id:=paste(start_id,ifelse(grepl('\\-',start_id),'acceptor','donor'))]
codingStart[,start_id:=NULL]

Coding_site=rbind(codingEnd,codingStart)
# security check. removes cases where crtytic acceptor/overlaps a constitutive donor
Coding_site=Coding_site[site_id%in%Const_coding_site,]

Coding_site[,Pct_cryptic_NS:=Nbread_cryptic_NS/Nbread_NS*100]
Coding_site[,Pct_cryptic_LPS:=Nbread_cryptic_LPS/Nbread_LPS*100]
Coding_site[,Pct_cryptic_PAM3CSK4:=Nbread_cryptic_PAM3CSK4/Nbread_PAM3CSK4*100]
Coding_site[,Pct_cryptic_R848:=Nbread_cryptic_R848/Nbread_R848*100]
Coding_site[,Pct_cryptic_IAV:=Nbread_cryptic_IAV/Nbread_IAV*100]



###### MELT the data (one line per metric/condition/splice site)
Coding_site_melted=melt(Coding_site,id.vars=c('site_id','Nbread_cryptic_NS','Nbread_NS','Pct_cryptic_NS'),variable.name='variable',
measure=c("Nbread_cryptic_NS","Nbread_cryptic_LPS","Nbread_cryptic_PAM3CSK4","Nbread_cryptic_R848","Nbread_cryptic_IAV","Nbread_NS","Nbread_LPS","Nbread_PAM3CSK4","Nbread_R848","Nbread_IAV","Pct_cryptic_NS","Pct_cryptic_LPS","Pct_cryptic_PAM3CSK4","Pct_cryptic_R848","Pct_cryptic_IAV"))
Coding_site_melted[,metric:=gsub('(Nbread_cryptic|Nbread|Pct_cryptic)_(NS|LPS|PAM3CSK4|R848|IAV)','\\1',Coding_site_melted$variable)]
Coding_site_melted[,condition:=gsub('(Nbread_cryptic|Nbread|Pct_cryptic)_(NS|LPS|PAM3CSK4|R848|IAV)','\\2',Coding_site_melted$variable)]
Coding_site_melted[,variable:=NULL]


###### CAST the data (one line per splice site/condition)
Coding_site_perCond= dcast(Coding_site_melted, site_id + condition + Nbread_cryptic_NS + Nbread_NS + Pct_cryptic_NS ~ metric,value.var='value')
Coding_site_perCond[,condition:=factor(condition,levels=c('NS','LPS','PAM3CSK4','R848','IAV'))]

# add infos on splice site 
SpliceSites_info=SpliceSites[,c('site_id','overlappingGene_withStrand','isImmune_overlappingGene_withStrand','overlappingSymbol_withStrand',
                                'maxlFC_overlappingGene_withStrand','log2FC_LPS_overlappingGene_withStrand','log2FC_PAM3CSK4_overlappingGene_withStrand','log2FC_R848_overlappingGene_withStrand','log2FC_IAV_overlappingGene_withStrand',
                                'maxFPKM_overlappingGene_withStrand','FPKM_LPS_overlappingGene_withStrand','FPKM_PAM3CSK4_overlappingGene_withStrand','FPKM_R848_overlappingGene_withStrand','FPKM_IAV_overlappingGene_withStrand','FPKM_NS_overlappingGene_withStrand')]
Coding_site_perCond=merge(Coding_site_perCond,SpliceSites_info)

Coding_site_perCond[,log2FC:=ifelse(condition=='LPS',log2FC_LPS_overlappingGene_withStrand,
                            ifelse(condition=='PAM3CSK4',log2FC_PAM3CSK4_overlappingGene_withStrand,
                            ifelse(condition=='R848',log2FC_R848_overlappingGene_withStrand,
                            ifelse(condition=='IAV',log2FC_IAV_overlappingGene_withStrand,NA))))]
Coding_site_perCond[,FPKM:=ifelse(condition=='NS',FPKM_NS_overlappingGene_withStrand,
                            ifelse(condition=='LPS',FPKM_LPS_overlappingGene_withStrand,
                            ifelse(condition=='PAM3CSK4',FPKM_PAM3CSK4_overlappingGene_withStrand,
                            ifelse(condition=='R848',FPKM_R848_overlappingGene_withStrand,
                            ifelse(condition=='IAV',FPKM_IAV_overlappingGene_withStrand,NA)))))]
Coding_site_perCond[,FoldChange:=cut(log2FC,c(-Inf,-1,1,Inf),labels=c('< -1','[-1,1]','>1'))]

######  Aggregate per gene 
NoisySplicing_perGene=Coding_site_perCond[,.(Rate_Noisy_Splicing_perIntron=mean(Pct_cryptic,na.rm=T),total_Pct_Noisy_Splicing=100-prod(1-Pct_cryptic/100,na.rm=T)*100,FPKM=max(FPKM,na.rm=T),log2FC=mean(log2FC,na.rm=T)),by=.(condition,overlappingGene_withStrand,overlappingSymbol_withStrand)]
NoisySplicing_perGene=NoisySplicing_perGene[!is.na(Rate_Noisy_Splicing_perIntron)]

###########################################################################################
#####       Annotate each gene for Number of exons, mean intron length, etc...        #####
###########################################################################################

introns_df=read.table(file=sprintf('%s/Annotation/IntronAnnotation_hg37_ens70_detectedERC.txt',HOME),header=T) # all introns of expressed transcript (FPKM>1, add estimated by cuffDiff)

TranscriptAnnot=read.table(paste(HOME,'/Annotation/TranscriptAnnotation_hg37_ens70.txt',sep=''),sep='\t',comment='',quote='',stringsAsFactors=FALSE)
TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='X']=23
TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='Y']=24
TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='MT']=26
TranscriptAnnot$Chromosome.Name=as.numeric(TranscriptAnnot$Chromosome.Name)

By=function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x}
mean_NbExons=By(TranscriptAnnot$Nb_exons,TranscriptAnnot$Ensembl.Gene.ID,mean)
max_NbExons=By(TranscriptAnnot$Nb_exons,TranscriptAnnot$Ensembl.Gene.ID,max)

TranscriptAnnot=TranscriptAnnot[order(TranscriptAnnot$Ensembl.Gene.ID,-TranscriptAnnot$Nb_exons),]
# For each gene, get the transcript with the largest number of exons
transcript_maxExon=TranscriptAnnot[!duplicated(TranscriptAnnot$Ensembl.Gene.ID),'Ensembl.Transcript.ID'] 
# introns are kept if they are present in the transcript with the largets number of exons
keep_intron=sapply(strsplit(introns_df$transcripts,'+',fixed=T),function(x){any(x%in%transcript_maxExon)}) 
MeanIntronLength_maxExon=By(introns_df$end[keep_intron]-introns_df$start[keep_intron],introns_df$gene_id[keep_intron],mean)

### add these informations to the Noisy Splicing datasetr
NoisySplicing_perGene$max_NbExons = max_NbExons[match(NoisySplicing_perGene$overlappingGene_withStrand,names(max_NbExons))]
NoisySplicing_perGene$MeanIntronLength=MeanIntronLength_maxExon[match(NoisySplicing_perGene$overlappingGene_withStrand,names(MeanIntronLength_maxExon))]
NoisySplicing_perGene=NoisySplicing_perGene[!is.na(NoisySplicing_perGene$max_NbExons),]
PATH=sprintf("%s/Maxime/Splicing/HISAT2/Results/",EVO_IMMUNO_POP)

# load information on intronic read counts at basal state
ByGene_Introns=fread(sprintf('%s/AllSamples_Intronic_Gene_count.txt',PATH))
ByGene_Introns=ByGene_Introns[-1,]
ByGene_Introns_cond=t(apply(as.data.frame(ByGene_Introns)[,grepl('RPKM',colnames(ByGene_Introns))],1,By,SampleAnnot$cond,mean,na.rm=T))
rownames(ByGene_Introns_cond)=ByGene_Introns$Gene
colnames(ByGene_Introns_cond)=condIndex

# Melt table (one line per gene/condition
ByGene_Introns_cond=melt(as.data.table(ByGene_Introns_cond,keep.rownames=T),id.vars='rn',variable.name='condition',value.name='Intron_RPKM')

# merge Intron read data with exon read data
NoisySplicing_perGene_withIntron=merge(NoisySplicing_perGene,ByGene_Introns_cond,by.x=c('overlappingGene_withStrand','condition'),by.y=c('rn','condition'))
NbCond=table(NoisySplicing_perGene_withIntron$overlappingGene_withStrand)
# filter out genes for which we could not compute Splicing noise in some conditions (eg. constitutive splice site is not active in at least one condition)
NoisySplicing_perGene_withIntron=NoisySplicing_perGene_withIntron[overlappingGene_withStrand%in%names(NbCond)[NbCond==5],]


#### cast table (one line per gene)
NoisySplicing_perGene_withIntron_cast=dcast(NoisySplicing_perGene_withIntron,overlappingGene_withStrand+overlappingSymbol_withStrand+max_NbExons+MeanIntronLength~condition,value.var=c('Rate_Noisy_Splicing_perIntron','FPKM','log2FC','Intron_RPKM'))
NoisySplicing_perGene_withIntron_cast[,log2FC_NS:=NULL]

# add infos on immune and max logFC
NoisySplicing_perGene_withIntron_cast[,is_immune:=ifelse(overlappingGene_withStrand%in%GoToTest_genes$immuneResponse,'yes','')]
NoisySplicing_perGene_withIntron_cast[,max_logFC:=pmax(log2FC_LPS,log2FC_PAM3CSK4,log2FC_R848,log2FC_IAV)]

write.table(NoisySplicing_perGene_withIntron,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/TableS2D.txt',HOME),row.name=F)
write.table(NoisySplicing_perGene_withIntron_cast,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/TableS2D_cast.txt',HOME),row.name=F)

########################################################
####     estimate pct of Noisy splicing Per sample  #### 
########################################################

Junc=fread(sprintf('%s/Maxime/Splicing/HISAT2/Results/All_junctions_bySample.txt',evo_immuno_pop))

Coding_siteByGene=list()
Pct_cryptic_sample=c()
for (samp in 1:970){
    cat(samp,'')
#    Junc_samp=do.call(rbind, Junc[paste(1:nrow(PosFile),samp)])
#    mm=match(paste(Junc_samp$start_id,Junc_samp$end_id),paste(Junc_annot$start_id,Junc_annot$end_id))
    Junc_samp=Junc[sample==SampleAnnot$sample_ID[samp],]
    Junc_samp[,gene:=Junc_annot$overlappingGene_withStrand[mm]]
    Junc_samp[,symbol:=Junc_annot$overlappingSymbol_withStrand[mm]]
    Junc_samp[,is_cryptic_end:=ifelse(end_id%in%cryptic_site,1,0)]
    Junc_samp[,is_cryptic_start:=ifelse(start_id%in%cryptic_site,1,0)]
    codingStart=Junc_samp[start_id%in%gsub(' (acceptor|donor)','',Const_coding_site),
                                     .(Nbread_cryptic=sum(Nb * is_cryptic_end), Nbread_total=sum(Nb)),by=.(start_id,gene,symbol,ind)]
    codingStart[is.na(codingStart)]=0
    codingStart[,site_id:=paste(start_id,ifelse(grepl('\\-',start_id),'acceptor','donor'))]
    codingStart[,start_id:=NULL]
    codingStart[,Pct_cryptic:=Nbread_cryptic/Nbread_total*100]

    codingEnd=Junc_samp[end_id%in%gsub(' (acceptor|donor)','',Const_coding_site),
                                     .(Nbread_cryptic=sum(Nb * is_cryptic_start), Nbread_total=sum(Nb)),by=.(end_id,gene,symbol,ind)]
    codingEnd[is.na(codingEnd)]=0
    codingEnd[,site_id:=paste(end_id,ifelse(grepl('\\+',end_id),'acceptor','donor'))]
    codingEnd[,end_id:=NULL]
    codingEnd[,Pct_cryptic:=Nbread_cryptic/Nbread_total*100]

    Coding_site=rbind(codingEnd,codingStart)
    # security check. removes cases where cryptic acceptor/overlaps a constitutive donor
    Coding_site=Coding_site[site_id%in%Const_coding_site,]
    
    Coding_siteByGene[[unique(Junc_samp$ind)]]=Coding_site[,.(Rate_Noisy_Splicing_perIntron=mean(Pct_cryptic,na.rm=T),total_Pct_Noisy_Splicing=100-prod(1-Pct_cryptic/100,na.rm=T)*100),by=.(gene,symbol,ind)]
    Pct_cryptic_sample[unique(Junc_samp$ind)]=mean(Coding_site$Pct_cryptic)
}
rm(Junc);gc()

Coding_siteByGene=rbindlist(Coding_siteByGene)

write.table(Coding_siteByGene,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/TableS2D_perInd.txt',HOME))
write.table(Pct_cryptic_sample,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/Pct_cryptic_sample_perInd.txt',HOME))


###################################################
##            Make Figures and number            ##
###################################################

NoisySplicing_perGene=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/TableS2D.txt',HOME))
TableS2D=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/TableS2D_cast.txt',HOME))

#############################################################
####         Analysis of mean rate of noisy splicing     ####
#############################################################

# melt (one line per gene/condition/metric)
TableS2D_melted=melt(TableS2D,id.vars=c('overlappingGene_withStrand','Rate_Noisy_Splicing_perIntron_NS','FPKM_NS','Intron_RPKM_NS','max_NbExons','MeanIntronLength','max_logFC','is_immune'),variable.name='variable',
measure=c("Rate_Noisy_Splicing_perIntron_NS","Rate_Noisy_Splicing_perIntron_LPS","Rate_Noisy_Splicing_perIntron_PAM3CSK4","Rate_Noisy_Splicing_perIntron_R848","Rate_Noisy_Splicing_perIntron_IAV","FPKM_NS","FPKM_LPS","FPKM_PAM3CSK4","FPKM_R848","FPKM_IAV","Intron_RPKM_NS","Intron_RPKM_LPS","Intron_RPKM_PAM3CSK4","Intron_RPKM_R848","Intron_RPKM_IAV",'log2FC_LPS','log2FC_PAM3CSK4','log2FC_R848','log2FC_IAV'))

TableS2D_melted[,metric:=gsub('(Rate_Noisy_Splicing_perIntron|FPKM|Intron_RPKM|log2FC)_(NS|LPS|PAM3CSK4|R848|IAV)','\\1',variable)]
TableS2D_melted[,condition:=gsub('(Rate_Noisy_Splicing_perIntron|FPKM|Intron_RPKM|log2FC)_(NS|LPS|PAM3CSK4|R848|IAV)','\\2',variable)]
TableS2D_melted[,variable:=NULL]
TableS2D_melted[,condition:=factor(condition,levels=c('NS','LPS','PAM3CSK4','R848','IAV'))]

# cast (one line per gene) keeping NS information for each gene 
TableS2D_diff= dcast(TableS2D_melted, overlappingGene_withStrand + condition + Rate_Noisy_Splicing_perIntron_NS + FPKM_NS + Intron_RPKM_NS + max_NbExons + MeanIntronLength + max_logFC + is_immune ~ metric ,value.var='value')
TableS2D_diff[,FoldChange:=cut(log2FC,c(-Inf,-1,1,Inf),labels=c('< -1','[-1,1]','>1'))]

############################
####      Figure 3B     ####
############################

hist(log10(pmax(1e-5,TableS2D$Rate_Noisy_Splicing_perIntron_NS)),col='grey',axes=F,main='',xlab='',br=50,xlim=c(-5,2))
axis(1,at=c(-4.95,-3,-2,-1,0,1,2),labels=c('0','0.001','0.01','0.1','1','10','100'))
axis(2,las=1,at=c(0,200,400,600))

############################
####      Figure 3C     ####
############################
# melt (one line per gene/condition)
melted=melt(TableS2D,id='overlappingGene_withStrand',measure=c('Rate_Noisy_Splicing_perIntron_NS','Rate_Noisy_Splicing_perIntron_LPS','Rate_Noisy_Splicing_perIntron_PAM3CSK4','Rate_Noisy_Splicing_perIntron_R848','Rate_Noisy_Splicing_perIntron_IAV'))
melted$variable=substr(melted$variable,31,40)
colnames(melted)=c('gene','condition','noisySplicing')

melted[,noisySplicing:=log10(pmax(1e-5,noisySplicing))]
melted[,condition:=factor(condition,c('NS','LPS','PAM3CSK4','R848','IAV'))]

boxplot(noisySplicing~condition,data=melted,notch=T,outcol='#12345600',col=colERC5,axes=F,ylim=c(-5,2))
axis(1,las=2,at=1:5,labels=condIndex)
axis(2,las=2,at=c(-5,-3,-2,-1,0,1,2),labels=rep('',7))


#####################################
###         Figure 3D           #####
#####################################


# TODO: add this figure 



#####################################
###         Figure 3E           #####
#####################################

library(pcaMethods)
NMD=-pca(scale(t(FPKM_gene[(GoToTest_genes$nonsense_mediated_decay)[grep('UPF|SMG', G2S(GoToTest_genes$nonsense_mediated_decay))],])))@scores[,1]
cor.test(NMD,Pct_cryptic_sample[names(NMD)])

#t = -42.893, df = 968, p-value < 4.513938e-226
#cor=-0.8094731 
#95 percent confidence interval: [-0.8301215 -0.7866076]

plot(NMD,bg=colERC[SampleAnnot$colSamp],Pct_cryptic_sample[names(NMD)] ,pch=21,cex=0.8,las=1,xlab='PC1 of expression of NMD genes\n (UPF* & SMG*)',ylab='Average % of splicing noise per site')


####################################
###     Figure S4A              ####
####################################

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/Fig_boxplot_Change_in_nb_of_cryptic_site_according_to_FC.pdf',HOME),height=3,width=9)
p <- ggplot(Nb_cryptic_SpliceSites_perCond,aes(x=FoldChange,y=DeltaCryptic,fill=FoldChange))+geom_boxplot(outlier.shape=NA,notch=TRUE)+ylim(c(-25,70))+theme_light()+facet_grid(~condition)+scale_fill_manual(values=FoldChange_col,name='Fold change')
p+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=.5)) + labs(y='Number of cryptic \nsplice site gained/lost',x='Fold change')
dev.off()


### significance analysis for Figure S2A
TestShift=Nb_cryptic_SpliceSites_perCond[condition!='NS' ,.(shift=median(DeltaCryptic,na.rm=T),shiftm=mean(DeltaCryptic,na.rm=T),P=wilcox.test(DeltaCryptic)$p.value),by=.(FoldChange,condition)]
TestShift[order(FoldChange,condition),]
#    FoldChange condition shift     shiftm             P
# 1:       < -1       LPS    -5 -7.0107383 1.114235e-100
# 2:       < -1  PAM3CSK4    -4 -5.4789311 1.059012e-108
# 3:       < -1      R848    -4 -4.8416988 1.466518e-211
# 4:       < -1       IAV    -2 -2.6989441 3.128249e-146
# 5:     [-1,1]       LPS     0 -0.1068534  1.066436e-15
# 6:     [-1,1]  PAM3CSK4     1  2.1976492 4.424230e-174
# 7:     [-1,1]      R848     1  2.7553139 3.453612e-164
# 8:     [-1,1]       IAV     2  4.9615140  0.000000e+00
# 9:         >1       LPS     9 12.8990291  3.152340e-79
#10:         >1  PAM3CSK4    10 15.9348199  1.327733e-90
#11:         >1      R848    15 20.1670673 1.636218e-131
#12:         >1       IAV    13 18.3256458 7.710660e-168



####################################
###     Figure S4B              ####
####################################
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/Fig_boxplot_change_in_rate_of_noisySplicing_according_to_FC_condition.pdf',HOME),height=3.5,width=10)
p <- ggplot(TableS2D_diff[condition!='NS',],aes(x=FoldChange,y=Rate_Noisy_Splicing_perIntron-Rate_Noisy_Splicing_perIntron_NS,fill=FoldChange))+geom_boxplot(outlier.shape=NA,notch=TRUE)+coord_cartesian(ylim=c(-0.4,0.7))+theme_light()+facet_grid(~condition)+scale_fill_manual(values=FoldChange_col,name='Fold Change')
p+ theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.5)) + labs(y='Change in % of Noisy splicing',x='Condition')+geom_hline(aes(yintercept=0),col='darkgreen')
dev.off()


### significance analysis for Figure S2B
TestShift=TableS2D_diff[condition!='NS' ,.(shift=median(Rate_Noisy_Splicing_perIntron-Rate_Noisy_Splicing_perIntron_NS,na.rm=T),P=wilcox.test(Rate_Noisy_Splicing_perIntron-Rate_Noisy_Splicing_perIntron_NS)$p.value),by=.(FoldChange,condition)]
TestShift[order(FoldChange,condition),]
#    FoldChange condition      shift             P
# 1:       < -1       LPS 0.03501537  5.191510e-63
# 2:       < -1  PAM3CSK4 0.04272885  1.101506e-95
# 3:       < -1      R848 0.05923776 3.405423e-236
# 4:       < -1       IAV 0.06239839  0.000000e+00
# 5:     [-1,1]       LPS 0.01838274  0.000000e+00
# 6:     [-1,1]  PAM3CSK4 0.02814597  0.000000e+00
# 7:     [-1,1]      R848 0.04337087  0.000000e+00
# 8:     [-1,1]       IAV 0.06350232  0.000000e+00
# 9:         >1       LPS 0.01744393  2.775732e-13
#10:         >1  PAM3CSK4 0.03123594  4.241589e-18
#11:         >1      R848 0.06411904  2.551466e-43
#12:         >1       IAV 0.08718479  9.385936e-44

##############################################################
############## Additional Figures For reviewers  #############
##############################################################


#### changes in rate of noisy splicing upon stimulation: measured at splice site level

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/Fig_boxplot_Nb_of_cryptic_splice_site_reads_according_to_condition.pdf',HOME),height=4.5,width=4.5)
p <- ggplot(Coding_site_perCond,aes(x=condition,y=Nbread_cryptic,fill=condition))+geom_boxplot(outlier.shape=NA,notch=TRUE)+coord_cartesian(ylim=c(0,2))+theme_light()+scale_fill_manual(values=colERC5,name='condition')
p+ theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.5)) + labs(y='% of Noisy splicing',x='Condition')
dev.off()


###################################################################################
###   correlation between gene expression, Intronic reads and Noisy splicing    ###
###################################################################################

w=which(TableS2D_diff$Intron_RPKM_NS/TableS2D_diff$FPKM_NS*100>0 & TableS2D_diff$cond=='NS')
layout(matrix(1:3,1))
COR=cor.test(TableS2D_diff$FPKM_NS[w]*100,TableS2D_diff$Rate_Noisy_Splicing_perIntron_NS[w],method='s')
smoothScatter(log10(TableS2D_diff$FPKM_NS[w]*100),log10(TableS2D_diff$Rate_Noisy_Splicing_perIntron_NS)[w],ylab='% splicing noise',xlab='FPKM',sub=paste('rho=',round(COR$est,2),', P=',COR$p.value))
COR=cor.test(TableS2D_diff$Intron_RPKM_NS[w]*100,TableS2D_diff$Rate_Noisy_Splicing_perIntron_NS[w],method='s')
smoothScatter(log10(TableS2D_diff$Intron_RPKM_NS[w]*100),log10(TableS2D_diff$Rate_Noisy_Splicing_perIntron_NS)[w],ylab='% splicing noise',xlab='Intron RPKM',sub=paste('rho=',round(COR$est,2),', P=',COR$p.value))
COR=cor.test(TableS2D_diff$Intron_RPKM_NS[w]/TableS2D_diff$FPKM_NS[w]*100,TableS2D_diff$Rate_Noisy_Splicing_perIntron_NS[w],method='s')
smoothScatter(log10(TableS2D_diff$Intron_RPKM_NS[w]/TableS2D_diff$FPKM_NS[w]*100),log10(TableS2D_diff$Rate_Noisy_Splicing_perIntron_NS)[w],ylab='% splicing noise',xlab='Intron/Exon Ratio',sub=paste('rho=',round(COR$est,2),', P=',COR$p.value))
