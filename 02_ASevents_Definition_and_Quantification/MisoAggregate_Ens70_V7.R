library(data.table)

# make PSI Matrix for Ens70 events (adding support)
# run SUPPA event_generator on gtf file to generate splicing events
# run MISO to quantify splicing events.
# aggregate output from SUPPA

# for each type of events, aggregate summary data from all samples with event annotation
for (type in c('A3','A5','SE','MX','RI','AF','AL')){
	cat('\n',type)
	DATA=read.table(sprintf('%s/Maxime/Splicing/MISO/hg19_ENS70/SplicingEvents_%s_strict.gtf',EVO_IMMUNO_POP,type),skip=1,sep='\t')
	# split last column
	colnames(DATA)=c('chrom','event_type','feature_type','start','end','score0','strand','score1','annot')
	y=sapply(strsplit(DATA$annot,';'),function(x){n=gsub('(.*) (.*);?','\\1',x);y=gsub('(.*) (.*);?','\\2',x);names(y)=n;y})
	if(is.matrix(y)){
		Y=t(y)
		}else{
		nn=unique(unlist(lapply(y,names)))
		Y=t(sapply(y,function(x){x[match(nn,names(x))]}))
		colnames(Y)=nn
		}
	DATA=cbind(DATA[-ncol(DATA)],Y)
	event_start=By(DATA$start,DATA$gene_id,min)
	event_end=By(DATA$end,DATA$gene_id,max)
	Annot_IOE=read.table(sprintf('%s/Maxime/Splicing/MISO/hg19_ENS70/SplicingEvents_%s_strict.ioe',EVO_IMMUNO_POP,type),skip=0,header=T)
	Annot_IOE$symbol=G2S(Annot_IOE$gene_id)
	Annot_IOE$ID=Make.unique(paste(Annot_IOE$symbol,type,sep=':'))
	Annot_IOE$start=event_start[match(substr(Annot_IOE$event_id,17,1000),names(event_start))]
	Annot_IOE$end=event_end[match(substr(Annot_IOE$event_id,17,1000),names(event_end))]

	for(ID in SampleAnnot$sample_ID){
		cat('.')
		samp1=read.table(sprintf('%s/Maxime/Splicing/MISO/summary/%s/Ens70_HISAT/summary/%s.miso_summary',EVO_IMMUNO_POP,ID,type),row.names=NULL)
		colnames(samp1)=c(cn(samp1)[2:5],'isoformA','isoformB',cn(samp1)[7:12])
		Annot_IOE[[paste(ID,'PSIpost',sep='_')]]=samp1$miso_posterior_mean[match(substr(Annot_IOE$event_id,17,1000),samp1$event_name)]
		Annot_IOE[[paste(ID,'CIrange',sep='_')]]=samp1$ci_high[match(substr(Annot_IOE$event_id,17,1000),samp1$event_name)]-samp1$ci_low[match(substr(Annot_IOE$event_id,17,1000),samp1$event_name)]
		NbAlt=gsub('0:([0-9]+)(,1:([0-9]+))?','\\3',samp1$assigned_counts)
		NbAlt=ifelse(NbAlt=='',0,as.numeric(NbAlt))
		NbRef=as.numeric(gsub('0:([0-9]+)(,1:([0-9]+))?','\\1',samp1$assigned_counts))
		Annot_IOE[[paste(ID,'Count',sep='_')]]=(NbRef+NbAlt)[match(substr(Annot_IOE$event_id,17,1000),samp1$event_name)]
		}
	write.table(Annot_IOE,file=sprintf('%s/Maxime/Splicing/MISO/aggregated/Ens70_HISAT/%s.summary.aggregated_V7.txt',EVO_IMMUNO_POP,type),quote=F,sep='\t',row.names=F)
	}

DATA_FULL=list()
for (type in c('A3','A5','SE','MX','RI','AF','AL')){
	cat('\n',type)
	DATA=read.table(sprintf('%s/Maxime/Splicing/MISO/hg19_ENS70/SplicingEvents_%s_strict.gtf',EVO_IMMUNO_POP,type),skip=1,sep='\t')
	# split last column
	colnames(DATA)=c('chrom','event_type','feature_type','start','end','score0','strand','score1','annot')
	y=sapply(strsplit(DATA$annot,';'),function(x){n=gsub('(.*) (.*);?','\\1',x);y=gsub('(.*) (.*);?','\\2',x);names(y)=n;y})
	if(is.matrix(y)){
		Y=t(y)
		}else{
		nn=unique(unlist(lapply(y,names)))
		Y=t(sapply(y,function(x){x[match(nn,names(x))]}))
		colnames(Y)=nn
		}
	DATA=cbind(DATA[-ncol(DATA)],Y)
	DATA_FULL[[type]]=DATA
	}
DATA_FULL=do.call(rbind,DATA_FULL)
save(DATA_FULL,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/Ens70_HISAT/PSI_events_ALL_coordinates.Rdata',sep=''))

# extract PSI, and annotations 
PSI_FULL=list()
PSI_Annot_FULL=list()
for (type in c('A3','A5','SE','MX','RI','AF','AL')){
	cat(type)
	SUMMARY=read.table(sprintf('%s/Maxime/Splicing/MISO/aggregated/Ens70_HISAT/%s.summary.aggregated_V7.txt',EVO_IMMUNO_POP,type),header=T,quote='',sep='\t')
	SUMMARY$strand=substr(SUMMARY$event_id,nchar(SUMMARY$event_id),nchar(SUMMARY$event_id))
	SUMMARY$event_type=type
	colnames(SUMMARY)[1]='chrom'
	PSI=as.matrix(SUMMARY[,grep('PSIpost',colnames(SUMMARY))])
	PctNA_cond=t(apply(is.na(PSI),1,By,substr(cn(PSI),8,8),mean))
	w=which(apply(PctNA_cond<0.05,1,any))
	MeanPSI_cond=t(apply(PSI,1,By,substr(colnames(PSI),8,8),mean,na.rm=T))
	colnames(PctNA_cond)=paste('PctNA',condIndex,sep='_')
	colnames(MeanPSI_cond)=paste('MeanPSI',condIndex,sep='_')

	annot_gene=SUMMARY[,c('event_id','event_type','chrom','start','end','strand','gene_id','symbol','ID')]

	# add support informations (Gene FPKM)
	Support_cond=GeneAnnot[match(annot_gene$gene_id,GeneAnnot$Ensembl.Gene.ID),grep('_mean',colnames(GeneAnnot))]
	Support_cond[is.na(Support_cond)]=0
	colnames(Support_cond)=paste('Support',condIndex,sep='_')

	PSI_Annot=cbind(annot_gene,PctNA_cond,MeanPSI_cond,Support_cond)
	PSI=PSI[,match(gsub('-','.',SampleAnnot$sample_ID,fixed=T),substr(colnames(PSI),1,8))]
	colnames(PSI)=SampleAnnot$sample_ID
	rownames(PSI)=PSI_Annot$ID
	rownames(PSI_Annot)=PSI_Annot$ID

	save(PSI,PSI_Annot,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/Ens70_HISAT/PSI_events_',type,'_V7.Rdata',sep=''))
	PSI_FULL[[type]]=PSI
	PSI_Annot_FULL[[type]]=PSI_Annot
	}
PSI=do.call(rbind,PSI_FULL)
PSI_Annot=do.call(rbind,PSI_Annot_FULL)


# define frequent alternative splicing events without junction reads
# Testable -> events that are frequent and will be considered for further analyses. 
Testable_cond=PSI_Annot[,grep('PctNA',colnames(PSI_Annot))]<0.05 & (PSI_Annot[,grep('MeanPSI',colnames(PSI_Annot))]>0.05 & PSI_Annot[,grep('MeanPSI',colnames(PSI_Annot))]<0.95) & PSI_Annot[,grep('Support',colnames(PSI_Annot))]> 10
colnames(Testable_cond)=paste('Testable',condIndex,sep='_')
PSI_Annot$NbTestableCond=apply(Testable_cond,1,sum)

########################################################################
# add junction read information to alternative splicing event definition (based on definition of all Junction reads Junction_annot_V7.2.txt
# TODO: recode this bit using data.table and related packages to make it more readable.

#################################################################
###    Annotation of Junctions associated to each AS event    ###
#################################################################


library(data.table)
Junc=as.data.frame(fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Junction_annot_V7.2.txt',HOME))

####  We split event ids to extract all relevant junctions 
event_info=strsplit(substr(PSI_Annot$event_id,17,1000),':')
### flatten event_info list into a vector
event_info_flat=unlist(event_info)
event_ID=rep(PSI_Annot$event_id,sapply(event_info,length)) # make a vector with matched event id

# and other infos: for each event...
gene=rep(sapply(event_info,function(x){x[1]}),sapply(event_info,length)) # ...1st element is gene 
chr=rep(sapply(event_info,function(x){x[2]}),sapply(event_info,length)) # ...2nd element is chr
strand=rep(sapply(event_info,function(x){x[length(x)]}),sapply(event_info,length)) # ...last element is strand

### extract junctions...
event_junc=event_info_flat[grep('[0-9]+-[0-9]+',event_info_flat)]
# ... and matched event ID
event_ID_junc=event_ID[grep('[0-9]+-[0-9]+',event_info_flat)]

# ... add other informations.
Junc_MISO=data.frame(event_id=event_ID_junc,chr=chr[grep('[0-9]+-[0-9]+',event_info_flat)],start=as.numeric(gsub('([0-9]+)-([0-9]+)','\\1',event_junc)),end=as.numeric(gsub('([0-9]+)-([0-9]+)','\\2',event_junc)),type=substr(event_ID_junc,17,18))

### For skipped exon events, we also need to extract the junction that skips the exon.
# ie: Event ids for SE are encoded as aaaa-bbbb:cccc-dddd. We want to extract junction aaaa-dddd
SE_event_info=event_info[sapply(event_info,function(y){y[1]=='SE'})]
SE_junc=sapply(SE_event_info,function(x){y=strsplit(x[3:4],'-');aaaa=y[[1]][1]; dddd= y[[2]][2]; paste(aaaa,dddd,sep='-')}) 
chr=sapply(SE_event_info,function(x){x[2]})
gene=sapply(SE_event_info,function(x){x[1]})
strand=sapply(SE_event_info,function(x){x[length(x)]})
Junc_MISO_SE=data.frame(event_id=PSI_Annot$event_id[sapply(event_info,function(y){y[1]=='SE'})],chr=chr,start=as.numeric(gsub('([0-9]+)-([0-9]+)','\\1',SE_junc)),end=as.numeric(gsub('([0-9]+)-([0-9]+)','\\2',SE_junc)),type="SE")

Junc_MISO=rbind(Junc_MISO,Junc_MISO_SE)
##### Junc_MISO now contains all junctions from known MISO event that we will annotate and analyse.

# match event junctions to observed junctions 
mm=match(paste(Junc_MISO$chr,Junc_MISO$start,Junc_MISO$end),paste(Junc$chrom_intron,Junc$start_intron-1,Junc$end_intron+1))

mean(is.na(mm)) # 0.223 
# 22% of MISO events have no read associated to their junction (when taking all reads),

# Let's annotate the Junctions from known AS events
Junc_MISO$Coverage=Junc$Total_count[mm]
Junc_MISO$GerpRS_start=Junc$GerpRS_start[mm]
Junc_MISO$GerpRS_end=Junc$GerpRS_end[mm]
Junc_MISO$GerpRS_start_min=pmin(Junc$GerpRS_start[mm],Junc$GerpRS_start2[mm])
Junc_MISO$GerpRS_end_min=pmin(Junc$GerpRS_end[mm],Junc$GerpRS_end2[mm])
Junc_MISO$Coverage_NS=Junc$Total_count_NS[mm]
Junc_MISO$Coverage_LPS=Junc$Total_count_LPS[mm]
Junc_MISO$Coverage_PAM3=Junc$Total_count_PAM3[mm]
Junc_MISO$Coverage_R848=Junc$Total_count_R848[mm]
Junc_MISO$Coverage_IAV=Junc$Total_count_IAV[mm]
# annotate GerpRS and coding status of start/end
Junc_MISO$GerpRS_start=(Junc$GerpRS_start[mm]+Junc$GerpRS_start2[mm])/2
Junc_MISO$GerpRS_end=(Junc$GerpRS_end[mm]+Junc$GerpRS_end2[mm])/2
Junc_MISO$coding_startSite=Junc$coding_startSite[mm]
Junc_MISO$coding_endSite=Junc$coding_endSite[mm]

# assign start/end id for each junction 
Junc_MISO$start_id=paste(Junc_MISO$chr,Junc_MISO$start+1,Junc_MISO$strand,ifelse(Junc_MISO$strand=='+','donor','acceptor'))
Junc_MISO$end_id=paste(Junc_MISO$chr,Junc_MISO$end-1,Junc_MISO$strand,ifelse(Junc_MISO$strand=='+','acceptor','donor'))

# For each junction, annotate the transcript in which they are present (1 if junction usage increases with PSI, 0 otherwise)
junc_ID=paste(Junc_MISO[,'start'],Junc_MISO[,'end'],sep='-')

# annotation for A3
first_A3=ifelse(Junc_MISO$strand[Junc_MISO$type=='A3']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='A3'],':'),function(x){x[4]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='A3'],':'),function(x){x[3]}))
strand_A3=Junc_MISO$strand[Junc_MISO$type=='A3']=='+'
Junc_MISO[Junc_MISO$type=='A3','transcript']=ifelse(junc_ID[Junc_MISO$type=='A3']==first_A3,0,1)*ifelse(strand_A3,1,-1)+ifelse(strand_A3,0,1)

# annotation for A5
first_A5=ifelse(Junc_MISO$strand[Junc_MISO$type=='A5']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='A5'],':'),function(x){x[4]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='A5'],':'),function(x){x[3]}))
strand_A5=Junc_MISO$strand[Junc_MISO$type=='A5']=='+'
Junc_MISO[Junc_MISO$type=='A5','transcript']=ifelse(junc_ID[Junc_MISO$type=='A5']==first_A5,0,1)*ifelse(strand_A5,1,-1)+ifelse(strand_A5,0,1)

# annotation for AF
first_AF=ifelse(Junc_MISO$strand[Junc_MISO$type=='AF']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='AF'],':'),function(x){x[4]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='AF'],':'),function(x){x[3]}))
strand_AF=Junc_MISO$strand[Junc_MISO$type=='AF']=='+'
Junc_MISO[Junc_MISO$type=='AF','transcript']=ifelse(junc_ID[Junc_MISO$type=='AF']==first_AF,1,0)*ifelse(strand_AF,1,-1)+ifelse(strand_AF,0,1)

# annotation for AL
first_AL=ifelse(Junc_MISO$strand[Junc_MISO$type=='AL']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='AL'],':'),function(x){x[3]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='AL'],':'),function(x){x[4]}))
strand_AL=Junc_MISO$strand[Junc_MISO$type=='AL']=='+'
Junc_MISO[Junc_MISO$type=='AL','transcript']=ifelse(junc_ID[Junc_MISO$type=='AL']==first_AL,0,1)*ifelse(strand_AL,1,-1)+ifelse(strand_AL,0,1)

# annotation for RI
Junc_MISO[Junc_MISO$type=='RI','transcript']=0

# annotation for MX
MX_id=paste(Junc_MISO[Junc_MISO$type=='MX','start'],Junc_MISO[Junc_MISO$type=='MX','end'],sep='-')
MX_1=ifelse(Junc_MISO$strand[Junc_MISO$type=='MX']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='MX'],':'),function(x){x[3]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='MX'],':'),function(x){x[5]}))
MX_2=ifelse(Junc_MISO$strand[Junc_MISO$type=='MX']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='MX'],':'),function(x){x[4]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='MX'],':'),function(x){x[6]}))
Junc_MISO[Junc_MISO$type=='MX','transcript']=ifelse(MX_id==MX_1 |MX_id==MX_2,1,0)
table(Junc_MISO$type, Junc_MISO$transcript)

# annotation for SE
Junc_MISO[Junc_MISO$type=='SE','transcript']=ifelse(junc_ID[Junc_MISO$type=='SE']%in%SE_junc,0,1)


###########################################################
###   filtering of low quality and redundant AS events  ###
###########################################################


######### extract coverage of each junction  (to )
Junc_coverage=as.matrix(Junc_MISO[,grep('Coverage_',colnames(Junc_MISO))])
# remove events where some junctions are covered by less than 30 leads 
cov_target=30
MISO_events_Junction_observed=apply(Junc_coverage,2,function(x){
    low_Coverage_junctions=x<cov_target | is.na(x)
    events_with_low_Coverage_junction=Junc_MISO[which(low_Coverage_junctions),'event_id']
    event_has_all_junctions_covered = !PSI_Annot$event_id%in%events_with_low_Coverage_junction
    event_has_all_junctions_covered })

By = function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x}
Total_junc_coverage_event=apply(Junc_coverage,2,By,Junc_MISO$event_id,sum,na.rm=T)
#sum(apply(MISO_events_Junction_observed,1,any))

###### redefine frequent alternative splicing events based on Junction read counts>30
Testable_cond=PSI_Annot[,grep('PctNA',colnames(PSI_Annot))]<0.05 & (PSI_Annot[,grep('MeanPSI',colnames(PSI_Annot))]>0.05 & PSI_Annot[,grep('MeanPSI',colnames(PSI_Annot))]<0.95) & PSI_Annot[,grep('Support',colnames(PSI_Annot))]> 10 & MISO_events_Junction_observed
colnames(Testable_cond)=paste('Testable',condIndex,sep='_')
PSI_Annot$NbTestableCond=apply(Testable_cond,1,sum)
for (i in 1:5){
	PSI_Annot[,grep('Testable',colnames(PSI_Annot))[i]]=Testable_cond[,i]
	}
Testable=apply(PSI_Annot[,grep('Testable',colnames(PSI_Annot))],1,any)

PSI_Annot$JuncCovered=apply(MISO_events_Junction_observed,1,any)
for (i in 1:5){
	PSI_Annot[,paste('JuncCovered',condIndex[i],sep='_')]=MISO_events_Junction_observed[,i]
	}
PSI_Annot$Total_junc_coverage=apply(Total_junc_coverage_event,1,sum,na.rm=T)[match(PSI_Annot$event_id,rn(Total_junc_coverage_event))]
for (i in 1:5){
	PSI_Annot[,paste('Total_junc_coverage',condIndex[i],sep='_')]=Total_junc_coverage_event[match(PSI_Annot$event_id,rn(Total_junc_coverage_event)),i]
	}

# extract sets of overlapping events 
PSI_Annot_testable=PSI_Annot[PSI_Annot$NbTestableCond>0,] # Junc covered is now included in the Testable set of events.
library(GenomicRanges)
GRPSI=makeGRangesFromDataFrame(PSI_Annot_testable)
oo=findOverlaps(GRPSI,GRPSI)
library(igraph)
gg=graph_from_edgelist(cbind(subjectHits(oo),queryHits(oo)), directed = FALSE)
clu <- components(gg)
#S2N=PSI_Annot_ok[,grep('Signal2Noise',colnames(PSI_Annot_ok))]
PSI_Annot_testable$EventCluster=clu$membership
PSI_Annot$OverlappingEvent_Cluster=PSI_Annot_testable$EventCluster[match(PSI_Annot$event_id,PSI_Annot_testable$event_id)]


# clustering algorithm walktrap
makeClusters=function(i,cor_th=0.8,abs=T){
	require(Hmisc)
	require(igraph)	
	w=which(PSI_Annot$OverlappingEvent_Cluster == i)
	if(length(w)>1){
		COR=rcorr(t(PSI[w,]))$r
		if(abs){
			cor_graph=graph_from_adjacency_matrix(abs(COR)>cor_th)
		}else{
			cor_graph=graph_from_adjacency_matrix(COR>cor_th)	
		}
		wc <- cluster_walktrap(cor_graph)
		res=paste(i,membership(wc))
		}else{
		if(length(w)==1){
			res=i
		}else{
			res=NULL
		}
	}
	cbind(w,res)
}


# compute correlations among each cluster of overlapping (testable event)
library(Hmisc)
mycor=c()
for (i in unique(PSI_Annot$OverlappingEvent_Cluster)){
	w=which(PSI_Annot$OverlappingEvent_Cluster == i)
	if(length(w)>1){
		COR=rcorr(t(PSI[w,]))$r
		mycor=c(mycor,as.numeric(COR[upper.tri(COR)]))
		}
	}



withinClustCOR=lapply(unique(PSI_Annot$OverlappingEvent_Cluster),makeClusters,cor_th=0.5)

# assign events to clusters
correlated_event_cluster=rep(NA,length(PSI_Annot$OverlappingEvent_Cluster))
for(i in 1:length(withinClustCOR)){
	correlated_event_cluster[as.numeric(withinClustCOR[[i]][,1])]=as.character(withinClustCOR[[i]][,2])
}

# choose representative event to keep for further analysis (highest Total Junction Count)

PSI_Annot$correlated_event_cluster=correlated_event_cluster
PSI_Annot=PSI_Annot[which(!is.na(PSI_Annot$correlated_event_cluster)),]
PSI_Annot=PSI_Annot[order(PSI_Annot$correlated_event_cluster,ifelse(is.na(PSI_Annot$Total_junc_coverage),1,0),-PSI_Annot$Total_junc_coverage,ifelse(PSI_Annot$gene_id %in% rownames(FPKM_gene),0,1)),]
toKeep=which(PSI_Annot$event_id%in%PSI_Annot$event_id[!duplicated(PSI_Annot$correlated_event_cluster) & PSI_Annot$NbTestableCond>0])

# add info on Correlation based clusters 
PSI_Annot$correlated_event_cluster=gsub(' ','_',PSI_Annot$correlated_event_cluster)
length(unique(PSI_Annot$correlated_event_cluster[toKeep])) # 16173
PSI_Annot$representative_event=PSI_Annot$event_id[toKeep][match(PSI_Annot$correlated_event_cluster,PSI_Annot$correlated_event_cluster[toKeep])]



###########################################################
###    Annotation of coding consequences of AS events   ###
###########################################################

Source='Ens70_HISAT'
library(data.table)

transcript_type_94=fread('/Volumes/@Home/Annotation/Ens94_Transcripts_GRCh37.txt')
Exons_94=fread('/Volumes/@Home/Annotation/Ens94_Exons_GRCh37.txt')
Exons_94=merge(Exons_94,transcript_type_94)
Exons_94=Exons_94[Exons_94$Chromosome%in%c(1:22,'X','Y','MT'),]
Exons_94$Strand=ifelse(Exons_94$Strand>0,'+','-')

# create SpliceSites_MISO containing all splice sites involved in frequent AS events (quantified with MISO)
SpliceSite_MISO_start=data.frame(site_id=Junc_MISO$start_id,event_id=Junc_MISO$event_id,Gerp_site=Junc_MISO$GerpRS_start,transcript=Junc_MISO$transcript,strand=Junc_MISO$strand,type=ifelse(Junc_MISO$strand=='+','donor','acceptor'),event_type=Junc_MISO$type,coding_site=Junc_MISO$coding_startSite)
SpliceSite_MISO_end=data.frame(site_id=Junc_MISO$end_id,event_id=Junc_MISO$event_id,Gerp_site=Junc_MISO$GerpRS_end,transcript=Junc_MISO$transcript,strand=Junc_MISO$strand,type=ifelse(Junc_MISO$strand=='-','donor','acceptor'),event_type=Junc_MISO$type,coding_site=Junc_MISO$coding_endSite)
SpliceSite_MISO=rbind(SpliceSite_MISO_start,SpliceSite_MISO_end)

SpliceSite_MISO$pos=gsub('([0-9XYMT]+ [0-9]+) [+-] (donor|acceptor)','\\1',SpliceSite_MISO$site_id)
Exons_94$pos_end=paste(Exons_94$'Chromosome/scaffold name',Exons_94$'Exon region end (bp)'+1)
Exons_94$pos_start=paste(Exons_94$'Chromosome/scaffold name',Exons_94$'Exon region start (bp)'-1)

Nb_splice_form=as.data.table(SpliceSite_MISO)[,.(Nb_splice_form_where_found=length(unique(transcript))),by=.(site_id, event_id)]

SpliceSite_MISO=as.data.table(merge(SpliceSite_MISO, Nb_splice_form))

trans_id_end=merge(SpliceSite_MISO,Exons_94[Exons_94$'Transcript type'=='protein_coding',c("pos_end","Transcript stable ID")],by.x='pos',by.y='pos_end')
trans_id_start=merge(SpliceSite_MISO,Exons_94[Exons_94$'Transcript type'=='protein_coding',c("pos_start","Transcript stable ID")],by.x='pos',by.y='pos_start')

SpliceSite_MISO_withCodingTranscript=rbind(trans_id_end,trans_id_start)
colnames(SpliceSite_MISO_withCodingTranscript)[colnames(SpliceSite_MISO_withCodingTranscript)=='Transcript stable ID']='transcript_id'

# for each possible splicing, count the number of sites of each splice form compatible with any given transcript, (splitting between shared and unique)
Nb_compatible_sites_between_AS_and_transcript=SpliceSite_MISO_withCodingTranscript[,.(Nb_sites=length(unique(site_id))),by=.(transcript,event_id,transcript_id,Nb_splice_form_where_found)]
Nb_compatible_sites_between_AS_and_transcript[,Nb_splice_form_where_found:=ifelse(Nb_splice_form_where_found==2,'shared','unique')]
Nb_compatible_site=dcast(Nb_compatible_sites_between_AS_and_transcript,event_id+transcript_id~transcript+Nb_splice_form_where_found,value.var=c('Nb_sites'),fill=0)
# safety check
if(! all(Nb_compatible_site$"1_shared"==Nb_compatible_site$"0_shared")){stop('the number of splice site of the splice form should be identical for sites that are shared')}
Nb_compatible_site[,'1_shared':=NULL]
setnames(Nb_compatible_site,"0_shared","Nb_sites_shared")
setnames(Nb_compatible_site,"0_unique","Nb_sites_0_unique")
setnames(Nb_compatible_site,"1_unique","Nb_sites_1_unique")

#Nb_compatible_sites_between_AS_and_transcript=merge(Nb_compatible_sites_between_AS_and_transcript,Nb_site_tot_transcript_split,by=c("transcript","event_id","Nb_splice_form_where_found"))
Nb_site_tot_transcript_split=SpliceSite_MISO[,.(Nb_sites_tot_transcript_shared_or_unique=sum(!duplicated(site_id))),by=.(transcript,event_id,Nb_splice_form_where_found)]
Nb_site_tot_transcript_split[,Nb_splice_form_where_found:=ifelse(Nb_splice_form_where_found==2,'shared','unique')]
Nb_expected_site=dcast(Nb_site_tot_transcript_split,event_id~transcript+Nb_splice_form_where_found,value.var=c('Nb_sites_tot_transcript_shared_or_unique'),fill=0)
# safety check
if(! all(Nb_expected_site$"1_shared"==Nb_expected_site$"0_shared")){stop('the number of splice site shared between splice forms should be the same for each splice form')}
Nb_expected_site[,"1_shared":=NULL]
setnames(Nb_expected_site,"0_shared","total_nb_sites_shared")
setnames(Nb_expected_site,"0_unique","total_nb_sites_0_unique")
setnames(Nb_expected_site,"1_unique","total_nb_sites_1_unique")

# annotate for each protein coding transcript the number of sites observed that are compatible with both forms, form 0 only, or form 1 only, + plus the number expected for shared/form 0/form 1
Nb_compatible_site=merge(Nb_compatible_site,Nb_expected_site)
compatibleTranscripts=Nb_compatible_site[Nb_sites_shared==total_nb_sites_shared & ((Nb_sites_0_unique==total_nb_sites_0_unique & Nb_sites_1_unique==0)|(Nb_sites_1_unique==total_nb_sites_1_unique & Nb_sites_0_unique==0)),]
compatibleTranscripts[,form:=ifelse(Nb_sites_0_unique==total_nb_sites_0_unique & Nb_sites_1_unique==0,0,1)]

compatibleTranscripts=compatibleTranscripts[,c('transcript_id','event_id','form')]
RI_GR=makeGRangesFromDataFrame(Junc_MISO[Junc_MISO$type=='RI',])
Exons_GR=makeGRangesFromDataFrame(Exons_94[Exons_94$'Transcript type'=='protein_coding',],seqnames='Chromosome/scaffold name',start.field='Exon region start (bp)',end.field='Exon region end (bp)',keep.extra=TRUE)
oo=findOverlaps(RI_GR,Exons_GR,type="within")

RI_compatibleTranscripts=data.table(transcript_id=Exons_GR$'Transcript stable ID'[subjectHits(oo)],event_id=Junc_MISO[which(Junc_MISO$type=='RI')[queryHits(oo)],'event_id'],form=1)

compatibleTranscripts=rbind(compatibleTranscripts,RI_compatibleTranscripts)

# add corresponding protein
prot_corresp=Exons_94[,c("Transcript stable ID","Protein stable ID")]
prot_corresp=prot_corresp[!duplicated(prot_corresp),]
compatibleTranscripts=merge(compatibleTranscripts,prot_corresp,by.x='transcript_id',by.y='Transcript stable ID')

NbCoding=dcast(compatibleTranscripts[!duplicated(compatibleTranscripts[,c('Protein stable ID','event_id','form')]),],event_id~form,value.var="Protein stable ID")
ID_Coding=dcast(compatibleTranscripts[!duplicated(compatibleTranscripts[,c('Protein stable ID','event_id','form')]),],event_id~form,value.var="Protein stable ID",fun=paste,collapse='//')
setnames(ID_Coding,c('0','1'),c('form0','form1'))
ID_Coding$type=ifelse(ID_Coding$form0=='' & ID_Coding$form1!='','gain of function',
                ifelse(ID_Coding$form1=='' & ID_Coding$form0!='','loss of function',
                ifelse(ID_Coding$form1==ID_Coding$form0 & ID_Coding$form1!='' & ID_Coding$form0!='','protein unchanged',
                ifelse(ID_Coding$form1!=ID_Coding$form0 & ID_Coding$form1!='' & ID_Coding$form0!='','modified protein',
                ifelse(ID_Coding$form1=='' & ID_Coding$form0=='','no protein',NA)))))

PSI_Annot$coding_type=ID_Coding$type[match(PSI_Annot$event_id,ID_Coding$event_id)]
PSI_Annot$coding_type[intersect(toKeep, which(is.na(PSI_Annot$coding_type)))]='no coding consequence'


###########################################################
###  Annotation of gain/loss of conserved splice sites  ###
###########################################################

coding_type_Event = SpliceSite_MISO[,.(conserved_0_any=any(Gerp_site>2 & Nb_splice_form_where_found!=2 & transcript==0),
                                      conserved_1_any=any(Gerp_site>2 & Nb_splice_form_where_found!=2 & transcript==1),
                                      conserved_0_all=(sum(Gerp_site>2 & Nb_splice_form_where_found!=2 & transcript==0)==sum(Nb_splice_form_where_found!=2 & transcript==0)),
                                      conserved_1_all=(sum(Gerp_site>2 & Nb_splice_form_where_found!=2 & transcript==1)==sum(Nb_splice_form_where_found!=2 & transcript==1)))
                                      , by=.(event_id)]

coding_type_Event$conserved_site_switch=coding_type_Event$conserved_0_any & coding_type_Event$conserved_1_any
coding_type_Event$conserved_site_gain=coding_type_Event$conserved_1_any & !coding_type_Event$conserved_0_any
coding_type_Event$conserved_site_lost=coding_type_Event$conserved_0_any & !coding_type_Event$conserved_1_any
coding_type_Event$conserved_site_switch_all=(coding_type_Event$conserved_0_any & coding_type_Event$conserved_0_all) & (coding_type_Event$conserved_1_any & coding_type_Event$conserved_1_all)
coding_type_Event$conserved_site_gain_all=(coding_type_Event$conserved_1_any & coding_type_Event$conserved_1_all) & !coding_type_Event$conserved_0_any
coding_type_Event$conserved_site_lost_all=(coding_type_Event$conserved_0_any & coding_type_Event$conserved_0_all) & !coding_type_Event$conserved_1_any
coding_type_Event$non_conserved_AS_site = !(coding_type_Event$conserved_1_any | coding_type_Event$conserved_0_any) 

mm=match(PSI_Annot$event_id,coding_type_Event$event_id)

PSI_Annot$conserved_site_GerpRS=ifelse(coding_type_Event$non_conserved_AS_site[mm], 'non conserved Alternative splice sites',
                        ifelse(coding_type_Event$conserved_site_switch_all[mm], 'change in conserved splice site (all)',
                        ifelse(coding_type_Event$conserved_site_switch[mm], 'change in conserved splice site',
                        ifelse(coding_type_Event$conserved_site_gain_all[mm], 'gain of conserved splice site (all)',
                        ifelse(coding_type_Event$conserved_site_gain[mm], 'gain of conserved splice site',
                        ifelse(coding_type_Event$conserved_site_lost_all[mm], 'loss of conserved splice site (all)',
                        ifelse(coding_type_Event$conserved_site_lost[mm], 'loss of conserved splice site',
                        "Other")))))))



##########################################################
### 		data cleaning, batch removal & PCA			##
##########################################################
# run everything ! 
library(pcaMethods)
library(lme4)
# assess all batches
PSI_prov0=PSI[toKeep,]
library(impute)
PSI_prov0=impute.knn(PSI_prov0)$data
tic=Sys.time()

coFactors=c("X260_280", "X260_230", "rin", "total_amount_axeq", 
			"RIN_axeq", "rRNA_ratio_axeq", "GC_pct_final", "Q30_pct",'sum_all_runs','library_mass_conc','FivePrimeBias')

tic=Sys.time()
VarBatch_notAdj=t(apply(logit(PSI_prov0,min=0.001,max=0.999),1,function(x){
		mod=lmer(x~population+condition+X260_280+X260_230+rin+total_amount_axeq+RIN_axeq+rRNA_ratio_axeq+GC_pct_final+Q30_pct+sum_all_runs+library_mass_conc+FivePrimeBias+(1|library_batch)+(1|batch)+(1|lane)+(1|experimentator)+(1|machine)+(1|index)+(1|date),data=SampleAnnot)
		c(sqrt(unlist(summary(mod)$varcor))	,summary(mod)$sigma,summary(mod)$coefficients[-1,3])}))
R2Batch_notAdj=VarBatch_notAdj[,1:7]^2/apply(VarBatch_notAdj[,1:8]^2,1,sum)
R2Cont_notAdj=P2R2(pnorm(abs(VarBatch_notAdj[,-(1:8)]),low=F)*2,ncol(FPKM_gene),1)
toc=Sys.time()
print(toc-tic)
save(VarBatch_notAdj,R2Batch_notAdj,R2Cont_notAdj,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/MISO_events_TechnicalCovariatesQuantification_V5.Rdata',HOME))

#load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/MISO_events_TechnicalCovariatesQuantification_V5.Rdata',HOME)))
#load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/MISO_cleaning/AdjustmentResults/PSI_%s_8221.Rdata',HOME,batch))
dir.create(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/',HOME))
mylogit=function(x,d=0){p=d+(1-2*d)*x; log(p/(1-p))}
PSI_prov=mylogit(PSI_prov0,d=0.001)
for(batch in c('library_batch','machine','date','batch')){
# batch='date'
	BATCH=SampleAnnot[,batch]
	if (any(table(SampleAnnot[,batch])<3)){
		condpop_PSI=apply(PSI_prov,1,By,SampleAnnot$CondPop,mean)
		PSI_scale=PSI_prov-t(condpop_PSI)[,as.character(SampleAnnot$CondPop)]
		PSI_batch=apply(PSI_scale,1,By,SampleAnnot[,batch],mean)
		inSmallBatch=SampleAnnot[,batch]%in%names(table(SampleAnnot[,batch])[table(SampleAnnot[,batch])<3])
		noBatch=SampleAnnot[,batch]==' ' 
		batchOK=setdiff(SampleAnnot[,batch],SampleAnnot[inSmallBatch| noBatch,batch])
		BATCH[inSmallBatch|noBatch]=batchOK[apply(cor(t(PSI_batch[batchOK,]),PSI_scale[,inSmallBatch|noBatch]),2,which.max)]
	}
	cat('Running ComBat\n')
	flush.console()
	library(sva)
	pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/Combat_priors_%s_V5.pdf',HOME,batch))
	tim=Sys.time()
	PSI_prov=ComBat(PSI_prov,BATCH, model.matrix(~condition*population,data=SampleAnnot),prior.plots=T,par.prior=T)
	print(Sys.time()-tim)  # 10s per 1000 lines ->15 min for 100 000 transcripts 
	dev.off()
	save(PSI_prov,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_%s_V5.Rdata',HOME,batch))
	gc()
}

PSI_prov=t(apply(PSI_prov,1,function(x){mod=lm(x~GC_pct_final+FivePrimeBias+condition*population,data=SampleAnnot);x-mod$coef[2]*scale(SampleAnnot$GC_pct_final,center=T,scale=F)-mod$coef[3]*scale(SampleAnnot$FivePrimeBias,center=T,scale=F)}))
colnames(PSI_prov)=colnames(FPKM_gene)

coFactors=c("X260_280", "X260_230", "rin", "total_amount_axeq", 
			"RIN_axeq", "rRNA_ratio_axeq", "GC_pct_final", "Q30_pct",'sum_all_runs','library_mass_conc','FivePrimeBias')

tic=Sys.time()
VarBatch_Adj=t(apply(PSI_prov,1,function(x){
		mod=lmer(x~population+condition+X260_280+X260_230+rin+total_amount_axeq+RIN_axeq+rRNA_ratio_axeq+GC_pct_final+Q30_pct+sum_all_runs+library_mass_conc+FivePrimeBias+(1|library_batch)+(1|batch)+(1|lane)+(1|experimentator)+(1|machine)+(1|index)+(1|date),data=SampleAnnot)
		c(sqrt(unlist(summary(mod)$varcor))	,summary(mod)$sigma,summary(mod)$coefficients[-1,3])}))
R2Batch_Adj=VarBatch_Adj[,1:7]^2/apply(VarBatch_Adj[,1:8]^2,1,sum)
R2Cont_Adj=P2R2(pnorm(abs(VarBatch_Adj[,-(1:8)]),low=F)*2,ncol(FPKM_gene),1)
toc=Sys.time()
print(toc-tic)
save(VarBatch_Adj,R2Batch_Adj,R2Cont_Adj,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/MISO_events_TechnicalCovariatesQuantification_V5_adj.Rdata',HOME))

PSI_prov=invlogit(PSI_prov)
save(PSI_prov,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_%s_V5.Rdata',HOME,length(toKeep)))

