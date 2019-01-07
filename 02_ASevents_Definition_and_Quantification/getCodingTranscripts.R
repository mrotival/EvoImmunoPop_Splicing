Source='Ens70_HISAT'
library(data.table)

transcript_type=fread('/Volumes/@Home/Annotation/Ens75_transcript_type.txt')
transcript_type_94=fread('/Volumes/@Home/Annotation/Ens94_Transcripts_GRCh37.txt')
Exons_94=fread('/Volumes/@Home/Annotation/Ens94_Exons_GRCh37.txt')
Exons_94=merge(Exons_94,transcript_type_94)
Exons_94=Exons_94[Exons_94$Chromosome%in%c(1:22,'X','Y','MT'),]
Exons_94$Strand=ifelse(Exons_94$Strand>0,'+','-')
# Annotate frequent AS events
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V5_withCounts.Rdata',sep=''))

Annot_IOE=list()
for (type in c('A3','A5','SE','MX','RI','AF','AL')){
	cat('\n',type)
	Annot_IOE[[type]]=read.table(sprintf('%s/Maxime/Splicing/MISO/hg19_ENS70/SplicingEvents_%s_strict.ioe',EVO_IMMUNO_POP,type),skip=0,header=T)
}
Annot_IOE=do.call(rbind,Annot_IOE)

#w=1:nrow(Annot_IOE)
#tim=Sys.time()
#total_transcripts_Nb=mapply(function(x,gene){isPresent=x%in%transcript_type[which(transcript_type$Gene==gene),"Transcript stable ID"]
#                                                               isCoding=x%in%transcript_type[which(transcript_type$'Transcript type'=='protein_coding' & transcript_type$Gene==gene),"Transcript stable ID"]
#                                                               c(length(isPresent),sum(isPresent),sum(isCoding))},strsplit(Annot_IOE[w,'total_transcripts'],','),Annot_IOE[w,'gene_id'])
#print(Sys.time()-tim)
#
#tim=Sys.time()
#alter_transcripts_Nb=mapply(function(x,gene){isPresent=x%in%transcript_type[which(transcript_type$Gene==gene),"Transcript stable ID"]
#                                                               isCoding=x%in%transcript_type[which(transcript_type$'Transcript type'=='protein_coding' & transcript_type$Gene==gene),"Transcript stable ID"]
#                                                               c(length(isPresent),sum(isPresent),sum(isCoding))},strsplit(Annot_IOE[w,'alternative_transcripts'],','),Annot_IOE[w,'gene_id'])
#print(Sys.time()-tim)
#
#Annot_IOE$alter_transcripts_Nb_Ens70=alter_transcripts_Nb[1,]
#Annot_IOE$alter_transcripts_Nb_Ens75=alter_transcripts_Nb[2,]
#Annot_IOE$alter_transcripts_Nb_coding_Ens75=alter_transcripts_Nb[3,]
#
#Annot_IOE$total_transcripts_Nb_Ens70=total_transcripts_Nb[1,]
#Annot_IOE$total_transcripts_Nb_Ens75=total_transcripts_Nb[2,]
#Annot_IOE$total_transcripts_Nb_coding_Ens75=total_transcripts_Nb[3,]
#
#Annot_IOE$event_lost_Ens75=(Annot_IOE$alter_transcripts_Nb_Ens75==Annot_IOE$total_transcripts_Nb_Ens75 | Annot_IOE$alter_transcripts_Nb_Ens75==0)
#Annot_IOE$codingLost=Annot_IOE$alter_transcripts_Nb_coding_Ens75==0 & Annot_IOE$total_transcripts_Nb_coding_Ens75>0
#Annot_IOE$codingGained=Annot_IOE$alter_transcripts_Nb_coding_Ens75==Annot_IOE$total_transcripts_Nb_coding_Ens75 & Annot_IOE$total_transcripts_Nb_coding_Ens75>0
#Annot_IOE$codingAny=Annot_IOE$total_transcripts_Nb_coding_Ens75>0
#
#Annot_IOE$codingBoth=Annot_IOE$alter_transcripts_Nb_coding_Ens75>0 & Annot_IOE$alter_transcripts_Nb_coding_Ens75<Annot_IOE$total_transcripts_Nb_coding_Ens75
#
#PSI_Annot$event_lost_Ens75=Annot_IOE$event_lost_Ens75[match(PSI_Annot$event_id,Annot_IOE$event_id)]
#PSI_Annot$codingBoth=Annot_IOE$codingBoth[match(PSI_Annot$event_id,Annot_IOE$event_id)]
#PSI_Annot$codingGained=Annot_IOE$codingGained[match(PSI_Annot$event_id,Annot_IOE$event_id)]
#PSI_Annot$codingLost=Annot_IOE$codingLost[match(PSI_Annot$event_id,Annot_IOE$event_id)]

#save(PSI_Annot,PSI,CIr,Count,toKeep,toKeep_S2N,toKeep_50,toKeep_80,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7_withCounts.Rdata',sep=''))


# annotate splice sites of frequent AS events  
#Source='Ens70_HISAT'
#load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7_withCounts.Rdata',sep=''))
#
#PSI_Annot$coding_type_Ens75=ifelse(PSI_Annot$event_lost_Ens75,'unknown',
#                        ifelse(PSI_Annot$codingBoth, 'change in protein-coding isoform',
#                        ifelse(PSI_Annot$codingGained, 'gain of protein-coding isoform',
#                        ifelse(PSI_Annot$codingLost, 'loss of protein-coding isoform',
#                        "Other"))))


library(data.table)
Junc_annot=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Junction_annot_V7.2.txt',HOME))

#load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup.Rdata',EVO_IMMUNO_POP))

x=strsplit(substr(PSI_Annot$event_id[toKeep],17,1000),':')
event_ID=rep(PSI_Annot$event_id[toKeep],sapply(x,length))
y=unlist(x)
gene=rep(sapply(x,function(x){x[1]}),sapply(x,length))
chr=rep(sapply(x,function(x){x[2]}),sapply(x,length))
strand=rep(sapply(x,function(x){x[length(x)]}),sapply(x,length))

event_ID_junc=event_ID[grep('[0-9]+-[0-9]+',y)]
event_junc=y[grep('[0-9]+-[0-9]+',y)]

# create Junc_MISO object containing all junctions from frequent AS events (quantified with MISO)
# for each junction, transcript is 0 if the junction is included when PSI=0, 1 if it is included when PSI=1
Junc_MISO=data.frame(event_id=event_ID_junc,chr=chr[grep('[0-9]+-[0-9]+',y)],start=as.numeric(gsub('([0-9]+)-([0-9]+)','\\1',event_junc)),end=as.numeric(gsub('([0-9]+)-([0-9]+)','\\2',event_junc)),type=substr(event_ID_junc,17,18),strand=strand[grep('[0-9]+-[0-9]+',y)],transcript=1)

# add splicing accross exon for SE events 
SEx=x[sapply(x,function(y){y[1]=='SE'})]
SE_junc=sapply(SEx,function(x){y=strsplit(x[3:4],'-');paste(y[[1]][1],y[[2]][2],sep='-')})
chr=sapply(SEx,function(x){x[2]})
gene=sapply(SEx,function(x){x[1]})
strand=sapply(SEx,function(x){x[length(x)]})
Junc_MISO_SE=data.frame(event_id=PSI_Annot$event_id[toKeep[sapply(x,function(y){y[1]=='SE'})]],chr=chr,start=as.numeric(gsub('([0-9]+)-([0-9]+)','\\1',SE_junc)),end=as.numeric(gsub('([0-9]+)-([0-9]+)','\\2',SE_junc)),type="SE",strand=strand,transcript=0)

Junc_MISO=rbind(Junc_MISO,Junc_MISO_SE)


RIx=PSI_Annot$event_id[toKeep[sapply(x,function(y){y[1]=='RI'})]]
Regx='(ENSG[0-9]+);RI:([0-9]+):([0-9]+):([0-9]+)-([0-9]+):([0-9]+):([+-])'
SpliceSites_MISO_RI_1=data.frame(event_id=RIx,chr=gsub(Regx,'\\2',RIx),start=as.numeric(gsub(Regx,'\\3',RIx))-1,strand=gsub(Regx,'\\7',RIx))
SpliceSites_MISO_RI_1$site_id=paste(SpliceSites_MISO_RI_1$chr,SpliceSites_MISO_RI_1$start,SpliceSites_MISO_RI_1$strand,ifelse(SpliceSites_MISO_RI_1$strand=='+','acceptor','donor'))
SpliceSites_MISO_RI_2=data.frame(event_id=RIx,chr=gsub(Regx,'\\2',RIx),start=as.numeric(gsub(Regx,'\\6',RIx))+1,strand=gsub(Regx,'\\7',RIx))
SpliceSites_MISO_RI_2$site_id=paste(SpliceSites_MISO_RI_2$chr,SpliceSites_MISO_RI_2$start,SpliceSites_MISO_RI_2$strand,ifelse(SpliceSites_MISO_RI_2$strand=='+','donor','acceptor'))


#mm=match(paste(Junc_MISO$chr,Junc_MISO$start,Junc_MISO$end),paste(Junc$chrom,Junc$start,Junc$end+1))

Junc_MISO$start_id=paste(Junc_MISO$chr,Junc_MISO$start+1,Junc_MISO$strand,ifelse(Junc_MISO$strand=='+','donor','acceptor'))
Junc_MISO$end_id=paste(Junc_MISO$chr,Junc_MISO$end-1,Junc_MISO$strand,ifelse(Junc_MISO$strand=='+','acceptor','donor'))

#mm=match(paste(Junc_MISO$chr,Junc_MISO$start,Junc_MISO$end),paste(Junc$chrom,Junc$start,Junc$end+1))
#mean(Junc_MISO$start_id%in%SpliceSites$site_id)
#mean(Junc_MISO$end_id%in%SpliceSites$site_id)


# Annotate transcript (for A3, A5, AF, AL)
first_A3=ifelse(Junc_MISO$strand[Junc_MISO$type=='A3']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='A3'],':'),function(x){x[4]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='A3'],':'),function(x){x[3]}))
first_A5=ifelse(Junc_MISO$strand[Junc_MISO$type=='A5']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='A5'],':'),function(x){x[4]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='A5'],':'),function(x){x[3]}))
first_AF=ifelse(Junc_MISO$strand[Junc_MISO$type=='AF']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='AF'],':'),function(x){x[4]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='AF'],':'),function(x){x[3]}))
first_AL=ifelse(Junc_MISO$strand[Junc_MISO$type=='AL']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='AL'],':'),function(x){x[3]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='AL'],':'),function(x){x[4]}))
strand_A3=Junc_MISO$strand[Junc_MISO$type=='A3']=='+'
strand_A5=Junc_MISO$strand[Junc_MISO$type=='A5']=='+'
strand_AF=Junc_MISO$strand[Junc_MISO$type=='AF']=='+'
strand_AL=Junc_MISO$strand[Junc_MISO$type=='AL']=='+'


junc_ID=paste(Junc_MISO[,'start'],Junc_MISO[,'end'],sep='-')
Junc_MISO[Junc_MISO$type=='A3','transcript']=ifelse(junc_ID[Junc_MISO$type=='A3']==first_A3,0,1)*ifelse(strand_A3,1,-1)+ifelse(strand_A3,0,1)
Junc_MISO[Junc_MISO$type=='A5','transcript']=ifelse(junc_ID[Junc_MISO$type=='A5']==first_A5,0,1)*ifelse(strand_A5,1,-1)+ifelse(strand_A5,0,1)
Junc_MISO[Junc_MISO$type=='AF','transcript']=ifelse(junc_ID[Junc_MISO$type=='AF']==first_AF,1,0)*ifelse(strand_AF,1,-1)+ifelse(strand_AF,0,1)
Junc_MISO[Junc_MISO$type=='AL','transcript']=ifelse(junc_ID[Junc_MISO$type=='AL']==first_AL,0,1)*ifelse(strand_AL,1,-1)+ifelse(strand_AL,0,1)

Junc_MISO[Junc_MISO$type=='RI','transcript']=0

# annotation for AL
Junc_MISO[Junc_MISO$type=='AL','transcript']=ifelse(junc_ID[Junc_MISO$type=='AL']==first_AL,0,1)*ifelse(strand_AL,1,-1)+ifelse(strand_AL,0,1)

# annotation for SE
Junc_MISO[Junc_MISO$type=='SE','transcript']=ifelse(junc_ID[Junc_MISO$type=='SE']%in%SE_junc,0,1)

# annotation for MX
MX_id=paste(Junc_MISO[Junc_MISO$type=='MX','start'],Junc_MISO[Junc_MISO$type=='MX','end'],sep='-')
MX_1=ifelse(Junc_MISO$strand[Junc_MISO$type=='MX']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='MX'],':'),function(x){x[3]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='MX'],':'),function(x){x[5]}))
MX_2=ifelse(Junc_MISO$strand[Junc_MISO$type=='MX']=='+',sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='MX'],':'),function(x){x[4]}),sapply(strsplit(Junc_MISO$event_id[Junc_MISO$type=='MX'],':'),function(x){x[6]}))
Junc_MISO[Junc_MISO$type=='MX','transcript']=ifelse(MX_id==MX_1 |MX_id==MX_2,1,0)
table(Junc_MISO$type, Junc_MISO$transcript)

# annotate Junctions from MISO events
mm=match(paste(Junc_MISO$chr,Junc_MISO$start,Junc_MISO$end),paste(Junc_annot$chrom_intron,Junc_annot$start_intron-1,Junc_annot$end_intron+1))
Junc_MISO$GerpRS_start=(Junc_annot$GerpRS_start[mm]+Junc_annot$GerpRS_start2[mm])/2
Junc_MISO$GerpRS_end=(Junc_annot$GerpRS_end[mm]+Junc_annot$GerpRS_end2[mm])/2
Junc_MISO$Coverage=Junc_annot$Total_count[mm]
Junc_MISO$Coverage_NS=Junc_annot$Total_count_NS[mm]
Junc_MISO$Coverage_LPS=Junc_annot$Total_count_LPS[mm]
Junc_MISO$Coverage_PAM3=Junc_annot$Total_count_PAM3[mm]
Junc_MISO$Coverage_R848=Junc_annot$Total_count_R848[mm]
Junc_MISO$Coverage_IAV=Junc_annot$Total_count_IAV[mm]
Junc_MISO$coding_startSite=Junc_annot$coding_startSite[mm]
Junc_MISO$coding_endSite=Junc_annot$coding_endSite[mm]

#SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.txt',HOME))

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

# We distinguish between :
# (a) AS events that do not alter the protein product (eg. Alternative UTR splicing) 
# (b) AS events that are associated to a change of protein isoform (type i), 
# (c) AS events that are associated to a gain/loss of function (eg nonsense splicing, type ii).

#SpliceSite_MISO$coding_site_Ens75 = SpliceSite_MISO$pos%in%paste(Exons_94$'Chromosome/scaffold name',Exons_94$'Genomic coding end'+1) | SpliceSite_MISO$pos%in%paste(Exons_94$'Chromosome/scaffold name',Exons_94$'Genomic coding start'-1)
#SpliceSite_MISO$coding_site_Ens75_withoutNMD = SpliceSite_MISO$pos %in% paste(Exons_94$'Chromosome/scaffold name',Exons_94$'Genomic coding end'+1)[Exons_94$'Transcript type'=='protein_coding'] | SpliceSite_MISO$pos%in%paste(Exons_94$'Chromosome/scaffold name',Exons_94$'Genomic coding start'-1)[Exons_94$'Transcript type'=='protein_coding']
#SpliceSite_MISO$site_of_coding_exon_Ens75_withoutNMD = SpliceSite_MISO$pos%in%Exons_94[!is.na(Exons_94$'Genomic coding end') & Exons_94$'Transcript type'=='protein_coding',get('pos_end')] | SpliceSite_MISO$pos%in%Exons_94[!is.na(Exons_94$'Genomic coding start') & Exons_94$'Transcript type'=='protein_coding',get('pos_start')]
#SpliceSite_MISO$site_of_coding_transcript_Ens75_withoutNMD = SpliceSite_MISO$pos%in%SpliceSite_MISO$pos%in%Exons_94[Exons_94$'Transcript type'=='protein_coding',get('pos_end')] | SpliceSite_MISO$pos%in%Exons_94[Exons_94$'Transcript type'=='protein_coding',get('pos_start')]
#SpliceSite_MISO$site_of_noncoding_transcript_Ens75_withoutNMD = SpliceSite_MISO$pos%in%SpliceSite_MISO$pos%in%Exons_94[Exons_94$'Transcript type'!='protein_coding',get('pos_end')] | SpliceSite_MISO$pos%in%Exons_94[Exons_94$'Transcript type'!='protein_coding',get('pos_start')]

coding_type_Event = SpliceSite_MISO[,.(
#                                      non_coding_Event=!any(coding_site_Ens75_withoutNMD),
#                                      non_coding_ASevent=!any(coding_site_Ens75_withoutNMD & Nb_splice_form_where_found!=2),
#
#                                      no_coding_site=!any(coding_site_Ens75_withoutNMD),
#                                      no_coding_AS_site=!any(coding_site_Ens75_withoutNMD & Nb_splice_form_where_found!=2),
#                                      coding_site_lost=any(coding_site_Ens75_withoutNMD & Nb_splice_form_where_found!=2 & transcript==0),
#                                      coding_site_gain=any(coding_site_Ens75_withoutNMD & Nb_splice_form_where_found!=2 & transcript==1),
#
#                                      no_coding_AS_exon=!any(site_of_coding_exon_Ens75_withoutNMD & Nb_splice_form_where_found!=2),
#                                      coding_exon_lost=any(site_of_coding_exon_Ens75_withoutNMD & Nb_splice_form_where_found!=2 & transcript==0),
#                                      coding_exon_gain=any(site_of_coding_exon_Ens75_withoutNMD & Nb_splice_form_where_found!=2 & transcript==1),
#
#                                      no_coding_AS_transcript=!any(site_of_coding_transcript_Ens75_withoutNMD & Nb_splice_form_where_found!=2),
#                                      coding_transcript_lost=any(site_of_coding_transcript_Ens75_withoutNMD & Nb_splice_form_where_found!=2 & transcript==0),
#                                      coding_transcript_gain=any(site_of_coding_transcript_Ens75_withoutNMD & Nb_splice_form_where_found!=2 & transcript==1),
                                      non_conserved=!any(Gerp_site>2),
                                      non_conserved_ASevent=!any(Gerp_site>2 & Nb_splice_form_where_found!=2),

                                      conserved_0_any=any(Gerp_site>2 & Nb_splice_form_where_found!=2 & transcript==0),
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
#
#coding_type_Event$coding_site_switch=coding_type_Event$coding_site_lost & coding_type_Event$coding_site_gain
#coding_type_Event$coding_site_gain=coding_type_Event$coding_site_gain & !coding_type_Event$coding_site_switch
#coding_type_Event$coding_site_lost=coding_type_Event$coding_site_lost & !coding_type_Event$coding_site_switch
#
#coding_type_Event$coding_exon_switch=coding_type_Event$coding_exon_gain & coding_type_Event$coding_exon_lost
#coding_type_Event$coding_exon_gain=coding_type_Event$coding_exon_gain & !coding_type_Event$coding_exon_switch
#coding_type_Event$coding_exon_lost=coding_type_Event$coding_exon_lost & !coding_type_Event$coding_exon_switch
#
#coding_type_Event$coding_transcript_switch=coding_type_Event$coding_transcript_gain & coding_type_Event$coding_transcript_lost
#coding_type_Event$coding_transcript_gain=coding_type_Event$coding_transcript_gain & !coding_type_Event$coding_transcript_switch
#coding_type_Event$coding_transcript_lost=coding_type_Event$coding_transcript_lost & !coding_type_Event$coding_transcript_switch


#coding_type_Event$coding_site_Ens94=ifelse(coding_type_Event$no_coding_AS_site, 'non coding Alternative splice sites',
#                        ifelse(coding_type_Event$coding_site_switch, 'change in protein-coding splice site',
#                        ifelse(coding_type_Event$coding_site_gain, 'gain of protein-coding splice site',
#                        ifelse(coding_type_Event$coding_site_lost, 'loss of protein-coding splice site',
#                        "Other"))))
#
#coding_type_Event$coding_exon_Ens94=ifelse(coding_type_Event$no_coding_AS_exon, 'non coding Alternative exon',
#                        ifelse(coding_type_Event$coding_exon_switch, 'change in protein-coding exon',
#                        ifelse(coding_type_Event$coding_exon_gain, 'gain of protein-coding exon',
#                        ifelse(coding_type_Event$coding_exon_lost, 'loss of protein-coding exon',
#                        "Other"))))
#
#coding_type_Event$coding_transcript_Ens94=ifelse(coding_type_Event$no_coding_AS_transcript, 'non coding Alternative isoform',
#                        ifelse(coding_type_Event$coding_transcript_switch, 'change in protein-coding isoform',
#                        ifelse(coding_type_Event$coding_transcript_gain, 'gain of protein-coding isoform',
#                        ifelse(coding_type_Event$coding_transcript_lost, 'loss of protein-coding isoform',
#                        "Other"))))

mm=match(PSI_Annot$event_id,coding_type_Event$event_id)

PSI_Annot$conserved_site_GerpRS=ifelse(coding_type_Event$non_conserved_AS_site[mm], 'non conserved Alternative splice sites',
                        ifelse(coding_type_Event$conserved_site_switch_all[mm], 'change in conserved splice site (all)',
                        ifelse(coding_type_Event$conserved_site_switch[mm], 'change in conserved splice site',
                        ifelse(coding_type_Event$conserved_site_gain_all[mm], 'gain of conserved splice site (all)',
                        ifelse(coding_type_Event$conserved_site_gain[mm], 'gain of conserved splice site',
                        ifelse(coding_type_Event$conserved_site_lost_all[mm], 'loss of conserved splice site (all)',
                        ifelse(coding_type_Event$conserved_site_lost[mm], 'loss of conserved splice site',
                        "Other")))))))

#PSI_Annot$coding_type_Ens70=ifelse(coding_type_Event$non_coding_ASevent[mm], 'non coding Alternative splice sites',
#                        ifelse(coding_type_Event$coding_switch[mm], 'change in protein-coding splice site',
#                        ifelse(coding_type_Event$coding_gain[mm], 'gain of protein-coding splice site',
#                        ifelse(coding_type_Event$coding_lost[mm], 'loss of protein-coding splice site',
#                        "Other"))))
#

#PSI_Annot$coding_site_Ens94=ifelse(coding_type_Event$no_coding_AS_site[mm], 'non coding Alternative splice sites',
#                        ifelse(coding_type_Event$coding_site_switch[mm], 'change in protein-coding splice site',
#                        ifelse(coding_type_Event$coding_site_gain[mm], 'gain of protein-coding splice site',
#                        ifelse(coding_type_Event$coding_site_lost[mm], 'loss of protein-coding splice site',
#                        "Other"))))
#
#PSI_Annot$coding_exon_Ens94=ifelse(coding_type_Event$no_coding_AS_exon[mm], 'non coding Alternative exon',
#                        ifelse(coding_type_Event$coding_exon_switch[mm], 'change in protein-coding exon',
#                        ifelse(coding_type_Event$coding_exon_gain[mm], 'gain of protein-coding exon',
#                        ifelse(coding_type_Event$coding_exon_lost[mm], 'loss of protein-coding exon',
#                        "Other"))))
#
#PSI_Annot$coding_transcript_Ens94=ifelse(coding_type_Event$no_coding_AS_transcript[mm], 'non coding Alternative isoform',
#                        ifelse(coding_type_Event$coding_transcript_switch[mm], 'change in protein-coding isoform',
#                        ifelse(coding_type_Event$coding_transcript_gain[mm], 'gain of protein-coding isoform',
#                        ifelse(coding_type_Event$coding_transcript_lost[mm], 'loss of protein-coding isoform',
#                        "Other"))))


save(PSI_Annot,PSI,Count,toKeep,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))

load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V5_nodup_withSelect_RBP.Rdata',EVO_IMMUNO_POP))
save(RESobs_nodup_1Mb,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/cis-psiQTL_MISO_counts_V7_nodup_withSelect_RBP_coding.Rdata',HOME)))


load(sprintf('%s/Maxime/Splicing/sQTL/MatrixEQTL-cis/cis-psiQTL_MISO_counts_V6.9_nodup.Rdata',EVO_IMMUNO_POP))
RESobs_nodup_1Mb$coding_type = PSI_Annot$coding_type[ match( RESobs_nodup_1Mb$event_id, PSI_Annot$event_id ) ]
save(RESobs_nodup_1Mb,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/cis-psiQTL_MISO_counts_V7.1_nodup_withSelect_RBP_coding_eQTL.Rdata',HOME))

color_coding=c("#225EA8","#FEE586","#E31A1C","#A1DAB4")
names(color_coding)=c("non coding","isoform switch","loss of function","gain of function")

barplot(table(PSI_Annot$coding_type, PSI_Annot$event_type),col=rev(color_coding))


#table(grepl('change in|gain of',PSI_Annot$conserved_site_GerpRS),PSI_Annot$coding_type)
#
########## debug
#
#PSI_Annot[which(PSI_Annot$coding_type=='Other')[1:5],]
#
#myevent=PSI_Annot[which(PSI_Annot$coding_type=='Other')[3],'event_id']
#myevent
#Annot_IOE[Annot_IOE$event_id==myevent,]
#mytranscripts=strsplit(Annot_IOE[Annot_IOE$event_id==myevent,'total_transcripts'],',')[[1]]
#mytranscripts
#coding_type_Event[coding_type_Event$event_id==myevent,]
#SpliceSite_MISO[SpliceSite_MISO$event_id==myevent,]
#PSI_Annot[which(PSI_Annot$event_id==myevent),]
#Exons_94[(Exons_94$'Exon region start (bp)'-1)%in%SpliceSite_MISO[SpliceSite_MISO$event_id==myevent,pos] | (Exons_94$'Exon region end (bp)'+1)%in%SpliceSite_MISO[SpliceSite_MISO$event_id==myevent,pos],]
#Exons_94[Exons_94$'Transcript stable ID'%in%mytranscripts,]

# case 1 : we stop at the first coding exon. so both splice sites are non coding but one of them mleads to a protein, the other doesn't.
# case 2 : one of the transcripts that use this splice site in Ens94 is coding, this transcript was absent from Ens70, and was thus not tested
# case 3 : 


as.data.table(SpliceSite_MISO)[order(event_type,Site_present_in_both_transcript_of_ASevent),.(Pct_conserved=mean(Gerp_site>2)),by=.(Site_present_in_both_transcript_of_ASevent, event_type)]
#    Site_present_in_both_transcript_of_ASevent event_type Pct_conserved
# 1:                                      FALSE         A3     0.6771429
# 2:                                       TRUE         A3     0.9020952
# 3:                                      FALSE         A5     0.5950533
# 4:                                       TRUE         A5     0.8913676
# 5:                                      FALSE         AF     0.5344400
# 6:                                       TRUE         AF     0.8919722
# 7:                                      FALSE         AL     0.6580977
# 8:                                       TRUE         AL     0.9108826
# 9:                                      FALSE         MX     0.6925795
#10:                                       TRUE         MX     0.9196113
#11:                                      FALSE         RI     0.8872792
#12:                                      FALSE         SE     0.6954114
#13:                                       TRUE         SE     0.9263449
