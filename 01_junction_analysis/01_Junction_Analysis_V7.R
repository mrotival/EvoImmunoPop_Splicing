


library(data.table)
library(GenomicRanges)

Junc_annot=fread(sprintf('%s/Maxime/Splicing/NoisySplicing/AllJunctions_V7.txt',EVO_IMMUNO_POP))
Junc_annot$start=Junc_annot$start+1 # by default in all previous scripts, start is the last bp of the previous exons and end is the last bp of the intron. we change it to firt/last bp of the intron

##### correct strand when Intro polis and HISAT disagree
wStrand=which(colnames(Junc_annot)=='strand')
colnames(Junc_annot)[wStrand[1]]='strand_HISAT'
colnames(Junc_annot)[wStrand[2]]='strand_Intropolis'
wStrandError=which(Junc_annot$strand_HISAT!=Junc_annot$strand_Intropolis & Junc_annot$strand_Intropolis!='') # cases where strand called y HISAT does not match strand in intropolis. We trust intropolis as a more reliable source of information.
# correct donor/acceptor status
Junc_annot$donor_site[wStrandError]=Junc_annot$acceptor[wStrandError] # this is inversed on purpose, because of the inversed strand error 
Junc_annot$acceptor_site[wStrandError]=Junc_annot$donor[wStrandError] # this is inversed on purpose, because of the inversed strand error 
Junc_annot$strand_correct=Junc_annot$strand_HISAT
Junc_annot$strand_correct[wStrandError]=Junc_annot$strand_Intropolis[wStrandError] # we correct the strand of HISAT in these cases

# annotate Donor/Acceptor status
Junc_annot$DA_type=ifelse(paste(Junc_annot$donor_site,Junc_annot$acceptor_site)%in%c('GT AG','GC AG','AT AC'),paste(Junc_annot$donor_site,Junc_annot$acceptor_site),'other')

##### Annotate junctions 
introns_df=as.data.frame(fread(file=sprintf('%s/Annotation/IntronAnnotation_hg37_ens70_detectedERC.txt',HOME)))
### Annotated junctions (based only on transcripts with FPKM > 1) #####
Junc_annot$Annotated_junc=paste(Junc_annot$chrom,Junc_annot$start,Junc_annot$end)%in%paste(introns_df$chrom,introns_df$start,introns_df$end)

Exons=read.table(paste(HOME,'/Annotation/ExonCoordinates_hg37_ens70.txt',sep=''),header=T,sep='\t',quote='',comment='')
rownames(Exons)=paste(Exons$Ensembl.Exon.ID,Exons$Ensembl.Transcript.ID)

### Annotate Gene associated to start/end site
mm_end=match(paste(Junc_annot$chr,as.N(Junc_annot$start)-1,sep=':'),paste(Exons$Chromosome.Name,as.N(Exons$Exon.Chr.End..bp.),sep=':'))
Junc_annot$Gene_start=Exons$Ensembl.Gene.ID[mm_end]
Junc_annot$symbol_start=Exons$Associated.Gene.Name[mm_end]

mm_start=match(paste(Junc_annot$chr,as.N(Junc_annot$end)+1,sep=':'),paste(Exons$Chromosome.Name,as.N(Exons$Exon.Chr.Start..bp.),sep=':'))
Junc_annot$Gene_end=Exons$Ensembl.Gene.ID[mm_start]
Junc_annot$symbol_end=Exons$Associated.Gene.Name[mm_start]

### Annotate gene(s) associated to junction 

# if any of the splice sites matches a known gene, junction is assigned to that gene. if both splice sites match a different gene donor/acceptor genes are reported.
Junc_annot$gene=ifelse(is.na(Junc_annot$Gene_start) & !is.na(Junc_annot$Gene_end),Junc_annot$Gene_end,
				ifelse(is.na(Junc_annot$Gene_end) & !is.na(Junc_annot$Gene_start),Junc_annot$Gene_start,
					ifelse(Junc_annot$Gene_start==Junc_annot$Gene_end,Junc_annot$Gene_start,
						ifelse(Junc_annot$strand_HISAT=='+',paste(Junc_annot$Gene_start,Junc_annot$Gene_end,sep='//'),
							paste(Junc_annot$Gene_end,Junc_annot$Gene_start,sep='//'))
								)
							)
						)

# if any of the splice sites matches a known gene, junction is assigned to that gene. if both splice sites match a different gene donor/acceptor genes are reported.
Junc_annot$symbol=ifelse(is.na(Junc_annot$symbol_start) & !is.na(Junc_annot$symbol_end),Junc_annot$symbol_end,
					ifelse(is.na(Junc_annot$symbol_end) & !is.na(Junc_annot$symbol_start),Junc_annot$symbol_start,
						ifelse(Junc_annot$symbol_start==Junc_annot$symbol_end,Junc_annot$symbol_start,
							ifelse(Junc_annot$strand_HISAT=='+',paste(Junc_annot$symbol_start,Junc_annot$symbol_end,sep='//'),
								paste(Junc_annot$symbol_end,Junc_annot$symbol_start,sep='//'))
									)
								)
							)

### Annotate coding status (v70) 
# of start/end sites 
Junc_annot$coding_start=paste(Junc_annot$chr,as.N(Junc_annot$start)-1,sep=':')%in%paste(Exons$Chromosome.Name,as.N(Exons$Genomic.coding.end),sep=':')
Junc_annot$coding_end=paste(Junc_annot$chr,as.N(Junc_annot$end)+1,sep=':')%in%paste(Exons$Chromosome.Name,as.N(Exons$Genomic.coding.start),sep=':')

# of junction
Junc_annot$coding=Junc_annot$coding_start & Junc_annot$coding_end

### Annotate GerpRS of end/start site (v70) 
Junc_annot$GerpRS_meanStart=(Junc_annot$GerpRS_start+Junc_annot$GerpRS_start2)/2
Junc_annot$GerpRS_meanEnd=(Junc_annot$GerpRS_end+Junc_annot$GerpRS_end2)/2

# clarify colnames
colnames(Junc_annot)=c("Total_count", "Total_count_NS", "Total_count_LPS", "Total_count_PAM3CSK4", 
"Total_count_R848", "Total_count_IAV", "chrom_intron", "start_intron", "end_intron", "strand_HISAT", 
"GerpRS_start", "GerpRS_start2", "GerpRS_end", "GerpRS_end2", "acceptor_site", "donor_site", "strand_Intropolis", 
"donor_Intropolis", "acceptor_Intropolis", "Nbread_Intropolis", "NbSample_Intropolis", "PctOfTotalReadsMappingtoStartSite",
"PctOfTotalReadsMappingtoEndSite", "TotalReadsMappingtoStartSite", "TotalReadsMappingtoEndSite","strand_correct",
"DA_type", "Annotated_junc", "Gene_startSite", "symbol_startSite",
"Gene_endSite", "symbol_endSite","gene_junction", "symbol_junction", "coding_startSite", "coding_endSite", "coding_junction", "GerpRS_meanStartSite", "GerpRS_meanEndSite")

### Annotate genes that overlap a junction.
# when a junction overlaps multiple genes we report both genes

GeneAnnot=read.table(paste(HOME,'/Annotation/GeneAnnotation_hg37_ens70.txt',sep=''),sep='\t',comment='',quote='',stringsAsFactors=FALSE,header=T,row.names=1)
colnames(GeneAnnot)[1]='Ensembl.Gene.ID'
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='X']=23
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='Y']=24
GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='MT']=26
GeneAnnot$Chromosome.Name=as.numeric(GeneAnnot$Chromosome.Name)
GeneAnnot=GeneAnnot[!is.na(GeneAnnot$Chromosome.Name),] # remove alternate chromosomes and unassigned scaffolds
GeneAnnot$Strand=ifelse(GeneAnnot$Strand>0,'+','-') # change strand from -1,1 to '+','-'

#### using strand information
Gene_GR=makeGRangesFromDataFrame(GeneAnnot,seqnames='Chromosome.Name',start.field='Gene.Start..bp.',end.field='Gene.End..bp.',strand.field='Strand')
Junc_GR=makeGRangesFromDataFrame(Junc_annot,seqnames='chrom_intron',start.field='start_intron',end.field='end_intron',strand.field='strand_correct')
oo=findOverlaps(Junc_GR,Gene_GR)

#find junctions that match a single gene and annotate them
is.unique=function(x){!x %in% x[duplicated(x)]}
uniq=is.unique(queryHits(oo))

Junc_annot$overlappingGene_withStrand=NA
Junc_annot$overlappingGene_withStrand[queryHits(oo)[uniq]]=GeneAnnot$Ensembl.Gene.ID[subjectHits(oo)[uniq]]
Junc_annot$overlappingSymbol_withStrand=NA
Junc_annot$overlappingSymbol_withStrand[queryHits(oo)[uniq]]=GeneAnnot$Associated.Gene.Name[subjectHits(oo)[uniq]]

#for junction mapping >1 gene, collapse genes and annotate junctions
collapseGene=function(x){paste(sort(unique(x)),collapse='//')}

junctions_overlappingMultiGene=as.data.table(cbind(juncLine=queryHits(oo)[!uniq],GeneAnnot[subjectHits(oo)[!uniq],c('Associated.Gene.Name','Ensembl.Gene.ID')]))
collapsed = junctions_overlappingMultiGene[, .(GeneCollapse = collapseGene(Ensembl.Gene.ID), 
                                                SymbolCollapse = collapseGene(Associated.Gene.Name)), by = (juncLine)]
Junc_annot$overlappingSymbol_withStrand[as.numeric(collapsed$juncLine)]=collapsed$SymbolCollapse
Junc_annot$overlappingGene_withStrand[as.numeric(collapsed$juncLine)]=collapsed$GeneCollapse


#### NOT using strand information
#find junctions that match a single gene and annotate them
strand(Gene_GR)='*'
oo=findOverlaps(Junc_GR,Gene_GR)
uniq=is.unique(queryHits(oo))

Junc_annot$overlappingGene_noStrand=NA
Junc_annot$overlappingGene_noStrand[queryHits(oo)[uniq]]=GeneAnnot$Ensembl.Gene.ID[subjectHits(oo)[uniq]]
Junc_annot$overlappingSymbol_noStrand=NA
Junc_annot$overlappingSymbol_noStrand[queryHits(oo)[uniq]]=GeneAnnot$Associated.Gene.Name[subjectHits(oo)[uniq]]

#for junction mapping >1 gene, collapse genes and annotate junctions
collapseGene=function(x){paste(sort(unique(x)),collapse='//')}

junctions_overlappingMultiGene=as.data.table(cbind(juncLine=queryHits(oo)[!uniq],GeneAnnot[subjectHits(oo)[!uniq],c('Associated.Gene.Name','Ensembl.Gene.ID')]))
collapsed = junctions_overlappingMultiGene[, .(GeneCollapse = collapseGene(Ensembl.Gene.ID), 
                                                SymbolCollapse = collapseGene(Associated.Gene.Name)), by = (juncLine)]
Junc_annot$overlappingGene_noStrand[as.numeric(collapsed$juncLine)]=collapsed$GeneCollapse
Junc_annot$overlappingSymbol_noStrand[as.numeric(collapsed$juncLine)]=collapsed$SymbolCollapse

########### Annotate start and end site
Junc_annot$start_id=paste(Junc_annot$chrom,Junc_annot$start_intron,Junc_annot$strand_correct)
Junc_annot$end_id=paste(Junc_annot$chrom,Junc_annot$end_intron,Junc_annot$strand_correct)

########### Annotate start and end site
colnames(Junc_annot)=c("Total_count", "Total_count_NS", "Total_count_LPS", "Total_count_PAM3CSK4", 
"Total_count_R848", "Total_count_IAV", "chrom_intron", "start_intron", "end_intron", "strand_HISAT", 
"GerpRS_start", "GerpRS_start2", "GerpRS_end", "GerpRS_end2", "acceptor_site", "donor_site", "strand_Intropolis", 
"donor_Intropolis", "acceptor_Intropolis", "Nbread_Intropolis", "NbSample_Intropolis", "PctOfTotalReadsMappingtoStartSite",
"PctOfTotalReadsMappingtoEndSite", "TotalReadsMappingtoStartSite", "TotalReadsMappingtoEndSite", "strand_correct","DA_type", 
"Annotated_junc", "Gene_startSite", "symbol_startSite","Gene_endSite", "symbol_endSite", "gene_junction","symbol_junction",
"coding_startSite", "coding_endSite", "coding_junction", "GerpRS_meanStartSite", "GerpRS_meanEndSite", 
"overlappingGene_withStrand", "overlappingSymbol_withStrand", "overlappingGene_noStrand", "overlappingSymbol_noStrand","start_id", "end_id")


write.table(Junc_annot,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Junction_annot_V7.2.txt',HOME),quote=F,sep='\t',row.names=F)
Junc_annot=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Junction_annot_V7.2.txt',HOME))

####### Start sites 
### + strand
StartSite=Junc_annot
#DonorSite_pos$start=DonorSite_pos$start+1 # added 1 to match UCSC coordinates
StartSite$type=ifelse(StartSite$strand_correct=='+','donor',ifelse(StartSite$strand_correct=='-','acceptor',NA))
StartSite$site_id=paste(StartSite$start_id,StartSite$type)

tim=Sys.time()
Tot_count_site = StartSite[, .(Total_count = sum(Total_count, na.rm = T), 
                    Total_count_NS = sum(Total_count_NS, na.rm = T),
                    Total_count_LPS = sum(Total_count_LPS, na.rm = T),
                    Total_count_PAM3CSK4 = sum(Total_count_PAM3CSK4, na.rm = T),
                    Total_count_R848 = sum(Total_count_R848, na.rm = T),
                    Total_count_IAV = sum(Total_count_IAV, na.rm = T)), by = site_id]
print(Sys.time()-tim)

# order junctions by starts_id, known/unknown status/ decreasing read count
StartSite=StartSite[order(StartSite$site_id,ifelse(is.na(StartSite$Nbread_Intropolis),1,0),-StartSite$Total_count),]
# keep only one junctiopn per start site, giving priority to known junctions, and junctions with highest read count
StartSite=StartSite[!duplicated(StartSite$site_id),]
StartSite$Junction_count=StartSite$Total_count
for(i in colnames(Tot_count_site)){
	StartSite[[i]]=Tot_count_site[match(StartSite$site_id,Tot_count_site$site_id),get(i)]
}
# annotate start site
StartSite$Gerp_site=StartSite$GerpRS_meanStartSite
StartSite$inEnsembl=!is.na(StartSite$Gene_startSite)
StartSite$coding=StartSite$coding_startSite
StartSite$gene_Mapped=StartSite$gene_junction
StartSite$gene=StartSite$Gene_startSite
StartSite$symbol=StartSite$symbol_startSite
StartSite$sequence=ifelse(StartSite$type=='acceptor',StartSite$acceptor_site,StartSite$donor_site)


StartConserv=as.data.frame(fread(sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteStartConserv_V2.txt',EVO_IMMUNO_POP)))
colnames(StartConserv)[1:12]=c("chrom_intron","start_intron","strand_HISAT","GerpRS_start","Total_count_site","PhylogeneticAge1","PhylogeneticAge2","PhylogeneticAge3","FirstAppeared_age1",'PhyloDist_FirstAppeared_age1','ancestralSequence','allSequence')
StartConserv$start_intron=StartConserv$start_intron+1 # add 1 to the start to match the 1st bp of the intron
#! test
	# age 1: first occurence of the current splice site in the ancestral sequences (not allowing missing data or loss and recovery of function)
	# age 2: first occurence of the current splice site in the ancestral sequences (allowing missing data)
	# age 3: first occurence of a splice site in the ancestral sequences (allowing loss and recovery of function, and missing data)
SummaryStats_start=t(sapply(strsplit(StartConserv$allSequence,':'),function(x){y=gsub('.*_([ATGC-]*)','\\1',x);c(FirstAlignedSeq=max(which(y!='--')),NbAlignedSeq=sum(y!='--'),NbAlignedSeqPrim=sum(y[1:11]!='--'),NbAlignedSeqMam=sum(y[1:18]!='--'),NbMatchingSeq=sum(y==y[1]),NbMatchingSeqPrim=sum(y[1:11]==y[1]),NbMatchingSeqMam=sum(y[1:18]==y[1]),NbObservedSeq=luq(y[y!='--']))}))
SummaryStats_start_age=sapply(strsplit(StartConserv$ancestralSequence,':'),function(x){y=gsub('.*_([ATGC-]*)','\\1',x);max(which(y==y[1]))})
StartConserv=cbind(StartConserv,SummaryStats_start,SummaryStats_start_age)
colnames(StartConserv)[13:21]=c("FirstAlignedSeq","NbAlignedSeq","NbAlignedSeqPrim","NbAlignedSeqMam","NbMatchingSeq","NbMatchingSeqPrim","NbMatchingSeqMam","NbObservedSeq","newAge")
write.table(StartConserv,file=sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteStartConserv_V7.txt',EVO_IMMUNO_POP),sep='\t',col.names=F,quote=F,row.names=F)

StartConserv=as.data.frame(fread(sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteStartConserv_V7.txt',EVO_IMMUNO_POP)))
colnames(StartConserv)[1:12]=c("chrom_intron","start_intron","strand_HISAT","GerpRS_start","Total_count_site","PhylogeneticAge1","PhylogeneticAge2","PhylogeneticAge3","FirstAppeared_age1",'PhyloDist_FirstAppeared_age1','ancestralSequence','allSequence')
colnames(StartConserv)[13:21]=c("FirstAlignedSeq","NbAlignedSeq","NbAlignedSeqPrim","NbAlignedSeqMam","NbMatchingSeq","NbMatchingSeqPrim","NbMatchingSeqMam","NbObservedSeq","newAge")

mm_Age=match(paste(StartSite$chrom_intron,StartSite$start_intron,StartSite$strand_HISAT),paste(StartConserv$chrom_intron,StartConserv$start_intron,StartConserv$strand_HISAT))
StartSite$FirstAppeared=StartConserv$FirstAppeared_age1[mm_Age]
StartSite$Time_FirstAppeared=StartConserv$PhyloDist_FirstAppeared_age1[mm_Age]
StartSite$NbAlignedPrimates=StartConserv$NbAlignedSeqPrim[mm_Age]
StartSite$NbMatchingSeqPrim=StartConserv$NbMatchingSeqPrim[mm_Age]
StartSite$allSequence=StartConserv$allSequence[mm_Age]
StartSite$ancestralSequence=StartConserv$ancestralSequence[mm_Age]
StartSite$newAge=StartConserv$newAge[mm_Age]

StartSite$start_site=StartSite$start_intron
StartSite$end_site=StartSite$start_intron+1



### - strand
EndSite=Junc_annot
#DonorSite_neg$start=DonorSite_neg$start+1 # added 1 to match UCSC coordinates
EndSite$type=ifelse(EndSite$strand_correct=='+','acceptor',ifelse(EndSite$strand_correct=='-','donor',NA))
EndSite$site_id=paste(EndSite$end_id,EndSite$type)

tim=Sys.time()
Tot_count_site = EndSite[, .(Total_count = sum(Total_count, na.rm = T), 
                    Total_count_NS = sum(Total_count_NS, na.rm = T),
                    Total_count_LPS = sum(Total_count_LPS, na.rm = T),
                    Total_count_PAM3CSK4 = sum(Total_count_PAM3CSK4, na.rm = T),
                    Total_count_R848 = sum(Total_count_R848, na.rm = T),
                    Total_count_IAV = sum(Total_count_IAV, na.rm = T)), by = site_id]
print(Sys.time()-tim)
EndSite=EndSite[order(EndSite$site_id,ifelse(is.na(EndSite$Nbread),1,0),-EndSite$Total_count),]
EndSite=EndSite[!duplicated(EndSite$site_id),]
EndSite$Junction_count=EndSite$Total_count
for(i in colnames(Tot_count_site)){
	EndSite[[i]]=Tot_count_site[match(EndSite$site_id,Tot_count_site$site_id),get(i)]
}
EndSite$Gerp_site=EndSite$GerpRS_meanEnd
EndSite$inEnsembl=!is.na(EndSite$Gene_end)
EndSite$gene_Mapped=EndSite$gene
EndSite$coding=EndSite$coding_end
EndSite$gene=EndSite$Gene_end
EndSite$symbol=EndSite$symbol_end

EndSite$sequence=ifelse(EndSite$type=='acceptor',EndSite$acceptor_site,EndSite$donor_site)

EndConserv=as.data.frame(fread(sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteEndConserv_V2.txt',EVO_IMMUNO_POP)))
colnames(EndConserv)[1:12]=c("chrom_intron","end_intron","strand_HISAT","GerpRS_end","Total_count_site","PhylogeneticAge1","PhylogeneticAge2","PhylogeneticAge3","FirstAppeared_age1",'PhyloDist_FirstAppeared_age1','ancestralSequence','allSequence')
	# age 1: first occurence of the current splice site in the ancestral sequences (not allowing missing data or loss and recovery of function)
	# age 2: first occurence of the current splice site in the ancestral sequences (allowing missing data)
	# age 3: first occurence of a splice site in the ancestral sequences (allowing loss and recovery of function, and missing data)-1
	
SummaryStats_end=t(sapply(strsplit(EndConserv$allSequence,':'),function(x){y=gsub('.*_([ATGC-]*)','\\1',x);c(FirstAlignedSeq=max(which(y!='--')),NbAlignedSeq=sum(y!='--'),NbAlignedSeqPrim=sum(y[1:11]!='--'),NbAlignedSeqMam=sum(y[1:18]!='--'),NbMatchingSeq=sum(y==y[1]),NbMatchingSeqPrim=sum(y[1:11]==y[1]),NbMatchingSeqMam=sum(y[1:18]==y[1]),NbObservedSeq=luq(y[y!='--']))}))
SummaryStats_end_age=sapply(strsplit(EndConserv$ancestralSequence,':'),function(x){y=gsub('.*_([ATGC-]*)','\\1',x);max(which(y==y[1]))})
EndConserv=cbind(EndConserv,SummaryStats_end,SummaryStats_end_age)
colnames(EndConserv)[13:21]=c("FirstAlignedSeq","NbAlignedSeq","NbAlignedSeqPrim","NbAlignedSeqMam","NbMatchingSeq","NbMatchingSeqPrim","NbMatchingSeqMam","NbObservedSeq","newAge")
write.table(EndConserv,file=sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteEndConserv_V7.txt',EVO_IMMUNO_POP),sep='\t',col.names=F,quote=F,row.names=F)
EndConserv=as.data.frame(fread(sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteEndConserv_V7.txt',EVO_IMMUNO_POP)))
colnames(EndConserv)[1:12]=c("chrom_intron","end_intron","strand_HISAT","GerpRS_end","Total_count_site","PhylogeneticAge1","PhylogeneticAge2","PhylogeneticAge3","FirstAppeared_age1",'PhyloDist_FirstAppeared_age1','ancestralSequence','allSequence')
colnames(EndConserv)[13:21]=c("FirstAlignedSeq","NbAlignedSeq","NbAlignedSeqPrim","NbAlignedSeqMam","NbMatchingSeq","NbMatchingSeqPrim","NbMatchingSeqMam","NbObservedSeq","newAge")

mm_Age=match(paste(EndSite$chrom,EndSite$end_intron,EndSite$strand_HISAT),paste(EndConserv$chrom_intron,EndConserv$end_intron,EndConserv$strand_HISAT))
EndSite$FirstAppeared=EndConserv$FirstAppeared_age1[mm_Age]
EndSite$Time_FirstAppeared=EndConserv$PhyloDist_FirstAppeared_age1[mm_Age]
EndSite$NbAlignedPrimates=EndConserv$NbAlignedSeqPrim[mm_Age]
EndSite$NbMatchingSeqPrim=EndConserv$NbMatchingSeqPrim[mm_Age]
EndSite$allSequence=EndConserv$allSequence[mm_Age]
EndSite$ancestralSequence=EndConserv$ancestralSequence[mm_Age]
EndSite$newAge=EndConserv$newAge[mm_Age]

EndSite$start_site=EndSite$end_intron-1
EndSite$end_site=EndSite$end_intron

SpliceSites=rbind(StartSite,EndSite)
write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.txt',HOME),quote=F,sep='\t',row.names=F)

############ Annotation is Immune gene 
Require('allGOterms')
GoToTest=list(immuneResponse='GO:0006955',
	defenseResponse='GO:0006952',
	InnateImmuneResponse='GO:0045087',
	AdaptiveImmuneResponse='GO:0002250',
	AntiviralResponse='GO:0051607',
	AntibacterialResponse='GO:0042742',
	TF_DNAbinding='GO:0003700',
	development='GO:0032502',
	olfactory_receptors='GO:0004984')
GoToTest_genes=lapply(GoToTest,function(x){allGOterms$gene[allGOterms$go==x]})
GoToTest_genes$all=unique(allGOterms$gene)
SpliceSites$isImmune=ifelse(SpliceSites$gene%in% GoToTest_genes$immuneResponse,'yes','')
SpliceSites$isImmune_overlappingGene_withStrand=ifelse(SpliceSites$overlappingGene_withStrand%in% GoToTest_genes$immuneResponse,'yes','')

SpliceSites$Quality_MultiZ46Align=FALSE
library(GenomicRanges)
hg19_align=list()
for (chr in 1:22){
cat(chr,'\n')
#Conservation_chr=fread(sprintf('%s/Annotation/Conservation/chr%s_mammalian_conservation_perBp.txt',HOME,chr))
hg19_align[[chr]]=fread(sprintf('%s/Annotation/Conservation/SEqAlign/Alignable_Human/chr%s.align.txt',HOME,chr))
hg19_align[[chr]]$chrom=gsub('hg19.chr','',hg19_align[[chr]]$chrom,fixed=T)
hg19_align[[chr]]=reduce(makeGRangesFromDataFrame(as.data.frame(hg19_align[[chr]]),seqnames='chrom',start.field='start',end.field='end'))
}
hg19_align=do.call(rbind,lapply(hg19_align,as.data.frame))
spliceSites_GR=makeGRangesFromDataFrame(as.data.frame(SpliceSites),seqnames='chrom_intron',start.field='start_site',end.field='end_site',keep.extra=FALSE)
hg19_align_GR=makeGRangesFromDataFrame(hg19_align)
oo=findOverlaps(spliceSites_GR,hg19_align_GR)
SpliceSites$Quality_MultiZ46Align[queryHits(oo)]=TRUE
write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.1.txt',HOME),quote=F,sep='\t',row.names=F)

wIntropolis=which(Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread_Intropolis))

Junc_GR=makeGRangesFromDataFrame(Junc_annot[wIntropolis,],seqnames='chrom_intron',start.field='start_intron',end.field='end_intron')
Sites_GR=makeGRangesFromDataFrame(SpliceSites,seqnames='chrom_intron',start.field='start_intron',end.field='end_intron')
oo=findOverlaps(Sites_GR,Junc_GR) # For each site what are the (known) junctions that overlap the most frequent junction using that splice site
mean(table(SpliceSites[queryHits(oo),'site_id'])==1) # 2.8%, => >97% of splice site have >1 overlapping junction in intropolis

#wTest=1:10000000
wTest=1:length(queryHits(oo))

tim=Sys.time()
Junc_overlapingSite=cbind(site_id=SpliceSites$site_id[queryHits(oo)[wTest]],Junc_annot[wIntropolis[subjectHits(oo)[wTest]],c('Total_count','Total_count_NS','Total_count_LPS','Total_count_PAM3CSK4','Total_count_R848','Total_count_IAV')])

# For each splice site, count the total number of reads overlapping the most frequent junction that uses the splice site
Tot_count_site = Junc_overlapingSite[, .(Total_reads_overlapping_MajorJunction = sum(Total_count, na.rm = T), 
                    Total_reads_overlapping_MajorJunction_NS = sum(Total_count_NS, na.rm = T),
                    Total_reads_overlapping_MajorJunction_LPS = sum(Total_count_LPS, na.rm = T),
                    Total_reads_overlapping_MajorJunction_PAM3CSK4 = sum(Total_count_PAM3CSK4, na.rm = T),
                    Total_reads_overlapping_MajorJunction_R848 = sum(Total_count_R848, na.rm = T),
                    Total_reads_overlapping_MajorJunction_IAV = sum(Total_count_IAV, na.rm = T)), by = site_id]
print(Sys.time()-tim)

SpliceSites[,'Total_reads_overlapping_MajorJunction']=Tot_count_site[match(SpliceSites$site_id,Tot_count_site$site_id),'Total_reads_overlapping_MajorJunction']
SpliceSites[,'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site']=SpliceSites[,'Total_count']/SpliceSites[,'Total_reads_overlapping_MajorJunction']

SpliceSites[,'Total_reads_overlapping_MajorJunction_NS']=Tot_count_site[match(SpliceSites$site_id,Tot_count_site$site_id),'Total_reads_overlapping_MajorJunction_NS']
SpliceSites[,'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_NS']=SpliceSites[,'Total_count_NS']/SpliceSites[,'Total_reads_overlapping_MajorJunction_NS']

SpliceSites[,'Total_reads_overlapping_MajorJunction_LPS']=Tot_count_site[match(SpliceSites$site_id,Tot_count_site$site_id),'Total_reads_overlapping_MajorJunction_LPS']
SpliceSites[,'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_LPS']=SpliceSites[,'Total_count_LPS']/SpliceSites[,'Total_reads_overlapping_MajorJunction_LPS']

SpliceSites[,'Total_reads_overlapping_MajorJunction_PAM3CSK4']=Tot_count_site[match(SpliceSites$site_id,Tot_count_site$site_id),'Total_reads_overlapping_MajorJunction_PAM3CSK4']
SpliceSites[,'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_PAM3CSK4']=SpliceSites[,'Total_count_PAM3CSK4']/SpliceSites[,'Total_reads_overlapping_MajorJunction_PAM3CSK4']

SpliceSites[,'Total_reads_overlapping_MajorJunction_R848']=Tot_count_site[match(SpliceSites$site_id,Tot_count_site$site_id),'Total_reads_overlapping_MajorJunction_R848']
SpliceSites[,'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_R848']=SpliceSites[,'Total_count_R848']/SpliceSites[,'Total_reads_overlapping_MajorJunction_R848']

SpliceSites[,'Total_reads_overlapping_MajorJunction_IAV']=Tot_count_site[match(SpliceSites$site_id,Tot_count_site$site_id),'Total_reads_overlapping_MajorJunction_IAV']
SpliceSites[,'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_IAV']=SpliceSites[,'Total_count_IAV']/SpliceSites[,'Total_reads_overlapping_MajorJunction_IAV']

write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.2.txt',HOME),quote=F,sep='\t',row.names=F)
SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.2.txt',HOME))
##########################################
##		Numbers in manuscript V7		##
##########################################

SpliceSites$isStart=ifelse(gsub(' (donor|acceptor)','',SpliceSites$site_id)==SpliceSites$start_id,TRUE,FALSE)
SpliceSites$isEnd=ifelse(gsub(' (donor|acceptor)','',SpliceSites$site_id)==SpliceSites$end_id,TRUE,FALSE)
SpliceSites$inEnsembl=ifelse(SpliceSites$isStart,!is.na(SpliceSites$Gene_start),!is.na(SpliceSites$Gene_end))
mean(SpliceSites$isStart & SpliceSites$isEnd)
mean(!SpliceSites$isStart & !SpliceSites$isEnd)
mean(SpliceSites$isStart | SpliceSites$isEnd)
SpliceSites$Zero_Gerp_site=FALSE
SpliceSites$Zero_Gerp_site[SpliceSites$isStart & (SpliceSites$GerpRS_start==0 | SpliceSites$GerpRS_start2==0)]=TRUE
SpliceSites$Zero_Gerp_site[SpliceSites$isEnd & (SpliceSites$GerpRS_end==0 | SpliceSites$GerpRS_end2==0)]=TRUE
SpliceSites$AssociatedGene=ifelse((SpliceSites$isStart & !is.na(SpliceSites$Gene_start)) | (SpliceSites$isEnd & is.na(SpliceSites$Gene_end)), SpliceSites$Gene_start, ifelse((SpliceSites$isStart & is.na(SpliceSites$Gene_start)) | (SpliceSites$isEnd & !is.na(SpliceSites$Gene_end)),SpliceSites$Gene_end,NA))
SpliceSites$AssociatedSymbol=ifelse((SpliceSites$isStart & !is.na(SpliceSites$symbol_start)) | (SpliceSites$isEnd & is.na(SpliceSites$symbol_end)), SpliceSites$symbol_start, ifelse((SpliceSites$isStart & is.na(SpliceSites$symbol_start)) | (SpliceSites$isEnd & !is.na(SpliceSites$symbol_end)),SpliceSites$symbol_end,NA))
SpliceSites$AnnotatedGene=ifelse(SpliceSites$isStart & !is.na(SpliceSites$Gene_start), SpliceSites$Gene_start, ifelse(SpliceSites$isEnd & !is.na(SpliceSites$Gene_end),SpliceSites$Gene_end,NA))
SpliceSites$AnnotatedSymbol=ifelse(SpliceSites$isStart & !is.na(SpliceSites$symbol_start), SpliceSites$symbol_start, ifelse(SpliceSites$isEnd & !is.na(SpliceSites$symbol_end),SpliceSites$symbol_end,NA))

#MeanExpr=GeneAnnot[grep('_mean',colnames(GeneAnnot))]
MeanExpr[is.na(MeanExpr)]=0
colnames(MeanExpr)=paste('FPKM',condIndex,sep='_')
lFC=log2(1+MeanExpr[,-1])-log2(1+MeanExpr[,1])%o%c(1,1,1,1)
colnames(lFC)=paste('log2FC',condIndex[-1],sep='_')

Expr_data=cbind(MeanExpr[match(SpliceSites$AssociatedGene,rn(MeanExpr)),],lFC[match(SpliceSites$AssociatedGene,rn(MeanExpr)),],maxFPKM=apply(MeanExpr,1,max)[match(SpliceSites$AssociatedGene,rn(MeanExpr))],maxlFC=apply(lFC,1,max)[match(SpliceSites$AssociatedGene,rn(MeanExpr))],max_abslFC=apply(abs(lFC),1,max)[match(SpliceSites$AssociatedGene,rn(MeanExpr))])
colnames(Expr_data)=paste(colnames(Expr_data),'AssociatedGene',sep='_')
SpliceSites=cbind(SpliceSites,Expr_data)

Expr_data=cbind(MeanExpr[match(SpliceSites$overlappingGene_withStrand,rn(MeanExpr)),],
                lFC[match(SpliceSites$overlappingGene_withStrand,rn(MeanExpr)),],
                maxFPKM=apply(MeanExpr,1,max)[match(SpliceSites$overlappingGene_withStrand,rn(MeanExpr))],
                maxlFC=apply(lFC,1,max)[match(SpliceSites$overlappingGene_withStrand,rn(MeanExpr))],
                max_abslFC=apply(abs(lFC),1,max)[match(SpliceSites$overlappingGene_withStrand,rn(MeanExpr))])
colnames(Expr_data)=paste(colnames(Expr_data),'overlappingGene_withStrand',sep='_')
SpliceSites=cbind(SpliceSites,Expr_data)

Expr_data=cbind(MeanExpr[match(SpliceSites$overlappingGene_noStrand,rn(MeanExpr)),],
                lFC[match(SpliceSites$overlappingGene_noStrand,rn(MeanExpr)),],
                maxFPKM=apply(MeanExpr,1,max)[match(SpliceSites$overlappingGene_noStrand,rn(MeanExpr))],
                maxlFC=apply(lFC,1,max)[match(SpliceSites$overlappingGene_noStrand,rn(MeanExpr))],
                max_abslFC=apply(abs(lFC),1,max)[match(SpliceSites$overlappingGene_noStrand,rn(MeanExpr))])
colnames(Expr_data)=paste(colnames(Expr_data),'overlappingGene_noStrand',sep='_')

SpliceSites=cbind(SpliceSites,Expr_data)

SpliceSites$gene_class_withStrand=ifelse(is.na(SpliceSites$overlappingGene_withStrand),'non genic',ifelse(grepl('//',SpliceSites$overlappingGene_withStrand),'multiple genes',ifelse(SpliceSites$maxFPKM_overlappingGene_withStrand>1,'expressed gene','non expressed gene')))
SpliceSites$gene_class_noStrand=ifelse(is.na(SpliceSites$overlappingGene_noStrand),'non genic',ifelse(grepl('//',SpliceSites$overlappingGene_noStrand),'multiple genes',ifelse(SpliceSites$maxFPKM_overlappingGene_noStrand>1,'expressed gene','non expressed gene')))

write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.3.txt',HOME),quote=F,sep='\t',row.names=F)

rk_time=data.frame(rank=rank(unique(SpliceSites$Time_FirstAppeared)),time=unique(SpliceSites$Time_FirstAppeared))
rk_time=rk_time[order(rk_time$time)[1:18],]
rk_time$Species=c('Humans','Chimpanzees','Gorillas','Orangutans','Gibbons','Monkeys','Tarsiers','Lemurs','Shrews','Rodents','Boroeutheria','Afroeutheria','Marsupials','Monotremes','Sauropsids','Amphibians','Vertebrates','Chordates')
rk_time$group=cut(rk_time$rank,c(-1,5,8,12,13,14,15,16,17,18))
levels(rk_time$group)=c('Apes','Primates','Placentals','Marsupials','Monotremes','Sauropsids','Tetrapods','Vertebrates','Chordates')

# Kumar et al, Nature, 1999 doi:10.1038/31927
rk_time$time_MA_kumar=c(0,5.5,6.7,8.2,14.6,23.3,NA,NA,85,90.8,92,105,173,NA,NA,360,450,564)
# Missing divergence times were inferred by calibration based on phylogenetic divergence.
rk_time$time_MA_kumar=round(approx(rk_time$time,rk_time$time_MA,rk_time$time)$y,1)


write.table(rk_time,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/rk_time_V7.txt',HOME),quote=F,sep='\t',row.names=F)

SpliceSites$oldAge=rk_time$rank[match(SpliceSites$FirstAppeared,rk_time$Species)]
SpliceSites$FirstAppeared_oldAge=SpliceSites$FirstAppeared
SpliceSites$FirstAppeared=NULL
SpliceSites$FirstAppeared_newAge=rk_time$Species[SpliceSites$newAge]
SpliceSites$Time_FirstAppeared_MA_NewAge=rk_time$time_MA_kumar[match(SpliceSites$newAge,rk_time$rank)]
SpliceSites$Time_FirstAppeared_MA_OldAge=rk_time$time_MA_kumar[match(SpliceSites$oldAge,rk_time$rank)]

SS_Activity=SpliceSites[,c("Total_count_NS", "Total_count_LPS", "Total_count_PAM3CSK4", "Total_count_R848", "Total_count_IAV")]/rep(1,nrow(SpliceSites))%o%table(SampleAnnot$cond)
colnames(SS_Activity)=gsub('Total_count','read_persample', colnames(SS_Activity))
SpliceSites=cbind(SpliceSites,SS_Activity)

SS_Pct_reads=SpliceSites[,c('Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_NS',
                            'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_LPS',
                            'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_PAM3CSK4',
                            'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_R848',
                            'Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site_IAV')]

SpliceSites$WeakAltConst=ifelse(apply(SS_Activity<1,1,all) | SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site<0.05,'weak',ifelse(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site >= 0.95,'constitutive','alternative'))

#SpliceSites$FirstAppeared=rk_time$rank[match(SpliceSites$newAge,rk_time$rank)]


write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.4.txt',HOME),quote=F,sep='\t',row.names=F)

