##################################################################################################
## INPUT:                                                                                       ##
## - ExonCoordinates_hg37_ens70.txt: file with coordinates of all exons in Ensemb v70           ##
##                                          with ioverlkapping exons                            ##
## OUTPUT:                                                                                      ##
## - All_Exons_merged_EnsV70.txt: file with coordinates of all exons in Ensemb v70              ##
##                                  when exons overlap, union is considered                     ##
## - All_Introns_merged_EnsV70.txt: file with coordinates of all introns in Ensemb v70          ##
##                                  when exons overlap, union is considered                     ##
##                                                                                              ##
## Calls R packages:                                                                            ##
##   GenomicRanges                                                                              ##
##   data.table                                                                                 ##
##################################################################################################

options(stringsAsFactors=FALSE,max.print=9999)


# local commands
HOME='/pasteur/homes/mrotival'
ExonCoordinates_hg37_ens70=paste(HOME,'/Annotation/ExonCoordinates_hg37_ens70.txt',sep='')
All_Exons_merged_EnsV70=sprintf('%s/Annotation/All_Exons_merged_EnsV70.gff',HOME)
All_Introns_merged_EnsV70=sprintf('%s/Annotation/All_Introns_merged_EnsV70.gff',HOME)

################################################################################################

Exons=as.data.frame(fread(ExonCoordinates_hg37_ens70,header=T))
colnames(Exons)=make.names(colnames(Exons))
rownames(Exons)=paste(Exons$Ensembl.Exon.ID,Exons$Ensembl.Transcript.ID)

library(GenomicRanges)

Exons$Strand=ifelse(Exons$Strand>0,'+',ifelse(Exons$Strand<0,'-','*'))
Exons=Exons[Exons$Chromosome.Name%in%c(1:22,'X','Y','MT'),]
Exon_GR=makeGRangesFromDataFrame(Exons,start.field='Exon.Chr.Start..bp.',end.field='Exon.Chr.End..bp.',seqnames='Chromosome.Name',keep.extra=T)
Exon_GR=Exon_GR[Exon_GR$Gene.Biotype=='protein_coding',]

w=1:length(Exon_GR)

tim=Sys.time()
EI=By(w,Exon_GR$Ensembl.Gene.ID[w],function(i){
    Exon_reduced_GR=reduce(Exon_GR[i,])
    intron_reduced_GR=gaps(Exon_reduced_GR)
    list(Exon_reduced_GR,intron_reduced_GR)
    })
print(Sys.time()-tim)
ExonRegion=sapply(EI,function(x){x[[1]]})
IntronRegion=sapply(EI,function(x){x[[2]][-1,]})

All_Introns=lapply(IntronRegion,as.data.frame)
All_Introns=do.call(rbind,All_Introns[sapply(IntronRegion,length)>0])
All_Introns$gene=rep(names(IntronRegion),sapply(IntronRegion,length))
All_Introns$ID=paste(All_Introns$gene,repseq(sapply(IntronRegion,length)[sapply(IntronRegion,length)>0]),sep='_')
All_Introns$symbol=GeneAnnot$Associated.Gene.Name[match(All_Introns$gene,GeneAnnot$Ensembl.Gene.ID)]
All_Introns$source='Ensemblv70'
All_Introns$feature='Intron'
All_Introns$score='.'
All_Introns$frame='.'
All_Introns$attribute=sprintf('gene_id=%s; intron_id=%s; symbol=%s',All_Introns$gene,All_Introns$ID,All_Introns$symbol)
colnames(All_Introns)[1]='seqname'
write.table(All_Introns[c("seqname",'source','feature','start','end','score','strand','frame','attribute')],file=All_Introns_merged_EnsV70,quote=F,row.names=F,sep='\t',col.names=F)

All_Exons=lapply(ExonRegion,as.data.frame)
All_Exons=do.call(rbind,All_Exons[sapply(ExonRegion,length)>0])
All_Exons$gene=rep(names(ExonRegion),sapply(ExonRegion,length))
All_Exons$ID=paste(All_Exons$gene,repseq(sapply(ExonRegion,length)[sapply(ExonRegion,length)>0]),sep='_')
All_Exons$symbol=GeneAnnot$Associated.Gene.Name[match(All_Exons$gene,GeneAnnot$Ensembl.Gene.ID)]
All_Exons$source='Ensemblv70'
All_Exons$feature='Exon'
All_Exons$score='.'
All_Exons$frame='.'
All_Exons$attribute=sprintf('gene_id=%s; exon_id=%s; symbol=%s',All_Exons$gene,All_Exons$ID,All_Exons$symbol)
colnames(All_Exons)[1]='seqname'
write.table(All_Exons[c("seqname",'source','feature','start','end','score','strand','frame','attribute')],file=All_Exons_merged_EnsV70,quote=F,row.names=F,sep='\t',col.names=F)

#    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#    source - name of the program that generated this feature, or the data source (database or project name)
#    feature - feature type name, e.g. Gene, Variation, Similarity
#    start - Start position of the feature, with sequence numbering starting at 1.
#    end - End position of the feature, with sequence numbering starting at 1.
#    score - A floating point value.
#    strand - defined as + (forward) or - (reverse).
#    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.


load("/Volumes/@Home/02_ERC_Data/SampleAnnot_withFlu_allSamples.Rdata")
SampleAnnot$library_Batch=paste(SampleAnnot$library_date,SampleAnnot$library_time)
SampleAnnot$experiment_date=SampleAnnot$date
SampleAnnot$experiment_RNAconcentration_ng_per_ul=SampleAnnot$ng_ul
SampleAnnot$experiment_rRNAratio_260_280=SampleAnnot$X260_280
SampleAnnot$experiment_rRNAratio_260_230=SampleAnnot$X260_230
SampleAnnot$experiment_rin=SampleAnnot$rin

SampleAnnot$condition=condIndex[SampleAnnot$condition]
SampleAnnot$population=substr(SampleAnnot$individu,1,3)


SampleAnnot=SampleAnnot[match(colnames(FPKM_gene),SampleAnnot$sample_ID),c( "sample_ID","individu", "condition", "population",
"experimentator","experiment_date", "experiment_RNAconcentration_ng_per_ul","experiment_rRNAratio_260_280", "experiment_rRNAratio_260_230", "experiment_rin", "quantity_sent", "shipment_date", "batch", 
"total_amount_axeq", "RIN_axeq", "rRNA_ratio_axeq", "results_axeq", "failed", 
"machine", "lane", "index", "index2", "total_bases", "read_count", "GC_pct", "Q20_pct", "Q30_pct", 
"machine_repeat", "lane_repeat", "index_repeat", "index2_repeat", "total_bases_repeat", "read_count_repeat", "GC_pct_repeat", "Q20_pct_repeat", "Q30_pct_repeat", 
"sum_all_runs", "rerun", "new_library", "library_Batch", "library_mass_conc", "library_molar_conc", "fragment_size", "GC_pct_final", "FivePrimeBias", "CoverageRegularity", 
"unmapped", "mapped_flu", "flu_pct_tot", "flu_FPKM")]
write.table(SampleAnnot, file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/SampleAnnot.txt',sep='\t',quote=F,row.names=F,col.names=T)

GeneAnnot=read.table(paste(HOME,'/Annotation/GeneAnnotation_hg37_ens70.txt',sep=''),sep='\t',comment='',quote='',stringsAsFactors=FALSE,header=T,row.names=1)
colnames(GeneAnnot)[1]='Ensembl.Gene.ID'
	GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='X']=23
	GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='Y']=24
	GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='MT']=26
	GeneAnnot$Chromosome.Name=as.numeric(GeneAnnot$Chromosome.Name)
    GeneAnnot$select[is.na(GeneAnnot$select)]=FALSE
    colnames(GeneAnnot)=c("Ensembl.Gene.ID", "Chromosome.Name", "Strand", "Gene Start (bp)","Gene End (bp)", "Associated Gene Name", "Description", "Gene Biotype","Expressed", "NS_mean", "LPS_mean", "PAM3CSK4_mean", "R848_mean", "IAV_mean","Nb_transcript")
    for(i in condIndex){
        GeneAnnot[[paste(i,'mean',sep='_')]][is.na(GeneAnnot[[paste(i,'mean',sep='_')]])]=0
        }
write.table(GeneAnnot, file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/GeneAnnotation_hg37_ens70.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(GeneAnnot[GeneAnnot$Expressed,], file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/GeneAnnot_expressed.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(FPKM_gene, file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt',sep='\t',quote=F,row.names=F,col.names=T)

TranscriptAnnot=read.table(paste(HOME,'/Annotation/TranscriptAnnotation_hg37_ens70.txt',sep=''),sep='\t',comment='',quote='',stringsAsFactors=FALSE)
TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='X']=23
TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='Y']=24
TranscriptAnnot$Chromosome.Name[TranscriptAnnot$Chromosome.Name=='MT']=26
TranscriptAnnot$Chromosome.Name=as.numeric(TranscriptAnnot$Chromosome.Name)

write.table(TranscriptAnnot, file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/TranscriptAnnotation_hg37_ens70.txt',sep='\t',quote=F,row.names=F,col.names=T)
isoAnnot=isoAnnot[,c("tracking_id", "gene_id", "gene_short_name", "length", "desc", "select", "NS_mean", "LPS_mean", "PAM3_mean", "R848_mean", "Flu_mean")]
colnames(isoAnnot)=c("transcript_id","gene_id", "symbol", "length", "description", "Expressed", "NS_mean", "LPS_mean", "PAM3CSK4_mean", "R848_mean", "IAV_mean")
write.table(isoAnnot, file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/02_TranscriptFPKM_cufflinks/TranscriptAnnot_expressed.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(FPKM, file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/02_TranscriptFPKM_cufflinks/FPKM_matrix.txt',sep='\t',quote=F,row.names=F,col.names=T)

Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V7.2_withCounts.Rdata',sep=''))
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME))
PSI=PSI_prov
PSI_Annot=PSI_Annot[toKeep,c("event_id","event_type","chrom","start","end","strand","gene_id","symbol","MeanPSI_NS","MeanPSI_LPS","MeanPSI_PAM3CSK4","MeanPSI_R848","MeanPSI_IAV","Support_NS","Support_LPS","Support_PAM3CSK4","Support_R848","Support_IAV","coding_type","conserved_site_GerpRS")]
colnames(PSI_Annot)=gsub('Support','FPKM',colnames(PSI_Annot))
write.table(PSI_Annot, file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/03_ASevents_MISO/PSI_Annot_frequent.txt',sep='\t',quote=F,row.names=F,col.names=T)
rownames(PSI)=PSI_Annot$event_id
write.table(PSI, file='/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/03_ASevents_MISO/PSI_matrix.txt',sep='\t',quote=F,row.names=F,col.names=T)

