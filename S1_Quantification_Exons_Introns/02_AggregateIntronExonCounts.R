##################################################################################################
## INPUT:                                                                                       ##
## - ${SAMPLE}_Intron_read_count.txt: file with HTseq counts of all introns in Ensemb v70       ##
## - ${SAMPLE}_Exon_read_count.txt: file with HTseq counts of all exons in Ensemb v70           ##
## - SampleAnnot.txt : file with information of RNA seq samples (ID, Nb of Reads)
## OUPUT                                                                                        ##
## - AllSample_Exon_read_count.txt: file with HTseq counts of all exons in Ensemb v70           ##
## - AllSample_Intron_read_count.txt: file with HTseq counts of all exons in Ensemb v70         ##
##                                                                                              ##
## RUN after ./align_and_counts_exon_intron_tars.sh                                             ##
##                                                                                              ##
## Calls R packages:                                                                            ##
##   data.table                                                                                 ##
##################################################################################################

library(data.table)
PATH=sprintf("%s/Maxime/Splicing/HISAT2/Results/",EVO_IMMUNO_POP)

All_Exons_merged_EnsV70=sprintf('%s/Annotation/All_Exons_merged_EnsV70.gff',HOME)
All_Introns_merged_EnsV70=sprintf('%s/Annotation/All_Introns_merged_EnsV70.gff',HOME)
#SampleAnnot=

All_Exons=as.data.frame(fread(All_Exons_merged_EnsV70))
colnames(All_Exons)=c("seqname",'source','feature','start','end','score','strand','frame','attribute')
All_Exons$exon_id=sapply(strsplit(All_Exons$attribute,'; '), function(x) gsub('exon_id=','',x[2]) )

All_Introns=as.data.frame(fread(All_Introns_merged_EnsV70))
colnames(All_Introns)=c("seqname",'source','feature','start','end','score','strand','frame','attribute')
All_Introns$intron_id=sapply(strsplit(All_Introns$attribute,'; '), function(x) gsub('intron_id=','',x[2]) )

i=0
for (samp in SampleAnnot$sample_ID){
    cat(i,'')
    introns=as.data.frame(fread(sprintf('%s/%s/%s_Intron_read_count.txt',PATH,samp,samp)))
    exons=as.data.frame(fread(sprintf('%s/%s/%s_Exon_read_count.txt',PATH,samp,samp)))
    if(i==0){
        Introns=introns
        colnames(Introns)[1:2]=c("intron_id",paste(samp,'count',sep='_'))
        Exons=exons
        colnames(Exons)[1:2]=c("exon_id",paste(samp,'count',sep='_'))
        Introns$length = (All_Introns$end-All_Introns$start+1)[match(Introns$intron_id,All_Introns$intron_id)]
        Exons$length = (All_Exons$end-All_Exons$start+1)[match(Exons$exon_id,All_Exons$exon_id)]
    }else{
        Introns[[paste(samp,'count',sep='_')]]=introns[[2]]
        Exons[[paste(samp,'count',sep='_')]]=exons[[2]]
    }
    Introns[[paste(samp,'RPKM',sep='_')]]=Introns[[paste(samp,'count',sep='_')]]/(Introns$length/1000)/(SampleAnnot$sum_all_runs[match(samp,SampleAnnot$sample_ID)]/1e6)
    Exons[[paste(samp,'RPKM',sep='_')]]=Exons[[paste(samp,'count',sep='_')]]/(Exons$length/1000)/(SampleAnnot$sum_all_runs[match(samp,SampleAnnot$sample_ID)]/1e6)
    i=i+1
}

par(mar=c(7,4,1,1))
layout(1:2)
hist(log2(Exons$'AFB011-1_RPKM'),col='grey',xlim=c(-15,15),breaks=seq(-15,15,1))
hist(log2(Introns$'AFB011-1_RPKM'),col='grey',xlim=c(-15,15),breaks=seq(-15,15,1))

Introns$Gene=sapply(strsplit(Introns$intron_id,'_'),function(x) x[1] )
Introns$previousExon=Introns$intron_id
Introns$nextExon=paste(Introns$Gene, sapply(strsplit(Introns$intron_id,'_'),function(x) as.numeric(x[2])+1 ),sep='_')

ByGene_Introns=Introns
ByGene_Introns=ByGene_Introns[order(ByGene_Introns$Gene),]
ByGene_Introns=ByGene_Introns[!duplicated(ByGene_Introns$Gene),]
ByGene_Introns$length=By(Introns$length,Introns$Gene,sum)
i=0
for (samp in SampleAnnot$sample_ID){
    cat(i,'')
    ByGene_Introns[[paste(samp,'count',sep='_')]]=By(Introns[[paste(samp,'count',sep='_')]],Introns$Gene,sum)
    ByGene_Introns[[paste(samp,'RPKM',sep='_')]]=ByGene_Introns[[paste(samp,'count',sep='_')]]/(ByGene_Introns$length/1000)/(SampleAnnot$sum_all_runs[match(samp,SampleAnnot$sample_ID)]/1e6)
    i=i+1
}

Exons$Gene=sapply(strsplit(Exons$exon_id,'_'),function(x) x[1] )
ByGene_Exons=Exons
ByGene_Exons=ByGene_Exons[order(ByGene_Exons$Gene),]
ByGene_Exons=ByGene_Exons[!duplicated(ByGene_Exons$Gene),]
ByGene_Exons$length=By(Exons$length,Exons$Gene,sum)
i=0
for (samp in SampleAnnot$sample_ID){
    cat(i,'')
    ByGene_Exons[[paste(samp,'count',sep='_')]]=By(Exons[[paste(samp,'count',sep='_')]],Exons$Gene,sum)
    ByGene_Exons[[paste(samp,'RPKM',sep='_')]]=ByGene_Exons[[paste(samp,'count',sep='_')]]/(ByGene_Exons$length/1000)/(SampleAnnot$sum_all_runs[match(samp,SampleAnnot$sample_ID)]/1e6)
    i=i+1
}
write.table(Introns,file=sprintf('%s/AllSamples_Intron_read_count.txt',PATH),quote=F,sep='\t',row.names=F)
write.table(Exons,file=sprintf('%s/AllSamples_Exon_read_count.txt',PATH),quote=F,sep='\t',row.names=F)
write.table(ByGene_Introns,file=sprintf('%s/AllSamples_Intronic_Gene_count.txt',PATH),quote=F,sep='\t',row.names=F)
write.table(ByGene_Exons,file=sprintf('%s/AllSamples_Exonic_Gene_count.txt',PATH),quote=F,sep='\t',row.names=F)
