evo_immuno_pop='/Volumes/evo_immuno_pop'
SampleFile=sprintf('%s/Maxime/Splicing/HISAT2/SampleFile.txt',evo_immuno_pop)

TotalCoverage=0 # used to be 200/ switched to 30 just now...
L=100000 # max intron length to consider a junction

library(data.table)
AllSamples=fread(SampleFile,col.names=c('samples','batch','merge'))

i=0
junctions=list()
for (sample in AllSamples[,get('samples')]){
    i=i+1
    cat(i,sample,'\n')
	juncFile=sprintf('%s/Maxime/Splicing/HISAT2/Results/%s/Junctions_%s.junc',evo_immuno_pop,sample,sample)
	junctions_samp=fread(juncFile,col.names=c('chrom','start','end','dot','count','strand'))
	##	add a sample name column 
	junctions_samp[,sample:=sample]
	##	filter to keep only relevant junctions 
	junctions[[i]] = junctions_samp
}
	
junctions=rbindlist(junctions)
junctions[,start:=start+1]
junctions[,start_id:=paste(chrom,start,strand)]
junctions[,end_id:=paste(chrom,end,strand)]
junctions[,junc_id:=paste(chrom,start,end,strand)]
junctions[,ind:=substr(sample,1,6)]
condIndex=c('NS','LPS','PAM3CSK4','R848','IAV')
junctions[,cond:=factor(condIndex[as.numeric(substr(sample,8,8))],levels=condIndex)]


outFile=sprintf('%s/Maxime/Splicing/HISAT2/Results/All_junctions_bySample.txt',evo_immuno_pop)
write.table(junctions,file=outFile,sep='\t',row.names=FALSE,col.names=T,quote=F)

# Compute total counts for each condition and add these infos.

junctions_ind=dcast(junctions,chrom+start+end+strand+start_id+end_id+junc_id+ind~cond,value.var='count',fill=0)
junctions_info=junctions_ind[!duplicated(junc_id),mget(c('chrom','start','end','strand','start_id','end_id','junc_id'))]

colnames(junctions_ind)[colnames(junctions_ind)%in%condIndex]=paste('count',condIndex,sep='_')
junctions_ind[,count_Total:=count_NS+count_LPS+count_PAM3CSK4+count_R848+count_IAV]
junctions_byCond=junctions_ind[,.(Total_count_NS=sum(count_NS),
                 Total_count_LPS=sum(count_LPS),
                 Total_count_PAM3CSK4=sum(count_PAM3CSK4),
                 Total_count_R848=sum(count_R848),
                 Total_count_IAV=sum(count_IAV),
                 Total_count=sum(count_Total))
                ,by=junc_id]
junctions_byCond=merge(junctions_info,junctions_byCond,by='junc_id')
outFile=sprintf('%s/Maxime/Splicing/HISAT2/Results/All_junctions_aggregated_bycond.txt',evo_immuno_pop)
write.table(junctions_byCond,file=outFile,sep='\t',row.names=FALSE,col.names=T,quote=F)

