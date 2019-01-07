i0=0
library(RcolorBrewer)
for (j in which(RESobs_nodup_1Mb$event_id%in%c("ENSG00000157601;RI:21:42793765:42793798-42793902:42794022:+","ENSG00000185499;A3:1:155162101-155162577:155162074-155162577:-"))){
i0=i0+1
cat(i0,'')
symbol=RESobs_nodup_1Mb$symbol[j]
type=RESobs_nodup_1Mb$event_type[j]
snp=RESobs_nodup_1Mb$snps[j]
map=getMapInfo(snp)
event_id=RESobs_nodup_1Mb$event_id[j]
gene=RESobs_nodup_1Mb$gene[j]
myExpr=2^FPKM_gene[gene,]-1
myPSI=PSI_prov[PSI_Annot[toKeep,"event_id"]==event_id,]
#ymax=ceiling(max(8,log2(1+myExpr)))
ymax=ceiling(max(log2(1+myExpr)))
ymin=min(log2(1+myExpr))

cat(ymax,ymin,'\n')
#ymax=14
#ymin=0

mySNP=2-getSNP(snp)
# we reverse the SNP to encode allele.1 as the major allele
mySNP_cond=mySNP[substr(names(myPSI),1,6)]
allele.1=map$allele.2 # minor allele 
allele.2=map$allele.1 # major allele
genotypes_alleles=paste(c(allele.1,allele.1,allele.2),c(allele.1,allele.2,allele.2),sep='')
MeanExprGenoCond=By(log2(1+myExpr),paste(SampleAnnot$condition,mySNP_cond),mean)
MeanExprGenoCond=MeanExprGenoCond[!grepl('NA',names(MeanExprGenoCond))]
myMain=paste(symbol,'-',type,'- NS: FPKM =',round(GeneAnnot[gene,'NS_mean']),', Pval=',Format(RESobs_nodup_1Mb$Pvalue_NS[j],2))
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/figures/sQTL_%s_PSI_%s_%s_%s.pdf',HOME,i0,symbol,type,snp),width=7,height=5)
plot(0:1,0:1,ylim=c(-0.15,1),xlim=c(0.5,20),col='#00000000',axes=F,xlab='',ylab='PSI',main=myMain)
axis(2,las=1)
axis(1,at=c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=4+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=8+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=12+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=16+c(1:3),labels=genotypes_alleles,las=2)
for(i in 1:5){
	w=which(SampleAnnot$condition==unique(SampleAnnot$condition)[i])
	boxplot(myPSI[w][mySNP_cond[w]==0],at=(i-1)*4+1,col=colERC5[i],notch=T,pch=16,cex=0.7,outcol='#00000088',add=T,axes=F,pars=list(boxwex=1.5))
	boxplot(myPSI[w][mySNP_cond[w]==1],at=(i-1)*4+2,col=colERC5[i],notch=T,pch=16,cex=0.7,outcol='#00000088',add=T,axes=F,pars=list(boxwex=1.5))
	boxplot(myPSI[w][mySNP_cond[w]==2],at=(i-1)*4+3,col=colERC5[i],notch=T,pch=16,cex=0.7,outcol='#00000088',add=T,axes=F,pars=list(boxwex=1.5))
}
   x=c(1:3,4+1:3,8+1:3,12+1:3,16+1:3);y=-0.1;dx=0.4;dy=0.04
   colors=grey(1-(MeanExprGenoCond-ymin)/(ymax-ymin))
   rect(x+dx,y-dy,x-dx,y+dy,col=colors)
axis(1,at=4+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=8+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=12+c(1:3),labels=genotypes_alleles,las=2)
axis(1,at=16+c(1:3),labels=genotypes_alleles,las=2)	

dev.off()
}