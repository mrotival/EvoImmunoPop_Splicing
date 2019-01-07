

TableS2=PSI_Annot[toKeep,]
# annotate coding consequence
TableS2$coding_type_AS=TableS2$coding_type
TableS2$coding_type_AS[TableS2$coding_type%in%c('gain of function','loss of function')]='gain/loss function'

table(TableS2$coding_type_AS)-
#gain/loss function     modified protein         no coding consequence 
#              6067               9098               1008 
table(TableS2$coding_type_AS)/nrow(TableS2)
#gain/loss function     isoform change         non coding 
#         0.3751314          0.5625425          0.0623261 


#########################################
########### Figure 2A 		#############
#########################################

colPSI=c(A3="#66C2A5", A5="#FC8D62", AF="#8DA0CB", AL="#E78AC3", MX="#A6D854",RI="#FFD92F", SE="#E5C494")
condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
colERC=c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
colERC5=c("#525252AA", "#E31A1CAA", "#33A02CAA", "#1F78B4AA", "#6A3D9AAA")

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/figures/Fig2A.pdf',HOME),height=4,width=4)
par(mar=c(4,4,1,1))
tab=table(PSI_Annot[toKeep,'event_type'])
axis(2)
x=as.numeric(barplot(tab,col=colPSI,ylim=c(0,max(tab)*1.2),las=1),axes=F)
text(x,tab+0.15*max(tab),labels=tab,srt=90)
dev.off()

####################################
########### Figure 2B ##############
####################################
library(pcaMethods)
PCs=pca(t(PSI_prov),scale='uv',nPcs=30)
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/figures/Fig2B_PCA_uv_16173_V5_adjusted.pdf',HOME),width=5.5,height=4.5)
par(mar=c(4,4,1,1))
plot(PCs@scores[,1:2],col=colERC[SampleAnnot$colSamp],pch=16,cex=0.7,las=1)
barplot(PCs@R2)
print(PCs@R2[1:2])
dev.off()


#########################################
########### Figure S2A & B ##############
#########################################

Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V8_withCounts.Rdata',sep=''))

#pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/sup_figures/FigS2A.pdf',HOME),width=5.5,height=6)
layout(1:2)

# for each condition, extract all PSI values of all quantified events (with >30 reads and <5% NA) 
PSI_cond=lapply(1:5,function(i){w=which(PSI_Annot[,paste('JuncCovered',condIndex[i],sep='_')] & PSI_Annot[,grep('PctNA',cn(PSI_Annot))[i]]<0.05)
                                # if PSI >0.5 in that condition, consider the percentage of REF isoform 
								MajorIso=PSI_Annot[w,grep('MeanPSI_',cn(PSI_Annot))[i]]>0.5
								ifelse(MajorIso,1,0)+ifelse(MajorIso,-1,1)*PSI[w,SampleAnnot$cond==i]})
# for each condition, extract all FPKM values (matched for each event)
Count_cond=lapply(1:5,function(i){  w=which(PSI_Annot[,paste('JuncCovered',condIndex[i],sep='_')] & PSI_Annot[,grep('PctNA',cn(PSI_Annot))[i]]<0.05)
									FPKM_gene[match(PSI_Annot$gene_id[w],rn(FPKM_gene)),SampleAnnot$cond==i]})

violinPlot_byGroup=function(x,group,colors=rep('grey',luq(group)),rectCols=rep('black',luq(group)),ylim=c(-max(abs(x),na.rm=T),max(abs(x),na.rm=T)),xlab='',ylab='',main='',splitSign=F,las=1,scaleWidth=T,h=bw.nrd0(x),...){
require(vioplot)
	group=as.factor(group)
	nlevels=length(levels(group))
	plot(1:nlevels,rep(0,nlevels),ylim=ylim,xlim=c(0.5,nlevels+0.5),axes=F,xlab=xlab,ylab=ylab,main=main,col='#00000000')
	if(splitSign){
		tab=table(group,x>0)
	}else{tab=table(group)
		scaleWex=max(sqrt(tab))
	}
	for (j in 1:nlevels){
		if(splitSign){
			wpos=which(group== levels(group)[j] & x>0)
			wneg=which(group== levels(group)[j] & x<0)
			if(scaleWidth){
				myWexPos=sqrt(length(wpos))/scaleWex
				myWexNeg=sqrt(length(wneg))/scaleWex
			}else{
				myWexPos=myWexNeg=1		
			}
			vioplot(x[wpos],at=j,add=T,wex=myWexPos,h=h,col=colors[j],...)
			vioplot(x[wneg],at=j,add=T,wex=myWexNeg,h=h,col=colors[j],...)
			abline(h=0,col='grey')
		}else{
			w=which(group== levels(group)[j] )
			if(scaleWidth){
				myWex=sqrt(length(w))/scaleWex
			}else{
				myWex=1		
			}
			vioplot(x[w],at=j,add=T,wex=myWex,h=h,col=colors[j],rectCol=rectCols[j],...)
			}
		}
	axis(2,las=2)
	axis(1,at=1:nlevels,labels=levels(group),las=las)
}

PSI_cond=unlist(PSI_cond)
Count_cond=unlist(Count_cond)
qCount=quantile(Count_cond,seq(0,1,0.25),na.rm=T)
w=which(!is.na(PSI_cond))
set.seed(0)
w=sample(w,100000)
quartile=cut(Count_cond[w],qCount)
# create violin plot by quartile 
vioPlot_byGroup(0.5-abs(PSI_cond[w]-0.5),quartile,ylim=c(0,0.5),las=3,h=0.05,ylab='Fraction of minor isoform',xlab='',axes=F)
ct=cor.test(0.5-abs(PSI_cond[w]-0.5),Count_cond[w])
ct #cor = -0.21

# for each condition, extract all PSI values of all testable events (with >30 reads and <5% NA, FPKM > 10 & Minor isoform > 5% )
PSI_cond=lapply(1:5,function(i){w=which(PSI_Annot[,grep('Testable_',cn(PSI_Annot))[i]])
								MajorIso=PSI_Annot[w,grep('MeanPSI_',cn(PSI_Annot))[i]]>0.5
#								MajorIso=rep(FALSE,length(w))
								ifelse(MajorIso,1,0)+ifelse(MajorIso,-1,1)*PSI[w,SampleAnnot$cond==i]})
# for each condition, extract all FPKM values (matched for each event)
Count_cond=lapply(1:5,function(i){FPKM_gene[match(PSI_Annot$gene_id[which(PSI_Annot[,grep('Testable_',cn(PSI_Annot))[i]])],rn(FPKM_gene)),SampleAnnot$cond==i]})
PSI_cond=unlist(PSI_cond)
Count_cond=unlist(Count_cond)
w=which(!is.na(PSI_cond))
set.seed(0)
w=sample(w,100000)
qCount=quantile(Count_cond[w],seq(0,1,0.25),na.rm=T)

quartile=cut(Count_cond[w],qCount)
# create violin plot by quartile 
vioPlot_byGroup(0.5-abs(PSI_cond[w]-0.5),quartile,ylim=c(0,0.5),las=3,h=0.05,ylab='Fraction of minor isoform',xlab='',varwidth=F)
ct=cor.test(abs(PSI_cond[w]-0.5),Count_cond[w]) #cor = -0.02883829 
dev.off()

load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/AdjustmentResults/PSI_adjust_all_16173_V5.Rdata',HOME))
range(PSI_prov)


#########################################
########### Figure S2B 		#############
#########################################
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/sup_figures/FigS2B_barplotTestableByCond_event.pdf',HOME),height=3.5*2/3,width=6.5*2/3)
testable_list=sapply(1:5,function(cond){unique(PSI_Annot[intersect(toKeep,which(PSI_Annot[,paste('Testable',condIndex[cond],sep='_')])),'event_id'])})
par(mar=c(4,9,1,1)*2/3)
x=sapply(testable_list,length)
names(x)=condIndex
barplot(rev(x),col=rev(colERC5),horiz=TRUE,las=1,axes=F)
axis(1,at=c(0,3,6,9,12)*1000,labels=c(0,3,6,9,12)*1000)
print(x)
dev.off()

#########################################
########### Figure S2C 		#############
#########################################

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/sup_figures/FigS2C_barplotTestableCondSpecific_gene.pdf',HOME),height=5*2/3,width=3.5*2/3)
par(mar=c(3,4,2,1))
par(xpd=T)
barplot(t(table(apply(PSI_Annot[toKeep,grep('Testable_',cn(PSI_Annot))],1,sum),PSI_Annot[toKeep,grep('Testable_',cn(PSI_Annot))[1]]))[2:1,],las=1,col=c(grey(0.8),colERC5[2]),ylim=c(0,8000),axes=F)
axis(2,at=seq(0,6000,1000),labels=seq(0,6000,1000),las=1)
legend(-2,9000,fill=c(colERC5[2],grey(0.8)),legend=c('stimulation-specific','other'),bty='n')
dev.off()

#########################################
########### Figure S2D 		#############
#########################################

library(VennDiagram)
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/sup_figures/FigS2D_vennTestable_event.pdf',HOME),height=2/3*6,width=2/3*6)
testable_list=sapply(1:5,function(cond){unique(PSI_Annot[intersect(toKeep,which(PSI_Annot[,paste('Testable',condIndex[cond],sep='_')])),'event_id'])})
names(testable_list)=rep('',5)
par(mar=c(1,1,1,1))
VD=venn.diagram(testable_list,col=colERC[2*0:4+2],fill=colERC[2*0:4+1],filename=NULL,margin=0.05,main='TestableGenes')
grid.newpage()
grid.draw(VD)
dev.off()

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/sup_figures/FigS2A_barplotTestableByCond_event.pdf',HOME),height=5.8,width=5.8)
barplot(sapply(testable_list,length),col=colERC5,ylab='Nb Testable',las=2)
dev.off()

#########################################
########### Figure S2E 		#############
#########################################

# /Volumes/@Home/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/PSI_basal_codingconsequences.pdf ---- 9x4 inches
DF=as.data.frame(PSI_Annot[toKeep,c('coding_type','MeanPSI_NS')])
DF$protein=DF$coding_type
DF$PSI=DF$MeanPSI_NS
DF$protein=factor(DF$protein,levels = unique(DF$protein)[c(2,1,4,3)])

color_coding=c("#225EA8","#FEE586","#E31A1C","#A1DAB4")
names(color_coding)=c("no coding consequence","modified protein","loss of function","gain of function")

ggplot(DF,aes(protein,PSI))+geom_violin(aes(fill=protein))+theme_minimal()+geom_boxplot(width=0.05,outlier.shape = NA)+scale_fill_manual(values=color_coding)


##########################################
###		prepare figure S2F          ######
##########################################
load("/Volumes/@home/03_Analysis/Splicing/Papier_Splicing/V5/data/MISO_events_TechnicalCovariatesQuantification_V5.Rdata")
load("/Volumes/@home/03_Analysis/Splicing/Papier_Splicing/V5/data/MISO_events_TechnicalCovariatesQuantification_V5_adj.Rdata")
	
par(mar=c(9,4,4,1))
R2_notAdj=cbind(R2Cont_notAdj,R2Batch_notAdj)
cutoffs=c(0,0.02,0.05,0.1,0.2,0.5,1)
Names=c("population", "LPS", "PAM3CSK4", "R848", "IAV", "total RNA", "RIN", "rRNA ratio", "GC (%)", "Q30 (%)", "library concentration", "5'/3' Bias", "index", "experiment date", "machine", "library batch", "lane", "Sequencing batch")
Levels=levels(cut(R2_notAdj[,1],cutoffs))
R2_group_NotAdj=apply(R2_notAdj,2,cut,cutoffs)
R2_group_NotAdj_counts=apply(R2_group_NotAdj,2,function(x){table(x)[Levels]})[6:1,]
R2_group_NotAdj_counts[is.na(R2_group_NotAdj_counts)]=0
R2_group_NotAdj_counts[6,]=0
R2_group_NotAdj_counts=R2_group_NotAdj_counts[,c(1:5,9:13,15:22)]
colnames(R2_group_NotAdj_counts)=Names

R2_Adj=cbind(R2Cont_Adj,R2Batch_Adj)
R2_group_Adj=apply(R2_Adj,2,cut,cutoffs)
R2_group_Adj_counts=apply(R2_group_Adj,2,function(x){table(x)[Levels]})[6:1,]
R2_group_Adj_counts[is.na(R2_group_Adj_counts)]=0
R2_group_Adj_counts[6,]=0
R2_group_Adj_counts=R2_group_Adj_counts[,c(1:5,9:13,15:22)]
colnames(R2_group_Adj_counts)=Names

##########################################
###		    Figure S2F 	    		######
##########################################
pdf('/Volumes/@home/03_Analysis/Splicing/Papier_Splicing/V5/sup_figures/FigureS2F.pdf',width=10.3,height=5.9)
par(mfrow=c(1,2))
par(mar=c(9,4,3,1))
BuPu=rev(colorRampPalette(brewer.pal(9, "BuPu"))(5))
barplot(R2_group_NotAdj_counts/nrow(R2_group_NotAdj)*100,main='Before Correction',las=2,col=BuPu,ylim=c(0,100),ylab="Percentage of splicing events")
# legend('topleft',fill=BuPu,legend=c('>50% R2','>20% R2','>10% R2','>5%  R2','>2%  R2'))
barplot(R2_group_Adj_counts/nrow(R2_group_Adj)*100,main='After Correction',las=2,col=BuPu,ylim=c(0,100),ylab="Percentage of splicing events")
legend('topright',fill=BuPu,legend=expression(paste(R^2>50,'%'),paste(R^2>20,'%'),paste(R^2>10,'%'),paste(R^2>5,'%'),paste(R^2>2,'%')),bty='n')
dev.off()
