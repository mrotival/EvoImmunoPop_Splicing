Junc_annot=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Junction_annot_V7.2.txt',HOME))
SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V7.4.txt',HOME))

# create supplementary figures S1A
wAnnot=which(Junc_annot$Annotated_junc)
wCanonical=which(Junc_annot$DA_type!='other')
wIntro_strand=which(Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread_Intropolis) & !is.na(Junc_annot$overlappingSymbol_withStrand))

PctReads=c(All=1,
Canonical=sum(as.numeric(Junc_annot$Total_count[wCanonical]))/sum(as.numeric(Junc_annot$Total_count)), #0.9991252
Intropolis=sum(as.numeric(Junc_annot$Total_count[wIntro_strand]))/sum(as.numeric(Junc_annot$Total_count)), #0.9984076
Ensembl=sum(as.numeric(Junc_annot$Total_count[wAnnot]))/sum(as.numeric(Junc_annot$Total_count))) #0.9772441

Nbjunc=c(All=nrow(Junc_annot),Canonical=length(wCanonical),Intropolis=length(wIntro_strand),Ensembl=length(wAnnot))

tab_total=table(SpliceSites$sequence[SpliceSites$sequence!='NN'])
Pct_nonCano_total=sum(tab_total[!names(tab_total)%in%c('AG','GT')])/sum(tab_total)

tab_intro=table(SpliceSites$sequence[!is.na(SpliceSites$Nbread_Intropolis) & SpliceSites$sequence!='NN' & !is.na(SpliceSites$overlappingGene_withStrand)])
Pct_nonCano_intro=sum(tab_intro[!names(tab_intro)%in%c('AG','GT')])/sum(tab_intro)

tab_annot=table(SpliceSites$sequence[SpliceSites$inEnsembl & !is.na(SpliceSites$Nbread_Intropolis) & SpliceSites$sequence!='NN' & !is.na(SpliceSites$overlappingGene_withStrand)])
Pct_nonCano_annot=sum(tab_annot[!names(tab_annot)%in%c('AG','GT')])/sum(tab_annot)

Pct_nonCanonical_SpliceSites_Intropolis=c()
for (j in 1:6){
    mycol=c('Total_count', 'Total_count_NS','Total_count_LPS','Total_count_PAM3CSK4','Total_count_R848','Total_count_IAV')[j]
    tab1=table(SpliceSites$sequence[SpliceSites$sequence!='NN' & !is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand) & SpliceSites[[mycol]]>0])
    Pct_nonCanonical_SpliceSites_Intropolis[j]=sum(tab1[!names(tab1)%in%c('AG','GT')])/sum(tab1)
    }
conservation_col=c(Neutral="#41B6C480", Weak="#FECC5CBB", Strong="#E31A1CBB")
splicing_col=c(Weak="#FEF9A9", Alternative="#A1DAB4", Constitutive="#225EA8")
condIndex=c("NS", "LPS", "PAM3CSK4", "R848", "IAV")
cond_color=c(NS="#525252AA", LPS="#E31A1CAA", PAM3CSK4="#33A02CAA", R848="#1F78B4AA", IAV="#6A3D9AAA")
database_col=c(All="lightgrey",Intropolis="#FECC5C",Ensembl="#225EA8")
PhyloAge_col <- colorRampPalette(brewer.pal(8,"BuPu"))(10)[2:10]

##########################################
### 	 Fig S1 impact of filters      ###
##########################################

##########################################
### 	 Fig S1A impact of filters      ## 5x 4 inches
##########################################
# Fig S1Abis_SpliceSite Sequence 4x5.5 inches
par(mar=c(6,4,1,4))
x=barplot(Nbjunc/1e6,las=2,ylab='discovered junctions (millions)')
y=(1-PctReads)*100
Factor=max(Nbjunc)/1e6/3
points(x,y*Factor,col='red',pch=16,type='b')
axis(4,at=seq(0,max(round(y)+0.5)*Factor,l=6),labels=seq(0,max(round(y))+0.5,l=6),las=2)
mtext('% of spliced reads removed',4,3)

##########################################
###   END  Fig S1A impact of filters    ###
##########################################

##########################################
### 	 Fig S1B impact of filters      ## 3.6x7.6 inches
##########################################
layout(matrix(rep(1:2,c(4,3)),nrow=1))
# Fig S1Bbis_SpliceSite Sequence 3.6x7.6 inches
#layout(matrix(rep(1:2,c(4,3)),nrow=1))
par(mar=c(6,6,1,0.1))
tab_total=rev(sort(tab_total))
barplot(tab_total/1000,col='lightgrey',las=2,cex.axis=1.3)
barplot(tab_intro[names(tab_total)]/1000,add=T,col=acol[3],axes=F,las=2,axisname=F)
barplot(tab_annot[names(tab_total)]/1000,add=T,col=acol[8],axes=F,las=2,axisname=F)
legend('topright',legend=c('All','Intropolis','Ensembl'),fill=database_col,bty='n',cex=1.2)

par(mar=c(6,6,1,1))
a=0.3;b=1.8
names(Pct_nonCanonical_SpliceSites_Intropolis)=c('All',condIndex)
y=c(Pct_nonCano_total,Pct_nonCano_intro,Pct_nonCano_annot,Pct_nonCanonical_SpliceSites_Intropolis[-1])*100
names(y)=c('All','Intropolis','Ensembl',condIndex)
barplot(y,col=c(database_col,cond_color),las=2,ylab='Percentage of non-canonical\n splice sites',cex.lab=1.2,space=c(a,a,a,b,a,a,a,a),cex.axis=1.3)

##########################################
###   END  Fig S1B impact of filters    ###
##########################################

sum(Junc_annot$Total_count/1e9) # 7.234256 -> 7.2 billion spliced reads
sum(Junc_annot$Total_count[is.na(Junc_annot$Nbread_Intropolis)]/1e9)/sum(Junc_annot$Total_count/1e9)*100 # 0.159241 -> 0.2 % reads are removed
sum(Junc_annot$Total_count[is.na(Junc_annot$Nbread_Intropolis) | is.na(Junc_annot$overlappingGene_noStrand)]/1e9)/sum(Junc_annot$Total_count/1e9)*100 # 0.2934486 -> 0.3 % reads are removed
sum(Junc_annot$Total_count[is.na(Junc_annot$Nbread_Intropolis) | is.na(Junc_annot$overlappingGene_withStrand)]/1e9)/sum(Junc_annot$Total_count/1e9)*100 # 0.4090983 -> 0.4 % reads are removed
# note there is a discordance because overlappingGene does not consider the strand for the junctions.

nrow(Junc_annot) # 4 472 078 junctions discovered
sum(!is.na(Junc_annot$Nbread_Intropolis))/1e6 # 2.099531 ->2.1 million junctions analysed
sum(!is.na(Junc_annot$Nbread_Intropolis) & !is.na(Junc_annot$overlappingGene_noStrand))/1e6 # 1.932976 ->  >1.9 million junctions analysed
sum(!is.na(Junc_annot$Nbread_Intropolis) & !is.na(Junc_annot$overlappingGene_withStrand))/1e6 # 1.734506 ->  >1.7 million junctions analysed
sum(!is.na(Junc_annot$Nbread_Intropolis) & !is.na(Junc_annot$overlappingGene_withStrand))/1e6 # 1.734506 ->  >1.7 million junctions analysed

nrow(SpliceSites) # 5 049 210 over 1.9 million
nrow(SpliceSites[!is.na(SpliceSites$Nbread_Intropolis)])  # 1983844 over 1.9 million splice sites
nrow(SpliceSites[!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_noStrand)])  # 1779666 over 1.7 million splice sites
nrow(SpliceSites[!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand)])  # 1517297 over 1.5 million splice sites

wIntro=which(!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand))
length(wIntro)
# 1517297

wIntro_AssignedUnambiguously=which(!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand) & !grepl('//',SpliceSites$overlappingGene_withStrand))
length(wIntro_AssignedUnambiguously)
#1389163

wIntro_ambigousMapped=!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand) & grepl('//',SpliceSites$overlappingGene_withStrand)
sum(wIntro_ambigousMapped)
#128134

wIntro_Expressed=!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene_withStrand) & !grepl('//',SpliceSites$overlappingGene_withStrand) & SpliceSites$maxFPKM_overlappingGene_withStrand>1
sum(wIntro_Expressed)
#1106965

#table(SpliceSites[,'WeakAltConst'],exclude='')
tabWeakAltConst=table(SpliceSites[wIntro_Expressed,'WeakAltConst'],exclude='')
tabWeakAltConst/sum(tabWeakAltConst)
# alternative constitutive         weak 
#  0.05063304   0.14611392   0.80325304 
sum(tabWeakAltConst[c('alternative','constitutive')])
SpliceSites=SpliceSites[wIntro_Expressed,]
wGerp=!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align
mean(SpliceSites$Gerp_site[which(wGerp)]<2)
# 0.6006095

####################################################################################
# Fig 1A: 6x 5.5 inches - Conservation of Splice sites and constitutive Splice sites
####################################################################################
library(data.table)
Gerp_mean=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_2bp_random.txt',HOME))$V1
DF=rbind(data.frame(Type=rep('Genome Wide',length(Gerp_mean[Gerp_mean!=0])),GerpRS=Gerp_mean[Gerp_mean!=0]),
         data.frame(Type=rep('All splice sites',sum(wGerp)),GerpRS=SpliceSites$Gerp_site[wGerp]),
         data.frame(Type=rep('Weak',sum(SpliceSites$WeakAltConst[wGerp]=='weak')),GerpRS=SpliceSites$Gerp_site[wGerp & SpliceSites$WeakAltConst=='weak']),
         data.frame(Type=rep('Alternative',sum(SpliceSites$WeakAltConst[wGerp]=='alternative')),GerpRS=SpliceSites$Gerp_site[wGerp & SpliceSites$WeakAltConst=='alternative']),
         data.frame(Type=rep('Constitutive',sum(SpliceSites$WeakAltConst[wGerp]=='constitutive')),GerpRS=SpliceSites$Gerp_site[wGerp & SpliceSites$WeakAltConst=='constitutive']))
DF$Type=factor(DF$Type,levels=c('Genome Wide','All splice sites','Weak','Alternative','Constitutive'))
library(ggplot2)

pdf('/Volumes/@Home/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/Fig1A_violin.pdf',width=10,height=3.5)
p<-ggplot(DF,aes(x=Type,y=GerpRS, fill=Type)) + geom_violin(width=1.3) + theme_minimal()+scale_fill_manual(values=c('Genome Wide'="#41B6C4",'All splice sites'="#FECC5C",splicing_col[1:3]))+ylim(low=-6, high=6)
p
dev.off()

pdf('/Volumes/@Home/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/Fig1A_density_small.pdf',width=6,height=3)
p <- ggplot(DF,aes(x=GerpRS, fill=Type,col=Type)) + geom_density() + theme_minimal() + xlim(low=-6, high=6)
p <- p + scale_color_manual(values=c('Genome Wide'="#41B6C4",'All splice sites'="#FECC5C",'Weak'="#BEB969",'Alternative'="#81BA94",'Constitutive'="#335EA8"))
p <- p + scale_fill_manual(values=c(alpha(c('Genome Wide'="#41B6C4",'All splice sites'="#FECC5C"),.5),alpha(splicing_col[1:3],c(.4,.3,.2))))
p+facet_grid(!Type%in%c('Genome Wide','All splice sites')~.)
dev.off()

pdf('/Volumes/@Home/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/NewFigures/Fig1A_density_small_horiz.pdf',width=7,height=2.5)
p+facet_grid(.~!Type%in%c('Genome Wide','All splice sites'))
dev.off()
##########################################
###   END  Fig 1A                      ###
##########################################

SpliceSites[wGerp,.(conserved=mean(Gerp_site>2)),by=WeakAltConst]
#  WeakAltConst conserved
#1:         weak 0.2547882
#2: constitutive 0.9772301
#3:  alternative 0.7194419

mean(!SpliceSites$inEnsembl[wGerp & SpliceSites$Gerp_site<2 & SpliceSites$WeakAltConst=='weak'])
# 0.9519104


#######################################################################################################
#### Fig 1B : 4x 5.5 inches - Conservation of weak constitutive and alternative Splice site  ######
#######################################################################################################

par(mar=c(4,5,1,0.5))
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align)
x=cut(SpliceSites$Gerp_site[w],c(-Inf,2,4,Inf))
levels(x)=c('<2','[2;4]','>4')
par(mar=c(4,5,1,0.5))
tab=table(x,SpliceSites$WeakAltConst[w])
tab=tab[,c('weak','alternative','constitutive')]
colnames(tab)=c('Weak','Alternative','Constitutive')
t(t(tab)/apply(tab,2,sum))*100
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),ylab='% of splice sites',las=1)
par(xpd=T)
legend(-2,-15,legend=levels(x),fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),bty='n',ncol=4)
par(xpd=F)
##########################################
###   END  Fig 1B                      ###
##########################################

#########################################################################################################################
#### Fig S1D : 6x 4.5 inches - coverage of conserved and non conserved splice sites by levels of gene expression   ######
#########################################################################################################################
lowExpressed=SpliceSites$maxFPKM_overlappingGene_withStrand>1 & SpliceSites$maxFPKM_overlappingGene_withStrand<10
midExpressed=SpliceSites$maxFPKM_overlappingGene_withStrand>10 & SpliceSites$maxFPKM_overlappingGene_withStrand<100
highExpressed=SpliceSites$maxFPKM_overlappingGene_withStrand>100
layout(1:3)
layout(matrix(1:3,1))
breaks=seq(-10,17,1)
SS_Activity=as.data.frame(SpliceSites[,c("read_persample_NS","read_persample_LPS","read_persample_PAM3CSK4","read_persample_R848","read_persample_IAV")])
w=which( wGerp & lowExpressed)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=breaks,plot=F)
hh$counts=hh$count/1000
w=which(wGerp & lowExpressed & SpliceSites$Gerp_site>2)
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(wGerp & lowExpressed & SpliceSites$Gerp_site>4 )
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=conservation_col['Neutral'],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=conservation_col['Weak'],add=T,border='#000000')
plot(hh2,col=conservation_col['Strong'],add=T,border='#000000')
#par(xpd=T);legend('topright',fill=conservation_col,legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(wGerp & midExpressed)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=breaks,plot=F)
hh$counts=hh$count/1000
w=which(wGerp & midExpressed & SpliceSites$Gerp_site>2 )
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(wGerp & midExpressed & SpliceSites$Gerp_site>4 )
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=conservation_col['Neutral'],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=conservation_col['Weak'],add=T,border='#000000')
plot(hh2,col=conservation_col['Strong'],add=T,border='#000000')
#par(xpd=T);legend('topright',fill=conservation_col,legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(wGerp & highExpressed)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=breaks,plot=F)
hh$counts=hh$count/1000
w=which(wGerp & highExpressed & SpliceSites$Gerp_site>2)
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(wGerp & highExpressed & SpliceSites$Gerp_site>4)
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=conservation_col['Neutral'],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=conservation_col['Weak'],add=T,border='#000000')
plot(hh2,col=conservation_col['Strong'],add=T,border='#000000')
#par(xpd=T);legend('topright',fill=conservation_col,legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)
par(xpd=T);legend('topright',fill=conservation_col,legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)
layout(1)
##########################################
###   END  Fig S1D                     ###
##########################################

####################################################################################
# Fig 1C: 4x 5.5 inches - age of weak constitutive and alternative Splice sites  ###
####################################################################################
w=which(wGerp)
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
par(mar=c(4,5,1,0.5))
tab=table(x,SpliceSites$WeakAltConst[w])
tab=tab[,c('weak','alternative','constitutive')]
colnames(tab)=c('weak','alt.','const.')
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(PhyloAge_col),ylab='% of splice sites',las=1)
#tab=table(x,SpliceSites$isImmune[w])
#barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

##########################################
###   END  Fig 1C                     ####
##########################################

###########################################################################################
## Fig 1D :  12x5.5 inches - age of non-immune,immune and stimulation induced splices sites ## 
###########################################################################################

wActive=which(wGerp & SpliceSites$WeakAltCons!='weak' )
wConst=which(wGerp & SpliceSites$WeakAltCons=='constitutive')
wAlt=which(wGerp & SpliceSites$WeakAltCons=='alternative')

#w=wConst
#NbNonConserved=By(SpliceSites$Gerp_site[w]<2,paste(SpliceSites$gene[w],SpliceSites$isImmune[w]),mean)
#NbConserved=By(SpliceSites$Gerp_site[w]>2,paste(SpliceSites$gene[w],SpliceSites$isImmune[w]),mean)

layout(matrix(c(1,1,1,2,2,2,3,3,3),nrow=1))
w=wActive
par(mar=c(4,5,1,0.5))
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC_overlappingGene_withStrand[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1])
tab3=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1 & SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(PhyloAge_col),ylab='% of splice sites',las=1)

par(mar=c(4,5,1,0.5))
w=wConst
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC_overlappingGene_withStrand[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1])
tab3=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1 & SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(PhyloAge_col),ylab='% of splice sites',las=1)

par(mar=c(4,5,1,0.5))
w=wAlt
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC_overlappingGene_withStrand[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1])
tab3=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1 & SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(PhyloAge_col),ylab='% of splice sites',las=1)

wActive=which(wGerp & SpliceSites$WeakAltCons!='weak')
wConst=which(wGerp & SpliceSites$WeakAltCons=='constitutive')
wAlt=which(wGerp & SpliceSites$WeakAltCons=='alternative')

layout(matrix(c(1,1,1,2,2,2,3,3,3),nrow=1))
w=wActive
par(mar=c(4,5,1,0.5))
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1])
tab3=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1 & SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

par(mar=c(4,5,1,0.5))
w=wConst
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC_overlappingGene_withStrand[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1])
tab3=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1 & SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

par(mar=c(4,5,1,0.5))
w=wAlt
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC_overlappingGene_withStrand[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1])
tab3=table(x[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1 & SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

w=wConst
GeneHasNonConservedConstSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$overlappingGene_withStrand[w],any)
NbConstSiteByGene=By(SpliceSites$Gerp_site[w]<2,SpliceSites$overlappingGene_withStrand[w],length)

mean(GeneHasNonConservedConstSite) # 20% of genes have a constitutive splice site that is non conserved
resGO_NcConst=GOSeq(names(GeneHasNonConservedConstSite)[GeneHasNonConservedConstSite], names(GeneHasNonConservedConstSite),addCI=T,FDR=0.05,bias.data=NbConstSiteByGene)
resGO_NcConst$type='constitutive'
resGO_NcConst=resGO_NcConst[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wAlt
GeneHasNonConservedAltSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$overlappingGene_withStrand[w],any)
mean(GeneHasNonConservedAltSite) # 59% of genes have a alternative splice site that is non conserved
NbAltSiteByGene=By(SpliceSites$Gerp_site[w]<2,SpliceSites$overlappingGene_withStrand[w],length)
resGO_NcAlt=GOSeq(names(GeneHasNonConservedAltSite)[GeneHasNonConservedAltSite], names(GeneHasNonConservedAltSite),addCI=T,FDR=0.05,bias.data=NbAltSiteByGene)
resGO_NcAlt$type='alternative'
resGO_NcAlt=resGO_NcAlt[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wActive
GeneHasNonConservedActSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$overlappingGene_withStrand[w],any)
NbActSiteByGene=By(SpliceSites$Gerp_site[w]<2,SpliceSites$overlappingGene_withStrand[w],length)

mean(GeneHasNonConservedActSite) # 54% of genes with an active splice site have one that is non conserved
resGO_NcAct=GOSeq(names(GeneHasNonConservedActSite)[GeneHasNonConservedActSite], names(GeneHasNonConservedActSite),addCI=T,FDR=0.05,bias.data=NbActSiteByGene)
resGO_NcAct$type='active'
resGO_NcAct=resGO_NcAct[,c(1:5,7:8,10:11,13:17,6,18,12)]

write.table(rbind(resGO_NcConst,resGO_NcAlt),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/TableS1A_GO_nonConserved.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(resGO_NcConst,resGO_NcAlt,resGO_NcAct),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/TableS1A_GO_nonConserved_withActive.txt',HOME),quote=F,sep='\t',row.names=FALSE)

odds.ratio.adj=function(y,x,z){
    z=as.data.frame(z)
    mod=summary(glm(y~x+.,data=z,bamily='binomial'))
    b=mod$coef[2,1]
    se=mod$coef[2,2]
    p=mod$coef[2,4]
    c(LowerCI=exp(b-1.96*se),OR=exp(b),UpperCI=exp(b+1.96*se),alpha=0.05,P=p)
}


w=wConst
GeneExpr=By(SpliceSites$maxFPKM_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
GeneExpr=cut(GeneExpr,c(1,5,10,50,100,500,1000,100000))

GenelFC=By(SpliceSites$maxlFC_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_all=odds.ratio.adj(GeneHasNonConservedConstSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_LPS_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_LPS=odds.ratio.adj(GeneHasNonConservedConstSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_PAM3CSK4_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_PAM3=odds.ratio.adj(GeneHasNonConservedConstSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_R848_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_R848=odds.ratio.adj(GeneHasNonConservedConstSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_IAV_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_IAV=odds.ratio.adj(GeneHasNonConservedConstSite,GenelFC>1,GeneExpr)
OR_LFC_Const=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='constitutive')

w=wAlt
GeneExpr=By(SpliceSites$maxFPKM_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
GeneExpr=cut(GeneExpr,c(1,5,10,50,100,500,1000,100000))

GenelFC=By(SpliceSites$maxlFC_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_all=odds.ratio.adj(GeneHasNonConservedAltSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_LPS_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_LPS=odds.ratio.adj(GeneHasNonConservedAltSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_PAM3CSK4_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_PAM3=odds.ratio.adj(GeneHasNonConservedAltSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_R848_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_R848=odds.ratio.adj(GeneHasNonConservedAltSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_IAV_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_IAV=odds.ratio.adj(GeneHasNonConservedAltSite,GenelFC>1,GeneExpr)

OR_LFC_Alt=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='alternative')

w=wActive
GeneExpr=By(SpliceSites$maxFPKM_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
GeneExpr=cut(GeneExpr,c(1,5,10,50,100,500,1000,100000))

GenelFC=By(SpliceSites$maxlFC_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_all=odds.ratio.adj(GeneHasNonConservedActSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_LPS_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_LPS=odds.ratio.adj(GeneHasNonConservedActSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_PAM3CSK4_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_PAM3=odds.ratio.adj(GeneHasNonConservedActSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_R848_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_R848=odds.ratio.adj(GeneHasNonConservedActSite,GenelFC>1,GeneExpr)
GenelFC=By(SpliceSites$log2FC_IAV_overlappingGene_withStrand[w],SpliceSites$overlappingGene_withStrand[w],max)
OR_IAV=odds.ratio.adj(GeneHasNonConservedActSite,GenelFC>1,GeneExpr)
OR_LFC_Act=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='active')

write.table(rbind(OR_LFC_Const,OR_LFC_Alt),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/TableS1B_OR_LFC.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(OR_LFC_Const,OR_LFC_Alt,OR_LFC_Act),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/TableS1B_OR_LFC_withActive.txt',HOME),quote=F,sep='\t',row.names=FALSE)

mean(SpliceSites[wConst,'Time_FirstAppeared_MA_NewAge']>=105)
# 0.9836494

# resample matching on gene age and expression
w=which(wGerp & SpliceSites$WeakAltCons!='weak')
GeneAge_overlappingGene_withStrand=By(SpliceSites$Time_FirstAppeared_MA_NewAge[w],SpliceSites$overlappingGene_withStrand[w],max,na.rm=T) # 
SpliceSites$GeneAge_overlappingGene_withStrand=GeneAge_overlappingGene_withStrand[match(SpliceSites$overlappingGene_withStrand,names(GeneAge_overlappingGene_withStrand))] #

FPKM_decile=quantile(SpliceSites$maxFPKM_overlappingGene_withStrand,seq(0,1,l=11),na.rm=T)
SpliceSites$FPKM_decile=cut(SpliceSites$maxFPKM_overlappingGene_withStrand,FPKM_decile)
SpliceSites$FPKM_group=cut(SpliceSites$maxFPKM_overlappingGene_withStrand,c(1,5,10,50,100,500,1000,100000))

w=which(wGerp & SpliceSites$WeakAltCons=='constitutive')
Age_all_const=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w]) #[1] 452.9037
Age_immune_const=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes']]) #[1] 417.5111
Age_up_const=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1]],na.rm=T) #[1] 397.9036
Age_all_const;Age_immune_const;Age_up_const;

100-Age_immune_const/Age_all_const*100
# 7.814593
100-Age_up_const/Age_all_const*100
# 12.14388


nsamp=1000
resamp_Age_maxFC_group=sapply(1:nsamp,function(i){toSample=By(SpliceSites$maxlFC_overlappingGene_withStrand[w]>1,paste(SpliceSites$GeneAge[w],SpliceSites$FPKM_group[w]),sum,na.rm=T)
	samp=unlist(lapply(names(toSample),function(i){sample(which(paste(SpliceSites$GeneAge[w],SpliceSites$FPKM_group[w])==i ),toSample[i])}))
	mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[samp]],na.rm=T)})
mean(resamp_Age_maxFC_group<=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1]],na.rm=T)) # P<10-4
# 0 -> p< 1e-3
2*pnorm(mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1]],na.rm=T),mean(resamp_Age_maxFC_group),sd(resamp_Age_maxFC_group),low=T)
# 4.049015e-81


resamp_Age_immune_group=sapply(1:nsamp,function(i){toSample=By(SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes',paste(SpliceSites$GeneAge_overlappingGene_withStrand[w],SpliceSites$FPKM_group[w]),sum,na.rm=T)
	samp=unlist(lapply(names(toSample),function(i){sample(which(paste(SpliceSites$GeneAge_overlappingGene_withStrand[w],SpliceSites$FPKM_group[w])==i),toSample[i])}))
	mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[samp]],na.rm=T)})
mean(resamp_Age_immune_group<=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes']],na.rm=T)) # P<10-4
#0.001
2*pnorm(mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes']],na.rm=T),mean(resamp_Age_immune_group),sd(resamp_Age_immune_group),low=T)
# 1.483473e-06

save(resamp_Age_immune_group,resamp_Age_maxFC_group,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/resamp_Age_spliceSites_FPKM_group_Constitutive_1000.Rdata',HOME))

w=which(wGerp & SpliceSites$WeakAltCons=='alternative')
Age_all_alt=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w]) # 412.808
Age_immune_alt=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes']]) # 377.2637
Age_up_alt=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1]],na.rm=T) # 333.7809
Age_all_alt;Age_immune_alt;Age_up_alt;
#[1] 292.5174
#[1] 247.7629
#[1] 221.5017
100-Age_immune_alt/Age_all_alt*100
#[1] 15.29976
100-Age_up_alt/Age_all_alt*100
#[1] 24.27741

nsamp=1000
resamp_Age_maxFC_group=sapply(1:nsamp,function(i){toSample=By(SpliceSites$maxlFC_overlappingGene_withStrand[w]>1,paste(SpliceSites$GeneAge[w],SpliceSites$FPKM_group[w]),sum,na.rm=T)
	samp=unlist(lapply(names(toSample),function(i){sample(which(paste(SpliceSites$GeneAge[w],SpliceSites$FPKM_group[w])==i ),toSample[i])}))
	mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[samp]],na.rm=T)})
mean(resamp_Age_maxFC_group<=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1]],na.rm=T)) # P<10-4
# 0 -> p< 1e-3
2*pnorm(mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1]],na.rm=T),mean(resamp_Age_maxFC_group),sd(resamp_Age_maxFC_group),low=T)
#[1] 1.49252e-44

resamp_Age_immune_group=sapply(1:nsamp,function(i){toSample=By(SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes',paste(SpliceSites$GeneAge_overlappingGene_withStrand[w],SpliceSites$FPKM_group[w]),sum,na.rm=T)
	samp=unlist(lapply(names(toSample),function(i){sample(which(paste(SpliceSites$GeneAge_overlappingGene_withStrand[w],SpliceSites$FPKM_group[w])==i),toSample[i])}))
	mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[samp]],na.rm=T)})
mean(resamp_Age_immune_group<=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes']],na.rm=T)) # P<10-4
#0.000
2*pnorm(mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes']],na.rm=T),mean(resamp_Age_immune_group),sd(resamp_Age_immune_group),low=T)
#[1] 0.0002147681

save(resamp_Age_immune_group,resamp_Age_maxFC_group,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/resamp_Age_spliceSites_FPKM_group_Alternative_1000.Rdata',HOME))


# Figure 1D
nrow(SpliceSites[wGerp & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 ,])
#[1] 193161
nrow(SpliceSites[wGerp & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & SpliceSites$NbMatchingSeqPrim==1,])
# 121
nrow(SpliceSites[wGerp & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & SpliceSites$coding ,])
#[1] 174138
nrow(SpliceSites[wGerp & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & SpliceSites$coding & SpliceSites$NbMatchingSeqPrim==1,])
# 28

reverseSeq=function(seq){paste(sapply(strsplit(seq,':'),function(x){nn=nchar(x);paste(substr(x,1,nn-3),reverse(substr(x,nn-1,nn)),sep='_')}),collapse=':')}

#plot(tree,type="fan")
strand='+'
reverse=function(Seq){
	Corresp=c('A','G','T','C','-')
	names(Corresp)=c('T','C','A','G','-')
	sapply(strsplit(Seq,''),function(x){paste(rev(Corresp[x]),collapse='')})
	}
	
HumanSpecificSpliceSites=SpliceSites[wGerp & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & SpliceSites$NbMatchingSeqPrim==1,]
# Chromosome	Splice site start	Splice site end	strand	type	Sequence of splice site	Gene associated to the splice site (Ensembl ID)	Gene associated to the splice site (symbol)	Present in over 95% of observed transcripts	Total read count NS	Gene FPKM NS	Total read count LPS	Gene FPKM LPS	Total read count Pam3Csk4	Gene FPKM Pam3Csk4	Total read count R848	Gene FPKM R848	Total read count IAV	Gene FPKM IAV	Mean GerpRS of splice site 	Fist apparition of the current splice site sequence in the phylogeny	Nb of primates with an alignment	Splice site sequence across common ancestors	Splice site sequence across 46 vertebrates
HumanSpecificSpliceSites$seq=ifelse(HumanSpecificSpliceSites$type=='donor',HumanSpecificSpliceSites$donor_site,HumanSpecificSpliceSites$acceptor_site)
HumanSpecificSpliceSites$allSequence[HumanSpecificSpliceSites$strand_correct=='-']=sapply(HumanSpecificSpliceSites$allSequence[HumanSpecificSpliceSites$strand_correct=='-'],reverseSeq)
HumanSpecificSpliceSites$ancestralSequence[HumanSpecificSpliceSites$strand_correct=='-']=sapply(HumanSpecificSpliceSites$ancestralSequence[HumanSpecificSpliceSites$strand_correct=='-'],reverseSeq)

HumanSpecificSpliceSites$ancestralSequence=gsub('treeShrew','treeshrew',HumanSpecificSpliceSites$ancestralSequence)
HumanSpecificSpliceSites$ancestralSequence=sapply(strsplit(HumanSpecificSpliceSites$ancestralSequence,':'),function(x){paste(toupper(substr(unlist(x),1,1)),substr(x,2,10000),sep='',collapse=':')})
HumanSpecificSpliceSites$ancestralSequence=gsub(':',' / ',HumanSpecificSpliceSites$ancestralSequence)
HumanSpecificSpliceSites$ancestralSequence=gsub('_',': ',HumanSpecificSpliceSites$ancestralSequence)

HumanSpecificSpliceSites$allSequence=gsub('X_tropicalis','X. tropicalis',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Rock_hyrax','Rock hyrax',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Kangaroo_rat','Kangaroo rat',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('TreeShrew','Treeshrew',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Zebra_finch','Zebra finch',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('lemur','Lemur',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Guinea_Pig','Guinea pig',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub(':',' / ',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('_',': ',HumanSpecificSpliceSites$allSequence)

HumanSpecificSpliceSites=HumanSpecificSpliceSites[,c("chrom_intron","start_site","end_site","strand_correct","seq","type","overlappingGene_withStrand","overlappingSymbol_withStrand","isImmune_overlappingGene_withStrand",'read_persample_NS',"FPKM_NS_overlappingGene_withStrand",'read_persample_LPS',"FPKM_LPS_overlappingGene_withStrand",'read_persample_PAM3CSK4',"FPKM_PAM3CSK4_overlappingGene_withStrand",'read_persample_R848',"FPKM_R848_overlappingGene_withStrand",'read_persample_IAV',"FPKM_IAV_overlappingGene_withStrand","Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site","WeakAltConst",'Gerp_site','NbAlignedPrimates',"allSequence","ancestralSequence","coding","inEnsembl")]
write.table(HumanSpecificSpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/TableS1C_HumanSpecificSpliceSites.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(HumanSpecificSpliceSites[HumanSpecificSpliceSites$coding,],file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/tables/TableS1C_HumanSpecificSpliceSites_coding.txt',HOME),quote=F,sep='\t',row.names=FALSE)

# Fig 1B left: 4x 5.5 inches - Conservation of weak constitutive and alternative Splice sites
par(mar=c(4,5,1,0.5))
w=which(wGerp)
x=cut(SpliceSites$Gerp_site[w],c(-Inf,2,4,Inf))
levels(x)=c('<2','[2;4]','>4')
par(mar=c(4,5,1,0.5))
tab=table(x,SpliceSites$WeakAltConst[w])
tab=tab[,c('weak','alternative','constitutive')]
colnames(tab)=c('weak','alt.','const.')
t(t(tab)/apply(tab,2,sum))*100
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),ylab='% of splice sites',las=1)
par(xpd=T)
legend(-2,-15,legend=levels(x),fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),bty='n',ncol=4)
par(xpd=F)




#=which(wGerp & SpliceSites$WeakAltCons!='weak')
#resamp_Age_maxFC_decile=sapply(1:nsamp,function(i){toSample=By(SpliceSites$maxlFC_overlappingGene_withStrand[w]>1,paste(SpliceSites$GeneAge[w],SpliceSites$FPKM_decile[w]),sum,na.rm=T)
#	samp=unlist(lapply(names(toSample),function(i){sample(which(paste(SpliceSites$GeneAge[w],SpliceSites$FPKM_decile[w])==i ),toSample[i])}))
#	mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[samp]],na.rm=T)})
#mean(resamp_Age_maxFC_decile<=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$maxlFC_overlappingGene_withStrand[w]>1]],na.rm=T)) # P<10-4
#0
#resamp_Age_immune_decile=sapply(1:nsamp,function(i){toSample=By(SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes',paste(SpliceSites$GeneAge[w],SpliceSites$FPKM_decile[w]),sum,na.rm=T)
#	samp=unlist(lapply(names(toSample),function(i){sample(which(paste(SpliceSites$GeneAge[w],SpliceSites$FPKM_decile[w])==i),toSample[i])}))
#	mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[samp]],na.rm=T)})
#mean(resamp_Age_immune_decile<=mean(SpliceSites$Time_FirstAppeared_MA_NewAge[w[SpliceSites$isImmune_overlappingGene_withStrand[w]=='yes']],na.rm=T)) # P<10-4
#0.001
#save(resamp_Age_immune_decile,resamp_Age_maxFC_decile,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V7/data/resamp_Age_spliceSites_FPKM_decile.Rdata',HOME))

#par(mar=c(4,5,5,0.5))
#hh=hist(c(-10,sample(Gerp_mean,length(wGerp),replace=T),Inf),br=50,plot=F)
#hh$counts=hh$count/1000
#plot(hh,xlim=c(-7,7),col=rcol[7],xlab='GerpRS',ylab='count',las=1,main='')
#hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[wGerp])),br=hh$br,plot=F);
#hh1$counts=hh1$count/1000;plot(hh1,col=gsub('80',"BB",rcol[3]),add=T,border='#00000022')
#hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='constitutive'])),br=hh$br,plot=F);
#hh1$counts=hh1$count/1000;plot(hh1,col=gsub('80',"BB",rcol[1]),add=T,border='#00000022')
#par(xpd=T);legend(-10,800,fill=c(rcol[7],gsub('80',"BB",rcol[c(3,1)])),legend=c('Genome Wide','splice sites','constitutive splice sites'),bty='n');par(xpd=F)

#wIntro=which(Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread_Intropolis))
#
#
#SpliceSitesOfknownjunction=unique(c(Junc_annot$start_id[wIntro],Junc_annot$end_id[wIntro]))
#luq(SpliceSitesOfknownjunction) # 1982792
#
#wIntroGenic=which(Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread_Intropolis) & !is.na(Junc_annot$overlappingSymbol))
#SpliceSitesOfknownGenicjunction=unique(c(Junc_annot$start_id[wIntroGenic],Junc_annot$end_id[wIntroGenic]))
#luq(SpliceSitesOfknownGenicjunction)
#
#wIntroNonGenic=which(Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread_Intropolis) & is.na(Junc_annot$overlappingSymbol)) 
#SpliceSitesOfknownNonGenicjunction=unique(c(Junc_annot$start_id[wIntroGenic],Junc_annot$end_id[wIntroGenic]))
#luq(SpliceSitesOfknownNonGenicjunction)
#
#knownGenicSpliceSites=SpliceSites[!is.na(SpliceSites$Nbread_Intropolis) & !is.na(SpliceSites$overlappingGene),get('site_id')]
#luq(knownGenicSpliceSites)
#knownNonGenicSpliceSites=SpliceSites[!is.na(SpliceSites$Nbread_Intropolis) & is.na(SpliceSites$overlappingGene),get('site_id')]
#luq(knownNonGenicSpliceSites)
#
#
#ToStudy=intersect(gsub(' acceptor| donor','',knownNonGenicSpliceSites),SpliceSitesOfknownGenicjunction)
#Junc_annot[Junc_annot$start_id%in%ToStudy[1:10] & Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread_Intropolis) & !is.na(Junc_annot$overlappingSymbol),]
#SpliceSites[gsub(' acceptor| donor','',SpliceSites$site_id)%in%ToStudy[1:10] & SpliceSites$strand_Intropolis!='' & !is.na(SpliceSites$Nbread_Intropolis),]


By(SpliceSites$Gerp_site[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align]>2,paste(ifelse(SpliceSites$inEnsembl,'annot','missing'),SpliceSites$WeakAltConst)[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align],mean)
#  annot alternative   annot constitutive     annot weak  missing alternative missing constitutive         missing weak 
#   0.7550677            0.9759993            0.6638520            0.1880780            0.1545101 			0.1972147


mean((SpliceSites$Gerp_site<2 & !SpliceSites$inEnsembl)[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='weak'])
# 0.702796
# While 70% of weak splice sites were cryptic splice sites that are absent from Ensembl annotation and present no sign of conservation, supporting the high prevalence of “noisy” splicing events previously reported in human cells {Melamud, 2009 #18;Pickrell, 2010 #16}, 56% of alternative and 97% of constitutive splice sites showed moderate to high conservation (GerpRS>2) 

mean((SpliceSites$Gerp_site>2)[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='weak'])
# 0.2552072

mean((!SpliceSites$inEnsembl)[SpliceSites$Gerp_site<2 & !SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='weak'])


mean((SpliceSites$Gerp_site<2 & !SpliceSites$inEnsembl)[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='weak'])


By(SpliceSites$newAge[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align]>12,SpliceSites$WeakAltConst[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align],mean,na.rm=T)
#alternative constitutive         weak 
#   0.4940938    0.9229044    0.2385476 
   
table(SpliceSites$WeakAltConst[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site < 2])/sum(!SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site < 2  & SpliceSites$Quality_MultiZ46Align)
#alternative constitutive         weak 
# 0.034131939  0.005453834  0.960414227 
# the vast majority of non conserved sites are weak

table(SpliceSites$WeakAltConst[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site > 2])/sum(!SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site > 2  & SpliceSites$Quality_MultiZ46Align)
# alternative constitutive         weak 
#  0.08406832   0.28880743   0.62712426 
# but weak sites are also the most frequent among conserved site

# Check sequence at recent but conserved sites
SpliceSites[which(SpliceSites$newAge<5 & SpliceSites$Gerp_site>4)[1:2],]
# Check sequence at old but non conserved sites
	
By(SpliceSites$newAge[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]<=12,SpliceSites$WeakAltConst[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site],mean,na.rm=T)
# alternative constitutive         weak 
#  0.50590621   0.07709557   0.76145239 
By(SpliceSites$newAge[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]==18,SpliceSites$WeakAltConst[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site],mean,na.rm=T)
# alternative constitutive         weak 
#  0.21282849   0.55860229   0.09426126 

mean(SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]!='weak',na.rm=T) # 0.0500166
mean(SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]=='constitutive',na.rm=T) # 0.0117122
sum(SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]!='weak',na.rm=T) # 57993
sum(SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]=='constitutive',na.rm=T) # 13580

mean(!SpliceSites$coding[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site & SpliceSites$WeakAltConst!='weak'],na.rm=T) # 
mean(!SpliceSites$coding[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site & SpliceSites$WeakAltConst=='constitutive'],na.rm=T) # 0.5478465
table(SpliceSites$coding[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site],SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site])
#        alternative constitutive    weak
#  FALSE       37189         6020 1080963
#  TRUE         7224         7560   20519


mean(SpliceSites$coding[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='constitutive'],na.rm=T) # 0.5478465


SpliceSites$[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='constitutive' & SpliceSites$coding]


#hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$inEnsembl])),br=hh$br,plot=F);
#hh1$counts=hh1$count/1000;plot(hh1,col=rcol[1],add=T,border='#00000022')
#hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Nbread>10000])),br=hh$br,plot=F);
#hh1$counts=hh1$count/1000;plot(hh1,density=30,add=T,border='#000000')




















#
## create supplementary figures S1A
#wAnnot=which(Junc_annot$Annotated_junc)
#wCanonical=which(Junc_annot$DA_type!='other')
#wIntro=which(Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread_Intropolis) & !is.na(Junc_annot$overlappingSymbol))
#
#PctReads=c(All=1,
#Canonical=sum(as.numeric(Junc_annot$Total_count[wCanonical]))/sum(as.numeric(Junc_annot$Total_count)), #0.9991252
#Intropolis=sum(as.numeric(Junc_annot$Total_count[wIntro]))/sum(as.numeric(Junc_annot$Total_count)), #0.9984076
#Ensembl=sum(as.numeric(Junc_annot$Total_count[wAnnot]))/sum(as.numeric(Junc_annot$Total_count))) #0.9772441
#
#Nbjunc=c(All=nrow(Junc_annot),Canonical=length(wCanonical),Intropolis=length(wIntro),Ensembl=length(wAnnot))
#
#tab_total=table(c(EndSite$sequence,StartSite$sequence)[c(EndSite$sequence,StartSite$sequence)!='NN'])
#Pct_nonCano_total=sum(tab_total[!names(tab_total)%in%c('AG','GT')])/sum(tab_total)
#
#tab_intro=table(c(EndSite$sequence,StartSite$sequence)[!is.na(c(EndSite$Nbread,StartSite$Nbread)) & c(EndSite$sequence,StartSite$sequence)!='NN' & c(!is.na(EndSite$overlappingGene),!is.na(StartSite$overlappingGene))])
#Pct_nonCano_intro=sum(tab_intro[!names(tab_intro)%in%c('AG','GT')])/sum(tab_intro)
#
#tab_annot=table(c(EndSite$sequence,StartSite$sequence)[c(EndSite$inEnsembl, StartSite$inEnsembl) & c(EndSite$sequence,StartSite$sequence)!='NN'])
#Pct_nonCano_annot=sum(tab_annot[!names(tab_annot)%in%c('AG','GT')])/sum(tab_annot)
#
## Fig S1Abis_SpliceSite Sequence 4x5.5 inches
#
###########################################
#### 	 Fig S1A impact of filters      ## 5x 4 inches
###########################################
#par(mar=c(6,4,1,4))
#x=barplot(Nbjunc/1e6,las=2,ylab='discovered junctions (millions)')
#y=(1-PctReads)*100
#Factor=max(Nbjunc)/1e6/3
#points(x,y*Factor,col='red',pch=16,type='b')
#axis(4,at=seq(0,max(round(y)+0.5)*Factor,l=6),labels=seq(0,max(round(y))+0.5,l=6),las=2)
#mtext('% of spliced reads removed',4,3)
#
###########################################
####   END  Fig S1A impact of filters    ###
###########################################
#layout(matrix(rep(1:2,c(4,3)),nrow=1))
#condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
## Fig S1Bbis_SpliceSite Sequence 3.6x7.6 inches
##layout(matrix(rep(1:2,c(4,3)),nrow=1))
#par(mar=c(6,6,1,0.1))
#tab_total=rev(sort(tab_total))
#barplot(tab_total/1000,col='lightgrey',las=2,cex.axis=1.3)
#barplot(tab_intro[names(tab_total)]/1000,add=T,col=acol[3],axes=F,las=2,axisname=F)
#barplot(tab_annot[names(tab_total)]/1000,add=T,col=acol[8],axes=F,las=2,axisname=F)
#legend('topright',legend=c('All','Intropolis','Ensembl'),fill=c('lightgrey',acol[c(3,8)]),bty='n',cex=1.2)
#
#par(mar=c(6,6,1,1))
#a=0.3;b=1.8
#names(Pct_nonCanonical_SpliceSites_Intropolis)=c('All',condIndex)
#y=c(Pct_nonCano_total,Pct_nonCano_intro,Pct_nonCano_annot,Pct_nonCanonical_SpliceSites_Intropolis[-1])*100
#names(y)=c('All','Intropolis','Ensembl',condIndex)
#barplot(y,col=c('lightgrey',acol[c(3,8)],colERC5),las=2,ylab='Percentage of non-canonical\n splice sites',cex.lab=1.2,space=c(a,a,a,b,a,a,a,a),cex.axis=1.3)
#
#
#tab_total=rev(sort(tab_total))
#barplot(tab_total/1000,col='lightgrey',las=2)
#barplot(tab_over1[names(tab_total)]/1000,add=T,col=acol[6],axes=F,las=2)
#barplot(tab_intro[names(tab_total)]/1000,add=T,col=acol[3],axes=F,las=2)
#barplot(tab_annot[names(tab_total)]/1000,add=T,col=acol[8],axes=F,las=2)
#legend('topright',legend=c('All',"Over 1 read",'Intropolis','Ensembl'),fill=c('lightgrey',acol[c(6,3,8)]),bty='n')
#barplot(c(Pct_nonCano_total,Pct_nonCano_over1,Pct_nonCano_intro,Pct_nonCano_annot)*100,ylab="Percentage of non canonical splice sites",col=c('lightgrey',acol[c(6,3,8)]),las=2)
#
#
## adding Genic regions (test)
#layout(matrix(c(1,1,2),nrow=1))
#tab_total=rev(sort(tab_total))
#barplot(tab_total/1000,col='lightgrey',las=2)
#barplot(tab_Genic[names(tab_total)]/1000,add=T,col=acol[8],axes=F,las=2,axisnames = FALSE)
#barplot(tab_intro[names(tab_total)]/1000,add=T,col=acol[6],axes=F,las=2,axisnames = FALSE)
#barplot(tab_annot[names(tab_total)]/1000,add=T,col=acol[3],axes=F,las=2,axisnames = FALSE)
#legend('topright',legend=c('All','Genic','Intropolis','Ensembl'),fill=c('lightgrey',acol[c(8,6,3)]),bty='n')
#barplot(c(Pct_nonCano_total,Pct_nonCano_Genic,Pct_nonCano_intro,Pct_nonCano_annot)*100,ylab="Percentage of non canonical splice sites",col=c('lightgrey',acol[c(8,6,3,1)]),las=2)
#
#
#Pct_nonCanonical_SpliceSites_Intropolis=c()
#for (j in 1:6){
#    mycol=c('Total_count', 'Total_count_NS','Total_count_LPS','Total_count_PAM3CSK4','Total_count_R848','Total_count_IAV')[j]
#    tab1=table(c(EndSite$sequence,StartSite$sequence)[c(EndSite$sequence,StartSite$sequence)!='NN' & !is.na(c(EndSite$Nbread,StartSite$Nbread)) & c(EndSite[[mycol]], StartSite[[mycol]])>0])
#    Pct_nonCanonical_SpliceSites_Intropolis[j]=sum(tab1[!names(tab1)%in%c('AG','GT')])/sum(tab1)
#    }

########################################
# 	 Fig S1C impact of filters      ###
########################################
#
#par(mar=c(7,4,1,4))
#x=barplot(Nbjunc/1e6,las=2,ylab='discovered junctions (millions)')
#y=(1-PctReads)*100
#Factor=max(Nbjunc)/1e6/3
#points(x,y*Factor,col='red',pch=16,type='b')
#axis(4,at=seq(0,max(round(y)+0.5)*Factor,l=6),labels=seq(0,max(round(y))+0.5,l=6),las=2)
#mtext('% of spliced reads removed',4,3)

########################################
#   END  Fig S1C impact of filters    ###
########################################


#Junc_annot_highCov=Junc_annot[Junc_annot$highCovAnyCond_junc & Junc_annot$strand_Intropolis!='',c("chrom","start","end","strand_Intropolis","donor","acceptor","gene","symbol","coding","Nbread","Total_count_NS","Total_count_LPS","Total_count_PAM3CSK4","Total_count_R848","Total_count_IAV")]
#colnames(Junc_annot_highCov)=c("chrom","start","end","strand","donor","acceptor","gene","symbol","is coding","Nb of Intropolis reads","Total read count - NS","Total read count - LPS","Total read count - PAM3CSK4","Total read count - R848","Total read count - IAV")
#Junc_annot_highCov$gene[is.na(Junc_annot_highCov$gene)]=''
#Junc_annot_highCov$symbol[is.na(Junc_annot_highCov$symbol)]=''
#Junc_annot_highCov$coding=ifelse(Junc_annot_highCov$coding,'yes','')
#write.table(Junc_annot_highCov,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_Junctions_HighCov.txt',HOME),sep='\t',quote=F,row.names=F)


#
#
#
#layout(1:3)
#barplot(table(as.factor(EndConserv$V13)),add=F,col=colERC5[1])
#barplot(table(EndConserv$V14),add=T,col=colERC5[2])
#barplot(table(EndConserv$V15),add=T,col=colERC5[3])
#barplot(table(EndConserv$V16),add=T,col=colERC5[4])
#
#barplot(table(as.factor(EndConserv$V13)[EndConserv$V4<2]),add=F,col=colERC5[1])
#barplot(table(EndConserv$V14[EndConserv$V4<2]),add=T,col=colERC5[2])
#barplot(table(EndConserv$V15[EndConserv$V4<2]),add=T,col=colERC5[3])
#barplot(table(EndConserv$V16[EndConserv$V4<2]),add=T,col=colERC5[4])
#
#barplot(table(as.factor(EndConserv$V13)[EndConserv$V4>2]),add=F,col=colERC5[1])
#barplot(table(as.factor(EndConserv$V14)[EndConserv$V4>2]),add=T,col=colERC5[2])
#barplot(table(as.factor(EndConserv$V15)[EndConserv$V4>2]),add=T,col=colERC5[3])
#barplot(table(as.factor(EndConserv$V16)[EndConserv$V4>2]),add=T,col=colERC5[4])
#barplot(table(as.factor(EndConserv$V13)[EndConserv$V4>4]),add=F,col=colERC5[1])
#
## sites that are alignable up to the lamprey are more frequent among conserved sites
#barplot(table(as.factor(EndConserv$V13)),add=F,col=colERC5[2])
#barplot(table(as.factor(EndConserv$V13)[EndConserv$V4<2 & EndConserv$V4> -2]),add=T,col=colERC5[1])
#hist(EndConserv$V4[EndConserv$V13==46]))
## hist(EndConserv$V4[EndConserv$V13<30],br=100,col=colERC[1])
#


#EndConserv=fread(sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteEndConserv.txt',EVO_IMMUNO_POP))
#StartConserv=fread(sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteStartConserv.txt',EVO_IMMUNO_POP))

#Donor sites 
#+ strand
#DonorSite_pos=Junc_annot[Junc_annot$strand_Intropolis=='+' & !is.na(Junc_annot$Nbread),]
#DonorSite_pos$start=DonorSite_pos$start+1 # added 1 to match UCSC coordinates
#DonorSite_pos$site_id=paste(DonorSite_pos$chrom,DonorSite_pos$start,DonorSite_pos$strand_Intropolis)
#Tot_count_site=apply(DonorSite_pos[,grep('Total_count',colnames(DonorSite_pos)),with=F],2,function(x){By(x,DonorSite_pos$site_id,sum,na.rm=T)})
#DonorSite_pos=DonorSite_pos[order(DonorSite_pos$site_id,-DonorSite_pos$Total_count),]
#DonorSite_pos=DonorSite_pos[!duplicated(DonorSite_pos$site_id),]
#DonorSite_pos$Junction_count=DonorSite_pos$Total_count
#for(i in colnames(Tot_count_site)){
#	DonorSite_pos[[i]]=Tot_count_site[match(DonorSite_pos$site_id,rownames(Tot_count_site)),i]
#}
#DonorSite_pos$matchingSite_MajorIsoform=paste(DonorSite_pos$chrom,DonorSite_pos$end,DonorSite_pos$strand_Intropolis)
#DonorSite_pos$Position_matchingSite_MajorIsoform=DonorSite_pos$end
#DonorSite_pos$end=DonorSite_pos$start+1
#DonorSite_pos$Gerp_site=DonorSite_pos$GerpRS_meanStart
#DonorSite_pos$inEnsembl=DonorSite_pos$Annotated_start
#DonorSite_pos$coding=DonorSite_pos$coding_start
#DonorSite_pos$constitutive=DonorSite_pos$constitutive_start
#DonorSite_pos$gene=DonorSite_pos$Gene_start
#DonorSite_pos$symbol=DonorSite_pos$symbol_start
#mm_Age=match(paste(DonorSite_pos$chrom,DonorSite_pos$start,DonorSite_pos$strand_HISAT),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))
#DonorSite_pos$FirstAppeared=StartConserv$V9[mm_Age]
#DonorSite_pos$Time_FirstAppeared=StartConserv$V10[mm_Age]
#DonorSite_pos$NbAlignedPrimates=StartConserv$NbAlignedSeqPrim[mm_Age]
#DonorSite_pos$NbMatchingSeqPrim=StartConserv$NbMatchingSeqPrim[mm_Age]
#DonorSite_pos$allSequence=StartConserv$V12[mm_Age]
#DonorSite_pos$ancestralSequence=StartConserv$V11[mm_Age]
#DonorSite_pos$newAge=StartConserv$newAge[mm_Age]
#
#
#- strand
#DonorSite_neg=Junc_annot[Junc_annot$strand_Intropolis=='-' & !is.na(Junc_annot$Nbread),]
#DonorSite_neg$start=DonorSite_neg$start+1 # added 1 to match UCSC coordinates
#DonorSite_neg$site_id=paste(DonorSite_neg$chrom,DonorSite_neg$end,DonorSite_neg$strand_Intropolis)
#Tot_count_site=apply(DonorSite_neg[,grep('Total_count',colnames(DonorSite_neg)),with=F],2,function(x){By(x,DonorSite_neg$site_id,sum,na.rm=T)})
#DonorSite_neg=DonorSite_neg[order(DonorSite_neg$site_id,-DonorSite_neg$Total_count),]
#DonorSite_neg=DonorSite_neg[!duplicated(DonorSite_neg$site_id),]
#DonorSite_neg$Junction_count=DonorSite_neg$Total_count
#for(i in colnames(Tot_count_site)){
#	DonorSite_neg[[i]]=Tot_count_site[match(DonorSite_neg$site_id,rownames(Tot_count_site)),i]
#}
#DonorSite_neg$matchingSite_MajorIsoform=paste(DonorSite_neg$chrom,DonorSite_neg$start,DonorSite_neg$strand_Intropolis)
#DonorSite_neg$Position_matchingSite_MajorIsoform=DonorSite_neg$start
#DonorSite_neg$start=DonorSite_neg$end-1
#DonorSite_neg$Gerp_site=DonorSite_neg$GerpRS_meanEnd
#DonorSite_neg$inEnsembl=DonorSite_neg$Annotated_end
#DonorSite_neg$coding=DonorSite_neg$coding_end
#DonorSite_neg$constitutive=DonorSite_neg$constitutive_end
#DonorSite_neg$gene=DonorSite_neg$Gene_end
#DonorSite_neg$symbol=DonorSite_neg$symbol_end
#mm_Age=match(paste(DonorSite_neg$chrom,DonorSite_neg$end,DonorSite_neg$strand_HISAT),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))
#DonorSite_neg$FirstAppeared=EndConserv$V9[mm_Age]
#DonorSite_neg$Time_FirstAppeared=EndConserv$V10[mm_Age]
#DonorSite_neg$NbAlignedPrimates=EndConserv$NbAlignedSeqPrim[mm_Age]
#DonorSite_neg$NbMatchingSeqPrim=EndConserv$NbMatchingSeqPrim[mm_Age]
#DonorSite_neg$allSequence=EndConserv$V12[mm_Age]
#DonorSite_neg$ancestralSequence=EndConserv$V11[mm_Age]
#DonorSite_neg$newAge=EndConserv$newAge[mm_Age]
#
#Acceptor sites 
#+ strand
#AcceptorSite_pos=Junc_annot[Junc_annot$strand_Intropolis=='+'  & !is.na(Junc_annot$Nbread),]
#AcceptorSite_pos$start=AcceptorSite_pos$start+1 # added 1 to match UCSC coordinates
#AcceptorSite_pos$site_id=paste(AcceptorSite_pos$chrom,AcceptorSite_pos$end,AcceptorSite_pos$strand_Intropolis)
#Tot_count_site=apply(AcceptorSite_pos[,grep('Total_count',colnames(AcceptorSite_pos)),with=F],2,function(x){By(x,AcceptorSite_pos$site_id,sum,na.rm=T)})
#AcceptorSite_pos=AcceptorSite_pos[order(AcceptorSite_pos$site_id,-AcceptorSite_pos$Total_count),]
#AcceptorSite_pos=AcceptorSite_pos[!duplicated(AcceptorSite_pos$site_id),]
#AcceptorSite_pos$Junction_count=AcceptorSite_pos$Total_count
#for(i in colnames(Tot_count_site)){
#	AcceptorSite_pos[[i]]=Tot_count_site[match(AcceptorSite_pos$site_id,rownames(Tot_count_site)),i]
#}
#AcceptorSite_pos$matchingSite_MajorIsoform=paste(AcceptorSite_pos$chrom,AcceptorSite_pos$start,AcceptorSite_pos$strand_Intropolis)
#AcceptorSite_pos$Position_matchingSite_MajorIsoform=AcceptorSite_pos$start+1 # added 1 to match UCSC coordinates
#AcceptorSite_pos$start=AcceptorSite_pos$end-1
#AcceptorSite_pos$Gerp_site=AcceptorSite_pos$GerpRS_meanEnd
#AcceptorSite_pos$inEnsembl=AcceptorSite_pos$Annotated_end
#AcceptorSite_pos$coding=AcceptorSite_pos$coding_end
#AcceptorSite_pos$constitutive=AcceptorSite_pos$constitutive_end
#AcceptorSite_pos$gene=AcceptorSite_pos$Gene_end
#AcceptorSite_pos$symbol=AcceptorSite_pos$symbol_end
#mm_Age=match(paste(AcceptorSite_pos$chrom,AcceptorSite_pos$end,AcceptorSite_pos$strand_HISAT),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))
#AcceptorSite_pos$FirstAppeared=EndConserv$V9[mm_Age]
#AcceptorSite_pos$Time_FirstAppeared=EndConserv$V10[mm_Age]
#AcceptorSite_pos$NbAlignedPrimates=EndConserv$NbAlignedSeqPrim[mm_Age]
#AcceptorSite_pos$NbMatchingSeqPrim=EndConserv$NbMatchingSeqPrim[mm_Age]
#AcceptorSite_pos$allSequence=EndConserv$V12[mm_Age]
#AcceptorSite_pos$ancestralSequence=EndConserv$V11[mm_Age]
#AcceptorSite_pos$newAge=EndConserv$newAge[mm_Age]
#
#- strand
#AcceptorSite_neg=Junc_annot[Junc_annot$strand_Intropolis=='-'  & !is.na(Junc_annot$Nbread),]
#AcceptorSite_neg$start=AcceptorSite_neg$start+1 # added 1 to match UCSC coordinates
#AcceptorSite_neg$site_id=paste(AcceptorSite_neg$chrom,AcceptorSite_neg$start,AcceptorSite_neg$strand_Intropolis)
#Tot_count_site=apply(AcceptorSite_neg[,grep('Total_count',colnames(AcceptorSite_neg)),with=F],2,function(x){By(x,AcceptorSite_neg$site_id,sum,na.rm=T)})
#AcceptorSite_neg=AcceptorSite_neg[order(AcceptorSite_neg$site_id,-AcceptorSite_neg$Total_count),]
#AcceptorSite_neg=AcceptorSite_neg[!duplicated(AcceptorSite_neg$site_id),]
#AcceptorSite_neg$Junction_count=AcceptorSite_neg$Total_count
#for(i in colnames(Tot_count_site)){
#	AcceptorSite_neg[[i]]=Tot_count_site[match(AcceptorSite_neg$site_id,rownames(Tot_count_site)),i]
#}
#AcceptorSite_neg$matchingSite_MajorIsoform=paste(AcceptorSite_neg$chrom,AcceptorSite_neg$end,AcceptorSite_neg$strand_Intropolis)
#AcceptorSite_neg$Position_matchingSite_MajorIsoform=AcceptorSite_neg$end
#AcceptorSite_neg$end=AcceptorSite_neg$start+1
#AcceptorSite_neg$Gerp_site=AcceptorSite_neg$GerpRS_meanStart
#AcceptorSite_neg$inEnsembl=AcceptorSite_neg$Annotated_start
#AcceptorSite_neg$coding=AcceptorSite_neg$coding_start
#AcceptorSite_neg$constitutive=AcceptorSite_neg$constitutive_start
#AcceptorSite_neg$gene=AcceptorSite_neg$Gene_start
#AcceptorSite_neg$symbol=AcceptorSite_neg$symbol_start
#mm_Age=match(paste(AcceptorSite_neg$chrom,AcceptorSite_neg$start,AcceptorSite_neg$strand_HISAT),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))
#AcceptorSite_neg$FirstAppeared=StartConserv$V9[mm_Age]
#AcceptorSite_neg$Time_FirstAppeared=StartConserv$V10[mm_Age]
#AcceptorSite_neg$NbAlignedPrimates=StartConserv$NbAlignedSeqPrim[mm_Age]
#AcceptorSite_neg$NbMatchingSeqPrim=StartConserv$NbMatchingSeqPrim[mm_Age]
#AcceptorSite_neg$allSequence=StartConserv$V12[mm_Age]
#AcceptorSite_neg$ancestralSequence=StartConserv$V11[mm_Age]
#AcceptorSite_neg$newAge=StartConserv$newAge[mm_Age]
#
#
#combine
#DonorSite=rbind(DonorSite_pos,DonorSite_neg)
#DonorSite$type='donor'
#combine
#AcceptorSite=rbind(AcceptorSite_pos,AcceptorSite_neg)
#AcceptorSite$type='acceptor'
## combine all
#SpliceSites=rbind(DonorSite,AcceptorSite)

############ Annotation is Immune gene 
#Require('allGOterms')
#GoToTest=list(immuneResponse='GO:0006955',
#	defenseResponse='GO:0006952',
#	InnateImmuneResponse='GO:0045087',
#	AdaptiveImmuneResponse='GO:0002250',
#	AntiviralResponse='GO:0051607',
#	AntibacterialResponse='GO:0042742',
#	TF_DNAbinding='GO:0003700',
#	development='GO:0032502',
#	olfactory_receptors='GO:0004984')
#GoToTest_genes=lapply(GoToTest,function(x){allGOterms$gene[allGOterms$go==x]})
#GoToTest_genes$all=unique(allGOterms$gene)
#SpliceSites$isImmune=ifelse(SpliceSites$gene%in% GoToTest_genes$immuneResponse,'yes','')
#
#write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V2.txt',HOME),quote=F,sep='\t',row.names=F)
#
#SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V2.txt',HOME))
########### Annotation Summary statistics
#library(data.table)
#StartConserv=fread(sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteStartConserv_V3.txt',EVO_IMMUNO_POP))
#StartConserv$V2=StartConserv$V2+1
#EndConserv=fread(sprintf('%s/Maxime/Splicing/Conservation/All_SpliceSiteEndConserv_V3.txt',EVO_IMMUNO_POP))
#
#wDonor_pos=which(SpliceSites$type=='donor' & SpliceSites$strand_Intropolis=='+')
#SpliceSites$NbAlignedPrimates[wDonor_pos]=StartConserv$NbAlignedSeqPrim[match(paste(SpliceSites$chrom[wDonor_pos], SpliceSites$start[wDonor_pos],SpliceSites$strand_HISAT[wDonor_pos]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#SpliceSites$NbMatchingSeqPrim[wDonor_pos]=StartConserv$NbMatchingSeqPrim[match(paste(SpliceSites$chrom[wDonor_pos], SpliceSites$start[wDonor_pos],SpliceSites$strand_HISAT[wDonor_pos]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#SpliceSites$allSequence[wDonor_pos]=StartConserv$V12[match(paste(SpliceSites$chrom[wDonor_pos], SpliceSites$start[wDonor_pos],SpliceSites$strand_HISAT[wDonor_pos]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#SpliceSites$ancestralSequence[wDonor_pos]=StartConserv$V11[match(paste(SpliceSites$chrom[wDonor_pos], SpliceSites$start[wDonor_pos],SpliceSites$strand_HISAT[wDonor_pos]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#SpliceSites$newAge[wDonor_pos]=StartConserv$newAge[match(paste(SpliceSites$chrom[wDonor_pos], SpliceSites$start[wDonor_pos],SpliceSites$strand_HISAT[wDonor_pos]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#
#wDonor_neg=which(SpliceSites$type=='donor' & SpliceSites$strand_Intropolis=='-')
#SpliceSites$NbAlignedPrimates[wDonor_neg]=EndConserv$NbAlignedSeqPrim[match(paste(SpliceSites$chrom[wDonor_neg], SpliceSites$end[wDonor_neg],SpliceSites$strand_HISAT[wDonor_neg]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#SpliceSites$NbMatchingSeqPrim[wDonor_neg]=EndConserv$NbMatchingSeqPrim[match(paste(SpliceSites$chrom[wDonor_neg], SpliceSites$end[wDonor_neg],SpliceSites$strand_HISAT[wDonor_neg]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#SpliceSites$allSequence[wDonor_neg]=EndConserv$V12[match(paste(SpliceSites$chrom[wDonor_neg], SpliceSites$end[wDonor_neg],SpliceSites$strand_HISAT[wDonor_neg]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#SpliceSites$ancestralSequence[wDonor_neg]=EndConserv$V11[match(paste(SpliceSites$chrom[wDonor_neg], SpliceSites$end[wDonor_neg],SpliceSites$strand_HISAT[wDonor_neg]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#SpliceSites$newAge[wDonor_neg]=EndConserv$newAge[match(paste(SpliceSites$chrom[wDonor_neg], SpliceSites$end[wDonor_neg],SpliceSites$strand_HISAT[wDonor_neg]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#
#wAcceptor_pos=which(SpliceSites$type=='acceptor' & SpliceSites$strand_Intropolis=='+')
#SpliceSites$NbAlignedPrimates[wAcceptor_pos]=EndConserv$NbAlignedSeqPrim[match(paste(SpliceSites$chrom[wAcceptor_pos], SpliceSites$end[wAcceptor_pos],SpliceSites$strand_HISAT[wAcceptor_pos]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#SpliceSites$NbMatchingSeqPrim[wAcceptor_pos]=EndConserv$NbMatchingSeqPrim[match(paste(SpliceSites$chrom[wAcceptor_pos], SpliceSites$end[wAcceptor_pos],SpliceSites$strand_HISAT[wAcceptor_pos]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#SpliceSites$allSequence[wAcceptor_pos]=EndConserv$V12[match(paste(SpliceSites$chrom[wAcceptor_pos], SpliceSites$end[wAcceptor_pos],SpliceSites$strand_HISAT[wAcceptor_pos]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#SpliceSites$ancestralSequence[wAcceptor_pos]=EndConserv$V11[match(paste(SpliceSites$chrom[wAcceptor_pos], SpliceSites$end[wAcceptor_pos],SpliceSites$strand_HISAT[wAcceptor_pos]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#SpliceSites$newAge[wAcceptor_pos]=EndConserv$newAge[match(paste(SpliceSites$chrom[wAcceptor_pos], SpliceSites$end[wAcceptor_pos],SpliceSites$strand_HISAT[wAcceptor_pos]),paste(EndConserv$V1,EndConserv$V2,EndConserv$V3))]
#
#wAcceptor_neg=which(SpliceSites$type=='acceptor' & SpliceSites$strand_Intropolis=='-')
#SpliceSites$NbAlignedPrimates[wAcceptor_neg]=StartConserv$NbAlignedSeqPrim[match(paste(SpliceSites$chrom[wAcceptor_neg], SpliceSites$start[wAcceptor_neg],SpliceSites$strand_HISAT[wAcceptor_neg]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#SpliceSites$NbMatchingSeqPrim[wAcceptor_neg]=StartConserv$NbMatchingSeqPrim[match(paste(SpliceSites$chrom[wAcceptor_neg], SpliceSites$start[wAcceptor_neg],SpliceSites$strand_HISAT[wAcceptor_neg]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#SpliceSites$allSequence[wAcceptor_neg]=StartConserv$V12[match(paste(SpliceSites$chrom[wAcceptor_neg], SpliceSites$start[wAcceptor_neg],SpliceSites$strand_HISAT[wAcceptor_neg]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#SpliceSites$ancestralSequence[wAcceptor_neg]=StartConserv$V11[match(paste(SpliceSites$chrom[wAcceptor_neg], SpliceSites$start[wAcceptor_neg],SpliceSites$strand_HISAT[wAcceptor_neg]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]
#SpliceSites$newAge[wAcceptor_neg]=StartConserv$newAge[match(paste(SpliceSites$chrom[wAcceptor_neg], SpliceSites$start[wAcceptor_neg],SpliceSites$strand_HISAT[wAcceptor_neg]),paste(StartConserv$V1,StartConserv$V2,StartConserv$V3))]

########### Annotation High Coverage 
#SpliceSites$highCovAnyCond=apply(sapply(1:5,function(i){SpliceSites[,grep('Total_count_',colnames(SpliceSites))[i],with=F]/table(SampleAnnot$cond)[i]>10}),1,any)
#SpliceSites$highCovAnyCond_junc=SpliceSites$site_id%in%c(Junc_annot$start_id[Junc_annot$highCovAnyCond_junc],Junc_annot$end_id[Junc_annot$highCovAnyCond_junc])  
#ConstitutiveJunc=Junc_annot[which(Junc_annot$highCovAnyCond_junc & Junc_annot$PctStart>0.95 & Junc_annot$PctEnd>0.95),]
#SpliceSites$hasConstitutiveJunc=SpliceSites$site_id%in%c(ConstitutiveJunc$start_id,ConstitutiveJunc$end_id)

#mean(SpliceSites$Gerp_site[SpliceSites$constitutive & SpliceSites$coding & SpliceSites$Total_count_NS>100])


#write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V2.txt',HOME),quote=F,sep='\t',row.names=F)
#SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites.txt',HOME))

########### Annotation Mean Expr and lFC
#
#MeanExpr=GeneAnnot[grep('_mean',cn(GeneAnnot))]
#colnames(MeanExpr)=paste('FPKM',condIndex,sep='_')
#lFC=log2(1+MeanExpr[,-1])-log2(1+MeanExpr[,1])%o%c(1,1,1,1)
#colnames(lFC)=paste('log2FC',condIndex[-1],sep='_')
#SpliceSites=cbind(SpliceSites,MeanExpr[match(SpliceSites$gene,rn(MeanExpr)),],lFC[match(SpliceSites$gene,rn(MeanExpr)),],maxFPKM=apply(MeanExpr,1,max)[match(SpliceSites$gene,rn(MeanExpr))],maxlFC=apply(lFC,1,max)[match(SpliceSites$gene,rn(MeanExpr))],max_abslFC=apply(abs(lFC),1,max)[match(SpliceSites$gene,rn(MeanExpr))])
#
#write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V2.txt',HOME),quote=F,sep='\t',row.names=F)
#
#any(duplicated(SpliceSites$site_id))
#
#SpliceSites$site_id=paste(SpliceSites$site_id,SpliceSites$type)
########### Annotation Constitutive clean

#library(GenomicRanges)
#SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V2.txt',HOME))
#Junc_GR=makeGRangesFromDataFrame(Junc_annot[Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread),])
#w=which(Junc_annot$strand_Intropolis!='' & !is.na(Junc_annot$Nbread))
#site_start=as.numeric(sapply(strsplit(SpliceSites$site_id,' '),function(x){x[2]}))
#site_end=SpliceSites$Position_matchingSite_MajorIsoform
#SpliceSites$junc_start=pmin(site_start,site_end)
#SpliceSites$junc_end=pmax(site_start,site_end)
#Sites_GR=makeGRangesFromDataFrame(SpliceSites,start.field='junc_start',end.field='junc_end')
#oo=findOverlaps(Sites_GR,Junc_GR)
#mean(table(SpliceSites[queryHits(oo),'site_id'])==1) # 3.4%, most splice site have several overlaping junctions mapping
#wTest=1:length(queryHits(oo))
#tic=Sys.time()
#x=Junc_annot$Total_count[w[subjectHits(oo)]] # much faster that way...
#y=SpliceSites$site_id[queryHits(oo)] # much faster that way...
#Total_reads_overlapping_MajorJunction=By(x,y,sum)
#toc=Sys.time()
#print(toc-tic)
#SpliceSites$Total_reads_overlapping_MajorJunction=Total_reads_overlapping_MajorJunction[match(SpliceSites$site_id,names(Total_reads_overlapping_MajorJunction))]
#SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site=SpliceSites$Total_count/SpliceSites$Total_reads_overlapping_MajorJunction

#n1=c(1e4,5e4,1e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,1e6)
#mytime=sapply(n1,function(n){tic=Sys.time()
#Total_reads_overlapping_MajorJunction=By(x[1:n],y[1:n],sum)
#toc=Sys.time()
#print(toc-tic)
#toc-tic
#})
#mytime=c(0.17,2.1,3.3,7.37,12.88,19.05,23.21,36.4,53.5,1.15*60,1.49*60)
#
#n1=5e6
#plot(n1,mytime)
#lines(n1,lm(mytime~ n1)$fitted,col='green')
#lines(n1,lm(mytime~ I(n1^2))$fitted,col='red')
#lines(n1,lm(mytime~ I(n1*log(n1)))$fitted,col='blue')
#lines(n1,lm(mytime~ I(n1^2)+n1+I(n1*log(n1)))$fitted,col='purple')

#Sum_OverlappingReads=Junc_annot[,sum(x),by=DonorSite_pos$site_id[queryHits(oo)]]
#Sum_OverlappingReads=apply(Junc_annot[,grep('Total_count',colnames(DonorSite_pos)),with=F],2,function(x){By(x,DonorSite_pos$site_id[queryHits(oo)],sum,na.rm=T)})
#for(i in colnames(Tot_count_site)){
#	DonorSite_pos[[paste('overlapping',i,sep='_')]]=Sum_OverlappingReads[match(DonorSite_pos$site_id,rownames(Tot_count_site)),i]
#}
##### change start and end to account for the shift in positions compared to UCSC
#dup_site_id=SpliceSites$site_id[duplicated(SpliceSites$site_id)]

#SpliceSites$isAlternative=ifelse(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site<0.95,'yes','')
#write.table(SpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V2.txt',HOME),quote=F,sep='\t',row.names=F)

##### end post-hoc change 



##########################################
##		Numbers in manuscript V6.2		##
##########################################

#SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V2.txt',HOME))
#SpliceSites$isStart=ifelse(gsub(' (donor|acceptor)','',SpliceSites$site_id)==SpliceSites$start_id,TRUE,FALSE)
#SpliceSites$isEnd=ifelse(gsub(' (donor|acceptor)','',SpliceSites$site_id)==SpliceSites$end_id,TRUE,FALSE)
#mean(SpliceSites$isStart & SpliceSites$isEnd)
#mean(!SpliceSites$isStart & !SpliceSites$isEnd)
#mean(SpliceSites$isStart | SpliceSites$isEnd)
#SpliceSites$Zero_Gerp_site=FALSE
#SpliceSites$Zero_Gerp_site[SpliceSites$isStart & (SpliceSites$GerpRS_start==0 | SpliceSites$GerpRS_start2==0)]=TRUE
#SpliceSites$Zero_Gerp_site[SpliceSites$isEnd & (SpliceSites$GerpRS_end==0 | SpliceSites$GerpRS_end2==0)]=TRUE
#SpliceSites$AssociatedGene=ifelse((SpliceSites$isStart & !is.na(SpliceSites$Gene_start)) | (SpliceSites$isEnd & is.na(SpliceSites$Gene_end)), SpliceSites$Gene_start, ifelse((SpliceSites$isStart & is.na(SpliceSites$Gene_start)) | (SpliceSites$isEnd & !is.na(SpliceSites$Gene_end)),SpliceSites$Gene_end,NA))
#SpliceSites$AssociatedSymbol=ifelse((SpliceSites$isStart & !is.na(SpliceSites$symbol_start)) | (SpliceSites$isEnd & is.na(SpliceSites$symbol_end)), SpliceSites$symbol_start, ifelse((SpliceSites$isStart & is.na(SpliceSites$symbol_start)) | (SpliceSites$isEnd & !is.na(SpliceSites$symbol_end)),SpliceSites$symbol_end,NA))
#SpliceSites$AnnotatedGene=ifelse(SpliceSites$isStart & !is.na(SpliceSites$Gene_start), SpliceSites$Gene_start, ifelse(SpliceSites$isEnd & !is.na(SpliceSites$Gene_end),SpliceSites$Gene_end,NA))
#SpliceSites$AnnotatedSymbol=ifelse(SpliceSites$isStart & !is.na(SpliceSites$symbol_start), SpliceSites$symbol_start, ifelse(SpliceSites$isEnd & !is.na(SpliceSites$symbol_end),SpliceSites$symbol_end,NA))
#
# GeneAnnot=read.table(paste(HOME,'/Annotation/GeneAnnotation_hg37_ens70.txt',sep=''),sep='\t',comment='',quote='',stringsAsFactors=FALSE,header=T,row.names=1)
# colnames(GeneAnnot)[1]='Ensembl.Gene.ID'
# GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='X']=23
# GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='Y']=24
# GeneAnnot$Chromosome.Name[GeneAnnot$Chromosome.Name=='MT']=26
# GeneAnnot$Chromosome.Name=as.numeric(GeneAnnot$Chromosome.Name)
# GeneAnnot=GeneAnnot[!is.na(GeneAnnot$Chromosome.Name),]
# GeneAnnot$Strand=ifelse(GeneAnnot$Strand>0,'+','-')
#Gene_GR=makeGRangesFromDataFrame(GeneAnnot,seqnames='Chromosome.Name',start.field='Gene.Start..bp.',end.field='Gene.End..bp.')
#Sites_GR=makeGRangesFromDataFrame(SpliceSites,start.field='start',end.field='end',strand='strand_HISAT')
#oo=findOverlaps(Sites_GR,Gene_GR)
#
#is.unique=function(x){!x %in% x[duplicated(x)]}
#By=function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x}
#
#uniq=is.unique(queryHits(oo))
#SpliceSites$overlappingGene=NA
#SpliceSites$overlappingGene[queryHits(oo)[uniq]]=GeneAnnot$Ensembl.Gene.ID[subjectHits(oo)[uniq]]
#GeneCollapse=By(GeneAnnot$Ensembl.Gene.ID[subjectHits(oo)[!uniq]],queryHits(oo)[!uniq],function(x){paste(sort(unique(x)),collapse='//')})
#SpliceSites$overlappingGene[as.numeric(names(GeneCollapse))]=GeneCollapse
#
#SpliceSites$overlappingSymbol=NA
#SpliceSites$overlappingSymbol[queryHits(oo)[uniq]]=GeneAnnot$Associated.Gene.Name[subjectHits(oo)[uniq]]
#SymbolCollapse=By(GeneAnnot$Associated.Gene.Name[subjectHits(oo)[!uniq]],queryHits(oo)[!uniq],function(x){paste(sort(unique(x)),collapse='//')})
#SpliceSites$overlappingSymbol[as.numeric(names(SymbolCollapse))]=SymbolCollapse
#sum(!is.na(SpliceSites$overlappingSymbol))
#
#MeanExpr=GeneAnnot[grep('_mean',cn(GeneAnnot))]
#MeanExpr[is.na(MeanExpr)]=0
#colnames(MeanExpr)=paste('FPKM',condIndex,sep='_')
#lFC=log2(1+MeanExpr[,-1])-log2(1+MeanExpr[,1])%o%c(1,1,1,1)
#colnames(lFC)=paste('log2FC',condIndex[-1],sep='_')
#
#Expr_data=cbind(MeanExpr[match(SpliceSites$AssociatedGene,rn(MeanExpr)),],lFC[match(SpliceSites$AssociatedGene,rn(MeanExpr)),],maxFPKM=apply(MeanExpr,1,max)[match(SpliceSites$AssociatedGene,rn(MeanExpr))],maxlFC=apply(lFC,1,max)[match(SpliceSites$AssociatedGene,rn(MeanExpr))],max_abslFC=apply(abs(lFC),1,max)[match(SpliceSites$AssociatedGene,rn(MeanExpr))])
#colnames(Expr_data)=paste(colnames(Expr_data),'AssociatedGene',sep='_')
#SpliceSites=cbind(SpliceSites,Expr_data)
#
#Expr_data=cbind(MeanExpr[match(SpliceSites$overlappingGene,rn(MeanExpr)),],lFC[match(SpliceSites$overlappingGene,rn(MeanExpr)),],maxFPKM=apply(MeanExpr,1,max)[match(SpliceSites$overlappingGene,rn(MeanExpr))],maxlFC=apply(lFC,1,max)[match(SpliceSites$overlappingGene,rn(MeanExpr))],max_abslFC=apply(abs(lFC),1,max)[match(SpliceSites$overlappingGene,rn(MeanExpr))])
#colnames(Expr_data)=paste(colnames(Expr_data),'overlappingGene',sep='_')
#SpliceSites=cbind(SpliceSites,Expr_data)
#
#SpliceSites$gene_class=ifelse(is.na(SpliceSites$overlappingGene),'non genic',ifelse(grepl('//',SpliceSites$overlappingGene),'multiple genes',ifelse(SpliceSites$maxFPKM_overlappingGene>1,'expressed gene','non expressed gene')))

#table(ifelse(!is.na(SpliceSites$AssociatedGene),'gene identified','no gene identified'), ifelse(SpliceSites$inEnsembl,'annotated','unannotated'))
#                    annotated unannotated
#  gene identified       435573      769399
#  no gene identified         0      778872


sum(Junc_annot$Total_count/1e9) # 7.234256 -> 7.2 billion spliced reads
sum(Junc_annot$Total_count[is.na(Junc_annot$Nbread]/1e9)/sum(Junc_annot$Total_count/1e9)*100 # 0.159241 -> 0.2 % reads are removed

sum(!is.na(Junc_annot$Nbread))/1e6 # 2.099531 ->2.1 million junctions analysed
nrow(SpliceSites) # 1983844 over 1.9 million
By(SpliceSites$Gerp_site[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align]>2,paste(ifelse(SpliceSites$inEnsembl,'annot','missing'),SpliceSites$WeakAltConst)[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align],mean)
#  annot alternative   annot constitutive     annot weak  missing alternative missing constitutive         missing weak 
#   0.7550677            0.9759993            0.6638520            0.1880780            0.1545101 			0.1972147

mean(SpliceSites$Gerp_site[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align]<2) 
# 0.6556069

By((!SpliceSites$inEnsembl & SpliceSites$WeakAltConst=='weak')[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align], SpliceSites$Gerp_site[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align]>2,mean)
#    FALSE      TRUE 
#0.9065251 0.4243913 

# 91 percent of these were cryptic splice sites absent from Ensembl annotations and presenting weak activity, supporting the high prevalence of “noisy” splicing events previously reported in human cells {Melamud, 2009 #18;Pickrell, 2010 #16}. 
# In contrast, 97% of constitutive and 56% of alternative splice sites while were conserved (GerpRS > 2),
# 96% of non-conserved splice sites exhibited weak activity (Fig 1B and S1D), conserved splice sites were enriched in active sites, with 56% of alternative splice sites and 97% of constitutive splice sites being conserved (GerpRS > 2, Fig 1C). 


mean((SpliceSites$Gerp_site<2 & !SpliceSites$inEnsembl)[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='weak'])
# 0.702796
# While 70% of weak splice sites were cryptic splice sites that are absent from Ensembl annotation and present no sign of conservation, supporting the high prevalence of “noisy” splicing events previously reported in human cells {Melamud, 2009 #18;Pickrell, 2010 #16}, 56% of alternative and 97% of constitutive splice sites showed moderate to high conservation (GerpRS>2) 

mean((SpliceSites$Gerp_site>2)[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='weak'])
# 0.2552072

mean((!SpliceSites$inEnsembl)[SpliceSites$Gerp_site<2 & !SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='weak'])


mean((SpliceSites$Gerp_site<2 & !SpliceSites$inEnsembl)[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='weak'])


By(SpliceSites$newAge[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align]>12,SpliceSites$WeakAltConst[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align],mean,na.rm=T)
#alternative constitutive         weak 
#   0.4940938    0.9229044    0.2385476 
   
table(SpliceSites$WeakAltConst[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site < 2])/sum(!SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site < 2  & SpliceSites$Quality_MultiZ46Align)
#alternative constitutive         weak 
# 0.034131939  0.005453834  0.960414227 
# the vast majority of non conserved sites are weak

table(SpliceSites$WeakAltConst[!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site > 2])/sum(!SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site > 2  & SpliceSites$Quality_MultiZ46Align)
# alternative constitutive         weak 
#  0.08406832   0.28880743   0.62712426 
# but weak sites are also the most frequent among conserved site

# Check sequence at recent but conserved sites
SpliceSites[which(SpliceSites$newAge<5 & SpliceSites$Gerp_site>4)[1:2],]
# Check sequence at old but non conserved sites
	
By(SpliceSites$newAge[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]<=12,SpliceSites$WeakAltConst[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site],mean,na.rm=T)
# alternative constitutive         weak 
#  0.50590621   0.07709557   0.76145239 
By(SpliceSites$newAge[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]==18,SpliceSites$WeakAltConst[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site],mean,na.rm=T)
# alternative constitutive         weak 
#  0.21282849   0.55860229   0.09426126 

mean(SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]!='weak',na.rm=T) # 0.0500166
mean(SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]=='constitutive',na.rm=T) # 0.0117122
sum(SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]!='weak',na.rm=T) # 57993
sum(SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site]=='constitutive',na.rm=T) # 13580

mean(!SpliceSites$coding[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site & SpliceSites$WeakAltConst!='weak'],na.rm=T) # 
mean(!SpliceSites$coding[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site & SpliceSites$WeakAltConst=='constitutive'],na.rm=T) # 0.5478465
table(SpliceSites$coding[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site],SpliceSites$WeakAltConst[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site])
#        alternative constitutive    weak
#  FALSE       37189         6020 1080963
#  TRUE         7224         7560   20519

mean(SpliceSites$coding[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='constitutive'],na.rm=T) # 0.5478465

SpliceSites$[SpliceSites$newAge<=12 & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='constitutive' & SpliceSites$coding]



# Fig 1A: 6x 5.5 inches - Conservation of Splice sites and constitutive Splice sites
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align)
library(data.table)
Gerp_mean=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_2bp_random.txt',HOME))$V1

par(mar=c(4,5,5,0.5))
hh=hist(c(-10,sample(Gerp_mean,sum(!SpliceSites$Zero_Gerp_site),replace=T),Inf),br=50,plot=F)
hh$counts=hh$count/1000
plot(hh,xlim=c(-7,7),col=rcol[7],xlab='GerpRS',ylab='count',las=1,main='')
hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[!SpliceSites$Zero_Gerp_site])),br=hh$br,plot=F);
hh1$counts=hh1$count/1000;plot(hh1,col=gsub('80',"BB",rcol[3]),add=T,border='#00000022')
#hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count>1000])),br=hh$br,col=rcol[1],add=T)
hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[!SpliceSites$Zero_Gerp_site & SpliceSites$WeakAltConst=='constitutive'])),br=hh$br,plot=F);
hh1$counts=hh1$count/1000;plot(hh1,col=gsub('80',"BB",rcol[1]),add=T,border='#00000022')
par(xpd=T);legend(-10,800,fill=c(rcol[7],gsub('80',"BB",rcol[c(3,1)])),legend=c('Genome Wide','splice sites','constitutive splice sites'),bty='n');par(xpd=F)
#hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$inEnsembl])),br=hh$br,plot=F);
#hh1$counts=hh1$count/1000;plot(hh1,col=rcol[1],add=T,border='#00000022')
#hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Nbread>10000])),br=hh$br,plot=F);
#hh1$counts=hh1$count/1000;plot(hh1,density=30,add=T,border='#000000')

# Fig S1B: 6x 4.5 inches - coverage of weak, alternative and constitutive splice sites
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak')
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst=='constitutive')
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=mergeCols(acol[3],acol[4],0.9),xlab='log2 Nb read per sample',ylab='count',las=1,main='')
plot(hh1,col=acol[6],add=T,border='#000000')
plot(hh2,col=acol[8],add=T,border='#000000')
par(xpd=T);legend('topright',fill=rev(c(acol[c(8,6)],mergeCols(acol[3],acol[4],0.9))),legend=c('Weak','Alternative','Constitutive'),bty='n');par(xpd=F)


# Fig S1B: 6x 4.5 inches - coverage of conserved and non conserved splice sites
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align )
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2)
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4)
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
#  Fig S1B: 6x 4.5 inches - 
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='')
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

# Fig S1B: 6x 4.5 inches - coverage of conserved and non conserved splice sites by levels of gene expression
layout(1:4)
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene<1)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene<1)
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene<1)
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

# Fig S1B: 6x 4.5 inches - coverage of conserved and non conserved splice sites by levels of gene expression
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10)
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10)
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)


w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100)
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100)
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>100)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>100)
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>100)
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

# Fig S1B: 6x 4.5 inches - coverage of conserved and non conserved splice sites by levels of gene expression among annotated splice sites
layout(1:4)
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene<1 & !is.na(SpliceSites$AnnotatedGene))
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene<1 & !is.na(SpliceSites$AnnotatedGene))
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene<1 & !is.na(SpliceSites$AnnotatedGene))
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10 & !is.na(SpliceSites$AnnotatedGene))
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10 & !is.na(SpliceSites$AnnotatedGene))
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10 & !is.na(SpliceSites$AnnotatedGene))
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100 & !is.na(SpliceSites$AnnotatedGene))
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100 & !is.na(SpliceSites$AnnotatedGene))
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100& !is.na(SpliceSites$AnnotatedGene))
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>100 & !is.na(SpliceSites$AnnotatedGene))
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>100 & !is.na(SpliceSites$AnnotatedGene))
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>100 & !is.na(SpliceSites$AnnotatedGene))
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

# Fig S1B: 6x 4.5 inches - coverage of conserved and non conserved splice sites by levels of gene expression among annotated splice sites
layout(1:4)
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene<1 & is.na(SpliceSites$AnnotatedGene))
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene<1 & is.na(SpliceSites$AnnotatedGene))
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene<1 & is.na(SpliceSites$AnnotatedGene))
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10 & is.na(SpliceSites$AnnotatedGene))
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10 & is.na(SpliceSites$AnnotatedGene))
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10 & is.na(SpliceSites$AnnotatedGene))
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100 & is.na(SpliceSites$AnnotatedGene))
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100 & is.na(SpliceSites$AnnotatedGene))
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100& is.na(SpliceSites$AnnotatedGene))
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>100 & is.na(SpliceSites$AnnotatedGene))
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>100 & is.na(SpliceSites$AnnotatedGene))
hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>100 & is.na(SpliceSites$AnnotatedGene))
hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='log2 Nb read per sample',ylab='count',las=1,main='',xlim=c(-9,17))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

# Fig S1B: 6x 4.5 inches - PctofReadCoveringSplice_LogitTransformed_AndConservation.pdf
layout(1:4)
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene<1)
hh=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene<1)
hh1=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene<1)
hh2=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='% of transcripts including splice site',ylab='count',las=1,main='',xlim=c(-8,8))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10)
hh=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10)
hh1=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>1 & SpliceSites$FPKM_NS_AssociatedGene<10)
hh2=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='% of transcripts including splice site',ylab='count',las=1,main='',xlim=c(-8,8))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)


w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100)
hh=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100)
hh1=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>10 & SpliceSites$FPKM_NS_AssociatedGene<100)
hh2=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='% of transcripts including splice site',ylab='count',las=1,main='',xlim=c(-8,8))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)

w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$FPKM_NS_AssociatedGene>100)
hh=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=30,plot=F)
hh$counts=hh$count/1000
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>2 & SpliceSites$FPKM_NS_AssociatedGene>100)
hh1=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=hh$br,plot=F)
hh1$counts=hh1$count/1000;
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$Gerp_site>4 & SpliceSites$FPKM_NS_AssociatedGene>100)
hh2=hist(logit(SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site[w],d=1e-3),br=hh$br,plot=F)
hh2$counts=hh2$count/1000;
par(mar=c(4,5,1,0.5))
plot(hh,col=rcol[7],xlab='% of transcripts including splice site',ylab='count',las=1,main='',xlim=c(-8,8))
plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
par(xpd=T);legend('topright',fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),legend=c(expression(GerpRS < 2),expression(paste(GerpRS%in%textstyle(' '),'[2,4]')),expression(GerpRS > 4)),bty='n');par(xpd=F)



w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align)
By(SpliceSites$Gerp_site[w]>2,sum)

##Fig S1C: 6x5 inches - coverage of recent and ancient  splice sites
#w=which(!SpliceSites$Zero_Gerp_site &SpliceSites$Quality_MultiZ46Align)
#hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
#hh$counts=hh$count/1000
#w=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$newAge<=12)
#hh0=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
#hh0$counts=hh0$count/1000;
#w=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$newAge<=9)
#hh1=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
#hh1$counts=hh1$count/1000;
#w=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$newAge<=5)
#hh2=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
#hh2$counts=hh2$count/1000;
#par(mar=c(4,5,1,0.5))
#plot(hh,col=acol[8],xlab='log2 Nb read per sample',ylab='count',las=1,main='')
#plot(hh0,col=mergeCols(acol[7],'white',0.2),add=T,border='#000000')
#plot(hh1,col=gsub('80',"FF",rcol[3]),add=T,border='#000000')
#plot(hh2,col=gsub('80',"BB",rcol[1]),add=T,border='#000000')
#par(xpd=T);legend('topright',fill=c(acol[8],mergeCols(acol[7],'white',0.2),acol[3],gsub('80',"BB",rcol[c(1)])),legend=c('vertebrates or older','placental mammals','primates','apes'),bty='n');par(xpd=F)

#  Fig S1C_V1: 6x4.5 inches - coverage of recent and ancient  splice sites
w=which(!SpliceSites$Zero_Gerp_site &SpliceSites$Quality_MultiZ46Align)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
plot(hh,col=acol[8],xlab='log2 Nb read per sample',ylab='count',las=1,main='')
j=0
for (i in rev(c(5,9,12,13,14,15,16,17,18))){
	j=j+1
	w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$newAge<=i)
	hh0=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
	hh0$counts=hh0$count/1000;
	plot(hh0,col=mycols[j],add=T,border='#00000000')
	}
plot(hh,col='#00000000',add=T)

#  Fig S1C_V2: 6x4.5 inches - coverage of recent and ancient  splice sites

w=which(!SpliceSites$Zero_Gerp_site &SpliceSites$Quality_MultiZ46Align)
hh=hist(log2(apply(SS_Activity,1,max)[w]),br=30,plot=F)
hh$counts=hh$count/1000
plot(hh,col=acol[8],xlab='log2 Nb read per sample',ylab='count',las=1,main='')
j=0
for (i in c(1,5,9,12,13,14,15,16,17)){
	j=j+1
	w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align & SpliceSites$newAge>=i)
	hh0=hist(log2(apply(SS_Activity,1,max)[w]),br=hh$br,plot=F)
	hh0$counts=hh0$count/1000;
	plot(hh0,col=rev(mycols)[j],add=T,border='#00000000')
	}
plot(hh,col='#00000000',add=T)

#par(xpd=T);legend('topright',fill=c(acol[8],mergeCols(acol[7],'white',0.2),acol[3],gsub('80',"BB",rcol[c(1)])),legend=c('vertebrates or older','placental mammals','primates','apes'),bty='n');par(xpd=F)

#hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count>1000])),br=hh$br,col=rcol[1],add=T)
hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[!SpliceSites$Zero_Gerp_site & SpliceSites$WeakAltConst=='constitutive'])),br=hh$br,plot=F);
hh1$counts=hh1$count/1000;plot(hh1,col=gsub('80',"BB",rcol[1]),add=T,border='#00000022')

# Fig 1B left: 4x 5.5 inches - Conservation of weak constitutive and alternative Splice sites
par(mar=c(4,5,1,0.5))
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align)
x=cut(SpliceSites$Gerp_site[w],c(-Inf,2,4,Inf))
levels(x)=c('<2','[2;4]','>4')
par(mar=c(4,5,1,0.5))
tab=table(x,SpliceSites$WeakAltConst[w])
tab=tab[,c('weak','alternative','constitutive')]
colnames(tab)=c('weak','alt.','const.')
t(t(tab)/apply(tab,2,sum))*100
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),ylab='% of splice sites',las=1)
par(xpd=T)
legend(-2,-15,legend=levels(x),fill=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)])),bty='n',ncol=4)
par(xpd=F)

# Fig 1B right: 4x 5.5 inches - Conservation of weak constitutive and alternative Splice sites
par(mar=c(4,5,1,0.5))
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align)
x=cut(SpliceSites$Gerp_site[w],c(-Inf,2,4,Inf))
tab=table(x,SpliceSites$WeakAltConst[w])
tab=tab[,c('constitutive','alternative','weak')]
rownames(tab)=c('< 2','[2-4]',' >4')
t(tab/apply(tab,1,sum))*100
barplot(t(tab/apply(tab,1,sum))*100,beside=F,border=NA,col=c(acol[c(8,6)],mergeCols(acol[3],acol[4],0.9)),ylab='% of splice sites',las=1)
par(xpd=T)
legend(-2.5,-15,legend=c('const.','alt.','weak'),fill=c(acol[c(8,6)],mergeCols(acol[3],acol[4],0.9)),bty='n',ncol=4)
par(xpd=F)

# Fig 1B bis: 6x 6 inches - Conservation of weak constitutive and alternative Splice sites
par(mar=c(4,5,1,0.5))
w=which(!SpliceSites$Zero_Gerp_site & SpliceSites$Quality_MultiZ46Align)
x=cut(SpliceSites$Gerp_site[w],c(-Inf,2,4,Inf))
levels(x)=c('<2','[2;4]','>4')
par(mar=c(4,5,1,0.5))
tab=table(x,SpliceSites$WeakAltConst[w])
tab=tab[,c('weak','alternative','constitutive')]

library(DescTools)
library(RColorBrewer)
makeTransparent=function(col,alphaAdd=1){RGB=col2rgb(col);alphaMax=1-alphaAdd*min(RGB)/255;rgb((RGB[1]-alphaAdd*min(RGB))/alphaMax/255,(RGB[2]-alphaAdd*min(RGB))/alphaMax/255,(RGB[3]-alphaAdd*min(RGB))/alphaMax/255,alphaMax*alphaAdd)}
source(sprintf('%s/01_scripts/Z_ProgUtils/PlotCirc2.R',HOME))
PlotCirc( tab,acol = c(acol[c(8,6)],mergeCols(acol[3],acol[4],0.9),rev(c('#91D3DB','#FCC34B','#E3444C'))),rcol = sapply(c('#91D3DB','#FCC34B','#E3444C'),makeTransparent,0.5),radius.in=8,gap=8)

}l, alpha = FALSE
# Fig 1C: 4x 5.5 inches - age of weak constitutive and alternative Splice sites

w=which(SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site)
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
par(mar=c(4,5,1,0.5))
tab=table(x,SpliceSites$WeakAltConst[w])
tab=tab[,c('weak','alternative','constitutive')]
colnames(tab)=c('weak','alt.','const.')
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)
#tab=table(x,SpliceSites$isImmune[w])
#barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)


# Fig 1D : 5x 8 inches - age of non-immune,immune and stimulation induced splices sites. 
wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & !SpliceSites$Zero_Gerp_site )
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & !SpliceSites$Zero_Gerp_site)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & !SpliceSites$Zero_Gerp_site)

wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site)
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site)

wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site )
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)

wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site & SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site > 0.05)
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site & SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site > 0.05)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site & SpliceSites$Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site > 0.05)



#w=wConst
#NbNonConserved=By(SpliceSites$Gerp_site[w]<2,paste(SpliceSites$gene[w],SpliceSites$isImmune[w]),mean)
#NbConserved=By(SpliceSites$Gerp_site[w]>2,paste(SpliceSites$gene[w],SpliceSites$isImmune[w]),mean)

x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')

layout(matrix(c(1,1,1,2,2,2,3,3,3),nrow=1))
w=wActive
par(mar=c(4,5,1,0.5))
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC[w]>1])
tab3=table(x[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

par(mar=c(4,5,1,0.5))
w=wConst
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC[w]>1])
tab3=table(x[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

par(mar=c(4,5,1,0.5))
w=wAlt
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
y=paste(ifelse(SpliceSites$isImmune[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC[w]>1])
tab3=table(x[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

# Fig 1D bis : 5x 8 inches - Conservation of  non-immune,immune and stimulation induced splices sites. 
colGerp=c(rcol[7],acol[3],gsub('80',"BB",rcol[c(1)]))
wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site)
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site)
#wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
#wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
#wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)


layout(matrix(c(1,1,1,2,2,2,3,3,3),nrow=1))
w=wActive
par(mar=c(4,5,1,0.5))
x=cut(SpliceSites$Gerp_site[w],c(-Inf,2,4,Inf))
levels(x)=c('<2','[2;4]','>4')
y=paste(ifelse(SpliceSites$isImmune[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC[w]>1])
tab3=table(x[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=colGerp,ylab='% of splice sites',las=1)

par(mar=c(4,5,1,0.5))
w=wConst
x=cut(SpliceSites$Gerp_site[w],c(-Inf,2,4,Inf))
levels(x)=c('<2','[2;4]','>4')

y=paste(ifelse(SpliceSites$isImmune[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC[w]>1])
tab3=table(x[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=colGerp,ylab='% of splice sites',las=1)

par(mar=c(4,5,1,0.5))
w=wAlt
x=cut(SpliceSites$Gerp_site[w],c(-Inf,2,4,Inf))
levels(x)=c('<2','[2;4]','>4')
y=paste(ifelse(SpliceSites$isImmune[w]=='yes','immune','non-immune'),ifelse(SpliceSites$maxlFC[w]>1,'up-regulated',''))
tab0=table(x)
tab1=table(x[SpliceSites$isImmune[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC[w]>1])
tab3=table(x[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=colGerp,ylab='% of splice sites',las=1)


# immunotype=ifelse(SpliceSites$isImmune=='yes',ifelse(SpliceSites$maxlFC>1,'immune upregulated','immune'),ifelse(SpliceSites$maxlFC>1,'upregulated','non immune'))
# By(SpliceSites$Time_FirstAppeared_MA[w], immunotype[w],mean,na.rm=T)
# table(immunotype)


wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)

wlist=list(wActive,wAlt,wConst)
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w])}) # 444.8737 371.3778 463.2273
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']])}) # 399.4866 292.1525 424.9786
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)}) # 389.7211 318.3872 414.5104


wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site & SpliceSites$maxFPKM>10)
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site & SpliceSites$maxFPKM>10)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site & SpliceSites$maxFPKM>10)

wlist=list(wActive,wAlt,wConst)
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w])}) # 424.2217 286.6537 464.9397
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']])}) # 367.0712 223.5247 411.9642
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)}) # 354.7741 239.3022 404.2407

wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )

wlist=list(wActive,wAlt,wConst)
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w])}) # 413.3139 299.4864 451.4007
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']])}) # 374.2398 241.2150 416.9097
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)}) # 349.9137 246.5937 400.4768




############# repeat with times from Kumar et al



wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)

wlist=list(wActive,wAlt,wConst)
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w])}) # 443.4848 372.7825 461.1408
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']])}) # 400.1325 296.1644 424.8250
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)}) # 390.9368 321.8874 414.9323


wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )

wlist=list(wActive,wAlt,wConst)
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w])}) # 412.8059 302.5242 449.7062
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']])}) # 375.6037 246.2769 417.0874
sapply(wlist,function(w){mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)}) # 351.9984 251.2187 401.3183




##########  nonConserved splice sites
wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site )

w=wConst
GeneHasNonConservedConstSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedConstSite) # 9% of genes have a constitutive splice site (>1 read per sample) that is non conserved
resGO_NcConst=GOSeq(names(GeneHasNonConservedConstSite)[GeneHasNonConservedConstSite], names(GeneHasNonConservedConstSite),addCI=T,FDR=0.2)
resGO_NcConst$type='constitutive'
resGO_NcConst=resGO_NcConst[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wActive
GeneHasNonConservedActSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedActSite) # 9% of genes have a constitutive splice site (>1 read per sample) that is non conserved
resGO_NcAct=GOSeq(names(GeneHasNonConservedActSite)[GeneHasNonConservedActSite], names(GeneHasNonConservedActSite),addCI=T,FDR=0.2)
resGO_NcAct$type='active'
resGO_NcAct=resGO_NcAct[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wAlt
GeneHasNonConservedAltSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedAltSite) # 19% of genes have a alternative splice site (>1 read per sample) that is non conserved
resGO_NcAlt=GOSeq(names(GeneHasNonConservedAltSite)[GeneHasNonConservedAltSite], names(GeneHasNonConservedAltSite),addCI=T,FDR=0.2)
resGO_NcAlt$type='alternative'
resGO_NcAlt=resGO_NcAlt[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wConst
GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
OR_all=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
GenelFC=By(SpliceSites$log2FC_LPS[w],SpliceSites$gene[w],max)
OR_LPS=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
GenelFC=By(SpliceSites$log2FC_PAM3[w],SpliceSites$gene[w],max)
OR_PAM3=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
GenelFC=By(SpliceSites$log2FC_R848[w],SpliceSites$gene[w],max)
OR_R848=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
GenelFC=By(SpliceSites$log2FC_IAV[w],SpliceSites$gene[w],max)
OR_IAV=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
OR_LFC_Const=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='constitutive')

#resGO_NcConst[order(-resGO_NcConst$lowerCI),1:16]

w=wAlt
GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
OR_all=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
GenelFC=By(SpliceSites$log2FC_LPS[w],SpliceSites$gene[w],max)
OR_LPS=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
GenelFC=By(SpliceSites$log2FC_PAM3[w],SpliceSites$gene[w],max)
OR_PAM3=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
GenelFC=By(SpliceSites$log2FC_R848[w],SpliceSites$gene[w],max)
OR_R848=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
GenelFC=By(SpliceSites$log2FC_IAV[w],SpliceSites$gene[w],max)
OR_IAV=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
OR_LFC_Alt=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='alternative')

#resGO_NcAlt[order(-resGO_NcAlt$lowerCI),1:16]
#resGO_NcAlt[,1:16]
w=wActive
GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
OR_all=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
GenelFC=By(SpliceSites$log2FC_LPS[w],SpliceSites$gene[w],max)
OR_LPS=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
GenelFC=By(SpliceSites$log2FC_PAM3[w],SpliceSites$gene[w],max)
OR_PAM3=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
GenelFC=By(SpliceSites$log2FC_R848[w],SpliceSites$gene[w],max)
OR_R848=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
GenelFC=By(SpliceSites$log2FC_IAV[w],SpliceSites$gene[w],max)
OR_IAV=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
OR_LFC_Act=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='active')


write.table(rbind(OR_LFC_Const,OR_LFC_Alt),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1B_OR_LFC.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(resGO_NcConst,resGO_NcAlt),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1A_GO_nonConserved.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(resGO_NcConst,resGO_NcAlt,resGO_NcAct),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1A_GO_nonConserved_withActive.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(OR_LFC_Const,OR_LFC_Alt,OR_LFC_Act),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1B_OR_LFC_withActive.txt',HOME),quote=F,sep='\t',row.names=FALSE)

# repeats with these definitions 
wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)

#wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site> -2)
#wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site> -2)
#wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site> -2)


w=wConst
GeneHasNonConservedConstSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedConstSite) # 9% of genes have a constitutive splice site (>1 read per sample) that is non conserved
resGO_NcCodingConst=GOSeq(names(GeneHasNonConservedConstSite)[GeneHasNonConservedConstSite], names(GeneHasNonConservedConstSite),addCI=T,FDR=0.2)
resGO_NcCodingConst$type='constitutive'
resGO_NcCodingConst=resGO_NcCodingConst[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wActive
GeneHasNonConservedActSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedActSite) # 9% of genes have a constitutive splice site (>1 read per sample) that is non conserved
resGO_NcCodingAct=GOSeq(names(GeneHasNonConservedActSite)[GeneHasNonConservedActSite], names(GeneHasNonConservedActSite),addCI=T,FDR=0.2)
resGO_NcCodingAct$type='active'
resGO_NcCodingAct=resGO_NcCodingAct[,c(1:5,7:8,10:11,13:17,6,18,12)]

mean(GeneHasNonConservedAltSite) # 19% of genes have a alternative splice site (>1 read per sample) that is non conserved
resGO_NcCodingAlt=GOSeq(names(GeneHasNonConservedAltSite)[GeneHasNonConservedAltSite], names(GeneHasNonConservedAltSite),addCI=T,FDR=0.2)
resGO_NcCodingAlt$type='alternative'
resGO_NcCodingAlt=resGO_NcCodingAlt[,c(1:5,7:8,10:11,13:17,6,18,12)]

# unbiased 
w=wAlt
GeneHasNonConservedAltSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
NbAltSiteByGene=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],length)
resGO_NcCodingAlt_unbiased=GOSeq(names(GeneHasNonConservedAltSite)[GeneHasNonConservedAltSite], names(GeneHasNonConservedAltSite),addCI=T,FDR=0.2,bias.data=NbAltSiteByGene)
resGO_NcCodingAlt_unbiased$type='alternative'
resGO_NcCodingAlt_unbiased=resGO_NcCodingAlt_unbiased[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wConst
GeneHasNonConservedConstSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
NbConstSiteByGene=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],length)
resGO_NcCodingConst_unbiased=GOSeq(names(GeneHasNonConservedConstSite)[GeneHasNonConservedConstSite], names(GeneHasNonConservedConstSite),addCI=T,FDR=0.2,bias.data=NbConstSiteByGene)
resGO_NcCodingConst_unbiased$type='constitutive'
resGO_NcCodingConst_unbiased=resGO_NcCodingConst_unbiased[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wActive
GeneHasNonConservedActSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
NbActSiteByGene=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],length)
resGO_NcCodingAct_unbiased=GOSeq(names(GeneHasNonConservedActSite)[GeneHasNonConservedActSite], names(GeneHasNonConservedActSite),addCI=T,FDR=0.2,bias.data=NbActSiteByGene)
resGO_NcCodingAct_unbiased$type='active'
resGO_NcCodingAct_unbiased=resGO_NcCodingAct_unbiased[,c(1:5,7:8,10:11,13:17,6,18,12)]

w=wConst
GeneHasNonConservedConstSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)

GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
OR_all=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
GenelFC=By(SpliceSites$log2FC_LPS[w],SpliceSites$gene[w],max)
OR_LPS=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
GenelFC=By(SpliceSites$log2FC_PAM3[w],SpliceSites$gene[w],max)
OR_PAM3=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
GenelFC=By(SpliceSites$log2FC_R848[w],SpliceSites$gene[w],max)
OR_R848=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
GenelFC=By(SpliceSites$log2FC_IAV[w],SpliceSites$gene[w],max)
OR_IAV=odds.ratio(table(GenelFC>1,GeneHasNonConservedConstSite))
OR_LFC_Coding_Const=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='constitutive')

#resGO_NcConst[order(-resGO_NcConst$lowerCI),1:16]

w=wAlt
GeneHasNonConservedAltSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)

GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
OR_all=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
GenelFC=By(SpliceSites$log2FC_LPS[w],SpliceSites$gene[w],max)
OR_LPS=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
GenelFC=By(SpliceSites$log2FC_PAM3[w],SpliceSites$gene[w],max)
OR_PAM3=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
GenelFC=By(SpliceSites$log2FC_R848[w],SpliceSites$gene[w],max)
OR_R848=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
GenelFC=By(SpliceSites$log2FC_IAV[w],SpliceSites$gene[w],max)
OR_IAV=odds.ratio(table(GenelFC>1,GeneHasNonConservedAltSite))
OR_LFC_Coding_Alt=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='alternative')

#resGO_NcAlt[order(-resGO_NcAlt$lowerCI),1:16]
#resGO_NcAlt[,1:16]
w=wActive
GeneHasNonConservedActSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
OR_all=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
GenelFC=By(SpliceSites$log2FC_LPS[w],SpliceSites$gene[w],max)
OR_LPS=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
GenelFC=By(SpliceSites$log2FC_PAM3[w],SpliceSites$gene[w],max)
OR_PAM3=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
GenelFC=By(SpliceSites$log2FC_R848[w],SpliceSites$gene[w],max)
OR_R848=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
GenelFC=By(SpliceSites$log2FC_IAV[w],SpliceSites$gene[w],max)
OR_IAV=odds.ratio(table(GenelFC>1,GeneHasNonConservedActSite))
OR_LFC_Coding_Act=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='active')



write.table(rbind(OR_LFC_Coding_Const,OR_LFC_Coding_Alt, OR_LFC_Coding_Act),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1B_OR_LFC_coding.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(resGO_NcCodingConst,resGO_NcCodingAlt,resGO_NcCodingAct),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1A_GO_nonConserved_Coding.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(resGO_NcCodingConst_unbiased,resGO_NcCodingAlt_unbiased,resGO_NcCodingAct_unbiased),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1A_GO_nonConserved_Coding_unbiased.txt',HOME),quote=F,sep='\t',row.names=FALSE)

wActive=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site & SpliceSites$maxFPKM>10)
wConst=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='constitutive' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site & SpliceSites$maxFPKM>10)
wAlt=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons=='alternative' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site & SpliceSites$maxFPKM>10)

write.table(rbind(OR_LFC_Const,OR_LFC_Alt),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1B_OR_LFC_FPKMover10.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(resGO_NcConst,resGO_NcAlt,resGO_NcAct),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1A_GO_nonConserved_FPKMover10.txt',HOME),quote=F,sep='\t',row.names=FALSE)


nrow(SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & !SpliceSites$Zero_Gerp_site,])
#[1] 235646
nrow(SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & SpliceSites$NbMatchingSeqPrim==1 & !SpliceSites$Zero_Gerp_site,])
# 293

nrow(SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & SpliceSites$coding & !SpliceSites$Zero_Gerp_site,])
#[1] 195774
nrow(SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & SpliceSites$coding & SpliceSites$NbMatchingSeqPrim==1 & !SpliceSites$Zero_Gerp_site,])
# 43

GeneAge=By(SpliceSites$newAge[w],SpliceSites$gene[w],max)
SpliceSites$GeneAge=GeneAge[match(SpliceSites$gene,names(GeneAge))]

reverseSeq=function(seq){paste(sapply(strsplit(seq,':'),function(x){nn=nchar(x);paste(substr(x,1,nn-3),reverse(substr(x,nn-1,nn)),sep='_')}),collapse=':')}

#plot(tree,type="fan")
strand='+'
reverse=function(Seq){
	Corresp=c('A','G','T','C','-')
	names(Corresp)=c('T','C','A','G','-')
	sapply(strsplit(Seq,''),function(x){paste(rev(Corresp[x]),collapse='')})
	}
	
HumanSpecificSpliceSites=SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=8 & SpliceSites$NbMatchingSeqPrim==1 & !SpliceSites$Zero_Gerp_site,]
# Chromosome	Splice site start	Splice site end	strand	type	Sequence of splice site	Gene associated to the splice site (Ensembl ID)	Gene associated to the splice site (symbol)	Present in over 95% of observed transcripts	Total read count NS	Gene FPKM NS	Total read count LPS	Gene FPKM LPS	Total read count Pam3Csk4	Gene FPKM Pam3Csk4	Total read count R848	Gene FPKM R848	Total read count IAV	Gene FPKM IAV	Mean GerpRS of splice site 	Fist apparition of the current splice site sequence in the phylogeny	Nb of primates with an alignment	Splice site sequence across common ancestors	Splice site sequence across 46 vertebrates
HumanSpecificSpliceSites$seq=ifelse(HumanSpecificSpliceSites$type=='donor',HumanSpecificSpliceSites$donor,HumanSpecificSpliceSites$acceptor)
HumanSpecificSpliceSites$allSequence[HumanSpecificSpliceSites$strand_Intropolis=='-']=sapply(HumanSpecificSpliceSites$allSequence[HumanSpecificSpliceSites$strand_Intropolis=='-'],reverseSeq)
HumanSpecificSpliceSites$ancestralSequence[HumanSpecificSpliceSites$strand_Intropolis=='-']=sapply(HumanSpecificSpliceSites$ancestralSequence[HumanSpecificSpliceSites$strand_Intropolis=='-'],reverseSeq)

HumanSpecificSpliceSites$ancestralSequence=gsub('treeShrew','treeshrew',HumanSpecificSpliceSites$ancestralSequence)
HumanSpecificSpliceSites$ancestralSequence=sapply(strsplit(HumanSpecificSpliceSites$ancestralSequence,':'),function(x){paste(toupper(substr(unlist(x),1,1)),substr(x,2,10000),sep='',collapse=':')})
HumanSpecificSpliceSites$ancestralSequence=gsub(':',' / ',HumanSpecificSpliceSites$ancestralSequence)
HumanSpecificSpliceSites$ancestralSequence=gsub('_',': ',HumanSpecificSpliceSites$ancestralSequence)

HumanSpecificSpliceSites$allSequence=gsub('X_tropicalis','X. tropicalis',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Rock_hyrax','Rock hyrax',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Kangaroo_rat','Kangaroo rat',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('TreeShrew','Treeshrew',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Zebra_finch','Zebra finch',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('lemur','Lemur',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Guinea_Pig','Guinea pig',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub(':',' / ',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('_',': ',HumanSpecificSpliceSites$allSequence)


HumanSpecificSpliceSites=HumanSpecificSpliceSites[,c("chrom","start","end","strand_Intropolis","seq","type","gene","symbol","isImmune",'read_persample_NS',"FPKM_NS",'read_persample_LPS',"FPKM_LPS",'read_persample_PAM3CSK4',"FPKM_PAM3CSK4",'read_persample_R848',"FPKM_R848",'read_persample_IAV',"FPKM_IAV","Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site","WeakAltConst",'Gerp_site','NbAlignedPrimates',"allSequence","ancestralSequence","coding","inEnsembl","AssociatedGene",'AssociatedSymbol')]
write.table(HumanSpecificSpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/TableS1C_V6_HumanSpecificSpliceSites',HOME),quote=F,sep='\t',row.names=FALSE)

SpliceSites$Time_FirstAppeared_MA=

w=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
mean(SpliceSites$Time_FirstAppeared_MA[w]) # 444.8737
mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']]) # 399.4866
mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T) # 389.7211

w=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltCons!='weak' & SpliceSites$inEnsembl & !SpliceSites$Zero_Gerp_site)
mean(SpliceSites$Time_FirstAppeared_MA[w]) # 413.3139
mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']]) # 374.2398
mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T) # 349.9137

resamp_Age_maxFC=sapply(1:1000,function(i){toSample=By(SpliceSites$maxlFC[w]>1,SpliceSites$GeneAge[w],sum,na.rm=T)
	samp=unlist(lapply(names(toSample),function(i){sample(which(SpliceSites$GeneAge[w]==i),toSample[i])}))
	mean(SpliceSites$Time_FirstAppeared_MA[w[samp]],na.rm=T)})
mean(resamp_Age_maxFC<=mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)) # P<10-4

resamp_Age_immune=sapply(1:1000,function(i){toSample=By(SpliceSites$isImmune[w]=='yes',SpliceSites$GeneAge[w],sum,na.rm=T)
	samp=unlist(lapply(names(toSample),function(i){sample(which(SpliceSites$GeneAge[w]==i),toSample[i])}))
	mean(SpliceSites$Time_FirstAppeared_MA[w[samp]],na.rm=T)})
mean(resamp_Age_maxFC<=mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)) # P<10-4

save(resamp_Age_immune,resamp_Age_maxFC,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/resamp_Age_spliceSites.Rdata',HOME))

shapiro.test(resamp_Age_immune) # P= 0.38
shapiro.test(resamp_Age_maxFC)	# P=0.1622

#tab=table(x,y)
#barplot(tab/apply(tab,2,sum)*100,beside=F,border=NA,col=rev(mycols),ylab='% of constitutive\n splice sites',las=1)


#tab=table(x,SpliceSites$isImmune[w])
#barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

x=cut(SpliceSites$newAge[wConst],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')

tab0=table(x)
tab1=table(x[SpliceSites$isImmune[wConst]=='yes'])
tab2=table(x[SpliceSites$maxlFC[wConst]>1])
tab3=table(x[SpliceSites$maxlFC[wConst]>1 & SpliceSites$isImmune[wConst]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of constitutive\n splice sites',las=1)

x=cut(SpliceSites$newAge[wAlt],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')

tab0=table(x)
tab1=table(x[SpliceSites$isImmune[wAlt]=='yes'])
tab2=table(x[SpliceSites$maxlFC[wAlt]>1])
tab3=table(x[SpliceSites$maxlFC[wAlt]>1 & SpliceSites$isImmune[wAlt]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of alternative\n splice sites',las=1)
#
#wilcox.test(as.numeric(x)[isImmune=='immune'],as.numeric(x))$p.value# p-value = 3.581088e-187
#wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1],as.numeric(x))$p.value# p-value = 4.79261e-281
#wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'],as.numeric(x)[SpliceSites$maxlFC[w]>1])$p.value# p-value = 2.220398e-115
#wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'],as.numeric(x)[SpliceSites$isImmune[w]=='yes'])$p.value# p-value = 4.688095e-113
#wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1],as.numeric(x)[SpliceSites$isImmune[w]=='yes'])$p.value# p-value = 0.02


w=which(SpliceSites$Quality_MultiZ46Align)
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
par(mar=c(4,5,1,0.5))
tab=table(x,SpliceSites$WeakAltConst[w])
tab=tab[,c('weak','alternative','constitutive')]
colnames(tab)=c('weak','alt.','const.')
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

########### study of recent active splice sites
placental_active=SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$newAge<=12 & SpliceSites$WeakAltConst!='weak',]	#// V5.2 
placental_active=SpliceSites[SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site & SpliceSites$newAge<=12 & SpliceSites$WeakAltConst!='weak',]	#// V5.3 added non ZeroGerpSite
primate_active=SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$newAge<=9 & SpliceSites$WeakAltConst!='weak',]
ape_active=SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$newAge<=5 & SpliceSites$WeakAltConst!='weak',]

 

table(primate_active$WeakAltConst)
# alternative constitutive
#       21017         3858
table(ape_active$WeakAltConst)
# alternative constitutive
#        6905         1448

mean(primate_active$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding']) # 0.2946332 => 29%
mean(ape_active$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding']) # 0.2645756	=> 26%
mean(placental_active$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding']) # 0.2645756	=> 26%


mean(placental_active$coding) # 0.2520302 => 25% // V5.3 # 0.2549273
mean(primate_active$coding) # 0.1420302 => 14%
mean(ape_active$coding) # 0.1318089 => 13%

#// V5.3
table(placental_active$WeakAltConst)
# alternative constitutive 
#       44413        13580 
table(placental_active$WeakAltConst[placental_active$coding])
# alternative constitutive 
#        7224         7560 

mean(primate_active$coding[primate_active$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding']],na.rm=T) # 0.4678674
mean(ape_active$coding[ape_active$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding']],na.rm=T) # 0.4800905


##########  Placental specific constitutive coding splice site 
w=which(SpliceSites$WeakAltConst=='constitutive'  & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site) # excluding weak splicing
GeneHasPlacentalSpecificConstSite=By(SpliceSites$newAge[w]<=12,SpliceSites$gene[w],any)
mean(GeneHasPlacentalSpecificConstSite) # 31% of genes have a constitutive splice site (>1 read per sample) that is Placental specific
resGO_PlaConst=GOSeq(names(GeneHasPlacentalSpecificConstSite)[GeneHasPlacentalSpecificConstSite], names(GeneHasPlacentalSpecificConstSite),addCI=T,FDR=0.2)
resGO_PlaConst$type='constitutive'
resGO_PlaConst=resGO_PlaConst[,c(1:5,7:8,10:11,13:17,6,18,12)]

GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
OR_all=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificConstSite))
GenelFC=By(SpliceSites$log2FC_LPS[w],SpliceSites$gene[w],max)
OR_LPS=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificConstSite))
GenelFC=By(SpliceSites$log2FC_PAM3[w],SpliceSites$gene[w],max)
OR_PAM3=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificConstSite))
GenelFC=By(SpliceSites$log2FC_R848[w],SpliceSites$gene[w],max)
OR_R848=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificConstSite))
GenelFC=By(SpliceSites$log2FC_IAV[w],SpliceSites$gene[w],max)
OR_IAV=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificConstSite))
OR_LFC_Const=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='constitutive')

resGO_PlaConst[1:10,1:16]
#    category over_represented_pvalue under_represented_pvalue numDEInCat numInCat ontology          FDR FoldEnrich       Pvalue nbInGrp nbInBckgd       OR        CI  lowerCI                                term         type
#3052  GO:0006955            2.399177e-18                        1        401      911       BP 3.913538e-14   1.400522 2.399177e-18    3263     10382 1.815658 [1.6-2.1] 1.577324                     immune response constitutive
#3049  GO:0006952            1.756449e-16                        1        396      917       BP 1.432560e-12   1.374010 1.756449e-16    3263     10382 1.749022   [1.5-2] 1.519702                    defense response constitutive
#869   GO:0002376            9.342017e-15                        1        584     1473       BP 5.079566e-11   1.261462 9.342017e-15    3263     10382 1.527616 [1.4-1.7] 1.360678               immune system process constitutive
#4052  GO:0009897            4.571986e-14                        1         79      124       CC 1.864456e-10   2.027073 4.571986e-14    3263     10382 3.899599 [2.7-5.8] 2.664011    external side of plasma membrane constitutive
#14875 GO:0098552            1.490399e-13                        1        104      182       CC 4.500689e-10   1.818134 1.490399e-13    3263     10382 2.971478 [2.2-4.1] 2.187221                    side of membrane constitutive
#1012  GO:0002682            1.655476e-13                        1        324      752       BP 4.500689e-10   1.370854 1.655476e-13    3263     10382 1.723255   [1.5-2] 1.477332 regulation of immune system process constitutive
#1947  GO:0004872            1.165584e-12                        1        195      414       MF 2.716143e-09   1.498643 1.165584e-12    3263     10382 2.002469 [1.6-2.5] 1.634550                   receptor activity constitutive
#1340  GO:0003677            1.721413e-12                        1        515     1309       MF 3.272151e-09   1.251790 1.721413e-12    3263     10382 1.492862 [1.3-1.7] 1.321716                         DNA binding constitutive
#1339  GO:0003676            1.805380e-12                        1        894     2439       MF 3.272151e-09   1.166245 1.805380e-12    3263     10382 1.361490 [1.2-1.5] 1.235977                nucleic acid binding constitutive
#4072  GO:0009986            2.661481e-12                        1        165      339       CC 4.341408e-09   1.548632 2.661481e-12    3263     10382 2.125689 [1.7-2.7] 1.699952                        cell surface constitutive

OR_LFC_Const
#             LowerCI       OR  UpperCI alpha            P condition         type
#odds ratio  0.9782934 1.220423 1.517693  0.05 7.250283e-02       LPS constitutive
#odds ratio1 0.9661507 1.195399 1.474613  0.05 9.455147e-02  PAM3CSK4 constitutive
#odds ratio2 1.2238014 1.457279 1.732901  0.05 1.868187e-05      R848 constitutive
#odds ratio3 1.6198897 1.893117 2.211286  0.05 6.794728e-16       IAV constitutive
#odds ratio4 1.4009288 1.590112 1.803907  0.05 6.231187e-13       all constitutive

##########  Placental specific alternative coding splice site
w=which(SpliceSites$WeakAltConst=='alternative'  & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & !SpliceSites$Zero_Gerp_site) # excluding weak splicing
GeneHasPlacentalSpecificAltSite=By(SpliceSites$newAge[w]<=12,SpliceSites$gene[w],any)
mean(GeneHasPlacentalSpecificAltSite) # 42% of genes have an Alternative splice site (>1 read per sample) that is Placental specific
resGO_PlaAlt=GOSeq(names(GeneHasPlacentalSpecificAltSite)[GeneHasPlacentalSpecificAltSite], names(GeneHasPlacentalSpecificAltSite),addCI=T,FDR=0.2)
resGO_PlaAlt$type='alternative'
resGO_PlaAlt=resGO_PlaAlt[,c(1:5,7:8,10:11,13:17,6,18,12)]


resGO_PlaAlt[1:10,1:16]

#        category over_represented_pvalue under_represented_pvalue numDEInCat numInCat ontology          FDR FoldEnrich       Pvalue nbInGrp nbInBckgd       OR        CI  lowerCI                               term        type
#2862  GO:0006952            1.686404e-17                        1        397      698       BP 2.553553e-13   1.353375 1.686404e-17    3539      8421 1.922878 [1.6-2.3] 1.639640                   defense response alternative
#823   GO:0002376            2.461042e-15                        1        591     1131       BP 1.863255e-11   1.243392 2.461042e-15    3539      8421 1.611824 [1.4-1.8] 1.418685              immune system process alternative
#2865  GO:0006955            9.814761e-14                        1        380      690       BP 4.953837e-10   1.310441 9.814761e-14    3539      8421 1.773926 [1.5-2.1] 1.512228                    immune response alternative
#1263  GO:0003677            3.618782e-11                        1        555     1114       MF 1.369890e-07   1.185471 3.618782e-11    3539      8421 1.438264 [1.3-1.6] 1.264917                        DNA binding alternative
#9509  GO:0045087            2.361412e-10                        1        253      456       BP 7.045844e-07   1.320197 2.361412e-10    3539      8421 1.774467 [1.5-2.2] 1.461272             innate immune response alternative
#13874 GO:0098542            2.791907e-10                        1        108      163       BP 7.045844e-07   1.576592 2.791907e-10    3539      8421 2.762270   [2-3.9] 1.972089 defense response to other organism alternative
#10852 GO:0050794            4.466251e-09                        1       1839     4163       BP 9.661139e-06   1.051135 4.466251e-09    3539      8421 1.190656 [1.1-1.3] 1.090803     regulation of cellular process alternative
#1262  GO:0003676            5.707229e-09                        1        901     1926       MF 1.080236e-05   1.113145 5.707229e-09    3539      8421 1.285159 [1.2-1.4] 1.158604               nucleic acid binding alternative
#10848 GO:0050789            6.636866e-09                        1       1922     4365       BP 1.116616e-05   1.047737 6.636866e-09    3539      8421 1.186646 [1.1-1.3] 1.087048   regulation of biological process alternative
#12281 GO:0065007            3.400432e-08                        1       1996     4562       BP 5.148934e-05   1.041090 3.400432e-08    3539      8421 1.167528 [1.1-1.3] 1.069255              biological regulation alternative

GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
OR_all=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificAltSite))
GenelFC=By(SpliceSites$log2FC_LPS[w],SpliceSites$gene[w],max)
OR_LPS=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificAltSite))
GenelFC=By(SpliceSites$log2FC_PAM3[w],SpliceSites$gene[w],max)
OR_PAM3=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificAltSite))
GenelFC=By(SpliceSites$log2FC_R848[w],SpliceSites$gene[w],max)
OR_R848=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificAltSite))
GenelFC=By(SpliceSites$log2FC_IAV[w],SpliceSites$gene[w],max)
OR_IAV=odds.ratio(table(GenelFC>1,GeneHasPlacentalSpecificAltSite))
OR_LFC_alt=cbind(rbind(OR_LPS,OR_PAM3,OR_R848,OR_IAV,OR_all),condition=c(condIndex[-1],'all'),type='alternative')

write.table(rbind(OR_LFC_Const,OR_LFC_alt),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/new/TableS1B_new_V5.3_OR_LFC.txt',HOME),quote=F,sep='\t',row.names=FALSE)
write.table(rbind(resGO_PlaConst,resGO_PlaAlt),file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/new/TableS1A_new_V5.3_GO_placentals.txt',HOME),quote=F,sep='\t',row.names=FALSE)

w=which(SpliceSites$WeakAltConst!='weak' & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align)
SpliceSites[w,'']

#nrow(SpliceSites[SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=5,])
##[1] 207387
#nrow(SpliceSites[SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=5 & SpliceSites$newAge<=12,])
##[1] 14059
#nrow(SpliceSites[SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=5 & SpliceSites$newAge<=12 & SpliceSites$NbMatchingSeqPrim==1,])
## 65

nrow(SpliceSites[SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=5 & !SpliceSites$Zero_Gerp_site,])
#[1] 207351
nrow(SpliceSites[SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=5 & SpliceSites$newAge<=12 & !SpliceSites$Zero_Gerp_site,])
#[1] 14024
nrow(SpliceSites[SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=5 & SpliceSites$newAge<=12 & SpliceSites$NbMatchingSeqPrim==1 & !SpliceSites$Zero_Gerp_site,])
# 55


reverseSeq=function(seq){paste(sapply(strsplit(seq,':'),function(x){nn=nchar(x);paste(substr(x,1,nn-3),reverse(substr(x,nn-1,nn)),sep='_')}),collapse=':')}

#plot(tree,type="fan")
strand='+'
reverse=function(Seq){
	Corresp=c('A','G','T','C','-')
	names(Corresp)=c('T','C','A','G','-')
	sapply(strsplit(Seq,''),function(x){paste(rev(Corresp[x]),collapse='')})
	}
	
HumanSpecificSpliceSites=SpliceSites[SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak' & SpliceSites$NbAlignedPrimates>=5 & SpliceSites$newAge<=12 & SpliceSites$NbMatchingSeqPrim==1 & !SpliceSites$Zero_Gerp_site,]
# Chromosome	Splice site start	Splice site end	strand	type	Sequence of splice site	Gene associated to the splice site (Ensembl ID)	Gene associated to the splice site (symbol)	Present in over 95% of observed transcripts	Total read count NS	Gene FPKM NS	Total read count LPS	Gene FPKM LPS	Total read count Pam3Csk4	Gene FPKM Pam3Csk4	Total read count R848	Gene FPKM R848	Total read count IAV	Gene FPKM IAV	Mean GerpRS of splice site 	Fist apparition of the current splice site sequence in the phylogeny	Nb of primates with an alignment	Splice site sequence across common ancestors	Splice site sequence across 46 vertebrates
HumanSpecificSpliceSites$seq=ifelse(HumanSpecificSpliceSites$type=='donor',HumanSpecificSpliceSites$donor,HumanSpecificSpliceSites$acceptor)
HumanSpecificSpliceSites$allSequence[HumanSpecificSpliceSites$strand_Intropolis=='-']=sapply(HumanSpecificSpliceSites$allSequence[HumanSpecificSpliceSites$strand_Intropolis=='-'],reverseSeq)
HumanSpecificSpliceSites$ancestralSequence[HumanSpecificSpliceSites$strand_Intropolis=='-']=sapply(HumanSpecificSpliceSites$ancestralSequence[HumanSpecificSpliceSites$strand_Intropolis=='-'],reverseSeq)

HumanSpecificSpliceSites$ancestralSequence=gsub('treeShrew','treeshrew',HumanSpecificSpliceSites$ancestralSequence)
HumanSpecificSpliceSites$ancestralSequence=sapply(strsplit(HumanSpecificSpliceSites$ancestralSequence,':'),function(x){paste(toupper(substr(unlist(x),1,1)),substr(x,2,10000),sep='',collapse=':')})
HumanSpecificSpliceSites$ancestralSequence=gsub(':',' / ',HumanSpecificSpliceSites$ancestralSequence)
HumanSpecificSpliceSites$ancestralSequence=gsub('_',': ',HumanSpecificSpliceSites$ancestralSequence)

HumanSpecificSpliceSites$allSequence=gsub('X_tropicalis','X. tropicalis',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Rock_hyrax','Rock hyrax',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Kangaroo_rat','Kangaroo rat',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('TreeShrew','Treeshrew',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Zebra_finch','Zebra finch',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('lemur','Lemur',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('Guinea_Pig','Guinea pig',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub(':',' / ',HumanSpecificSpliceSites$allSequence)
HumanSpecificSpliceSites$allSequence=gsub('_',': ',HumanSpecificSpliceSites$allSequence)


HumanSpecificSpliceSites=HumanSpecificSpliceSites[,c("chrom","start","end","strand_Intropolis","seq","type","gene","symbol","isImmune",'read_persample_NS',"FPKM_NS",'read_persample_LPS',"FPKM_LPS",'read_persample_PAM3CSK4',"FPKM_PAM3CSK4",'read_persample_R848',"FPKM_R848",'read_persample_IAV',"FPKM_IAV","Pct_of_reads_overlapping_MajorJunction_thatMap_to_splice_site","WeakAltConst",'Gerp_site','NbAlignedPrimates',"allSequence","ancestralSequence")]
write.table(HumanSpecificSpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/tables/new/TableS1C_V5.3_HumanSpecificSpliceSites',HOME),quote=F,sep='\t',row.names=FALSE)

##############################################################
##					END HERE								##
##############################################################

##########  Primate specific active coding splice site
w=which(SpliceSites$WeakAltConst!='weak'  & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align) # excluding weak splicing
GeneHasPrimateSpecificSite=By(SpliceSites$newAge[w]<=9,SpliceSites$gene[w],any)
GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
odds.ratio(table(GenelFC>1,GeneHasPrimateSpecificSite))
#           LowerCI       OR UpperCI alpha            P
#odds ratio 1.322276 1.529977 1.76686  0.05 1.211268e-08

mean(GeneHasPrimateSpecificSite) # 16% of genes have an active splice site (>1 read per sample) that is Human specific
resGO=GOSeq(names(GeneHasPrimateSpecificSite)[GeneHasPrimateSpecificSite], names(GeneHasPrimateSpecificSite),addCI=T,FDR=0.2)
#       category over_represented_pvalue under_represented_pvalue numDEInCat numInCat                                                                       term ontology          FDR
#9608  GO:0043169            8.746142e-09                1.0000000        473     2465                                                             cation binding       MF 9.360146e-05
#11017 GO:0046872            1.133395e-08                1.0000000        468     2440                                                          metal ion binding       MF 9.360146e-05
#3538  GO:0008270            3.078667e-07                1.0000000        172      784                                                           zinc ion binding       MF 1.695012e-03
#7185  GO:0032196            3.467086e-06                0.9999999          8        9                                                              transposition       BP 1.431646e-02
#9606  GO:0043167            1.249964e-05                0.9999904        654     3695                                                                ion binding       MF 3.714443e-02
#10262 GO:0045087            2.148242e-05                0.9999867        136      622                                                     innate immune response       BP 3.714443e-02
#8119  GO:0034340            2.253558e-05                0.9999936         23       60                                              response to type I interferon       BP 3.714443e-02
#12756 GO:0060337            2.253558e-05                0.9999936         23       60                                        type I interferon signaling pathway       BP 3.714443e-02
#13883 GO:0071357            2.253558e-05                0.9999936         23       60                                     cellular response to type I interferon       BP 3.714443e-02
#3089  GO:0006952            2.395921e-05                0.9999839        200      974                                                           defense response       BP 3.714443e-02
#8120  GO:0034341            2.473747e-05                0.9999911         32       97                                               response to interferon-gamma       BP 3.714443e-02
#11043 GO:0046914            4.882295e-05                0.9999667        185      909                                               transition metal ion binding       MF 6.720072e-02
#3217  GO:0007175            9.618586e-05                1.0000000          5        5 negative regulation of epidermal growth factor-activated receptor activity       BP 1.141739e-01
#4235  GO:0010528            1.096440e-04                0.9999971          6        7                                                regulation of transposition       BP 1.141739e-01
#4236  GO:0010529            1.096440e-04                0.9999971          6        7                                       negative regulation of transposition       BP 1.141739e-01
#3092  GO:0006955            1.106002e-04                0.9999221        195      969                                                            immune response       BP 1.141739e-01
#4030  GO:0009615            1.698812e-04                0.9999088         60      241                                                          response to virus       BP 1.650546e-01
#                                                                            Term FoldEnrich       Pvalue nbInGrp nbInBckgd        OR           CI  lowerCI
#9608                                                              cation binding   1.185817 8.746142e-09    1791     11068  1.312426    [1.2-1.5] 1.165812
#11017                                                          metal ion binding   1.185304 1.133395e-08    1791     11068  1.310352    [1.2-1.5] 1.163527
#3538                                                            zinc ion binding   1.355770 3.078667e-07    1791     11068  1.504112    [1.3-1.8] 1.252180
#7185                                                               transposition   5.493145 3.467086e-06    1791     11068 41.614064 [5.6-1828.2] 5.571582
#9606                                                                 ion binding   1.093797 1.249964e-05    1791     11068  1.179494    [1.1-1.3] 1.059587
#10262                                                     innate immune response   1.351208 2.148242e-05    1791     11068  1.486363    [1.2-1.8] 1.211026
#8119                                               response to type I interferon   2.368919 2.253558e-05    1791     11068  3.248265    [1.8-5.6] 1.837810
#12756                                        type I interferon signaling pathway   2.368919 2.253558e-05    1791     11068  3.248265    [1.8-5.6] 1.837810
#13883                                     cellular response to type I interferon   2.368919 2.253558e-05    1791     11068  3.248265    [1.8-5.6] 1.837810
#3089                                                            defense response   1.268950 2.395921e-05    1791     11068  1.380947    [1.2-1.6] 1.165484
#8120                                                response to interferon-gamma   2.038693 2.473747e-05    1791     11068  2.577949      [1.6-4] 1.628109
#11043                                               transition metal ion binding   1.257713 4.882295e-05    1791     11068  1.360795    [1.1-1.6] 1.141586
#3217  negative regulation of epidermal growth factor-activated receptor activity   6.179788 9.618586e-05    1791     11068       Inf    [4.8-Inf] 4.753680
#4235                                                 regulation of transposition   5.296961 1.096440e-04    1791     11068 31.140859 [3.8-1422.8] 3.777691
#4236                                        negative regulation of transposition   5.296961 1.096440e-04    1791     11068 31.140859 [3.8-1422.8] 3.777691
#3092                                                             immune response   1.243611 1.106002e-04    1791     11068  1.342209    [1.1-1.6] 1.131010
#4030                                                           response to virus   1.538536 1.698812e-04    1791     11068  1.741805    [1.3-2.4] 1.272514
##########  Primate specific constitutive coding splice site 
w=which(SpliceSites$WeakAltConst=='constitutive'  & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align) # excluding weak splicing
GeneHasPrimateSpecificConstSite=By(SpliceSites$newAge[w]<=9,SpliceSites$gene[w],any)
mean(GeneHasPrimateSpecificConstSite) # 5.8% of genes have a constitutive splice site (>1 read per sample) that is Primate specific
resGO_PrimConst=GOSeq(names(GeneHasPrimateSpecificConstSite)[GeneHasPrimateSpecificConstSite], names(GeneHasPrimateSpecificConstSite),addCI=T,FDR=0.2)
GenelFC=By(SpliceSites$maxlFC[w],SpliceSites$gene[w],max)
odds.ratio(table(GenelFC>1,GeneHasPrimateSpecificConstSite))
           LowerCI       OR  UpperCI alpha            P
odds ratio 1.336126 1.683441 2.106468  0.05 9.598833e-06

#  category  over_represented_pvalue under_represented_pvalue numDEInCat numInCat                            term ontology          FDR                            Term FoldEnrich       Pvalue nbInGrp nbInBckgd             OR        CI  lowerCI
# GO:0003676            1.790855e-11                1.0000000        200     2439            nucleic acid binding       MF 2.921243e-07            nucleic acid binding   1.390135 1.790855e-11     613     10392		1.630701 [1.4-1.9] 1.361321
# GO:0003677            8.028764e-11                1.0000000        124     1309                     DNA binding       MF 6.548260e-07                     DNA binding   1.605908 8.028764e-11     613     10392		1.838824 [1.5-2.3] 1.483722
# GO:0043169            9.404602e-10                1.0000000        186     2322                  cation binding       MF 5.113595e-06                  cation binding   1.357968 9.404602e-10     613     10392		1.558572 [1.3-1.9] 1.296192
# GO:0046872            2.224324e-09                1.0000000        183     2298               metal ion binding       MF 9.070793e-06               metal ion binding   1.350019 2.224324e-09     613     10392		1.542082 [1.3-1.9] 1.281307
# GO:0097159            7.781534e-06                0.9999952        245     3605 organic cyclic compound binding       MF 2.538647e-02 organic cyclic compound binding   1.152125 7.781534e-06     613     10392		1.271820 [1.1-1.5] 1.071762
# GO:1901363            1.759020e-05                0.9999885        242     3585   heterocyclic compound binding       MF 4.782188e-02   heterocyclic compound binding   1.144366 1.759020e-05     613     10392		1.255746 [1.1-1.5] 1.057735

##########  Primate specific alternative coding splice site
w=which(SpliceSites$WeakAltConst=='alternative'  & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align) # excluding weak splicing
GeneHasPrimateSpecificAltSite=By(SpliceSites$newAge[w]<=9,SpliceSites$gene[w],any)
mean(GeneHasPrimateSpecificAltSite) # 16% of genes have an Alternative splice site (>1 read per sample) that is Primate specific
resGO=GOSeq(names(GeneHasPrimateSpecificAltSite)[GeneHasPrimateSpecificAltSite], names(GeneHasPrimateSpecificAltSite),addCI=T,FDR=0.2)

#      category over_represented_pvalue under_represented_pvalue numDEInCat numInCat                term ontology        FDR                   Term FoldEnrich       Pvalue nbInGrp nbInBckgd       OR        CI      lowerCI
# GO:0006952            1.027381e-07                1.0000000        162      699       defense response       BP 0.00155566       defense response   1.410456 1.027381e-07    1386      8435 1.604889 [1.3-1.9]		1.323911
# GO:0008270            3.536121e-06                0.9999982        142      641       zinc ion binding       MF 0.02064124       zinc ion binding   1.348193 3.536121e-06    1386      8435 1.498252 [1.2-1.8]	 	1.222509
# GO:0006955            4.089534e-06                0.9999977        155      691        immune response       BP 0.02064124        immune response   1.365135 4.089534e-06    1386      8435 1.529906 [1.3-1.9]		1.258118
# GO:0045087            7.617110e-06                0.9999958        108      457 innate immune response       BP 0.02883457 innate immune response   1.438234 7.617110e-06    1386      8435 1.622230   [1.3-2]		1.284267

w=which(SpliceSites$WeakAltConst!='weak'  & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align) # excluding weak splicing
GeneHasPrimateSpecificSite=By(SpliceSites$newAge[w]<=9,SpliceSites$gene[w],any)
mean(GeneHasPrimateSpecificSite) # 16% of genes have an active splice site (>1 read per sample) that is Human specific
resGO=GOSeq(names(GeneHasPrimateSpecificSite)[GeneHasPrimateSpecificSite], names(GeneHasPrimateSpecificSite),addCI=T)

w=which(SpliceSites$WeakAltConst!='weak' & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align) # excluding weak splicing
GeneHasApeSpecificSite=By(SpliceSites$newAge[w]<=5,SpliceSites$gene[w],any)
mean(GeneHasApeSpecificSite) # 5.7% of genes have an active splice site (>1 read per sample) that is Human specific
resGO=GOSeq(names(GeneHasApeSpecificSite)[GeneHasApeSpecificSite], names(GeneHasApeSpecificSite),addCI=T)
# only enrichment: transition metal/ zinc ion binding

########### study of consitutive non conserved sites
constit_ncons=SpliceSites[!SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site<2 & SpliceSites$WeakAltConst=='constitutive',]

FullGeneAnnot=fread(sprintf('%s/Annotation/GeneAnnotation_hg37_ens70.txt',HOME),drop=1)
constit_ncons$gene[is.na(constit_ncons$gene)]=''
mean(constit_ncons$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding'])
mean(constit_ncons$coding[constit_ncons$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding']])

table(SpliceSites[!SpliceSites$Zero_Gerp_site & SpliceSites$Gerp_site<2 & SpliceSites$inEnsembl ,'coding'])

constit_ncons[constit_ncons$newAge==1 & constit_ncons$coding & constit_ncons$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding'],]
SpliceSites[SpliceSites$newAge==1 & !SpliceSites$Zero_Gerp_site & SpliceSites$WeakAltConst=='constitutive' & SpliceSites$coding & SpliceSites$gene%in%FullGeneAnnot$Ensembl.Gene.ID[FullGeneAnnot$Gene.Biotype=='protein_coding'],]
w=which(!isInactive_allCond & !SpliceSites$Zero_Gerp_site) # including alternative splicing


#isActive_anyCond=apply(SS_Activity>10,1,any)
#isActive_allCond=apply(SS_Activity>10,1,all)
#mean(isActive_anyCond[!isInactive_allCond]) # 55.8%

#isCryptic=SpliceSites$Total_count<
#isCryptic_all=apply(SpliceSites[,2:6]<30,1, all)
#mean(isCryptic) # 0.7089847 -> 70.9% of all splice sites

#sum(isActive_anyCond) # 155231
#sum(isActive_allCond) # 89290
#
#mean(isActive_anyCond) # 7.8 % of all splice sites
#mean(isActive_allCond) # 4.5% of all splice sites

sum(!is.na(Junc_annot$Nbread))/1e6 # 2.099531 ->2.1 million

#mean(SpliceSites$Gerp_site[isCryptic & SpliceSites$Gerp_site!=0]>2) # 0.2088718
mean(SpliceSites$Gerp_site[isInactive_allCond & !SpliceSites$Zero_Gerp_site]>2) # 0.9405729
conserved=SpliceSites$Gerp_site>2

summary(glm(conserved[!SpliceSites$Zero_Gerp_site]~ isInactive_allCond[!SpliceSites$Zero_Gerp_site] + isActive_anyCond[!SpliceSites$Zero_Gerp_site] + isActive_allCond[!SpliceSites$Zero_Gerp_site]+ isAlternative + coding +isImmune+I(maxlFC>1),data=SpliceSites[!SpliceSites$Zero_Gerp_site,], family=binomial))

                                                    Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                          0.99564    0.02013   49.46   <2e-16 ***
#isInactive_allCond[!SpliceSites$Zero_Gerp_site]TRUE -0.86782    0.01587  -54.70   <2e-16 ***
#isActive_anyCond[!SpliceSites$Zero_Gerp_site]TRUE    0.73842    0.02344   31.50   <2e-16 ***
#isActive_allCond[!SpliceSites$Zero_Gerp_site]TRUE    0.46464    0.03092   15.03   <2e-16 ***
#isAlternativeyes                                    -1.68143    0.01825  -92.12   <2e-16 ***
#codingTRUE                                           2.43223    0.01393  174.61   <2e-16 ***
#isImmuneyes                                         -0.27104    0.02551  -10.63   <2e-16 ***
#I(maxlFC > 1)TRUE                                   -0.26787    0.01918  -13.96   <2e-16 ***

#mean(SpliceSites$Gerp_site[isActive_anyCond & !SpliceSites$Zero_Gerp_site]>2) # 0.9405729
#mean(SpliceSites$Gerp_site[isActive_allCond & !SpliceSites$Zero_Gerp_site]>2) # 0.970662
#mean(SpliceSites$Gerp_site[isActive_allCond & !SpliceSites$Zero_Gerp_site & SpliceSites$constitutive]>2) # 0.9405729

w=which(SpliceSites$WeakAltConst=='constitutive' & !SpliceSites$Zero_Gerp_site & SpliceSites$coding) # including alternative splicing
GeneHasNonConservedSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedSite) # 56% of genes have an active splice site (>1 read per sample) that is not conserved (including alternative sites)
Gene_FPKM_over10=names(GeneHasNonConservedSite)%in%SpliceSites$gene[SpliceSites$maxFPKM>=10]
resGO=GOSeq(names(GeneHasNonConservedSite)[GeneHasNonConservedSite], names(GeneHasNonConservedSite),addCI=T)
resGO[order(-resGO$lowerCI),-12]
# strongest FE/OR/lowerCI: negative regulation of viral replication
resGO=GOSeq(names(GeneHasNonConservedSite)[GeneHasNonConservedSite & Gene_FPKM_over10], names(GeneHasNonConservedSite)[Gene_FPKM_over10],addCI=T)
resGO[order(-resGO$lowerCI),-12]
# not much to see (plasma membrane...)


mean(rk_time$time_MA[SpliceSites[SpliceSites$Quality_MultiZ46Align & SpliceSites$WeakAltConst!='weak',]$newAge])


# enrichment of Human specific active sites
w=which(!isInactive_allCond & SpliceSites$Quality_MultiZ46Align) # excluding weak splicing
GeneHasHumanSpecificSite=By(SpliceSites$newAge[w]==1,SpliceSites$gene[w],any)
mean(GeneHasHumanSpecificSite) # 3% of genes have an active splice site (>1 read per sample) that is Human specific
resGO=GOSeq(names(GeneHasHumanSpecificSite)[GeneHasHumanSpecificSite], names(GeneHasHumanSpecificSite),addCI=T)
# no Enrichment

# enrichment of primate specific active sites
w=which(!isInactive_allCond & SpliceSites$Quality_MultiZ46Align) # excluding weak splicing
GeneHasPrimateSpecificSite=By(SpliceSites$newAge[w]<=9,SpliceSites$gene[w],any)
mean(GeneHasPrimateSpecificSite) # 34% of genes have an active splice site (>1 read per sample) that is Human specific
resGO=GOSeq(names(GeneHasPrimateSpecificSite)[GeneHasPrimateSpecificSite], names(GeneHasPrimateSpecificSite),addCI=T)
# strongest FE/OR/lowerCI: defense response to other organism



GeneHasNonConservedSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
Gene_FPKM_over10=names(GeneHasNonConservedSite)%in%SpliceSites$gene[SpliceSites$maxFPKM>=10]
mean(GeneHasNonConservedSite) # 56% of genes have an active splice site (>1 read per sample),
resGO=GOSeq(names(GeneHasNonConservedSite)[GeneHasNonConservedSite], names(GeneHasNonConservedSite),addCI=T)
# strongest FE: positive regulation of type I interferon production
resGO=GOSeq(names(GeneHasNonConservedSite)[GeneHasNonConservedSite & Gene_FPKM_over10], names(GeneHasNonConservedSite)[Gene_FPKM_over10],addCI=T)
resGO[order(-resGO$lowerCI),-12]
# not much to see (plasma membrane...)

w=which(!isInactive_allCond & SpliceSites$Gerp_site!=0 & SpliceSites$isAlternative=='yes') # excluding alternative splicing
mean(SpliceSites$Gerp_site[w]>2)
GeneHasNonConservedSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedSite) # 56% of genes have an active splice site (>1 read per sample),
resGO=GOSeq(names(GeneHasNonConservedSite)[GeneHasNonConservedSite], names(GeneHasNonConservedSite),addCI=T)
resGO[order(-resGO$lowerCI),-12]
# strongest FE: positive regulation of type I interferon production


w=which(!isInactive_allCond_butFlu & !SpliceSites$Zero_Gerp_site) # excluding alternative splicing
GeneHasNonConservedSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedSite) # 56% of genes have an active splice site (>1 read per sample),
resGO=GOSeq(names(GeneHasNonConservedSite)[GeneHasNonConservedSite], names(GeneHasNonConservedSite),addCI=T)
resGO[order(-resGO$lowerCI),-12]


# Genes with old non Conserved splicing 
w=which(SpliceSites$WeakAltConst=='constitutive' & !SpliceSites$Zero_Gerp_site & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align) # excluding alternative splicing
mean(SpliceSites$Gerp_site[w]>2)
GeneHasOldNonConservedSite=By(SpliceSites$newAge[w]<=9,SpliceSites$gene[w],any)
mean(GeneHasOldNonConservedSite) # 66% of genes have an old non conserved consitutive splice site (>1 read per sample),
resGO=GOSeq(names(GeneHasOldNonConservedSite)[GeneHasOldNonConservedSite], names(GeneHasOldNonConservedSite),addCI=T)
resGO[order(-resGO$lowerCI),-12]

# strongest FE: positive regulation of type I interferon production



SpliceSites=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4A_SpliceSites_V2.txt',HOME))
SpliceSites$isStart=ifelse(gsub(' (donor|acceptor)','',SpliceSites$site_id)==SpliceSites$start_id,TRUE,FALSE)
SpliceSites$isEnd=ifelse(gsub(' (donor|acceptor)','',SpliceSites$site_id)==SpliceSites$end_id,TRUE,FALSE)
mean(SpliceSites$isStart & SpliceSites$isEnd)
mean(!SpliceSites$isStart & !SpliceSites$isEnd)
mean(SpliceSites$isStart | SpliceSites$isEnd)
SpliceSites$Zero_Gerp_site=FALSE
SpliceSites$Zero_Gerp_site[SpliceSites$isStart & (SpliceSites$GerpRS_start==0 | SpliceSites$GerpRS_start2==0)]=TRUE
SpliceSites$Zero_Gerp_site[SpliceSites$isEnd & (SpliceSites$GerpRS_end==0 | SpliceSites$GerpRS_end2==0)]=TRUE



plotGO=function(resGO,mar=18,...){
	splitname=function(x,sep=' ',nmax=40){y=strsplit(x,sep)
										y=sapply(y,function(z){
														countchar=cumsum(nchar(z));
														group=countchar%/%nmax;
														paste(By(z,group,paste,collapse=' '),collapse='\n')
														})
										x[nchar(x)>(nmax+10)]=y[nchar(x)>(nmax+10)]
										x
			}														
	params=par()
	par(mar=c(4,mar,1,1))
	resGO=resGO[nrow(resGO):1,]
	barplot(-log10(resGO$FDR), main="", horiz=TRUE, names.arg=splitname(resGO$Term),col=ifelse(resGO$ontology=='BP',colPSI[1],ifelse(resGO$ontology=='MF',colPSI[2],colPSI[3])),xlab=expression(-log[10](P[adj])),las=1,...)
	par(mar=params$mar)
	legend("bottomright",fill=colPSI[1:3],legend=c('BP','MF','CC'),bty='n')
	}
pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/SpliceTime_vs_GerpRS.pdf',HOME))
plotGO(TableS4B,mar=20)
plotGO(TableS4B_coding,mar=20)
dev.off()

w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)
GeneHasNonConservedSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
mean(GeneHasNonConservedSite);length(GeneHasNonConservedSite)
TableS4B_highCov=GOSeq(names(GeneHasNonConservedSite)[GeneHasNonConservedSite], names(GeneHasNonConservedSite),addCI=T)
write.table(TableS4B_highCov,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4B_EnrichmentNonConserved.txt',HOME),quote=F,row.names=F,sep='\t')
plotGO(TableS4B_highCov[1:15,],mar=20)


##############################################
# 	FINAL VERSION FOR MANUSCRIPT			##
# constitutive and alternative # Fig4C		##
##############################################
w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align)
mean(SpliceSites$Gerp_site[w]>2) #98.5
w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)
GeneHasNonConservedSite=By(SpliceSites$Gerp_site[w]<2,SpliceSites$gene[w],any)
GeneAge=By(SpliceSites$newAge[w],SpliceSites$gene[w],max) #16%
mean(GeneHasNonConservedSite);length(GeneHasNonConservedSite)

lFC=log2(1+MeanExpr[,-1])-log2(1+MeanExpr[,1])%o%c(1,1,1,1)
colnames(lFC)=paste('log2FC',condIndex[-1],sep='_')
maxlFC=apply(lFC[match(names(GeneHasNonConservedSite),rn(lFC)),],1,max)
odds.ratio(table(GeneHasNonConservedSite,maxlFC>1)) 
#           LowerCI       OR  UpperCI alpha            P
#odds ratio 1.987867 2.369055 2.818087  0.05 2.064338e-21
table(GeneAge,maxlFC>1)
resamp_nbNonConserved=sapply(1:10000,function(i){toSample=By(maxlFC>1,GeneAge,sum)
	samp=unlist(lapply(names(toSample),function(i){sample(which(GeneAge==i),toSample[i])}))
	sum(GeneHasNonConservedSite[samp])})
mean(resamp_nbNonConserved>sum(GeneHasNonConservedSite[maxlFC>1])) # P<10-4

table(GeneAge,maxlFC>1)

resamp_nbNonConserved=sapply(1:10000,function(i){toSample=By(names(GeneHasNonConservedSite)%in%,GeneAge,sum)
	samp=unlist(lapply(names(toSample),function(i){sample(which(GeneAge==i),toSample[i])}))
	sum(GeneHasNonConservedSite[samp])})
mean(resamp_nbNonConserved>sum(GeneHasNonConservedSite[maxlFC>1])) # P<10-4

hist(resamp_nbNonConserved,col='grey');points(sum(GeneHasNonConservedSite[maxlFC>1]),pch=16,col'=red')

TableS4B_highCov_coding=GOSeq(names(GeneHasNonConservedSite)[GeneHasNonConservedSite], names(GeneHasNonConservedSite),addCI=T)
write.table(TableS4B_highCov_coding,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4B_EnrichmentNonConserved_coding.txt',HOME),quote=F,row.names=F,sep='\t')


SpliceSites$isAlternative=''

TableS4C=SpliceSites[highCovAnyCond_junc & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10 & SpliceSites$Gerp_site<2,c('site_id','chrom','start','end','strand_Intropolis','type',"donor",'acceptor','gene','symbol',"Total_count_NS","Total_count_LPS","Total_count_PAM3CSK4","Total_count_R848","Total_count_IAV",'Gerp_site','Time_FirstAppeared_MA','newAge',"FPKM_NS","FPKM_LPS","FPKM_PAM3CSK4","FPKM_R848","FPKM_IAV",'isAlternative','isImmune')]
TableS4C$Seq=ifelse(TableS4C$type=='donor',TableS4C$donor,TableS4C$acceptor)
TableS4C$FistAppearedInNewAge=rk_time$Species[TableS4C$newAge]
TableS4C=TableS4C[,c('site_id','chrom','start','end','strand_Intropolis','type','Seq','gene','symbol','isImmune',"Total_count_NS","FPKM_NS","Total_count_LPS","FPKM_LPS","Total_count_PAM3CSK4","FPKM_PAM3CSK4","Total_count_R848","FPKM_R848","Total_count_IAV","FPKM_IAV",'Gerp_site','isAlternative','FistAppearedInNewAge','Time_FirstAppeared_MA')]
write.table(TableS4C,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4C_NonConserved_coding_HighCoverage.txt',HOME),quote=F,row.names=F,sep='\t')

w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10 & SpliceSites$NbAlignedPrimates-1>=5 & SpliceSites$coding & !SpliceSites$Zero_Gerp_site)
HumanSpecificSpliceSites=SpliceSites[w[which(SpliceSites$NbAlignedPrimates[w]-1>=5 & SpliceSites$FirstAppeared[w]=='Humans' & SpliceSites$NbMatchingSeqPrim[w]==1)],c('site_id','chrom','start','end','strand_Intropolis','type',"donor",'acceptor','gene','symbol','coding',"Total_count_NS","Total_count_LPS","Total_count_PAM3CSK4","Total_count_R848","Total_count_IAV",'Gerp_site','Time_FirstAppeared_MA','newAge',"FPKM_NS","FPKM_LPS","FPKM_PAM3CSK4","FPKM_R848","FPKM_IAV",'NbAlignedPrimates','allSequence','ancestralSequence','isImmune','isAlternative'),]
HumanSpecificSpliceSites$Seq=ifelse(HumanSpecificSpliceSites$type=='donor',HumanSpecificSpliceSites$donor,HumanSpecificSpliceSites$acceptor)
HumanSpecificSpliceSites$FistAppearedInNewAge=rk_time$Species[HumanSpecificSpliceSites$newAge]
HumanSpecificSpliceSites$coding=ifelse(HumanSpecificSpliceSites$coding,'yes','')
HumanSpecificSpliceSites$allSequence[HumanSpecificSpliceSites$strand=='-']=sapply(HumanSpecificSpliceSites$allSequence[HumanSpecificSpliceSites$strand=='-'],reverseSeq)
HumanSpecificSpliceSites$ancestralSequence[HumanSpecificSpliceSites$strand=='-']=sapply(HumanSpecificSpliceSites$ancestralSequence[HumanSpecificSpliceSites$strand=='-'],reverseSeq)
HumanSpecificSpliceSites=HumanSpecificSpliceSites[HumanSpecificSpliceSites$coding=='yes' & HumanSpecificSpliceSites$FistAppearedInNewAge=='Humans',c('site_id','chrom','start','end','strand_Intropolis','type','Seq','gene','symbol','isImmune',"Total_count_NS","FPKM_NS","Total_count_LPS","FPKM_LPS","Total_count_PAM3CSK4","FPKM_PAM3CSK4","Total_count_R848","FPKM_R848","Total_count_IAV","FPKM_IAV",'Gerp_site','isAlternative','FistAppearedInNewAge','NbAlignedPrimates','ancestralSequence','allSequence')]

reverseSeq=function(seq){paste(sapply(strsplit(seq,':'),function(x){nn=nchar(x);paste(substr(x,1,nn-3),reverse(substr(x,nn-1,nn)),sep='_')}),collapse=':')}

#plot(tree,type="fan")
strand='+'
reverse=function(Seq){
	Corresp=c('A','G','T','C','-')
	names(Corresp)=c('T','C','A','G','-')
	sapply(strsplit(Seq,''),function(x){paste(rev(Corresp[x]),collapse='')})
	}
tree$tip.label=sapply(leafLabels,function(x){x1=substr(x,1,nchar(x)-3);x2=substr(x,nchar(x)-1,nchar(x));paste(gsub('-$','- ',gsub('^-',' -',ifelse(strand=='+',x2,reverse(x2)))),x1,sep=' : ')})

overlap_snp=mapply(function(chrom,start,end){x=getMapInfo_locus(chrom,start,end);if(nrow(x)==0){c('','')}else{c(x$snp.name,x$SNPfreq)}},HumanSpecificSpliceSites$chrom,HumanSpecificSpliceSites$start,HumanSpecificSpliceSites$end)
write.table(HumanSpecificSpliceSites,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Table_S4D_HumanSpecificSpliceSites.txt',HOME),quote=F,sep='\t',row.names=F)


w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)
boxplot(SpliceSites$Gerp_site[w]~ifelse(SpliceSites$isImmune[w]=='yes','immune',''),notch=T,col=mycols[c(3,7)],axes=T,pch=16,cex=0.3,ylab='GerpRS',ylim=c(-2,6.5),las=1)
wilcox.test(SpliceSites$Gerp_site[w]~SpliceSites$isImmune[w]) # P<10E-16
#x0=cut(SpliceSites$rank_FirstAppeared,c(-1,10,11,12,13,14,15,16,17,18))
w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)
mean(SpliceSites$Time_FirstAppeared_MA[w]) # 157.395
mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']],na.rm=T)
mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='']],na.rm=T)
mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)
GeneAge=By(SpliceSites$newAge[w],SpliceSites$gene[w],max)
SpliceSites$GeneAge=GeneAge[match(SpliceSites$gene,names(GeneAge))]

resamp_Age_maxFC=sapply(1:10000,function(i){toSample=By(SpliceSites$maxlFC[w]>1,SpliceSites$GeneAge[w],sum)
	samp=unlist(lapply(names(toSample),function(i){sample(which(SpliceSites$GeneAge[w]==i),toSample[i])}))
	mean(SpliceSites$Time_FirstAppeared_MA[w[samp]],na.rm=T)})
mean(resamp_Age_maxFC<=mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T)) # P<10-4

resamp_Age_immune=sapply(1:10000,function(i){toSample=By(SpliceSites$isImmune[w]=='yes',SpliceSites$GeneAge[w],sum)
	samp=unlist(lapply(names(toSample),function(i){sample(which(SpliceSites$GeneAge[w]==i),toSample[i])}))
	mean(SpliceSites$Time_FirstAppeared_MA[w[samp]],na.rm=T)})
mean(resamp_Age_immune<=mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']],na.rm=T)) 

save(resamp_Age_immune,resamp_Age_maxFC,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/resamp_Age_spliceSites.Rdata',HOME))

shapiro.test(resamp_Age_immune) # P= 0.38
shapiro.test(resamp_Age_maxFC)	# P=0.1622
GeneAge=By(SpliceSites$newAge[w],SpliceSites$gene[w],max)

pnorm(mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$isImmune[w]=='yes']],na.rm=T),mean(resamp_Age_immune),sd(resamp_Age_immune))	# 1.928401e-18
pnorm(mean(SpliceSites$Time_FirstAppeared_MA[w[SpliceSites$maxlFC[w]>1]],na.rm=T),mean(resamp_Age_maxFC),sd(resamp_Age_maxFC))	# 1.928401e-99



##############################################
# 	FINAL VERSION FOR MANUSCRIPT			##
# constitutive and alternative # Fig4C		##
##############################################

#  Figure S6 C and D
w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)
wConst=which(SpliceSites$highCovAnyCond_junc & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10 & SpliceSites$hasConstitutiveJunc)
wAlt=which(SpliceSites$highCovAnyCond_junc & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10 & !SpliceSites$hasConstitutiveJunc)

# constitutive and alternative # Fig4C
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')

layout(matrix(c(1,1,2,2,2,3,3,3),nrow=1))
par(mar=c(4,5,1,0.5))
tab=table(x,ifelse(SpliceSites$hasConstitutiveJunc[w],'const','alt'))
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)
#tab=table(x,SpliceSites$isImmune[w])
#barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

x=cut(SpliceSites$newAge[wConst],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')

tab0=table(x)
tab1=table(x[SpliceSites$isImmune[wConst]=='yes'])
tab2=table(x[SpliceSites$maxlFC[wConst]>1])
tab3=table(x[SpliceSites$maxlFC[wConst]>1 & SpliceSites$isImmune[wConst]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of constitutive\n splice sites',las=1)

x=cut(SpliceSites$newAge[wAlt],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')

tab0=table(x)
tab1=table(x[SpliceSites$isImmune[wAlt]=='yes'])
tab2=table(x[SpliceSites$maxlFC[wAlt]>1])
tab3=table(x[SpliceSites$maxlFC[wAlt]>1 & SpliceSites$isImmune[wAlt]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of alternative\n splice sites',las=1)

wilcox.test(as.numeric(x)[isImmune=='immune'],as.numeric(x))$p.value# p-value = 3.581088e-187
wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1],as.numeric(x))$p.value# p-value = 4.79261e-281
wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'],as.numeric(x)[SpliceSites$maxlFC[w]>1])$p.value# p-value = 2.220398e-115
wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'],as.numeric(x)[SpliceSites$isImmune[w]=='yes'])$p.value# p-value = 4.688095e-113
wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1],as.numeric(x)[SpliceSites$isImmune[w]=='yes'])$p.value# p-value = 0.02





plotGO(TableS4B_highCov_coding[1:15,],mar=20)

rk_time=data.frame(rank=rank(unique(SpliceSites$Time_FirstAppeared)),time=unique(SpliceSites$Time_FirstAppeared))
rk_time$Species=SpliceSites$FirstAppeared[match(rk_time$time,SpliceSites$Time_FirstAppeared)]
rk_time=rk_time[order(rk_time$time),]
rk_time$group=cut(rk_time$rank,c(-1,5,9,12,13,14,15,16,17,18))
levels(rk_time$group)=c('apes','primates','Placentals','Marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
rk_time$time_MA=rep(NA,nrow(rk_time))
rk_time$time_MA[1:4]=c(0,5.5,8,14)
rk_time$time_MA[10:14]=c(75,85,90,173,185)
rk_time$time_MA[10:18]=c(75,85,90,173,185,250,360,460,570)
rk_time$time_MA=round(approx(rk_time$time,rk_time$time_MA,rk_time$time)$y)

rk_time
#   rank    time      Species       group time_MA
#8     1 0.00000       Humans        apes       0
#17    2 0.00659  Chimpanzees        apes       6
#13    3 0.00877     Gorillas        apes       8
#4     4 0.01871   Orangutans        apes      14
#5     5 0.03297      Gibbons        apes      22
#2     6 0.04297      Monkeys    primates      27
#15    7 0.09988     Tarsiers    primates      58
#7     8 0.11119       Lemurs    primates      64
#18    9 0.12649       Shrews    primates      72
#1    10 0.13118      Rodents  Placentals      75
#11   11 0.15185 Boroeutheria  Placentals      85
#3    12 0.17513 Afroeutheria  Placentals      90
#16   13 0.30788   Marsupials  Marsupials     173
#12   14 0.38031   Monotremes  monotremes     185
#6    15 0.49021   Sauropsids  sauropsids     250
#14   16 0.65636   Amphibians   tetrapods     360
#9    17 0.95676  Vertebrates vertebrates     460
#10   18 1.46805    Chordates   chordates     570

SpliceSites$rank_FirstAppeared=rk_time$rank[match(SpliceSites$Time_FirstAppeared,rk_time$time)]
SpliceSites$Time_FirstAppeared_MA=rk_time$time_MA[match(SpliceSites$newAge,rk_time$rank)]
mycols <- colorRampPalette(brewer.pal(8,"BuPu"))(10)[2:10]



###### Fig S7B Sup: GerpRS by Age Group and 
w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$coding & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
GeneAge=By(SpliceSites$newAge[w],SpliceSites$gene[w],max)
Gerp=SpliceSites$Gerp_site[w]
isImmune=ifelse(SpliceSites$isImmune[w]=='yes','immune','non immune')
layout(matrix(c(1,1,rep(2,6)),nrow=1))
par(mar=c(7,4,2,1))
boxplot(Gerp~ifelse(isImmune=='non immune','a','b'),notch=T,col=mycols[c(3,7)],axes=F,pch=16,cex=0.3,ylab='GerpRS',las=1)
axis(2,las=1); 
axis(1,at=1:2,labels=c('non immune', 'immune'),las=2); 

boxplot(Gerp~paste(as.numeric(x),ifelse(isImmune=='non immune','a','b')),notch=T,col=mycols[c(3,7)],axes=F,pch=16,cex=0.3,ylab='GerpRS',las=1)
axis(2,las=1); 
par(xpd=T)
text(9.5,-9,labels='estimated age')
legend(11,-3,c('non immune','immune'),fill=mycols[c(3,7)],bty='n');axis(1,at=1:2,labels=c('',''))
for (i in 1:9){P=wilcox.test(Gerp[isImmune=='immune' & as.numeric(x)==i],Gerp[isImmune=='non immune'& as.numeric(x)==i])$p.value; 
				text(2*(i-1)+1.5,6.5,labels=ifelse(P<0.001,"***",ifelse(P<0.01,'**',ifelse(P<0.05,'*','ns'))),col='grey')}
for (i in 1:9){axis(1,at=2*(i-1)+1:2,labels=c('',''))}
for (i in 1:9){axis(1,at=2*(i-1)+1.5,tick=F,labels=levels(x)[i],las=2)}


###### Fig 4D 

CommonAncestorName=c('Humans','Chimpanzees','Gorillas','Orangutans','Gibbons','Monkeys','Tarsiers','Lemurs','Shrews','Rodents','Boroeutheria','Afroeutheria','Marsupials','Monotremes','Sauropsids','Amphibians','Vertebrates','Chordates')
library(ape)

leafLabels_CSF3=strsplit(as.data.frame(HumanSpecificSpliceSites[HumanSpecificSpliceSites$symbol=='CSF3','allSequence'])[[1]],':')[[1]]
leafLabels_CSF3R=strsplit(as.data.frame(HumanSpecificSpliceSites[HumanSpecificSpliceSites$symbol=='CSF3R','allSequence'])[[1]],':')[[1]]

#plot(tree,type="fan")
strand='+'
reverse=function(Seq){
	Corresp=c('A','G','T','C','-')
	names(Corresp)=c('T','C','A','G','-')
	sapply(strsplit(Seq,''),function(x){paste(rev(Corresp[x]),collapse='')})
	}
leafLabels_CSF3=leafLabels_CSF3[1:11]
leafLabels_CSF3R=leafLabels_CSF3R[1:11]

tree=read.tree(file=sprintf('%s/Annotation/Conservation/SEqAlign/tree46vertebrates.nwk',HOME))

layout(matrix(1:2,nrow=1))
par(mar=c(4,1,1,3),xpd=T)
tree$tip.label=sapply(leafLabels_CSF3,function(x){x1=substr(x,1,nchar(x)-3);x2=substr(x,nchar(x)-1,nchar(x));paste(gsub('-$','- ',gsub('^-',' -',ifelse(strand=='+',x2,reverse(x2)))),x1,sep=' : ')})
plot(subtrees(tree)[[10]],use.edge.length=F,tip.color=ifelse(substr(leafLabels_CSF3,nchar(leafLabels_CSF3)-1,nchar(leafLabels_CSF3))==substr(leafLabels_CSF3[1],nchar(leafLabels_CSF3[1])-1,nchar(leafLabels_CSF3[1])),'#DC143C','black'))
tree$tip.label=sapply(leafLabels_CSF3R,function(x){x1=substr(x,1,nchar(x)-3);x2=substr(x,nchar(x)-1,nchar(x));paste(gsub('-$','- ',gsub('^-',' -',ifelse(strand=='+',x2,reverse(x2)))),x1,sep=' : ')})
plot(subtrees(tree)[[10]],use.edge.length=F,tip.color=ifelse(substr(leafLabels_CSF3R,nchar(leafLabels_CSF3R)-1,nchar(leafLabels_CSF3R))==substr(leafLabels_CSF3R[1],nchar(leafLabels_CSF3R[1])-1,nchar(leafLabels_CSF3R[1])),'#DC143C','black'))

treePrim=tree
treePrim$tip.labels=treePrim$tip.labels[1:11]
subtreeplot(tree)

#age=c(,Boroeuteria=90e6,
#tab=table(x[ConstitutiveJunc$GerpRS_start<2],ConstitutiveJunc$isImmune[ConstitutiveJunc$GerpRS_start<2])
#colnames(tab)=c('non immune','immune')
#chisq.test(tab[12:18,])
#wilcox.test(ConstitutiveJunc$T_FirstAppeared.x[ConstitutiveJunc$GerpRS_start<2  & ConstitutiveJunc$isImmune],ConstitutiveJunc$T_FirstAppeared.x[ConstitutiveJunc$GerpRS_start<2 & !ConstitutiveJunc$isImmune])# p-value = 0.8429
#By(ConstitutiveJunc$T_FirstAppeared.x[ConstitutiveJunc$GerpRS_start<3 ],ConstitutiveJunc$isImmune[ConstitutiveJunc$GerpRS_start<3],mean) # NI 10.851852  IM: 9.607306 
#By(ConstitutiveJunc$T_FirstAppeared.x[ConstitutiveJunc$GerpRS_start<2 ],ConstitutiveJunc$isImmune[ConstitutiveJunc$GerpRS_start<2],mean) # NI 7.720149 7.708333  IM: 7.708333 
#By(ConstitutiveJunc$T_FirstAppeared.x[ConstitutiveJunc$GerpRS_start>4 ],ConstitutiveJunc$isImmune[ConstitutiveJunc$GerpRS_start>4],mean) # 16.12420 15.60438 

# constitutive # Fig4C 
w=which(SpliceSites$hasConstitutiveJunc & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)

# constitutive and alternative # Fig4C
w=which(SpliceSites$highCovAnyCond_junc & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')

tab=table(x,ifelse(SpliceSites$hasConstitutiveJunc[w],'const','alt'))
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)
tab=table(x,SpliceSites$isImmune[w])
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)

tab0=table(x)
tab1=table(x[SpliceSites$isImmune[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC[w]>1])
tab3=table(x[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of splice sites',las=1)
wilcox.test(as.numeric(x)[isImmune=='immune'],as.numeric(x))$p.value# p-value = 3.581088e-187
wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1],as.numeric(x))$p.value# p-value = 4.79261e-281
wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'],as.numeric(x)[SpliceSites$maxlFC[w]>1])$p.value# p-value = 2.220398e-115
wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'],as.numeric(x)[SpliceSites$isImmune[w]=='yes'])$p.value# p-value = 4.688095e-113
wilcox.test(as.numeric(x)[SpliceSites$maxlFC[w]>1],as.numeric(x)[SpliceSites$isImmune[w]=='yes'])$p.value# p-value = 0.02

SpliceSites[w[which(SpliceSites$NbAlignedPrimates[w]-1>=4 & SpliceSites$FirstAppeared[w]=='Humans' & SpliceSites$NbMatchingSeqPrim[w]==1)],]
# alternative only
w=which(!SpliceSites$hasConstitutiveJunc & SpliceSites$highCovAnyCond_junc & SpliceSites$Quality_MultiZ46Align & SpliceSites$maxFPKM>10)
x=cut(SpliceSites$newAge[w],c(-1,5,9,12,13,14,15,16,17,18))
levels(x)=c('apes','primates','placentals','marsupials','monotremes','sauropsids','tetrapods','vertebrates','chordates')
tab0=table(x)
tab1=table(x[SpliceSites$isImmune[w]=='yes'])
tab2=table(x[SpliceSites$maxlFC[w]>1])
tab3=table(x[SpliceSites$maxlFC[w]>1 & SpliceSites$isImmune[w]=='yes'])
barplot(t(rbind(tab0/sum(tab0),tab1/sum(tab1),tab2/sum(tab2),tab3/sum(tab3)))*100,beside=F,border=NA,col=rev(mycols),ylab='% of alternative splice sites',las=1)

# legend
par(mar=c(7,4,2,1))
X=rbind(1:9,1:9)-5
colnames(X)=levels(x)
Image(X,col=rev(mycols))
barplot(t(t(tab)/apply(tab,2,sum))*100,beside=F,border=NA,col=mycols,ylab='% of splice sites',las=1)
wilcox.test(as.numeric(x)[isImmune=='immune'],as.numeric(x)[isImmune=='non immune'])# p-value = 3.885689e-234

pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/SpliceSiteAge_Constitutive_immuneVSnonimmune.pdf',HOME))
barplot(tab,beside=F,border=NA,ylab='Nb of constitutive, non conserved splice sites',las=1)
dev.off()



mean(SpliceSites$Gerp_site[SpliceSites$coding & SpliceSites$Quality_MultiZ46Align]>2) # 99.1% ok
luq(SpliceSites$symbol[SpliceSites$highCovAnyCond_junc & SpliceSites$Quality_MultiZ46Align]) # 8187 gene symbols

#mean(paste(ConstitutiveJunc$chrom,ConstitutiveJunc$end+1)%in%paste(Exons[,'Chromosome Name'],Exons[,'Genomic coding start']))
#mm=match(paste(ConstitutiveJunc$chrom,ConstitutiveJunc$end+1),paste(Exons[,'Chromosome Name'],Exons[,'Genomic coding start']))
#ConstitutiveJunc$gene_end=Exons[mm,'Ensembl Gene ID']
#ConstitutiveJunc$symbol_end=Exons[mm,'Associated Gene Name']
#mean(paste(ConstitutiveJunc$chr,ConstitutiveJunc$start)%in%paste(Exons[,'Chromosome Name'],Exons[,'Genomic coding end']))
#mm=match(paste(ConstitutiveJunc$chr,ConstitutiveJunc$start),paste(Exons[,'Chromosome Name'],Exons[,'Genomic coding end']))
#ConstitutiveJunc$gene_start=Exons[mm,'Ensembl Gene ID']
#ConstitutiveJunc$symbol_start=Exons[mm,'Associated Gene Name']
#ConstitutiveJunc=ConstitutiveJunc[which(!is.na(ConstitutiveJunc$gene_end) & !is.na(ConstitutiveJunc$gene_start)),]
#mean(ConstitutiveJunc$GerpRS_start<2)
#mean(ConstitutiveJunc$GerpRS_end<2)
#1-mean(ConstitutiveJunc$GerpRS_end<2 | ConstitutiveJunc$GerpRS_start<2) # 0.9814851



w=which(SpliceSites$Quality_MultiZ46Align)
mean(SpliceSites$Gerp_site[w]>2)
mean(SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_NS>100]>2)
mean(SpliceSites$Gerp_site[SpliceSites$isImmune=='yes' & SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_NS>100]>2)
mean(SpliceSites$Gerp_site[SpliceSites$isImmune=='yes' & SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_IAV>100]>2)
mean(SpliceSites$Gerp_site[SpliceSites$isImmune=='yes' & SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_R848>100]>2)
mean(SpliceSites$Gerp_site[SpliceSites$isImmune=='yes' & SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_NS>100]>2)

#tic=Sys.time()
#SummaryStats_start=t(sapply(strsplit(as.data.frame(StartConserv[1:1000,])$V12,':'),function(x){y=gsub('.*_([ATGC-]*)','\\1',x);c(FirstAlignedSeq=max(which(y!='--')),NbAlignedSeq=sum(y!='--'),NbMatchingSeq=sum(y==y[1]),NbObservedSeq=luq(y[y!='--']))}))
#toc=Sys.time()
#print(toc-tic)
#SummaryStats_end=t(sapply(strsplit(as.data.frame(EndConserv)$V12,':'),function(x){y=gsub('.*_([ATGC-]*)','\\1',x);c(FirstAlignedSeq=max(which(y!='--')),NbAlignedSeq=sum(y!='--'),NbMatchingSeq=sum(y==y[1]),NbObservedSeq=luq(y[y!='--']))}))

# sapply(strsplit(as.data.frame(StartConserv[V2%in%x$start,])$V12,':'),function(x){y=gsub('.*_([ATGC-]*)','\\1',x);c(FirstAlignedSeq=max(which(y!='--')),NbAlignedSeq=sum(y!='--'),NbMatchingSeq=sum(y==y[1]),NbObservedSeq=luq(y[y!='--']))})

wNS=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_NS>100)
odds.ratio(table(SpliceSites$isImmune[wNS],SpliceSites$Gerp_site[wNS]<2))  # immune genes are enriched in non conserved even among NS splice site (not the effect of stim specific)
wR848=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_R848>100)
odds.ratio(table(SpliceSites$isImmune[wR848],SpliceSites$Gerp_site[wR848]<2))# immune genes are enriched in non conserved even among R848 splice site (not the effect of stim specific)
wIAV=which(SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_IAV>100)
odds.ratio(table(SpliceSites$isImmune[wIAV],SpliceSites$Gerp_site[wIAV]<2)) # immune genes are not enriched in non conserved among IAV splice site

odds.ratio(table(SpliceSites$Quality_MultiZ46Align,SpliceSites$isImmune))

Gerp_mean=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_2bp_random.txt',HOME))$V1
hh=hist(Gerp_mean,br=100,xlim=c(-7,7))
hist(SpliceSites$Gerp_site)

Junc_annot$highCovAnyCond_junc=apply(sapply(1:5,function(i){Junc_annot[,grep('Total_count_',colnames(Junc_annot))[i],with=F]/table(SampleAnnot$cond)[i]>10}),1,any)
SpliceSites$highCovAnyCond=apply(sapply(1:5,function(i){SpliceSites[,grep('Total_count_',colnames(SpliceSites))[i],with=F]/table(SampleAnnot$cond)[i]>10}),1,any)


hh=hist(c(-10,sample(Gerp_mean,sum(SpliceSites$Quality_MultiZ46Align),replace=T),Inf),br=50,plot=F)
hh$counts=hh$count/1000
plot(hh,xlim=c(-7,7),col=rcol[6],xlab='GerpRS',ylab='count',las=1,main='')
hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align])),br=hh$br,plot=F);
hh1$counts=hh1$count/1000;plot(hh1,col=rcol[3],add=T,border='#00000022')
#hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count>1000])),br=hh$br,col=rcol[1],add=T)
hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & highCovAnyCond])),br=hh$br,plot=F);
hh1$counts=hh1$count/1000;plot(hh1,col=rcol[1],add=T,border='#00000022')
par(xpd=T);legend(-10,700,fill=rcol[c(6,3,1)],legend=c('Genome Wide','splice sites','splice sites, >10 reads per sample'),bty='n');par(xpd=F)
#hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$inEnsembl])),br=hh$br,plot=F);
#hh1$counts=hh1$count/1000;plot(hh1,col=rcol[1],add=T,border='#00000022')
#hh1=hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Nbread>10000])),br=hh$br,plot=F);
#hh1$counts=hh1$count/1000;plot(hh1,density=30,add=T,border='#000000')

hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_NS<10 & SpliceSites$Total_count>1000])),br=hh$br,col=rcol[3])
hh=hist(c(-10,sample(Gerp_mean,sum(SpliceSites$Quality_MultiZ46Align),replace=T),Inf),br=100,xlim=c(-7,7))
hist(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_NS<10 & SpliceSites$Total_count>100])),br=hh$br,col=rcol[3])
mean(pmax(-10,pmin(7,SpliceSites$Gerp_site[SpliceSites$Quality_MultiZ46Align & SpliceSites$Total_count_NS<10 & SpliceSites$Total_count>100]))>2)

# % of Novel splice site by condition, for various coverage thresholds
barplot(sapply(1:5,function(i){mean(is.na(Junc_annot$Nbread[Junc_annot[[i+1]]>20]))}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[3]}),add=F)
barplot(sapply(1:5,function(i){mean(is.na(Junc_annot$Nbread[Junc_annot[[i+1]]>50]))}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[2]}),add=T)
barplot(sapply(1:5,function(i){mean(is.na(Junc_annot$Nbread[Junc_annot[[i+1]]>100]))}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[1]}),add=T)

# % of Novel splice site by condition, for various coverage thresholds, only AG GT AG GC or AT AC
par(mar=c(6,4,4,1))
barplot(sapply(1:5,function(i){mean(is.na(Junc_annot$Nbread[Junc_annot[[i+1]]>20 & Junc_annot$DA_type!='other']))}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[3]}),add=F,ylab='% of novel Junctions',las=2,names.arg=condIndex)
barplot(sapply(1:5,function(i){mean(is.na(Junc_annot$Nbread[Junc_annot[[i+1]]>50 & Junc_annot$DA_type!='other']))}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[2]}),add=T,axes=F)
barplot(sapply(1:5,function(i){mean(is.na(Junc_annot$Nbread[Junc_annot[[i+1]]>100 & Junc_annot$DA_type!='other']))}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[1]}),add=T,axes=F)
par(xpd=T)
text(0,0.06,'read coverage',adj=c(0,1))
legend(-3.5,0.053,'top',fill=colorRampPalette(c(colERC[2],'white'))(5)[1:3],legend=c('>100', '>50', '>20'),bty='n',ncol=3,x.intersp=0.1)

# % of non canonical splice site by condition, for various coverage thresholds, only AG GT AG GC or AT AC
par(mar=c(7,4,7,1))
barplot(sapply(1:5,function(i){mean(Junc_annot$DA_type[Junc_annot[[i+1]]>20 ]=='other')}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[3]}),add=F,ylab='% of non canonical Junctions',las=2,names.arg=condIndex)
barplot(sapply(1:5,function(i){mean(Junc_annot$DA_type[Junc_annot[[i+1]]>50 ]=='other')}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[2]}),add=T,axes=F)
barplot(sapply(1:5,function(i){mean(Junc_annot$DA_type[Junc_annot[[i+1]]>100]=='other')}),col=sapply(colERC[2+0:4*2],function(x){colorRampPalette(c(x,'white'))(5)[1]}),add=T,axes=F)
par(xpd=T)

barplot(sapply(1:5,function(i){table(Junc_annot$DA_type[Junc_annot[[i+1]]>20 ]=='other')})/1000,col=colERC,add=F,ylab='Nb of junction (x1000, canon vs non canon)',las=2,names.arg=condIndex,beside=T)
barplot(sapply(1:5,function(i){table(is.na(Junc_annot$Nbread[Junc_annot[[i+1]]>20 & Junc_annot$DA_type!='other']))})/1000,col=colERC,add=F,ylab='Nb of Junctions (x1000, novel/known)',las=2,names.arg=condIndex,beside=T)


By(Junc_annot$DA_type!='other',paste(Junc_annot$coverage_level,ifelse(is.na(Junc_annot$Nbread),'new','known')),mean)
sapply(2:5,function(i){odds.ratio(cbind(table(is.na(Junc_annot$Nbread[Junc_annot$DA_type!='other' & Junc_annot[[2]]>20])),table(is.na(Junc_annot$Nbread[Junc_annot$DA_type!='other' & Junc_annot[[i+1]]>20]))))})
# Enrichment in novel junctions, compared to NS 
#        [,1]         [,2]         [,3]         [,4]    
#LowerCI 1.049402     1.16371      1.535575     4.526615
#OR      1.112438     1.23104      1.619346     4.737682
#UpperCI 1.179258     1.302382     1.707996     4.960714
#alpha   0.05         0.05         0.05         0.05    
#P       0.0003169576 2.667785e-13 1.629805e-73 0       

par(mar=c(12,4,1,1))
barplot(1-By(Junc_annot$DA_type!='other',paste(Junc_annot$coverage_level,ifelse(is.na(Junc_annot$Nbread),'new','known')),mean),las=2,ylab='% of false positives')





hist(c(Junc_annot$GerpRS_start[wIntro[!duplicated(paste(Junc_annot$chrom,Junc_annot$start)[wIntro])]],Junc_annot$GerpRS_end[wIntro[!duplicated(paste(Junc_annot$chrom,Junc_annot$end)[wIntro])]]),add=F)
hist(c(Junc_annot$GerpRS_start[wIntro[!duplicated(paste(Junc_annot$chrom,Junc_annot$start)[wGTAG])]],Junc_annot$GerpRS_end[wIntro[!duplicated(paste(Junc_annot$chrom,Junc_annot$end)[wGTAG])]]),add=F)

sapply(2:5,function(i){odds.ratio(cbind(table(Junc_annot$GerpRS_meanEnd[wIntro][Junc_annot$Total_count_NS[wIntro]>100]<2),table(Junc_annot$GerpRS_meanEnd[wIntro][Junc_annot[wIntro,grep('Total_count',cn(Junc_annot))[i+1]]>100]<2)))})
sapply(2:5,function(i){odds.ratio(cbind(table(Junc_annot$GerpRS_meanEnd[wIntro][Junc_annot$Total_count_NS[wIntro]>100]<2),table(Junc_annot$GerpRS_meanEnd[wIntro][Junc_annot[wIntro,grep('Total_count',cn(Junc_annot))[i+1]]>100]<2)))})

mean(Junc_annot$GerpRS_meanEnd[wIntro][which(abs(log2(Junc_annot$Total_count_R848[wIntro]/Junc_annot$Total_count_NS[wIntro]))< 0.2)]>2) # 0.5460719
mean(Junc_annot$GerpRS_meanEnd[wIntro][which(log2(Junc_annot$Total_count_R848[wIntro]/Junc_annot$Total_count_NS[wIntro])>1)]>2) # 0.4006715



# less conserved site active after stimulation.

######## load sequence hg19 in R
#library("Biostrings")
#hg19 <- readDNAStringSet("my.fasta")
#seq_name = names(hg19)
#sequence = paste(hg19)
#df <- data.frame(seq_name, sequence)

library(data.table)
library(GenomicRanges)
Gerp_removed=list()
Gerp_chr=list()
Map_select=as.data.frame(fread(sprintf('cut -f 1-10 %s/Annotation/Select/Map_Select_LD_allChr_V2.txt',HOME)))
for (chr in 1:22){
	cat(chr,'\n')
#	Conservation_chr=fread(sprintf('%s/Annotation/Conservation/chr%s_mammalian_conservation_perBp.txt',HOME,chr))
	hg19_align=fread(sprintf('%s/Annotation/Conservation/SEqAlign/Alignable_Human/chr%s.align.txt',HOME,chr))
	hg19_align$chrom=gsub('hg19.chr','',hg19_align$chrom,fixed=T)
	hg19_align_GR=reduce(makeGRangesFromDataFrame(as.data.frame(hg19_align),seqnames='chrom',start.field='start',end.field='end'))
	Gerp_chr[[chr]]=fread(sprintf('%s/Annotation/Conservation/hg19.GERP_scores/chr%s.maf.rates',HOME,chr))
	hg19_nalign_GR=setdiff(GRanges(seqnames=chr,range=IRanges(start=1,end=nrow(Gerp_chr[[chr]]))),hg19_align_GR)
	colnames(Gerp_chr[[chr]])=c('GerpExp','GerpRS')
	Gerp_GR=GRanges(seqnames=chr,range=IRanges(start=1:nrow(Gerp_chr[[chr]]),end=1:nrow(Gerp_chr[[chr]])))
	oo=findOverlaps(Gerp_GR,hg19_nalign_GR)
	rm(Gerp_GR)
	Gerp_removed[[chr]]=Gerp_chr[[chr]]$GerpRS[queryHits(oo)]
	Gerp_chr[[chr]]$GerpRS[queryHits(oo)]=NA
	rm(oo)
}
Gerp_all=unlist(lapply(Gerp_chr,function(x){x$GerpRS}))
Gerp_removed=unlist(Gerp_removed)
samp=sample(1:length(Gerp_all),1000000)
Gerp_mean=(Gerp_all[samp]+Gerp_all[samp+1])/2
Gerp_mean=Gerp_mean[!is.na(Gerp_mean)]

Gerp_around_snp=sapply(seq(-30,30),function(x){mean(unlist(unlist(sapply(1:22,function(chr){Gerp_chr[[chr]]$GerpRS[as.integer(Map_select$position[which(Map_select$chrom==chr)]+x)]}))),na.rm=T)})
#[1] -0.10195841 -0.09530038 -0.09535411 -0.10235445 -0.09561574 -0.09672540
# [7] -0.10482750 -0.09597736 -0.09633529 -0.10379186 -0.09720675 -0.09718729
#[13] -0.10570854 -0.09857206 -0.10051198 -0.10656766 -0.09839888 -0.09528731
#[19] -0.10350920 -0.09714013 -0.10131286 -0.10699713 -0.09969746 -0.09841873
#[25] -0.10956851 -0.10201386 -0.10499130 -0.11109377 -0.12045880 -0.20233655
#[31] -0.39967812 -0.20308892 -0.12064903 -0.11029945 -0.10491373 -0.10241877
#[37] -0.11215068 -0.10112411 -0.10323354 -0.11043238 -0.10177825 -0.10022454
#[43] -0.10722151 -0.09731001 -0.10003002 -0.10771311 -0.10189674 -0.10054210
#[49] -0.10570594 -0.09819324 -0.09731862 -0.10392663 -0.09608125 -0.09779310
#[55] -0.10516188 -0.09734554 -0.09665850 -0.10245400 -0.09739701 -0.09615181
#[61] -0.10358595

Gerp_snp=unlist(sapply(1:22,function(chr){Gerp_chr[[chr]]$GerpRS[as.integer(Map_select$position[which(Map_select$chrom==chr)])]}))

Map_select$GerpRS=Gerp_snp
write.table(Map_select[,c('snp.name','chromosome','position','SNPfreq','GerpRS')],file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Map_Gerp_snp.txt',HOME),quote=F,row.names=F,col.names=F,sep='\t')
write.table(Gerp_snp,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_snp.txt',HOME),quote=F,row.names=F,col.names=F,sep='\t')
write.table(Gerp_all[samp],file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_1bp_random.txt',HOME),quote=F,row.names=F,col.names=F,sep='\t')

write.table(Gerp_mean,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_2bp_random.txt',HOME),quote=F,row.names=F,col.names=F,sep='\t')
write.table(Gerp_removed,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_removed.txt',HOME),quote=F,row.names=F,col.names=F,sep='\t')

Gerp_mean=fread(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/NoisySplicing/Gerp/Gerp_2bp_random.txt',HOME))$V1
#pdf(sprintf('%s/03_Analysis/Splicing/NoisySplicing/Conservation_GerpRS_DA_site_VS_Genome_%s.pdf',HOME,sample))
hh=hist(Gerp_mean,br=100,freq=F,xlim=c(-6,6),xlab='GerpRS',main='Conservation of Junctions')
hist(Gerp_all[intersect(wGerp,Map_select$position)],add=T,col=rcol[8],freq=F,br=hh$breaks,border='#00000000')

Junc_annot$is_goodAlignment=FALSE
for (chr in 1:22){
	cat(chr,'\n')
#	Conservation_chr=fread(sprintf('%s/Annotation/Conservation/chr%s_mammalian_conservation_perBp.txt',HOME,chr))
	hg19_align=fread(sprintf('%s/Annotation/Conservation/SEqAlign/Alignable_Human/chr%s.align.txt',HOME,chr))
	hg19_align$chrom=gsub('hg19.chr','',hg19_align$chrom,fixed=T)
	hg19_align_GR=reduce(makeGRangesFromDataFrame(as.data.frame(hg19_align),seqnames='chrom',start.field='start',end.field='end'))
	Junc_GR=makeGRangesFromDataFrame(Junc_annot,seqnames='chrom',start.field='start',end.field=
	Junc_annot$is_goodAlignment[Junc_annot$chrom==chr]=

hist(c(Junc_annot$GerP_acceptor,Junc_annot$GerP_donor)[Junc_annot$count>100],add=T,col=SetAlpha(acol[1],0.3),freq=F,br=hh$breaks,border='#00000000')
hist(c(Junc_annot$GerP_acceptor,Junc_annot$GerP_donor)[Junc_annot$count==1],add=T,col=SetAlpha(acol[3],0.5),freq=F,br=hh$breaks,border='#00000000')
#hist(Conservation_chr22[,'PhyloP'],add=T,freq=F,br=hh$breaks)
legend('topleft',fill=c('#00000000',rcol[8],SetAlpha(acol[3],0.5),SetAlpha(acol[1],0.3)),legend=c('Genome Wide','SNPs','junctions 1 read','junctions >100 read'))
#dev.off()
mean(Junc_annot$GerpRS_meanEnd[wIntro]>2)

