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



