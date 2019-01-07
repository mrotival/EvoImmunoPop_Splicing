par(mar=c(7,4,3,1))
tab=table(RESobs_nodup_1Mb[which(RESobs_nodup_1Mb$r2_sQTL_eQTL>0.8),'BestModel'], is_rsQTL_any[which(RESobs_nodup_1Mb$r2_sQTL_eQTL>0.8)])
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)
colnames(tab)=c('basal sQTL','rsQTL')
x=barplot(tab,col=getColors(11)[9:11],las=3)
par(xpd=T);text(x,1.1,labels=n_sQTL);par(xpd=F,srt=90)
par(xpd=T);#legend(-2.5,1.5,legend=rownames(tab),fill=getColors(11)[9:11],bty='n',cex=0.85);par(xpd=F)


par(mar=c(7,4,3,1))
tab=table(RESobs_nodup_1Mb[which(RESobs_nodup_1Mb$r2_sQTL_eQTL>0.8),'BestModel'], is_rsQTL_any[which(RESobs_nodup_1Mb$r2_sQTL_eQTL>0.8)])
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)
colnames(tab)=c('basal sQTL','rsQTL')
x=barplot(tab,col=getColors(11)[9:11],las=3)
par(xpd=T);text(x,1.1,labels=n_sQTL);par(xpd=F,srt=90)
par(xpd=T);#legend(-2.5,1.5,legend=rownames(tab),fill=getColors(11)[9:11],bty='n',cex=0.85);par(xpd=F)

################# causality VS overlapping an eQTL
#pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/figures/Fig2_bestModel_rsQTL_horiz.pdf',HOME),height=5.6,width=3)
w=which(RESobs_nodup_1Mb$P_sQTL_geneExpr<0.05 & RESobs_nodup_1Mb$BestModel_Prob>0.8)
tab=table(RESobs_nodup_1Mb[w,'BestModel'],RESobs_nodup_1Mb$r2_sQTL_eQTL[w]>0.8)
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)*100
colnames(tab)=c('no','yes')
x=barplot(tab,col=getColors(11)[9:11],las=3)
par(xpd=T);
text(x,1.1,labels=n_sQTL);par(xpd=F,srt=90)

#pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/figures/Fig2_bestModel_legend.pdf',HOME),height=5.6,width=3)
plot.new()
legend('top',fill=getColors(11)[9:11],legend=rownames(tab),bty='n')
#tab
#                            no        yes
#  SNP->Expr->Splice 0.03672316 0.39887640
#  SNP->Splice->Expr 0.26553672 0.07303371
#  Splice<-SNP->Expr 0.69774011 0.52808989


pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/figures/Fig2_bestModel_sQTL_P_expr_horiz.pdf',HOME),height=5.6,width=3)
w=which(RESobs_nodup_1Mb$P_sQTL_geneExpr<0.05)
tab=table(RESobs_nodup_1Mb[w,'BestModel'],ifelse(RESobs_nodup_1Mb$P_sQTL_geneExpr[w]<1e-7,'yes','no'))
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)*100
x=barplot(tab,col=getColors(11)[9:11],las=3)
par(xpd=T);
text(x,1.1*max(apply(tab,2,sum)),labels=n_sQTL);par(xpd=F,srt=90)

w=which(RESobs_nodup_1Mb$P_sQTL_geneExpr < 0.05)
tab=table(RESobs_nodup_1Mb[w,'BestModel'],ifelse(RESobs_nodup_1Mb$P_sQTL_geneExpr[w]<1e-7,'yes','no'))
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)*100

x=barplot(tab,col=getColors(11)[9:11],las=3)
par(xpd=T);
text(x,1.1*max(apply(tab,2,sum)),labels=n_sQTL);par(xpd=F,srt=90)

# P_expr<0.05
# P_expr<1e-7
#                            no        yes
#  SNP->Expr->Splice  0.6578947 40.1913876
#  SNP->Splice->Expr 32.0175439  6.6985646
#  Splice<-SNP->Expr 67.3245614 53.1100478
#pdf(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/figures/Fig2_bestModel_rsQTL_noProbLimit_horiz.pdf',HOME),height=5.6,width=3)
w=which(RESobs_nodup_1Mb$P_sQTL_geneExpr<0.05)
tab=table(RESobs_nodup_1Mb[w,'BestModel'],RESobs_nodup_1Mb$r2_sQTL_eQTL[w]>0.8)
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)*100
colnames(tab)=c('no','yes')
x=barplot(tab[,2:1],col=getColors(11)[9:11],las=3,yaxt='n')
par(xpd=T);text(x,1.2*max(apply(tab,2,sum)),labels=paste('N =',rev(n_sQTL)));par(xpd=F,srt=90)
axis(4)
#dev.off()
#                          no        yes
#  SNP->Expr->Splice 0.03151261 0.38095238
#  SNP->Splice->Expr 0.30252101 0.08465608
#  Splice<-SNP->Expr 0.66596639 0.53439153

# % of sQTL affecting expression (10^-3)
barplot(by(RESobs_nodup_1Mb$P_sQTL_geneExpr<0.05,cut(RESobs_nodup_1Mb$R2,c(0,0.2,0.5,1)),mean),las=2)

w0=which(RESobs_nodup_1Mb$BestModel_Prob>0.8)
w1=which(RESobs_nodup_1Mb$BestModel_Prob>0.8 & RESobs_nodup_1Mb$P_sQTL_geneExpr<0.05)
w2=which(RESobs_nodup_1Mb$BestModel_Prob>0.8 & RESobs_nodup_1Mb$r2_sQTL_eQTL>0.8)

tab=cbind(all=table(RESobs_nodup_1Mb[w0,'BestModel']),expr=table(RESobs_nodup_1Mb[w1,'BestModel']),eQTL=table(RESobs_nodup_1Mb[w2,'BestModel']))
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)
colnames(tab)=c('basal sQTL','rsQTL')
x=barplot(tab,col=getColors(11)[9:11],las=3)
par(xpd=T);text(x,1.1,labels=n_sQTL);par(xpd=F,srt=90)
par(xpd=T);#legend(-2.5,1.5,legend=rownames(tab),fill=getColors(11)[9:11],bty='n',cex=0.85);par(xpd=F)

w1=which(RESobs_nodup_1Mb$P_sQTL_geneExpr<0.05)
w2=which(RESobs_nodup_1Mb$r2_sQTL_eQTL>0.8)

tab=cbind(expr=table(RESobs_nodup_1Mb[w1,'BestModel']),eQTL=table(RESobs_nodup_1Mb[w2,'BestModel']))
n_sQTL=apply(tab,2,sum)
tab=t(t(tab)/n_sQTL)
colnames(tab)=c('basal sQTL','rsQTL')
x=barplot(tab,col=getColors(11)[9:11],las=3)
par(xpd=T);text(x,1.1,labels=n_sQTL);par(xpd=F,srt=90)
par(xpd=T);#legend(-2.5,1.5,legend=rownames(tab),fill=getColors(11)[9:11],bty='n',cex=0.85);par(xpd=F)


tab=table(RESobs_nodup_1Mb[,'BestModel'])
colors=getColors(11)[9:11]
names(colors)=names(tab)
plot(-log10(RESobs_nodup_1Mb$pvalue),-log10(RESobs_nodup_1Mb$P_sQTL_geneExpr),pch=16,col=colors[RESobs_nodup_1Mb$BestModel])
