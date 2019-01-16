library(RcolorBrewer)
PheWas=read.csv2('Table_PheWAS_CYP3A5.csv',stringsAsFactor=F)

############################################################################################################################################
############ Table Table_PheWAS_CYP3A5.csv contain Phewas Pvalues for rs776746 retrieved from http://geneatlas.roslin.ed.ac.uk/ ############
############################################################################################################################################

######### PHEWAS phenotypes have been manually assigned to the following categories: (type variable)
# "Anthropometric measurements", "Bone and joints disorders", "Cancer", "Cardiovascular disorders", "Digestive and metabolic disorders", 
# "Heamatological measurements", "Immune disorders", "Lifestyle", "Nervous and psychiatric disorders", "Pulmonary disease", "Other"

PheWas$type=as.factor(PheWas$type)
PheWas$type=relevel(PheWas$type,'Other')
levels(PheWas$type)=c("Other", "Anthropometric measurements", "Bone and joints disorders", 
"Cancer", "Cardiovascular disorders", "Digestive and metabolic disorders", 
"Heamatological measurements", "Immune disorders", "Lifestyle", 
"Nervous and psychiatric disorders", "Pulmonary disease")

colors=rev(c(brewer.pal(8, 'Set1'),brewer.pal(8, 'Set2')))
colors=colors[c(1:8,10,13,16)]
PheWas$color=as.numeric(as.factor(PheWas$type))
PheWas$p.value=as.numeric(PheWas$p.value)

oo=order(-PheWas$color,PheWas$p.value)
PheWas$fdr=p.adjust(PheWas$p.value)

#####################################################
###                  Figure 6F                    ###
#####################################################
#4 x 7.75 inches

par(mar=c(5.1,4.1,0.5,2.1))
plot(-log10(PheWas[oo,'p.value']),col=ifelse(PheWas$fdr[oo]<0.05,'black',colors[PheWas[oo,'color']]),
                                  pch=ifelse(PheWas$fdr[oo]<0.05,ifelse(PheWas$Beta[oo]<0,24,25),16),
                                  cex=ifelse(PheWas$fdr[oo]<0.05,0.8,0.4),
                                  bg=colors[PheWas[oo,'color']],
                                  ylab=expression(-log[10](p)), 
                                  xlab='',
                                  axes=F)
axis(2, las=1) 
axis(1, at=c(0,cumsum(rev(table(PheWas$type)))), labels=FALSE) 
abline(h = -log10(0.05/779), lty = 2, col = 'grey')
par(xpd=T)
legend(-120, -1.5, ncol=2, bty = 'n', fill = rev(colors), legend = rev(levels(PheWas$type)), cex = 0.7)
par(xpd=F)

