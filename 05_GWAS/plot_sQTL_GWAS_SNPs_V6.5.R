#GWAS_sQTL=read.table('/Users/mrotival/WORK/Splicing/GWAS_sQTL_cyto_2.txt',header=T,sep='\t',stringsAsFactors=F)
Source='Ens70_HISAT'
load(paste(EVO_IMMUNO_POP,'/Maxime/Splicing/MISO/aggregated/',Source,'/PSI_events_ALL_V5_withCounts.Rdata',sep=''))

library(igraph)
GWAS_sQTL=read.table(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V6/data/GWAS_sQTL_forGraph.txt',HOME),header=T,sep='\t',stringsAsFactors=F,quote='')
# GWAS_sQTL=read.table('/Users/mrotival/WORK/06_papers/PapierSplicing_20180126/V6.5/GWAS_sQTL_forGraph.txt',header=T,sep='\t',stringsAsFactors=F,quote='')
GWAS_sQTL$LPS_lFC=log2(GWAS_sQTL$LPS_Expr)-log2(GWAS_sQTL$NS_Expr)
GWAS_sQTL$PAM3CSK4_lFC=log2(GWAS_sQTL$PAM3CSK4_Expr)-log2(GWAS_sQTL$NS_Expr)
GWAS_sQTL$R848_lFC=log2(GWAS_sQTL$R848_Expr)-log2(GWAS_sQTL$NS_Expr)
GWAS_sQTL$IAV_lFC=log2(GWAS_sQTL$IAV_Expr)-log2(GWAS_sQTL$NS_Expr)
GWAS_sQTL$Max_lFC=pmax(GWAS_sQTL$LPS_lFC,GWAS_sQTL$PAM3CSK4_lFC,GWAS_sQTL$R848_lFC,GWAS_sQTL$IAV_lFC)

#GWAS_sQTL=GWAS_sQTL[GWAS_sQTL$dup_event==0,]
luq=function(x){length(unique(x))}

#for(i in grep('Testable_',colnames(PSI_Annot))){
#	GWAS_sQTL[[colnames(PSI_Annot)[i]]]=PSI_Annot[match(GWAS_sQTL$event_id,PSI_Annot$event_id),i]
#	}
#GWAS_sQTL$nsSpecific_event=GWAS_sQTL$Testable_NS & !apply(GWAS_sQTL[,grep('Testable_',cn(GWAS_sQTL))[-1]],1,any)
#GWAS_sQTL$stimSpecific_event= !GWAS_sQTL$Testable_NS & apply(GWAS_sQTL[,grep('Testable_',cn(GWAS_sQTL))[-1]],1,any)


colPSI=structure(c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
"#FFD92F", "#E5C494"), .Names = c("A3", "A5", "AF", "AL", "MX","RI", "SE"))

library(igraph)

nodes=data.frame(id=c(unique(GWAS_sQTL$GWAS_Trait_R2_1E5),unique(GWAS_sQTL$event_id),unique(GWAS_sQTL$haplo)),type=rep(c('GWAS','Splice','SNP'),c(luq(GWAS_sQTL$GWAS_Trait_R2_1E5),luq(GWAS_sQTL$event_id),luq(GWAS_sQTL$haplo))))
nodes$color=colPSI[GWAS_sQTL$event_type[match(nodes$id,GWAS_sQTL$event_id)]]
nodes$color[is.na(nodes$color)]='grey'
nodes$color[nodes$type=='SNP']='darkgrey'
nodes$size=ifelse(nodes$type=='GWAS',14,18)
nodes$size[nodes$type=='SNP']=5
nodes$vertex.label.color=ifelse(nodes$type=='GWAS',"black","black")
nodes$vertex.label.color[nodes$type=='SNP']='#00000000'
nodes$vertex.label.family='arial'
nodes$vertex.label.cex=ifelse(nodes$type=='GWAS',0.8,0.6)
nodes$vertex.frame.width= ifelse(GWAS_sQTL$Max_lFC[match(nodes$id,GWAS_sQTL$event_id)]>1, 2, 1)
nodes$vertex.frame.width[is.na(nodes$vertex.frame.width)]=0.01
nodes$vertex.frame.color=substr(colERC5[1],1,7)
#nodes$vertex.frame.color[which(GWAS_sQTL$is_stim_sQTL[match(nodes$id,GWAS_sQTL$symbol)]==1)]=substr(colERC5[4],1,7)
nodes$vertex.frame.color[which(GWAS_sQTL$Max_lFC[match(nodes$id,GWAS_sQTL$event_id)]>1)]='black'
nodes$vertex.frame.color[which(GWAS_sQTL$Max_lFC[match(nodes$id,GWAS_sQTL$event_id)]>1 & GWAS_sQTL$NS_Expr[match(nodes$id,GWAS_sQTL$event_id)]<10)]=colERC[4]
nodes$vertex.label=nodes$id
nodes$vertex.label[nodes$type=='SNP']=''
nodes$vertex.label[nodes$type=='Splice']=GWAS_sQTL$symbol[match(nodes$id[nodes$type=='Splice'],GWAS_sQTL$event_id)]



mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

add.vertex.shape("fcircle", clip=igraph.shape.noclip,
		plot=mycircle, parameters=list(vertex.frame.color=1,
                                  vertex.frame.width=1))

X=GWAS_sQTL[,c('event_id','haplo','event_type','isSignif_8','is_stim_sQTL','Pvalue_NS')]
X$isSignif_8=FALSE
X=X[!duplicated(X),-3]
colnames(X)[1]='id'
X$width=ifelse(X$Pvalue_NS>0.001,3,1)
X$color=colERC[8]
#X$color[which(X$is_stim_sQTL)]=substr(colERC[3],1,7)
X$color[which(X$Pvalue_NS>0.001 )]=substr(colERC[4],1,7)

Y=GWAS_sQTL[,c('GWAS_Trait_R2_1E5','haplo',"isSignif_8",'is_stim_sQTL','Pvalue_NS')]
Y[,4:5]=FALSE
Y=Y[!duplicated(Y),]
colnames(Y)[1]='id'
Y$width=ifelse(Y$isSignif_8,1,1)
Y$color=ifelse(Y$isSignif_8,grey(0.2),grey(0.6))

edges=rbind(X,Y)
#edges$width=ifelse(edges$isSignif_8,1,1)
#edges$color=ifelse(edges$isSignif_8,grey(0.2),grey(0.6))

net <- graph_from_data_frame(d=edges,vertices=nodes,directed=F)
load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/GWAS_sQTL_plot_coords_withSNP_v2.Rdata',HOME))
#l3[,2]=(l3[,2]-mean(range(l3[,2])))/diff(range(l3[,2]))*2
#l3[,1]=(l3[,1]-mean(range(l3[,1])))/diff(range(l3[,1]))*2
l3[,1]=(l3[,1]-mean(range(l3[,1])))
l3[,2]=(l3[,2]-mean(range(l3[,2])))
par()$mar
# 13x9 inches
plot(net,layout=l3,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color=grey(0.2),
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label,rescale=F)

# black circle : genes is up regulated after stimulated
# red circle :gene is expressed only upon stimulation
# red: line sQTL is detected only upon stimulation (P<1e-3)


plot(net,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color='black',
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label)


load(sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/GWAS_sQTL_plot_coords_withSNP_v2.Rdata',HOME))
#l3[,2]=(l3[,2]-mean(range(l3[,2])))/diff(range(l3[,2]))*2
#l3[,1]=(l3[,1]-mean(range(l3[,1])))/diff(range(l3[,1]))*2
l3[,1]=(l3[,1]-mean(range(l3[,1])))
par()$mar
plot(net,layout=l3,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color=grey(0.2),
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label,rescale=F)


tkid <- tkplot(net) #tkid is the id of the tkplot that will open
l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
l2=l
l2[,2]=(l2[,2]-500)/600
l2[,1]=(l2[,1]-750)/600
#save(l2,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/GWAS_sQTL_plot_coords_withSNP.Rdata',HOME))

nodes$vertex.label.color=ifelse(nodes$type=='GWAS',"black","black")
nodes$vertex.label.family='arial'
nodes$vertex.label.cex=ifelse(nodes$type=='GWAS',0.8,0.6)

nodes$vertex.frame.width= ifelse(GWAS_sQTL$is_rsQTL[match(nodes$id,GWAS_sQTL$event_id)]==1, 3, 1)
nodes$vertex.frame.width[is.na(nodes$vertex.frame.width)]=0.01
nodes$vertex.frame.color=colERC[2]
nodes$vertex.frame.color[which(GWAS_sQTL$is_nsSpecific[match(nodes$id,GWAS_sQTL$event_id)]==1)]=substr(colERC[6],1,7)
nodes$vertex.frame.color[which(GWAS_sQTL$is_stimSpecific[match(nodes$id,GWAS_sQTL$event_id)]==1)]=substr(colERC[4],1,7)

nodes$vertex.frame.color="#00000000"
nodes$vertex.frame.width= 1
nodes$vertex.frame.color[which(GWAS_sQTL$nsSpecific_event[match(nodes$id,GWAS_sQTL$event_id)])]=substr(colERC[6],1,7)
nodes$vertex.frame.color[which(GWAS_sQTL$stimSpecific_event[match(nodes$id,GWAS_sQTL$event_id)])]=substr(colERC[4],1,7)

net <- graph_from_data_frame(d=edges,vertices=nodes,directed=F)

plot(net,layout=l2,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color=grey(0.2),
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label,rescale=F)
		
tk_set_coords(tkid,l2*250+500)
l <- tkplot.getcoords(tkid)
l3=l
l3[,2]=(l3[,2]-250)/250
l3[,1]=(l3[,1]-750)/250


net <- graph_from_data_frame(d=edges,vertices=nodes,directed=F)
plot(net,layout=l3,vertex.shape="fcircle",
		vertex.frame.color=nodes$vertex.frame.color,
		vertex.frame.width=nodes$vertex.frame.width,
		vertex.label.cex=nodes$vertex.label.cex,
		vertex.label.color=grey(0.2),
		vertex.label.family='Arial',
		vertex.label=nodes$vertex.label,rescale=F)
save(l3,file=sprintf('%s/03_Analysis/Splicing/Papier_Splicing/V5/data/GWAS_sQTL_plot_coords_withSNP_v2.Rdata',HOME))


plot(net, layout=l3,mark.groups=list(c(1,4,5,8), c(15:17)), mark.col=c("#C5E5E7","#ECD89A"), mark.border=NA)

		