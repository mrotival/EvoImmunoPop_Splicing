
# plot output from GOSeq

plotGO = function(resGO,mar=18,...){
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
	legend("bottomright",fill=colPSI[1:3],legend=c('BP','MF','CC'),bty='n')
	par(mar=params$mar)
	}


# Fisher exact test with easy top use output
odds.ratio = function(tab, alpha = 0.05){	
    test=fisher.test(tab,conf.level=1-alpha)
    oframe <- data.frame(LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], alpha = alpha,P=test$p.value)
    oframe
}

# goseq with annotation of genes that contribute to a GO, corresponding Odds ratio and its 95% confidence interval
function(geneList,background,bias.data=NULL,FDR=0.05,biasCorrect=TRUE,overOnly=T,addGenes=T,addCI=F,allGOTerms){
	require(GO.db)
	require(goseq)
	
    if(addGenes & missing(allGOTerms)){
        allGOTerms=as.data.frame(fread('allGOterms_EnsGRC37_13042017.txt'))
    }
    if(is.numeric(geneList) | is.logical(geneList)){
        geneList=background[geneList]
        }
    DE=background%in%geneList
    names(DE)=background
    if(is.null(bias.data)){
        nullP=nullp(DE, 'hg19', 'ensGene',plot.fit=TRUE)
    }else{
            nullP=nullp(DE, 'hg19', 'ensGene', bias.data=bias.data,plot.fit=TRUE)
    }
    if(!biasCorrect){
        nullP$pwf=mean(DE)
        }
    res=goseq(nullP, 'hg19', 'ensGene')
    resUp=res
    resUp$FDR=p.adjust(res$over_represented_pvalue,'fdr')
    resUp=resUp[resUp$FDR<FDR,]
    resUp$Term=sapply(mget(resUp[,1],GOTERM),function(x){x@Term})
    resUp$FoldEnrich=resUp[,4]/resUp[,5]/mean(DE)
    if(!overOnly){
        resDn=res
        resDn$FDR=p.adjust(res$under_represented_pvalue,'fdr')
        resDn=resDn[resDn$FDR<FDR,]
resDn$Term=sapply(mget(resDn[,1],GOTERM),function(x){x@Term})
        resDn$FoldEnrich=resDn[,4]/resDn[,5]/mean(DE)
        res=rbind(resUp,resDn)
res$Pvalue=ifelse(res$FoldEnrich>1,res$over_represented_pvalue,res$under_represented_pvalue)
    }else{
        res=resUp
        res$Pvalue=res$over_represented_pvalue
    }
    if(addGenes){
        GeneInGO=allGOterms[allGOterms$gene%in%geneList & allGOterms$go%in%res[,1] ,]
GeneInGO=By(GeneInGO$gene,GeneInGO$go,function(x){paste(sort(unique(G2S(x))),collapse=' // ')})
        res$genes=GeneInGO[match(res[,1],names(GeneInGO))]
    }
    if(addCI){
    	if(nrow(res)>0){
    	res$nbInGrp=rep(length(geneList),nrow(res))
    	res$nbInBckgd=rep(length(background),nrow(res))
        counts2tab = function(xy,x,y,tot){matrix(c(tot-x-y+xy,x-xy,y-xy,xy),2)}
    	OR=sapply(as.data.frame(t(mapply(function(xy,x,y,tot){odds.ratio(counts2tab(xy,x,y,tot))},res$numDEInCat,res$numInCat,length(geneList),res$nbInBckgd))),as.numeric)
    	if(is.matrix(OR)){
    		OR=as.data.frame(OR)
    	}else{
			nn=names(OR)
    		OR=as.data.frame(matrix(OR,1))
    		colnames(OR)=nn
    	}
    	
    	OR$CI=paste('[',round(OR[,'LowerCI'],1),'-',round(OR[,'UpperCI'],1),']',sep='')
		res$OR=OR$OR
		res$CI=OR$CI
		res$lowerCI=OR$LowerCI
    		}
		}    
    res
    }