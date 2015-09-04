## version dedicated to mouse / Affy Gene ST arrays
## Author: Fathi Elloumi
## updated by Bong-Hyun Kim
## version 2

#try( install.packages( 'getopt' ) )
#library( getopt )

args <- commandArgs(TRUE)
outdir <- args[3]
if ( !file.exists(outdir) ) {
	dir.create(outdir) ;
}

project <- args[1]
inputdir <- normalizePath(args[2])
outdir <- normalizePath(outdir)
phenofile <- normalizePath(args[4])
contrasfile <- normalizePath(args[5])
species <- "mo" #args[6] #"mo" or "hu" 

spl1 <- paste0( "pd.",species,"gene.2.0.st" ) 
spl2 <- paste0( species,"gene20sttranscriptcluster.db" )
#library(spl2, character.only=TRUE)

if (!file.exists(inputdir) ){
	exit(-1) ;
}

#install & load libraries
libs = c( "limma", "gplots", "oligo", "geneplotter", "topGO", "Rgraphviz", spl1, spl2  )
for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    cat( paste( c( "installing", i) )  ) ;
    cat( "\n"  )
    install.packages(i, repos='http://watson.nci.nih.gov/cran_mirror/')
    #cat( "install should be run... Not implemented yet!\n" ) ;
  }
  library(i,character.only = TRUE)
}

marray_pipeline_mm_affy2st= function(project,inputdir,outdir,phenofile,contrasfile) {
 
## load the packages for Affy Gene ST arrays
# 
#library(limma)
#library(gplots)
# load the oligo library
#library(oligo)
#library(geneplotter)

#speicies specific library
#library(pd.mogene.2.0.st)
#library(st_library, character.only=TRUE)


## parameters

# project="CCBR402"
# file for phenoData, 1 group of intersest SampleGroup, SampleName, SampleID
# file for contrasts


## take input ====================================================================================

setwd(inputdir)

## create a GeneFeatureset object ========================================

# Read in the CEL files in the directory
cels <- list.celfiles()
print( cels )
pd<-read.AnnotatedDataFrame(phenofile,header=TRUE,row.name="SampleName" ,sep="\t")
celfiles <- read.celfiles(cels, phenoData=pd)
##
	facs <- factor(pData(celfiles)$SampleGroup)
	labfacs=levels(facs)
	nbfacs=length(labfacs)
	
# read contrasts
# need to test that the groups belong to levels !!!

	
	contra=read.delim(contrasfile)
	nb=dim(contra)[1]
	cons=c()
	for (k in 1:nb) {
		if ((contra[k,1] %in% labfacs) & (contra[k,2] %in% labfacs) )
		{ 
			cons=c(cons,paste(contra[k,1],"-",contra[k,2],sep="")) 
		}
		else 
		{
			cat("One of the groups in contrasts file at line :",k+1,"does not match a group in phenotype file..Quitting!!!\n")
			print( contra )
			return (-1)
		}
    }

## outputs ====================================================================================
cat("Raw data QC Plots ......\n ")	
setwd(outdir)
	
## raw data plots


	
for (i in 1:nbfacs)

{
  #  MAplot
  # MA plots are then used to visualize intensity-dependent ratio for each group
  igp=which(pData(celfiles)$SampleGroup==labfacs[i])
  png(paste(project,"-MAplot-",labfacs[i],"-rawdata.png",sep=""),width=1000,height=1000)
  MAplot(celfiles[,igp],pairs=TRUE,plotFun=smoothScatter,main="MVA plot") # 
  dev.off()
  cat("group : ",i, "plot done\n")
}



png(paste(project,"-histog-rawdata.png",sep=""))
hist(celfiles,which="all", main =" Samples distribution")
dev.off()


# boxplot(celfiles[,1:2], col=c(3,2)) # default both
png(paste(project,"-boxplot-rawdata.png",sep=""),width=1000,height=700)
par(mar=c(8, 4, 4, 8) + 0.1)
boxplot(celfiles, col=c(3,2),which="all", main="Boxplots before normalization",las=2)
dev.off()



## quality control =================================================

cat("  RLE and NUSE Plots ...\n ")
# RLE and NUSE plots / GOOD for Oligo
celfiles.qc <- fitProbeLevelModel(celfiles)

png(paste(project,"-rle-plot.png",sep=""),width=1000,height=700)
par(mar=c(8, 4, 4, 8) + 0.1,las=2)
RLE(celfiles.qc, main="RLE")
dev.off()

png(paste(project,"-nuse-plot.png",sep=""),width=1000,height=700)
par(mar=c(8, 4, 4, 8) + 0.1,las=2)
NUSE(celfiles.qc, main="NUSE")
dev.off()

cat("  Probset presence Plot ...\n ")

# present Absent Calls
# PAp = paCalls(celfiles[, 1:2])
PApset = paCalls(celfiles,"PSDABG")

#head(PApset)
 
PApsetNC <- dim(PApset)[2]
PApsetNR <- dim(PApset)[1]
# 269751     43
# present probes set

pres=c()
for ( i in 1:PApsetNC ) 
  pres=c(pres,length(which(PApset[,i]<0.05)))
names(pres)=colnames(PApset)

pres=round((pres/PApsetNR)*100,2)

# hist(pres)
# plot(density(pres))

png(paste(project,"-probesetcalls.png",sep=""),width=1000,height=1000)
par(mar=c(8, 4, 4, 8) + 0.1)
plot(pres, type="l",ylab="Presence percentage", main ="Probsets presence",xaxt="n", xlab="")
axis(1,at=c(1:PApsetNC),labels=names(pres),las=2,cex.lab=0.4)
dev.off()

## QC using arrayQualityMetrics ==========================================
library (arrayQualityMetrics)

# for R3.1.2 need gridSVG 1.4-3

	cat("  ArrayQualityMetrics for raw data...\n ")
#converting raw expression matrix to ExpressionSet object
esetR=ExpressionSet(assayData=exprs(celfiles),phenoData=pd,annotation="pd.mogene.2.0.st")

arrayQualityMetrics(esetR, intgroup="SampleGroup", outdir = "QC_raw", force =TRUE, do.logtransform = TRUE)

cat("Normalization and post-Normalization plots ......\n ")
	
# QC after normalization ===================================================

# new
# genePS <- rma(celfiles, target = "probeset") NO 
# geneCore <- rma(celfiles, target = "core")

celfiles.rma =rma(celfiles, background=TRUE, normalize=TRUE, subset=NULL, target="core")

# object.size(celfiles.rma)
# [1] 2680904
# object.size(celfiles)
# [1] 145553672

## some plots
png(paste(project,"-boxplot-rma.png",sep=""),width=1000,height=1000)
par(mar=c(8, 4, 4, 8) + 0.1)
boxplot(celfiles.rma,col=c(3,2), main="Boxplots after RMA normalization",las=2)
dev.off()

# title is not changing ??
# MA plots for groups after RMA
for (i in 1:nbfacs)

{
#  MAplot
	igp=which(pData(celfiles.rma)$SampleGroup==labfacs[i])
	png(paste(project,"-MAplot-",labfacs[i],"-rma.png",sep=""),width=1000,height=1000)
	MAplot(celfiles.rma[,igp],pairs=TRUE,plotFun=smoothScatter,main="MVA plot") # 
	dev.off()
	cat("group : ",i, "plot done\n")
}


# trad.scatter.plot(exprs(celfiles.rma)[,1],exprs(celfiles.rma)[,2],fc.line.col="lightblue",col="blue")
cat("  ArrayQualityMetrics after normalization ...\n ")
arrayQualityMetrics(expressionset = celfiles.rma, intgroup=c("SampleGroup"), outdir = "QC_normalized", force =TRUE)
# dev.off()

cat("  PCA and clustering plots ...\n ")
## pca , clustering and heatmap ====================================
# histo
	
	png(paste(project,"-Histo-rma.png",sep=""),width=1000,height=1000)
	hist(celfiles.rma, main="Distribution after Normalization")
	dev.off()
#pca
# 
pr1=prcomp(t(exprs(celfiles.rma)),scale=T)
summary(pr1)
pr1$x # scores for the samples
(pr1$sdev)**2 ## variances for all , equiv to eigen values
head(pr1$rotation[,1:2]) # loading for pc1 and pc2 (pc vs genes)
# pc are the eigen vectors

igp1=which(pData(celfiles.rma)$SampleGroup==labfacs[1])
xrange=range(pr1$x[,1])
yrange=range(pr1$x[,2])

pc1.var=100*round(((pr1$sdev)**2)[1]/sum((pr1$sdev)**2),digits=2) # %var pc1 
pc2.var=100*round(((pr1$sdev)**2)[2]/sum((pr1$sdev)**2),digits=2) # % var pc2

xlab=paste("PC1 - ",pc1.var," % of variation",sep="")
ylab=paste("PC2 - ",pc2.var," % of variation",sep="")

png(paste(project,"-pca-rma.png",sep=""),width=1000,height=1000)
plot(pr1$x[igp1,1],pr1$x[igp1,2],col=1,main="prcomp after scaling",ylab=ylab,xlab=xlab,xlim=xrange,ylim=yrange)
for (i in 2:nbfacs)
{
	igp=which(pData(celfiles.rma)$SampleGroup==labfacs[i])
	points(pr1$x[igp,1],pr1$x[igp,2],col=i)

}
legend('bottom', labfacs,lty=1, col=1:nbfacs, bty='n', cex=.6)
dev.off()

#hclust

eset <- exprs(celfiles.rma) # now a matrix
colnames(eset)=paste(pData(celfiles)$SampleID,"_",pData(celfiles)$SampleGroup,sep="")
# colnames(eset)=pData(celfiles.rma)$barcode.s.

# standardize before clustering
eset.s=scale(eset)
distance.s <- dist(t(eset.s))
clusters.s <- hclust(distance.s)
d2=as.dist(1-cor(eset.s))

png(paste(project,"-HC-EuclPears-Scale-rma.png",sep=""),width=1000,height=1000)
par(mfrow=c(1,2))
plot(clusters.s,main="Euclidean with scale")
# rect.hclust(clusters.s, k=2, border="red")
plot(hclust(d2), main="Pearson Correlation with scale")
# rect.hclust(hclust(d2), k=2, border="red")
dev.off()

png(paste(project,"-HC-Eucl-Scale-rma.png",sep=""),width=1000,height=1000)
plot(clusters.s,main="Euclidean with scale")
dev.off()
	
png(paste(project,"-HC-Pearson-Scale-rma.png",sep=""),width=1000,height=1000)
plot(hclust(d2), main="Pearson Correlation with scale")
dev.off()
	
## heatmap
#
distylog2=dist(t(eset))
mat = as.matrix(distylog2)
heatmap.2(mat, trace="none", margin=c(10,10))
# 
dev.copy(png,paste(project,"-heatmap_samplebysample-noscale.png",sep=""),width=1000,height=1000)
dev.off()

cat("Limma, DEG and MAplot ......\n ")

## DE probesets ==============================================
# ------------------------------------------------------------

#library(limma)
# limma uses an empirical Bayes method to moderate the standard errors of the estimated log-fold changes. This results in more stable
# inference and improved power, especially for experiments with small numbers of arrays




myfactor <- factor(pData(celfiles.rma)$SampleGroup)
design1 <- model.matrix(~0+myfactor)
colnames(design1) <- levels(myfactor)

fit1 <- lmFit(celfiles.rma,design1)
contrast.matrix <- makeContrasts(contrasts=cons,levels=design1)

fit2 <- contrasts.fit(fit1, contrast.matrix)
ebayes.fit2=eBayes(fit2) # smooths the std error

# #EXTRACTING ALL GENES FOR EACH CONTRAST

##ANNOTATE PROBESET IDS FROM ANNOTATION PACKAGE FROM BIOCONDUCTOR
## load libraries as sources of annotation

#library(mogene20sttranscriptcluster.db)
Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))

res=matrix(0,nb,3)
colnames(res)=c("adjPvalue<0.05","+abs(fc)>=1","+abs(fc)>=2")
rownames(res)=cons

for (i in 1:nb)
{

   all.genes.con = topTable(ebayes.fit2, coef = i, number=nrow(ebayes.fit2))

  ## generate plotMA from geneplotter
	
  dataf=data.frame("m"=all.genes.con$AveExpr,"fc"=all.genes.con$logFC,"sig"=all.genes.con$adj.P.Val<0.05)
  png(paste(project,"-MAplot-",cons[i],".png",sep=""),width=1000,height=700)
  plotMA(dataf,main=cons[i])
  dev.off()
	
  # name=paste("all.genes.con",i,sep="")
  res[i,1]=length(which(all.genes.con$adj.P.Val<0.05))
  res[i,2]=length(which(all.genes.con$adj.P.Val<0.05 & abs(all.genes.con$logFC)>=1))
  res[i,3]=length(which(all.genes.con$adj.P.Val<0.05 & abs(all.genes.con$logFC)>=2))

  # Merge data frames together (like a database table join)

  all <- merge(all.genes.con, Annot,by.x=0, by.y=0, all.x=T)
  all=all[order(all$P.Value),]
  colnames(all)[1]="probsetID"
	

  # Write out to a file:
  write.table(all,file=paste(project,"-",cons[i],"_all_genes.txt",sep=""),sep="\t",row.names=F)
  cat("Contrast: ",i," done \n")

  #GO BP Term Analysis
  sampleGOdata <- new ( "topGOdata", description=paste0(project," GO BP Term Analysis"), 
                        ontology="BP", allGenes=geneList, geneSel=topDiffGenes, 
                        nodeSize=10, annot=annFUN.db, affyLib="mogene20sttranscriptcluster.db")
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
  topBP <- GenTable(sampleGOdata, classicFisher = resultFisher, 
                     classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", 
                     ranksOf = "classicFisher", topNodes = 30)
  write.table( topBP, paste0(project,"_",cons[i],"_GO_BP_top30.txt"))
  pdf( paste0(project,"_",cons[i],"_GO_BP_top.pdf") )
  showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo ="all")
  dev.off()
  
  #GO MF Term Analysis
  sampleGOdata <- new ( "topGOdata", description=paste0(project," GO MF Term Analysis"), 
                        ontology="MF", allGenes=geneList, geneSel=topDiffGenes, 
                        nodeSize=10, annot=annFUN.db, affyLib="mogene20sttranscriptcluster.db")
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
  topMF <- GenTable(sampleGOdata, classicFisher = resultFisher, 
                    classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", 
                    ranksOf = "classicFisher", topNodes = 30)
  write.table( topMF, paste0(project,"_",cons[i],"_GO_MF_top30.txt"))
  pdf( paste0(project,"_",cons[i],"_GO_MF_top.pdf") )
  showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo ="all")
  dev.off()
  
  #GO CC Term Analysis
  sampleGOdata <- new ( "topGOdata", description=paste0(project," GO CC Term Analysis"), 
                        ontology="CC", allGenes=geneList, geneSel=topDiffGenes, 
                        nodeSize=10, annot=annFUN.db, affyLib="mogene20sttranscriptcluster.db")
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
  topCC <- GenTable(sampleGOdata, classicFisher = resultFisher, 
                    classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", 
                    ranksOf = "classicFisher", topNodes = 30)
  write.table( topCC, paste0(project,"_",cons[i],"_GO_CC_top30.txt"))
  pdf( paste0(project,"_",cons[i],"_GO_CC_top.pdf") )
  showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo ="all")
  dev.off()
## end for
}
write.table(res,paste(project,"_deg_summary.txt",sep=""),sep="\t",col.names=NA)

return (0)

}


#call the main function
marray_pipeline_mm_affy2st(project,inputdir,outdir,phenofile,contrasfile)

#Fathi's example
#marray_pipeline_mm_affy2st("CCBR-402","/Users/elloumif/Documents/project402","/Users/elloumif/Documents/project402/results","ccbr402Sampledata.txt","contrastsList.txt")

