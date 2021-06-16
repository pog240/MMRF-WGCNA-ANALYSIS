setwd("D:/Anoval_Analysis")
#load("/Users/oadebayo/Documents/Olayinka_New_Analysis/saved.image.BM.noOutliers.lmm.regressedAgeSex_StageProtected21.Rdata")
load("D:/Anoval_Analysis/ola_21Modules_april_22_21.RData")
#traits<-read.csv("traits270x16_rvd.csv", header = TRUE, row.names = 1)
#cleanDat<-read.csv("/Users/oadebayo/Documents/Olayinka_New_Analysis/cleadat35187x732.csv", header = TRUE, row.names = 1, check.names = FALSE)
#numericMeta<-traits
#head(cleanDat)
#head(traits)
dim(cleanDat)
dim(numericMeta)
rownames(traits)==colnames(cleanDat) #sanity check -- are sample names in same order?
# If the sanity test is TRUE no need to run the nest line of code
#cleanDat<-cleanDat[,na.omit(match(rownames(numericMeta),colnames(cleanDat)))] #cull to keep only samples we have traits for, and match the column (sample) order of cleanDat to the row order of numericMeta

#numericMeta <- numericMeta[match(colnames(cleanDat),rownames(numericMeta)),] #use this line instead of above if you have more samples in your traits file than you do in abundance data; trait sample (row) order will be matched to column names of cleanDat; or use both lines if matching needs to be enforced due to different missing samples in both traits and abundance data
#rownames(numericMeta)==colnames(cleanDat) #sanity check -- are sample names in same order?
dim(cleanDat) # check diamension of cleansat
## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
sdout=3
install.packages("WGCNA")
install.packages("doParallel")
install.packages("biomaRt")
install.packages("NMF")
install.packages("plotly")
install.packages("stringr")
install.packages("cluster")
install.packages("boot")
install.packages("Cairo")
install.packages("callr")

library(WGCNA)
library(NMF)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(Cairo)
library(doParallel)
library(biomaRt)
library(NMF)
library(plotly)
library(stringr)
library(cluster)
library(boot)
library("callr")

#powers <- seq(2,12,by=1)
#sft <- pickSoftThreshold(t(cleanDat),blockSize=nrow(cleanDat)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         #powerVector=powers,
                         #corFnc="bicor",networkType="signed")
#jpeg(file="oosavingsah_plottablesft12.jpeg")
#tableSFT<-sft[[2]]
#plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R?")
#dev.off()

#powers <- seq(6.5,20,by=0.5)
#sft2 <- pickSoftThreshold(t(cleanDat),blockSize=nrow(cleanDat)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         #powerVector=powers,
                         #corFnc="bicor",networkType="signed")
# power is set to 7 base on the 
#power=7
#net <- blockwiseModules(t(cleanDat),power=power,deepSplit=1.5,minModuleSize=180,
                        #mergeCutHeight=0.15,TOMdenom="mean", #detectCutHeight=0.9999,                        #TOMdenom="mean" may get more small modules here.
                        #corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        #..........verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDat)+1000,reassignThresh=0.05)
nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))
net.ds2<-net
FileBaseName=paste0("MyNetworkDescription_power_",power,"_PAMstageTRUE")
CairoPDF(file="05_6_2021.GlobalNetworkPlots-FileBaseName_21.pdf",width=16,height=12)
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

numericIndices<-unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))
#Warnings OK; This determines which traits are numeric and if forced to numeric values, non-NA values do not sum to 0

geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDat),use="pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]

plotDendroAndColors(dendro=net$dendrograms[[1]],
                    colors=t(rbind(net$colors,geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=TRUE,
                    dendroLabels=FALSE,
                    groupLabels=c("Module Colors",colnames(numericMeta)[numericIndices]))
head(MEs)
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL
plotEigengeneNetworks(tmpMEs, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))
head(numericMeta)
MM_Vital_Status<-numericMeta$MM_Vital_Status
MM_Vital_Status[numericMeta$MM_Vital_Status==1]<-"death"  #only necessary if Group was numerically encoded; does nothing if the Grouping vector has no numeric values
MM_Vital_Status[numericMeta$MM_Vital_Status==0]<-"alive"

#ethnicity<-numericMeta$ethnicity
#ethnicity[numericMeta$ethnicity==1]<-"European American"
#ethnicity[numericMeta$ethnicity==2]<-"African American"
#ethnicity[numericMeta$ethnicity==3]<-"other"

#OS<-numericMeta$OS
#PFS<-numericMeta$PFS
#TUMOR_TYPE[numericMeta$TUMOR_TYPE==1]<-"PAAT"  #there is one of these lines for each numeric value set in the Group column of the traits.csv file
#TUMOR_TYPE[numericMeta$TUMOR_TYPE==2]<-"PAOS"
head(numericMeta)
regvars <- data.frame(as.factor( numericMeta$MM_Vital_Status )) #, as.numeric(numericMeta$ethnicity), as.numeric(numericMeta$OS),as.numeric(numericMeta$PFS))
colnames(regvars) <- c("MM_Vital_Status")
lm1 <- lm(data.matrix(MEs)~MM_Vital_Status,data=regvars)
pvec <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec) <- colnames(MEs)
kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")
library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs,numericMeta[,numericIndices])
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p
textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));
cexy <- if(nModules>75) { 0.8 } else { 1 }
colvec <- rep("white",1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
labeledHeatmap(Matrix = apply(moduleTraitPvalue,2,as.numeric),
               xLabels = colnames(numericMeta)[numericIndices],
               yLabels = paste0("ME",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colvec,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               cex.lab.y= cexy,
               zlim = c(0,0.15),
               main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Patients correlation p value"),
               cex.main=0.8)
numericMetaCustom<-numericMeta[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p
moduleTraitPvalue<-signif(moduleTraitPvalue, 1)
moduleTraitPvalue[moduleTraitPvalue > as.numeric(0.05)]<-as.character("")
textMatrix = moduleTraitPvalue; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
                                #textMatrix = gsub("()", "", textMatrix,fixed=TRUE)
labelMat<-matrix(nrow=(length(names(MEs))), ncol=2,data=c(rep(1:(length(names(MEs)))),labels2colors(1:(length(names(MEs))))))
labelMat<-labelMat[match(names(MEs),labelMat[,2]),]
for (i in 1:(length(names(MEs)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
for (i in 1:length(names(MEs))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }
xlabAngle <- if(nModules>75) { 90 } else { 45 }
par(mar=c(16, 12, 3, 3) )
par(mfrow=c(1,1))
bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))
colvec<-c(bw(50),wr(50))
labeledHeatmap(Matrix = t(moduleTraitCor)[,],
               yLabels = colnames(numericMetaCustom),
               xLabels = labelMat[,2],
               xSymbols = labelMat[,1],
               xColorLabels=TRUE,
               colors = colvec,
               textMatrix = t(textMatrix)[,],
               setStdMargins = FALSE,
               cex.text = 1.2,
               cex.lab.x = cexy,
               xLabelsAngle = xlabAngle,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(-1,1),
               main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
               cex.main=0.8)
toplot <- MEs
colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot)
pvec <- pvec[match(names(pvec),rownames(toplot))]
#rownames(toplot) <- paste(rownames(toplot),"\np = ",signif(pvec,2),sep="")
rownames(toplot) <- paste(orderedModules[match(colnames(MEs),orderedModules[,2]),1]," ",rownames(toplot),"  |  t-test p=",signif(pvec,2),sep="")
# add any traits of interest you want to be in the legend
#ethnicity<-numericMeta$ethnicity
#ethnicity[numericMeta$ethnicity==1]<-"European American"
#ethnicity[numericMeta$ethnicity==2]<-"African American"
#ethnicity[numericMeta$ethnicity==3]<-"other"

#OS<-numericMeta$OS
#PFS<-numericMeta$PFS
metdat=data.frame(MM_Vital_Status=MM_Vital_Status)#,ethnicity=ethnicity, OS=OS,PFS=PFS)

# set colors for the traits in the legend
heatmapLegendColors=list('MM_Vital_Status'=c("midnightblue","red"), #,"seagreen3","hotpink","purple"),
                         
                         #'ethnicity'=c("pink","dodgerblue","seagreen3","hotpink"),
                         #'Age'=c("white","darkgreen"), #young to old
                         #'Gender'=c("pink","dodgerblue"), #F, M
                         'Modules'=sort(colnames(MEs)))
library(NMF)
par(mfrow=c(1,1))
aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES IN ORIGINAL, e.g. BATCH OR REGION ORDER",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=1.2, ## Character sizes
         cexCol=1.2,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=NA) ## Do not cluster columns - keep given order
aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE,Colv=TRUE,) ## Cluster columns
colnames(numericMeta)[numericIndices]
table(MM_Vital_Status)
par(mfrow=c(4,6)) #rows,columns on each page (try to choose to keep all plots for each module on one page or row(s), without squeezing too much in)
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(MM_Vital_Status,c("death","alive")),col=colnames(MEs)[i],ylab="Eigengene Value",main=rownames(toplot)[i],xlab=NULL,las=2)
  #boxplot(toplot[i,]~factor(ethnicity,c("white","black","other")),col=colnames(MEs)[i],ylab="Eigengene Value",main=rownames(toplot)[i],xlab=NULL,las=2)#you choose the order of groups for boxplots
  #verboseScatterplot(x=numericMeta[,"ethnicity"],y=toplot[i,],xlab="race (1=white; 2=blac; 6=other;)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  verboseScatterplot(x=numericMeta[,"MM_Vital_Status"],y=toplot[i,],xlab="progrssion (0=no; 1=yes)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"days.to.death"],y=toplot[i,],xlab="Days to Death",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"gender"],y=toplot[i,],xlab="Gender (1=male; 2=female)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"TUMOR_TYPE"],y=toplot[i,],xlab="TUMOR_TYPE at Sampling (1-2)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"in_treatment"],y=toplot[i,],xlab="Number of Days in Treatment",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"NoTreatment"],y=toplot[i,],xlab="Treatment: None (--)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  for (i in 1:7) frame()
}
dev.off()
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(cleanDat)),rownames(cleanDat),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
write.csv(kMEtable,file="ola21Module_2021.ModuleAssignments-.csv")
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
save.image('ola_21Modules_april_22_21.RData')




save.image("ola27modules.RData")
