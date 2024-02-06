rm(list=ls())
library('WGCNA')
setwd('/Users/lichaoran/dataprocess/2016-nature-HE/WGCNA')
##exp matrix
data=read.csv("/Users/lichaoran/dataprocess/2016-nature-HE/WGCNA/he201655_tpm_5_2fc_14_wgcna.csv",
              row.names = 1,check.names = F)
data=log2(data+1)
write.csv(data,"/Users/lichaoran/dataprocess/2016-nature-HE/WGCNA/log-he201655_tpm_5_2fc_14_wgcna.csv")
options(stringsAsFactors = FALSE)
dim(data)
names(data)
head(data)
##extract
datExpr0=as.data.frame(t(data))
dim(datExpr0)
head(datExpr0)
#check missing value
gsg=goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
#clustering
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main=2)

#delete outliers
abline(h = 15,col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 15,minSize = 10)
table( clust)



##read trait data
traitData=read.csv("29-sam-traitdata.csv",
                   row.names = 1,check.names = F)
row.names(traitData)=rownames(datExpr0)
sameSample=intersect(rownames(datExpr0), rownames(traitData))
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]
dev.off()
sampleTree2 = hclust(dist(datExpr0), method = "average")
plot(sampleTree2)

traitColors = numbers2colors(datTraits, signed = FALSE)

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datExpr0, datTraits, file = "/Users/lichaoran/dataprocess/2016-nature-HE/wgcna-new/01-HE-55-29-dataInput-chatrait.RData")

##build network, identify modules
rm(list = ls())
load(file = "01-HE-55-29-dataInput-chatrait.RData")
disableWGCNAThreads() 
#datExpr0=as.data.frame(t(datExpr0))
#datExpr1 <- datExpr0[order(apply(datExpr0,1,mad), decreasing = T)[1:5000], ]
#datExpr1=as.data.frame(t(datExpr1))
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers,verbose = 5,blockSize = 5000)
if (is.na(power)){ 
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18), 
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16), 
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14), 
                               ifelse(type == "unsigned", 6, 12))        
                 ) 
  ) 
}
sft$powerEstimate
sizeGrWindow(9, 5)
par(mfrow = c(1,2)) 
cex1 = 0.9 #字符大小
##Scale-free Topological Fit Index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab=" Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2" , type="n" ,
     main = paste(" Scale independence" ));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1 , col="red");
abline(h=0.85, col="red") #查看位于0.9以上的点

##mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)" ,ylab="Mean Connectivity", type="n" ,
     main = paste( "Mean connectivity" ))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


cor =WGCNA::cor 
net =blockwiseModules(datExpr0,  maxBlockSize = 5000,power = sft$powerEstimate,
                      TOMType = "unsigned", minModuleSize = 30,
                      reassignThreshold = 0,mergeCutHeight=0.2,
                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                      saveTOMs = TRUE,
                      deepSplit=4,
                      saveTOMFileBase = "55-29samTOM",
                      verbose = 3)
cor = stats::cor #改回去



#view the number of modules and the number of genes
table(net$colors)


sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang=0.03,
                    addGuide = TRUE,guideHang = 0.05)

##Export the gene list for each module
mog=data.frame(Mog=net[["colors"]])
mog=cbind(rownames(mog),mog)
colnames(mog)=c('geneid','md')
b=split(mog$geneid,mog$md)

library(tidyr)
mc=cbind(rownames(mog),data.frame(moduleColors))
colnames(mc)=c('geneid','mc')
mc1=split(mc$geneid,mc$mc)
##Cbind
mc2 <- data.table::rbindlist(lapply(mc1, function(x) data.table::data.table(t(x))),fill = TRUE) %>% t() %>% 
  data.frame(row.names = seq(1:max(lengths(mc1)) )) %>% 
  purrr::set_names(names(mc1))
write.csv(mc2,"he201629_difstages_genelist.csv")

# Save the assigned modules and information about genes contained in each
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms [[1]];
save(MEs,moduleLabels, moduleColors, geneTree,
     file = "02-5529-ME.RData")


lnames = load(file = "01-HE-55-29-dataInput-chatrait.RData")
lnames 
lnames = load(file = "02-5529-ME.RData")
lnames
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p"); #Computed the correlation for each module and each genes
moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nSamples);

sizeGrWindow(8,10)
#Display the correlation coefficients and P-values between module
textMatrix = paste(signif(moduleTraitCor, 2),"\n(",
                   signif(moduleTraitPvalue, 1), ")",sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(7, 8.5, 4, 4));
#Show the correlation coefficients using a heatmap.
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs) ,
               colorLabels = FALSE ,
               colors = blueWhiteRed(50),
               xColorWidth=0.1,
               textMatrix = textMatrix,
               setStdMargins = FALSE ,
               cex.text = 0.5,
               cex.lab.x = 0.5,
               cex.lab.y = 1,
               cex.lab = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

library(eoffice)
library(ggpubr)
#topptx(wgheatmap,"module-trait.pptx")
unique(moduleColors)
which.module="red"
names(datExpr0)[moduleColors=="red"]
df=t(scale(datExpr0[,moduleColors==which.module ]))

transition = as.data.frame(datTraits[5]);
names(transition) = "lbph" ;
modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
#Analysis of the correlation between genes and each module.
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples));
names(geneModuleMembership) = paste( "MM" , modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, transition, use = "p")) ;#和中期过度关联
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(transition), sep="");
names(GSPvalue)= paste("p.GS.", names(transition), sep="");
#########
ADJ=abs(cor(datExpr0,use="p"))^6
Alldegrees =intramodularConnectivity(ADJ, moduleColors)
write.csv(Alldegrees, file = "intramodularConnectivity-1.csv")
##################
colorlevels=unique(moduleColors)
pdf("GS vs. degree_weight-1.pdf",width = 14,height = 8)
par(mfrow=c(4,5))
par(mar = c(5,5,3,3))
for (i in c(1:length(colorlevels))) 
{
  whichmodule=colorlevels[[i]]; 
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees$kWithin[restrict1], 
                     geneTraitSignificance[restrict1,1], 
                     col=moduleColors[restrict1],
                     main=whichmodule, 
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
dev.off()
# KME
datKME=signedKME(datExpr0, MEs, outputColumnName="kME_MM.")
write.csv(datKME, "kME_MM.csv")
# grey60
 datKME=signedKME(datExpr0, MEs, outputColumnName="kME_MM.")
 module = "grey60"
 column = match(module, modNames)
 moduleGenes = moduleColors==module
 grey60_module<-as.data.frame(dimnames(data.frame(datExpr0))[[2]][moduleGenes])
 names(grey60_module)="genename"
 grey60_KME<-as.data.frame(datKME[moduleGenes,column]) 
 names(grey60_KME)="KME"
 rownames(grey60_KME)=grey60_module$genename
 FilterGenes = abs(grey60_KME$KME) > 0.8
 table(FilterGenes)
FilterGenes
#FALSE  TRUE 
#24    39 
 module = "grey60"
 column = match(module, modNames)
 moduleGenes = moduleColors==module
 grey60_module<-as.data.frame(dimnames(data.frame(datExpr0))[[2]][moduleGenes])
 names(grey60_module)="genename"
 MM<-abs(geneModuleMembership[moduleGenes,column])
 GS<-abs(geneTraitSignificance[moduleGenes, 1])
 grey60_MMGS<-as.data.frame(cbind(MM,GS))
 rownames(grey60_MMGS)=grey60_module$genename
 hub_b<-abs(grey60_MMGS$MM)>0.8&abs(grey60_MMGS$GS)>0.2
 table(hub_b)
 grey60_hub_b<-subset(grey60_MMGS, abs(grey60_MMGS$MM)>0.8&abs(grey60_MMGS$GS)>0.2)
 write.csv(grey60_hub_b, "hubgene_MMGS_grey60.csv")
 #########hub
HubGenes <- chooseTopHubInEachModule(datExpr0,moduleColors)
write.csv (HubGenes,file = "TopHubGenes_of_each_module.csv",quote=F)
#GS和MM(grey60)
module = "grey60"
column = match(module,modNames);
moduleGenes= moduleColors == module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for transition",
                   main = paste("Module membership vs. gene significance\n"),
                   ceх.main=1.2,cex.lab = 1.2,cex.axis =1.2 ,col = module)
