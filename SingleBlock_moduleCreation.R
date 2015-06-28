

### This is soft-thresholding and single block module creation 
setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_lnc")

### Load packages
library(WGCNA)
library(permute)
library(lattice)
library(vegan)
library(scatterplot3d)
library(MASS)
library(nlme)
library(mgcv)
library(cluster)
library(rgl)
library(ape)
library(picante)
library(gdata)
library(WriteXLS)
library(plyr)
library(caret)
library(BiodiversityR)
library(gtools)
library(AppliedPredictiveModeling)
library(limma)
library(vegetarian)
library(survival)
library(randomForest)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(flashClust)
transparentTheme(trans = 0.4)

#### Following on from file prep 

#load(file="lncRNA_total_input.RData")
#load(file = "lncRNA_trimmed_input.RData")
load(file = "lncRNA_top10K_input.RData")

#### Set Pam50 factor again if need

Pam50 <- as.factor(datClin$PAM50Call_RNAseq)

################
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to=18, by=2))
# Call the network topology analysis function

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 3, 
                        blockSize = 5000, networkType = "signed", moreNetworkConcepts = TRUE)


# Plot the results:
sizeGrWindow(12,9)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red");

# this line corresponds to using an R^2 cut-off of 0.9
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


##################################################

##To check the efficiency of the chosen soft power, put it in here and plot
### Might try values as indicated by above (11 & then 12)

k=softConnectivity(datE=datExpr,power=12, type = "signed")

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

### After checking - use value of x for softthresholding
####  Proceed to module creation using automatic, dynamic cutoffs. 

softPower = 12;
adjacency = adjacency(datExpr, power = softPower, type = "signed");
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM

######## Make the tree
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

###### Now we're going to go ahead with module creation
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.3
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


##### Now we'll make a dataframe to store the moduleColours for genes, 
### We'll keep adding to this baby later
### Match the names in the blocks to the colours

#### 
Modules <- data.frame(colnames(datExpr), moduleColors)
lookfor <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "EZH2", "DNMT3L", "MTHFR")

Interest_mods <- Modules[Modules$colnames.datExpr %in% lookfor,] 


# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "Dynamic_merged_lncRNA.RData")

staticCut <- as.character(cutreeStaticColor(geneTree, cutHeight=.99, minSize=50))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(geneTree, 
                    colors = data.frame( moduleColors,staticCut),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")
table(staticCut)




##################

# Rename to moduleColors
moduleColors = staticCut
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;



MEList = moduleEigengenes(datExpr, colors = moduleColors)
MEs = MEList$eigengenes
AEs = MEList$averageExpr

##Bind MEs and AEs together

ModExp <- cbind(MEs, AEs)
## add rownames of samples
rownames(ModExp) <- rownames(datExpr)
### Write out file
write.table(ModExp, "ModuleInfo_static.txt", sep ="\t")


save(MEs, moduleLabels, moduleColors, geneTree, file = "Static_lncRNA.RData")



############## Create a data frame of just the modules of interest, to see other members

BlackMod <- subset(Modules,moduleColours == "black")

GreenMod <- subset(Modules, moduleColours == "green")

MidBlue <- subset(Modules, moduleColours == "midnightblue")

##################################################


#### Add the mergedColor info into Modules frame

Modules$mergedColours <- mergedColors

Interest_mods <- Modules[Modules$colnames.datExpr %in% lookfor,] 
LncMods <- Modules[Modules$colnames.datExpr %in% lookmore,]


Modules <- data.frame(colnames(datExpr), staticCut)
Modules$dynamicCut <- moduleColors
lookfor <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "EZH2", "DNMT3L", "MTHFR")

Interest_mods$Static <- Modules[Modules$colnames.datExpr %in% lookfor,] 


####Write out the interest modules file and the total module assignation file

write.table(Modules, "Module_assignation.txt", sep = "\t")
write.table(Interest_mods, "Interest_modules.txt", sep = "\t")

### They have not changed module - if they had, find new members of module via:

BlueMod <- subset(Modules,moduleColours == "blue")
BlackMod <-subset(Modules, moduleColours == "black")
#MageMod <- subset(Modules, moduleColours == "magenta")

#MidBlue <- subset(Modules, moduleColours == "midnightblue")