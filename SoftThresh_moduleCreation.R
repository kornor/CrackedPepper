

### This is soft-thresholding and blockwise module creation 
setwd("~/Bioinformatics Work/TCGA initial project/WGCNA_lnc")

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

load(file="lncRNA_total_input.RData")


####



################
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to=18, by=2))
# Call the network topology analysis function

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, 
                        blockSize = 5000, networkType = "signed")


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
k=softConnectivity(datE=datExpr,power=12, type = "signed")

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")



bwnet = blockwiseModules(datExpr, power = 12, networkType = "signed",
                         TOMType = "signed", minModuleSize = 100,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "LncTOM_signed",deepsplit = 2,
                         verbose = 5)
################
## Create a frame of just datExp for genes in each block

block1 <- datExpr[bwnet$blockGenes[[1]]]
block2 <- datExpr[bwnet$blockGenes[[2]]]
block3 <- datExpr[bwnet$blockGenes[[3]]]
block4 <- datExpr[bwnet$blockGenes[[4]]]
block5 <- datExpr[bwnet$blockGenes[[5]]]


### See if you can find a gene in a block
match("FOXM1", names(block1))
match("FOXM1", names(block2))
match("FOXM1", names(block3))
match("FOXM1", names(block4))
match("FOXM1", names(block5))

match("DNMT1", names(block1))
match("DNMT1", names(block2))
match("DNMT1", names(block3))
match("DNMT1", names(block4))
match("DNMT1", names(block5))



# Convert labels to colors for plotting

# Relabel blockwise modules
moduleLabels = bwnet$colors
moduleColours = labels2colors(bwnet$colors)

## Have a look at a table of the modules
table(moduleColours)

##### Now we'll make a dataframe to store the moduleColours for genes, 
### We'll keep adding to this baby later
### Match the names in the blocks to the colours


Modules <- data.frame(colnames(datExpr), moduleColours)
lookfor <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "EZH2", "DNMT3L")

Interest_mods <- Modules[Modules$colnames.datExpr %in% lookfor,] 

# open a graphics window
sizeGrWindow(12, 9)
#### Manual Blockwise
##open a graphics window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], moduleColours[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colours in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], moduleColours[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colours in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 3
plotDendroAndColors(bwnet$dendrograms[[3]], moduleColours[bwnet$blockGenes[[3]]],
                    "Module colors", main = "Gene dendrogram and module colours in block 3",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 4
plotDendroAndColors(bwnet$dendrograms[[4]], moduleColours[bwnet$blockGenes[[4]]],
                    "Module colors", main = "Gene dendrogram and module colours in block 4",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 5
plotDendroAndColors(bwnet$dendrograms[[5]], moduleColours[bwnet$blockGenes[[5]]],
                    "Module colors", main = "Gene dendrogram and module colours in block 5",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


##############################################################
####  Now let's try the static cut

static1 <- as.character(cutreeStaticColor(bwnet$dendrograms[[1]], cutHeight=.98, minSize=100))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    colors = data.frame( moduleColours[bwnet$blockGenes[[1]]],static1),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")

static2 <- as.character(cutreeStaticColor(bwnet$dendrograms[[2]], cutHeight=.98, minSize=100))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(bwnet$dendrograms[[2]], 
                    colors = data.frame( moduleColours[bwnet$blockGenes[[2]]],static2),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")


static3 <- as.character(cutreeStaticColor(bwnet$dendrograms[[3]], cutHeight=.98, minSize=100))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(bwnet$dendrograms[[3]], 
                    colors = data.frame( moduleColours[bwnet$blockGenes[[3]]],static3),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")


static4 <- as.character(cutreeStaticColor(bwnet$dendrograms[[4]], cutHeight=.98, minSize=100))

# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(bwnet$dendrograms[[4]], 
                    colors = data.frame(moduleColours[bwnet$blockGenes[[4]]],static4),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")



static5 <- as.character(cutreeStaticColor(bwnet$dendrograms[[5]], cutHeight=.98, minSize=100))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(bwnet$dendrograms[[5]], 
                    colors = data.frame( moduleColours[bwnet$blockGenes[[5]]],static5),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")


#### Now concatonate all the static colours
staticColours <- c(static1, static2, static3, static4, static5)
table(staticColours)

### Add to the Modules file

Modules$staticColours <- staticColours
lookmore <- c("LINC00511", "LINC00665","LINC00680", "DANCR", "Z83851.1", "RP9P-002", "RP11.132A1.4",
              "AC073871.2", "MFI2.AS1")
Interest_mods <- Modules[Modules$colnames.datExpr %in% lookfor,] 
LncMods <- Modules[Modules$colnames.datExpr %in% lookmore,]
### you can see that the DNMTs pass out of the module in this case.


################ Make a frame of the module eigengenes for the modules in auto cut

MEList = moduleEigengenes(datExpr, colors = moduleColours)
MEs = MEList$eigengenes
AEs = MEList$averageExpr

##Bind MEs and AEs together

ModExp <- cbind(MEs, AEs)
## add rownames of samples
rownames(ModExp) <- rownames(datExpr)
### Write out file
write.table(ModExp, "ModuleInfo_dynamic_AutoCut.txt", sep ="\t")


############## Create a data frame of just the modules of interest, to see other members

BlueMod <- subset(Modules,moduleColours == "blue")

MageMod <- subset(Modules, moduleColours == "magenta")

MidBlue <- subset(Modules, moduleColours == "midnightblue")

##############  Save out the autos
geneTree = c(bwnet$dendrograms[[1]],bwnet$dendrograms[[2]],
             bwnet$dendrograms[[3]],bwnet$dendrograms[[4]],bwnet$dendrograms[[5]]);
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "WGCNA_module_labels_allgenes.RData")


######## Merging like modules
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, moduleColours, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

table(mergedColors)

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(moduleColours, mergedColors),
                    c("Auto", "Merged auto"),
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)

#dev.off()

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(mergedColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, mergedColors, geneTree, 
     file = "Merged_autoModules_allGenes.RData")
##################################################


#### Add the mergedColor info into Modules frame

Modules$mergedColours <- mergedColors

Interest_mods <- Modules[Modules$colnames.datExpr %in% lookfor,] 
LncMods <- Modules[Modules$colnames.datExpr %in% lookmore,]

### They have not changed module - if they had, find new members of module via:

#BlueMod <- subset(Modules,moduleColours == "blue")

#MageMod <- subset(Modules, moduleColours == "magenta")

#MidBlue <- subset(Modules, moduleColours == "midnightblue")