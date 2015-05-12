### This is soft-thresholding and blockwise module creation 

#### Following on from file prep 

load(file="lncRNA_input.RData")


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


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting

moduleLabels = bwnet$colors
moduleColours = labels2colors(bwnet$colors)

## Have a look at a table of the modules
table(moduleColours)

##### Now we'll make a dataframe to store the moduleColours for genes, 
### We'll keep adding to this baby later
Modules <- data.frame(colnames(exp), moduleColours)
lookfor <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "EZH2", "DNMT3L")

Interest_mods <- Modules[Modules$colnames.exp %in% lookfor,] 


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

MEs = bwnet$MEs;
geneTree = bwnet$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "WGCNA_module_labels_allgenes.RData")
