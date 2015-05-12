setwd("~/Bioinformatics Work/TCGA initial project/WGCNA_lnc")

### This analysis includes all the genes for mRNA analysis of samples, + lncs
### NOT trimmed to only the methylation related ones

total.exp <- read.table( "Complete_exp.txt", sep = "\t", header = TRUE, row.names = 1)
total.exp <- as.data.frame(total.exp)
lnc <- read.table("lnc_for_merge.txt", sep = "\t", header = TRUE, row.names = 1)
lnc <- as.data.frame(t(lnc))
merge.exp <- rbind(total.exp, lnc)

exp <- merge.exp
## look for genes with missing values
gsg = goodSamplesGenes(merge.exp, verbose = 3);
gsg$allOK


## if return is true, all good to go
### Otherwise

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(exp)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(exp)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exp = exp[gsg$goodSamples, gsg$goodGenes]
}

########### Trying again


exp <- read.table("Total_com_merged.txt", sep = "\t",header = TRUE, row.names = 1)

clincom <- read.table( "Clincom_merged.txt", sep = "\t",header = TRUE, row.names = 1)

methTraits.com <- read.table( "MethTraits_com_merged.txt", sep = "\t",header = TRUE, row.names = 1)




## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(exp), method = "average");

# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2)

## if an outlier there - cut the tree to remove it
## Plot a line to show the cut

abline(h = 3e+05, col = "red");

# Determine clusters under the line (use h from above)

clust = cutreeStatic(sampleTree, cutHeight = 3e+05, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.

keepSamples = (clust==1)
datExpr = exp#[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)



clincom <- read.table( "Clincom_merged.txt", sep = "\t",header = TRUE, row.names = 1)

methTraits.com <- read.table( "MethTraits_com_merged.txt", sep = "\t",header = TRUE, row.names = 1)


##Traits

Samples = rownames(datExpr);
traitRows = match(Samples, rownames(methTraits.com));

datTraits = methTraits.com[traitRows,];
#rownames(datTraits) = methTraits[traitRows, 1]


#rownames(datTraits) = clin[traitRows,];
collectGarbage();

save(datExpr, datTraits, file = "lncRNA_input.RData")

load(file="lncRNA_input.RData")
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

##Now trim to top genes and search for DNMT1
k=softConnectivity(datE=datExpr,power=9, type = "signed")

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

datExpr=datExpr[, rank(-k,ties.method="first" )<=10000]

match("DNMT1", colnames(datExpr))

softPower = 9;
adjacency = adjacency(datExpr, power = softPower, type = "signed");
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM


match("DNMT1", colnames(datExpr))

######## Make the tree
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

###### Now we're going to go ahead with module creation
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
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

MEDissThres = 0.25
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
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "Dynamic_merged_lncRNA_sub.RData")

staticCut <- as.character(cutreeStaticColor(geneTree, cutHeight=.99, minSize=30))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(geneTree, 
                    colors = data.frame( moduleColors,staticCut),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")
table(staticCut)

##################
MEList = moduleEigengenes(datExpr, colors = staticCut)
MEs = MEList$eigengenes
AEs = MEList$averageExpr

##Bind MEs and AEs together

ModExp <- cbind(MEs, AEs)
## add rownames of samples
rownames(ModExp) <- rownames(datExpr)
### Write out file
write.table(ModExp, "ModuleInfo_static_sub.txt", sep ="\t")


save(MEs, staticCut, geneTree, file = "Static_sub_lncRNA.RData")


########################

################################Define numbers of genes and samples to do traits comparison
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, Ave_meth, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

##########Will only display correlations
textMatrix = paste(signif(moduleTraitCor, 2),sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(Ave_meth),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))
##############

Modules <- data.frame(colnames(datExpr), moduleColors)
lookfor <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1.1", "TET1")

DNMT_mods <- Modules[Modules$colnames.datExpr %in% lookfor,] 

