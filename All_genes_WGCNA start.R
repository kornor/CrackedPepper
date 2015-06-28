setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_lnc")

### This analysis includes all the genes for mRNA analysis of samples, + lncs
### NOT trimmed to only the methylation related ones

### Load packages
library(WGCNA)
library(flashClust)
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
library(sparcl)


library(plyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)

transparentTheme(trans = 0.4)


## #########################################  This can be skipped - move down to load prepped file
### Load exp file & lnc expression file
total.exp <- read.table( "Complete_exp.txt", sep = "\t", header = TRUE, row.names = 1)
total.exp <- as.data.frame(t(total.exp))

topK.exp <- read.table("Top10K_exp.txt", sep = "\t", header = TRUE, row.names = 1)
topK.exp <- as.data.frame(t(topK.exp))


lnc <- read.table("lnc_for_merge.txt", sep = "\t", header = TRUE, row.names = 1)

merge.exp <- cbind(total.exp, lnc)
merge.exp <- cbind(topK.exp, lnc)

exp <- merge.exp

## look for genes with missing values
##Exclude genes with results for less than 400 samples; exclude samples with 
#results for less than 20 000 genes. 
gsg = goodSamplesGenes(exp, minNSamples = 413, minNGenes = 11221, verbose = 5);
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

exp <- as.data.frame(t(exp))
write.table(exp, "Top10K_exp_lnc.txt", sep = "\t")

exp <- read.table("Top10K_exp_lnc.txt", sep = "\t", header = TRUE, row.names = 1)
exp <- as.data.frame(t(exp))




############ Now need to match and trim the rows (samples) for all sets


list <- intersect(rownames(methTraits), rownames(clincom))

clincom <- clincom[list,]
exp <- exp[list,]
methTraits <- methTraits[list,]

exp <- as.data.frame(t(exp))

write.table(exp, "Final_top10K_exp_lnc.txt", sep = "\t")
write.table(clincom, "Final_clin_lnc.txt", sep = "\t")
write.table(methTraits, "Final_traits_lnc.txt", sep = "\t")

###########  **************   START HERE NOW   ****************   #############
########### Add in the exp & clinical files and the methylation traits of interest
exp <- read.table("Final_top10K_exp_lnc.txt", sep = "\t", header = TRUE, row.names = 1)
exp <- as.data.frame(t(exp))

clincom <- read.table( "Final_clin_lnc.txt", sep = "\t",header = TRUE, row.names = 1)

methTraits <- read.table("Final_traits_lnc.txt", sep = "\t",header = TRUE, row.names = 1)
methTraits1 <- read.table("methTraits_fang200i.txt", sep = "\t", header = TRUE, row.names = 1)
methTraits2 <- read.table("methTraits2_Stir200i.txt", sep = "\t", header = TRUE, row.names = 1)
allTraits <- read.table("Traits_ALL.txt", sep = "\t", header = TRUE, row.names = 1)


## look for genes with missing values
##Exclude genes with results for less than 400 samples; exclude samples with 
#results for less than 20 000 genes. 
gsg = goodSamplesGenes(exp, minNSamples = 408, minNGenes = 11221, verbose = 5);
gsg$allOK

## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(exp), method = "average");

#### Make factor labels for the subtypes
Pam50 <- as.factor(clincom$PAM50Call_RNAseq)
# Relabel blockwise modules


PamColours = labels2colors(Pam50)
count(PamColours)

# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)

ColorDendrogram(sampleTree, y = PamColours, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     labels = FALSE, branchlength = 35)

legend("topright",legend=levels(factor(Pam50)),
       fill = (c("turquoise","blue", "brown", "yellow", "green" )), cex = 1)
##########   Different way?

plotDendroAndColors(sampleTree, PamColours,
                    groupLabels = "Pam50 subtype",
                    main = "Sample dendrogram and subtype", dendroLabels = FALSE)
legend("topright",legend=levels(factor(Pam50)),
       fill = (c("turquoise","blue", "brown", "yellow", "green" )), cex = 1)


## Might keep the top basal set, just because
## Plot a line to show the cut

abline(h = 210, col = "red");

# Determine clusters under the line (use h from above)

#clust = cutreeStatic(sampleTree, cutHeight = 300, minSize = 10)
clust = cutreeStatic(sampleTree, cutHeight = 210, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.

keepSamples = (clust==1)
datExpr = exp[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


##Traits

Samples = rownames(datExpr);
traitRows = match(Samples, rownames(methTraits));
datTraits = methTraits[traitRows,];
datTraits1 = methTraits1[traitRows,];
datTraits2 = methTraits2[traitRows,];
datTraits3 = allTraits[traitRows,];

datClin = clincom[traitRows,]

collectGarbage();



###################

# Re-cluster samples to look for underlying patterns

sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, 
#grey means missing entry
# need to consider how to do this for discrete variables?

traitColors = numbers2colors(datTraits, 
                             signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap", dendroLabels = FALSE)


#### Note that the side panel - the "basal" group, is hypomethylated across all (comparatively)


### Save out the data and proceed to soft thresholding script ("SingleBlock_moduleCreation.R")

#save(datExpr, datTraits, datClin, file = "lncRNA_trimmed_input.RData")
#save(datExpr, datTraits, datClin, file = "lncRNA_total_input.RData")
save(datExpr, datTraits, datTraits1, datTraits2, datTraits3, datClin, file = "lncRNA_top10K_input.RData")




