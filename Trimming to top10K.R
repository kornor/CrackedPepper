setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_lnc")

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


########### Load in the exp file

exp <- read.table( "Complete_exp.txt", sep = "\t", header = TRUE, row.names = 1)
exp <- as.data.frame(t(exp))

## look for genes with missing values
##Exclude genes with results for less than 400 samples; exclude samples with 
#results for less than 20 000 genes. 
gsg = goodSamplesGenes(exp, minNSamples = 413, minNGenes = 20000, verbose = 5);
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

datExpr <- exp 
#### Get the soft threshold

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



##To check the efficiency of the chosen soft power, put it in here and plot
### Might try 8 and 9

k=softConnectivity(datE=datExpr,power=9, type = "signed")

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

###To trim to top 10 000 genes (which will probably exclude most lncs?)
datExpr=datExpr[, rank(-k,ties.method="first" )<=10000]

match("DNMT1", colnames(datExpr))


datExpr <- t(datExpr)
write.table(datExpr, "Top10K_exp.txt", sep = "\t")
