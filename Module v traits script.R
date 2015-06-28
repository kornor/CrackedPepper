

#### Next step - look at expression of module vs trait



################################Define numbers of genes and samples to do traits comparison
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColours)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p");
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
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))
######################################

###### Set the trait of most interest in determining the modules

y = datTraits$MethScore

#### Calculate the modules' correlation to this trait

datME=moduleEigengenes(datExpr,moduleColours)$eigengenes
signif(cor(datME, use="p"), 2)

dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")


sizeGrWindow(8,9)
plotMEpairs(datME,y=y)

#######Heatmaps of module expression - two most interconnected modules

sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="blue";
plotMat(t(scale(datExpr[,moduleColours==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )

# for the other modules... we use
which.module="black";
plotMat(t(scale(datExpr[,moduleColours==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="turquoise";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )



#######
GS1=as.numeric(cor(y,datExpr, use="p"))

GeneSig <- data.frame(GS1) 
rownames(GeneSig)  <- colnames(datExpr)
list <- GeneSig[BlueMod$colnames.datExpr.,]
BlueSig <- as.data.frame(GeneSig[BlueMod$colnames.datExpr.,])
rownames(BlueSig) <- BlueMod$colnames.datExpr.
 

BlackSig <- as.data.frame(GeneSig[BlackMod$colnames.datExpr.,])
rownames(BlackSig) <- BlackMod$colnames.datExpr. 

BlackSig$DNMT1

  
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, moduleColours, mean, na.rm=T)

sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,moduleColours)

p.values = corPvalueStudent(cor(y,datME, use="p"), nSamples = length(y))
write.table(p.values, "Module correlation Trait.txt", sep ="\t")

ADJ1=abs(cor(datExpr,use="p"))^12
Alldegrees1=intramodularConnectivity(ADJ1, moduleColours)
head(Alldegrees1)


colorlevels=unique(moduleColours)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColours==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColours[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

######Just one module
whichmodule="blue";
restrict1 = (moduleColours==whichmodule);
verboseScatterplot(Alldegrees1$kWithin[restrict1],
                   GeneSignificance[restrict1], col=moduleColours[restrict1],
                   main=whichmodule,
                   xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)





##### Output the measures of gene signifigance

# Output of eigengene information:
datMEy = data.frame(y, datME)
eigengeneSignificance = cor(datMEy, y);
eigengeneSignificance[1,1] = (1+max(eigengeneSignificance[-1, 1]))/2
eigengeneSignificance.pvalue = corPvalueStudent(eigengeneSignificance, nSamples = length(y))
namesME=names(datMEy)
# Form a summary data frame
out1 = data.frame(t(data.frame(eigengeneSignificance,eigengeneSignificance.pvalue, namesME, t(datMEy))))
# Set appropriate row names
dimnames(out1)[[1]][1]="EigengeneSignificance"
dimnames(out1)[[1]][2]="EigengeneSignificancePvalue"
dimnames(out1)[[1]][3]="ModuleEigengeneName"
dimnames(out1)[[1]][-c(1:3)]=dimnames(datExpr)[[1]]
# Write the data frame into a file
write.table(out1, file="MEResultsNetworkScreening.csv", row.names=TRUE, col.names = TRUE, sep=",")
# Display the first few rows:
head(out1)

#####


ConnectBlue <- Alldegrees1[moduleColours == "blue",]
ConnectBrown <- Alldegrees1[moduleColors == "brown",]
ConnectTurq <- Alldegrees1[moduleColors == "turquoise",]


write.table(ConnectBlue, "Blue_intramodular_connect.txt", sep ="\t")
write.table(ConnectBrown, "Brown_intramodular_connect.txt", sep ="\t")
write.table(ConnectTurq, "Turquoise_intramodular_connect.txt", sep ="\t")


Her2 <- as.factor(clin$)
Group <- as.factor(clin$PAM)


aov <- aov(MEs$MEbrown ~ Her2)
summary(aov)
TukeyHSD(aov)
aov <- aov(MEs$MEturquoise ~ Her2)
summary(aov)
aov <- aov(MEs$MEred ~ ER)
summary(aov)

#############################################################


source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
biocLite("org.Hs.eg.db")
biocLite("AnnotationDbi")
library(Go.db)
library(AnnotationDbi)
library(org.Hs.eg.db)


source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
#load library
library(org.Hs.eg.db)

GOenr = GOenrichmentAnalysis(dynamicColors, allLLIDs, organism = "human", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R’s output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

###################### MDS scaling

cmd1=cmdscale(as.dist(dissTOM),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(dynamicColors), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")


________


#### Cytoscape

### First reduce the network to the top genes only
list <- row.names(Top30Black)
TOMblack <- datExpr[,list]

TOMBlue <- subset(datExpr, moduleColours == "blue")

# Recalculate topological overlap on only the new small expression set
TOM = TOMsimilarityFromExpr(TOMblue, power = 12);

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = names(TOMblack));





