setwd("~/Bioinformatics Work/Meth & RNA/Methylation files_modified")
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
library(sm)

library(plyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)

######  This is code to check the status of the distribution of methylation in the meth sets
#### Clinical file and make factor
clincom <- read.table( "Final_clin_lnc.txt", sep = "\t",header = TRUE, row.names = 1)
Pam50 <- as.factor(clincom$PAM50Call_RNAseq)


### Load in the fang 200i file first

fang200i <- read.table("Fang_TSS200_island_final.txt", sep = "\t", header = TRUE, row.names = 1)

hist(fang200i$ARHGEF7)
hist(fang200i$RASGRF2)
hist(fang200i$SFRP2)
hist(fang200i$AOX1, main = "Histogram of Beta values \n for AOX1 (TSS200, islands)", 
     xlab = "Methylation Beta Value",
     col = "cornflowerblue")

hist(fang200i$RASGRF2, main = "Histogram of Beta values \n for RASGRF2 (TSS200, islands)", 
     xlab = "Methylation Beta Value",
     col = "cornflowerblue")

hist(fang200i$SFRP2, main = "Histogram of Beta values \n for SFRP2 (TSS200, islands)", 
     xlab = "Methylation Beta Value",
     col = "tomato1")


fang200a <- read.table("Fang_TSS200_all_final.txt", sep = "\t", header = TRUE, row.names = 1)

hist(fang200a$AOX1)
hist(fang200a$ACTA1)
hist(fang200a$RASGRF2)
hist(fang200a$ARHGEF7)

fang1500a <- read.table("Fang_TSS1500_all_final.txt", sep = "\t", header = TRUE, row.names = 1)

hist(fang1500a$AOX1)
hist(fang1500a$ACTA1)
hist(fang1500a$RASGRF2)
hist(fang1500a$ARHGEF7)


stir200a <- read.table("Stir_TSS200_all_final.txt", sep = "\t", header = TRUE, row.names = 1)

hist(stir200a$SLC6A3)



############ Looking at dis by subtype?

attach(fang200i)

# create value labels 
cyl.f <- factor(cyl, levels= c(4,6,8),
                labels = c("4 cylinder", "6 cylinder", "8 cylinder")) 

# plot densities 
sm.density.compare(RASGRF2, Pam50, xlab="Methylation Beta value")
title(main="Methylation beta value \n by subtype")

# add legend via mouse click
colfill<-c(2:(2+length(levels(Pam50)))) 
legend(locator(1), levels(Pam50), fill=colfill)


# Kernel density plots for mpg
# grouped by number of gears (indicated by color)
qplot(RASGRF2, data=fang200i, geom="density", fill=Pam50, alpha=I(.5), 
      main="Distribution of Methylation Beta Value \n for RASGRF2", xlab="Beta Value", 
      ylab="Density")


qplot(RASGRF2, data=fang200i, geom="density", colour=Pam50, alpha=I(.5), 
      main="Distribution of Methylation Beta Value \n for RASGRF2 (TSS200, islands)", xlab="Beta Value", 
      ylab="Density")


qplot(AOX1, data=fang200i, geom="density", colour=Pam50, alpha=I(.5), 
      main="Distribution of Methylation Beta Value \n for AOX1 (TSS200, islands)", xlab="Beta Value", 
      ylab="Density")


qplot(SFRP2, data=fang200i, geom="density", colour=Pam50, alpha=I(.6), 
      main="Distribution of Methylation Beta Value \n for SFRP2 (TSS200, islands)", xlab="Beta Value", 
      ylab="Density")

fang200i$Ave <- rowMeans(fang200i[,1:21])
qplot(Ave, data=fang200i, geom="density", colour=Pam50, alpha=I(.6), 
      main="Distribution of 'MethScore' \n (using Fang, TSS200, islands)", xlab="MethScore", 
      ylab="Density")


fang <- read.table("Fang_TSS200_island_exp_set.txt", sep = "\t", header = TRUE, row.names = 1)

qplot(MethScore, data=fang, geom="density", colour=Pam50, alpha=I(.6), 
      main="Distribution of 'MethScore' \n (using Fang, TSS200, islands)", xlab="MethScore", 
      ylab="Density")
