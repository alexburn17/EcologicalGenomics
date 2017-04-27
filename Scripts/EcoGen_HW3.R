################################################################################################
# P. Alexander Burnham
#April 4, 2017
# Eco Gen
# Homework 3:
################################################################################################

################################################################################################
# Prelims and data loading and prepping:


# set working directory
setwd("~/EcologicalGenomics/RawData")

# List the files in this directory -- you should see your results output from VCFTools if the download was successful
list.files()

# ...and load the libraries
library(vcfR)
library(adegenet)

#Read the vcf SNP data into R
minMaxDat <- read.vcfR("SSW_by24inds.minMaxAllele.recode.vcf")
MAFdat <- read.vcfR("SSW_by24inds.MAF.recode.vcf")


# create genlight object for minMaxDat
gl1 <- vcfR2genlight(minMaxDat)
print(gl1) # Looks good! Right # of SNPs and individuals!


# create genlight object for lessMiss
gl2 <- vcfR2genlight(MAFdat)
print(gl2) # Looks good! Right # of SNPs and individuals!


# For info, try:
gl1$ind.names
gl1$loc.names[1:10]
gl1$chromosome[1:3]

# Notice there's nothing in the field that says "pop"? Let's fix that...
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort by Individual ID, just like the VCF file

# assigning location and meta data to both lessMiss and mafDat

#mafDat
gl1$pop <- ssw_meta$Location # assign locality info
gl1$other <- as.list(ssw_meta$Trajectory) # assign disease status

# lessMiss 
gl2$pop <- ssw_meta$Location # assign locality info
gl2$other <- as.list(ssw_meta$Trajectory) # assign disease status

###############################################################################################
# END OF DATA LOADING ETC.



################################################################################################
# DATA ANALYSIS

# WE can explore the structure of our SNP data using the glPlot function, which gives us a sample x SNP view of the VCF file
glPlot(gl1, posi="bottomleft")
glPlot(gl2, posi="bottomleft")


# row, let's compute the PCA on the SNP genotypes and plot it:
pca1 <- glPca(gl1, nf=4) # nf = number of PC axes to retain (here, 4)
pca1 # prints summary


# row, let's compute the PCA on the SNP genotypes and plot it:
pca2 <- glPca(gl2, nf=4) # nf = number of PC axes to retain (here, 4)
pca2 # prints summary



################################################################################################
# PLOTTING



#PLOTTING FOR min and max allele freq of (2)
# Plot the individuals in SNP-PCA space, with locality labels:
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=gl1$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2",
     main="PCA on MAF w/no missing data (3375 SNPs)")


legend("topleft", 
       legend=unique(gl1$pop), 
       pch=20, 
       col=c("black", "red"))






# perhaps we want to show disease status instead of locality:
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=as.factor(unlist(gl1$other)), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2",
     main="PCA on MAF w/no missing data (3375 SNPs)"
)

legend("topleft", 
       legend=unique(as.factor(unlist(gl1$other))), 
       pch=20, 
       col=as.factor(unique(unlist(gl1$other))))



#--------------------------------------------------------------------------------------------



#PLOTTING FOR NO MISSING DATA
# Plot the individuals in SNP-PCA space, with locality labels:
plot(pca2$scores[,1], pca2$scores[,2], 
     cex=2, pch=20, col=gl2$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2",
     main ="PCA on Min & Max alleles = 2 (275 SNPs)")

legend("topleft", 
       legend=unique(gl2$pop), 
       pch=20, 
       col=c("black", "red")
       )






# perhaps we want to show disease status instead of locality:
plot(pca2$scores[,1], pca2$scores[,2], 
     cex=2, pch=20, col=as.factor(unlist(gl2$other)), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2",
     main="PCA on Min & Max alleles = 2 (275 SNPs)")

legend("topleft", 
       legend=unique(as.factor(unlist(gl2$other))), 
       pch=20, 
       col=as.factor(unique(unlist(gl2$other))))





################################################################################################
# EXPLORING RESULTS

# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(pca1$loadings[,1]),
            threshold=quantile(abs(pca1$loadings), 0.999))

# Get their locus names
threshold <- quantile(abs(pca1$loadings), 0.999)
gl1$loc.names[which(abs(pca1$loadings)>threshold)]



