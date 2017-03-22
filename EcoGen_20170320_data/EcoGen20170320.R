################################################################################################
# P. Alexander Burnham
# March 20, 2017
# Eco Gen
# working with frq and fst files:
################################################################################################
################################################################################################
################################################################################################

# set working directory
setwd("~/EcologicalGenomics/EcoGen_20170320_data")

# List the files in this directory -- you should see your results output from VCFTools if the download was successful
list.files()

# Let's do the allele freq comparisons first:
H_freq <- read.table("H_AlleleFreqs.frq", header=T)
S_freq <- read.table("S_AlleleFreqs.frq", header=T)

# Since these files have identical numbers of SNPs in the exact same order, we can concatenate them together into one large dataframe:
All_freq <- merge(H_freq, S_freq, by=c("CHROM", "POS"))

# Check the results of your merge to make sure things look OK
str(All_freq) # shows the structure of the data
head(All_freq)

# Looks good, now let's calculate the difference in minor allele frequency at each SNP and plot as a histogram
All_freq$diff <- (All_freq$H_ALT - All_freq$S_ALT)

hist(All_freq$diff, breaks=50, col="red", main="Allele frequency difference (H-S)")

# Looks like most loci show little difference (i.e., likely drift), but perhaps a few show very large differences between healthy and sick (drift or selection?)

# How do these highly divergent frequenices compare to Fst at the same SNPs?
fst <- read.table("HvS_Fst.weir.fst", header=T)

All_freq.fst <- merge(All_freq, fst, by=c("CHROM", "POS"))

# showing raw difference by variation in each of S and H pops. (as increase allele freq diff. in either direction you increase Fst)

# allele freq is higher in one that the other group (Fst=proportion of allele freq variation that is distributed between groups)
plot(All_freq.fst$diff, All_freq.fst$WEIR_AND_COCKERHAM_FST, xlab="Allele frequency difference (H-S)", ylab="Fst", main="Healthy vs. Sick SNP divergence")

# Which are the genes that are showing the highest divergence between Healthy and Sick?
All_freq.fst[which(All_freq.fst$WEIR_AND_COCKERHAM_FST>0.2),]
