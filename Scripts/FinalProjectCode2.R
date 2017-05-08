###########################################################################################
# P. Alexander Burnham and Lauren Ash
# April 23, 2017
# Ecological Genomics
# Team Sherlock Data Analysis
###########################################################################################

# Preliminaries:
# Clear memory of characters:
ls()
rm(list=ls())

# set working directory: (files are in /RawData)
setwd("~/EcologicalGenomics")

###########################################################################################
# load packages
library(plyr)
library(ggplot2)

###########################################################################################
# Lauren's null model code:
# subsettting data:

library(DESeq2)

### Subsetting data
conds <- read.delim("RawData/cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
colData1 <- as.data.frame(conds)
colDataDay6<-subset(colData1, colData1$day=="day06")
colDataDay6
colData<-subset(colDataDay6, colDataDay6$indiv != "37")
colData # 7 sick # 11 healthy
# Day 15 H: 10, 24, 27, 31, 33, 34
# Day 15 S: 8, 9, 15, 19, 20, 23
## no 7, 22 (sick) in day 6
#sample log2foldchange of 279 genes with replacement - calculate average 10000 times to generate distribution; which proportion are sig on their own


colDataDay15<-subset(colData1, colData1$day=="day15")
colDataDay15
str(countsTable)
countsTable <- read.delim('RawData/countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
countData<-countData[, which(colnames(countData) %in% row.names(colData))]

## DESeq Model: Day 6, no MMs, 18 individuals (7 sick, 11 healthy); model with location and health 
library(DESeq2)
ddsFULL <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)

ddsFULL <- ddsFULL[ rowSums(counts(ddsFULL)) > 100, ]

colData(ddsFULL)$health <- factor(colData(ddsFULL)$health, levels=c("H","S"))

ddsFULL <- DESeq(ddsFULL)
resFULL <- results(ddsFULL)
resFULL <- resFULL[order(resFULL$padj),] #sorts according to pvalue
head(resFULL)
summary(resFULL)
names(resFULL)

###########################################################################################
# Reading in Data for each group of genes:


### reading in list of immune related genes (259)
immuneNames <- read.delim("RawData/immune.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(immuneNames) <- immuneNames$trinID
immuneRes<-resFULL[which(row.names(resFULL) %in% row.names(immuneNames)),]
dim(immuneRes) 

#--------------------------------------------------------------------------------------------

### Reading in stress related genes (276)
stressNames <- read.delim("RawData/Pisaster_stress.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(stressNames) <- stressNames$trinID
stressRes<-resFULL[which(row.names(resFULL) %in% row.names(stressNames)),]
dim(stressRes) 

#--------------------------------------------------------------------------------------------

### Reading in virus related genes (325)
virusNames <- read.delim("RawData/virus.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(virusNames) <- virusNames$trinID
virusRes<-resFULL[which(row.names(resFULL) %in% row.names(virusNames)),]
dim(virusRes) 

#--------------------------------------------------------------------------------------------

### Reading in bacteria related genes (1320)
bacteriaNames <- read.delim("RawData/bacteria.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(bacteriaNames) <- bacteriaNames$trinID
bacteriaRes<-resFULL[which(row.names(resFULL) %in% row.names(bacteriaNames)),]
dim(bacteriaRes)


###########################################################################################
# Random distribution of full data set 
# set parameters to original
oPar <- par(no.readonly=TRUE)

par(mfrow=c(2,2))
par(mar=c(4.5,4.5,2,1))

### IMMUNE

# Immune Random null distribution
logfoldrandoMean <- vector(mode="numeric")
for(i in 1:10000){
  logfoldrando<-vector(mode='numeric')
  logfoldrando <- sample(x=resFULL$log2FoldChange, size=dim(immuneRes)[1], replace=T)
  logfoldrandoMean[i]<-mean(logfoldrando)
}

## Histogram immune
hist(logfoldrandoMean, 
     breaks=50, 
     main="IMMUNE GENES",
     xlab=NULL,
     ylim=c(0,800),
     xlim=c(-.3, .15),
     cex.axis=1,
     cex.lab=1.2,
     font.lab=2 
     )

Interval <- quantile(x=logfoldrandoMean, probs=c(0.025,0.975))
abline(v = Interval,col="black",lwd=2,lty="dotted")
abline(v=mean(immuneRes$log2FoldChange), col="red",lwd=2)

#--------------------------------------------------------------------------------------------

### STRESS

# Stress Random null distribution
logfoldrandoMean1 <- vector(mode="numeric")
for(i in 1:10000){
  logfoldrando1<-vector(mode='numeric')
  logfoldrando1 <- sample(x=resFULL$log2FoldChange, size=dim(stressRes)[1], replace=T)
  logfoldrandoMean1[i]<-mean(logfoldrando1)
}

## Histogram - stress
hist(logfoldrandoMean1, 
     breaks=50, 
     main="STRESS GENES",
     xlab=NULL,
     ylim=c(0,800),
     xlim=c(-.3, .15),
     cex.axis=1,
     cex.lab=1.2,
     font.lab=2 
)

Interval <- quantile(x=logfoldrandoMean1, probs=c(0.025,0.975))
abline(v = Interval,col="black",lwd=2,lty="dotted")
abline(v=mean(stressRes$log2FoldChange), col="red",lwd=2)

#--------------------------------------------------------------------------------------------

### VIRUS

# Stress Random null distribution
logfoldrandoMean2 <- vector(mode="numeric")
for(i in 1:10000){
  logfoldrando2<-vector(mode='numeric')
  logfoldrando2 <- sample(x=resFULL$log2FoldChange, size=dim(virusRes)[1], replace=T)
  logfoldrandoMean2[i]<-mean(logfoldrando2)
}

## Histogram - stress
hist(logfoldrandoMean2,
     breaks=50,
     main="VIRUS GENES",
     xlab="log(fold change)",
     ylim=c(0,800),
     xlim=c(-.3, .15),
     cex.axis=1,
     cex.lab=1.2,
     font.lab=2 
)

Interval <- quantile(x=logfoldrandoMean2, probs=c(0.025,0.975))
abline(v = Interval,col="black",lwd=2,lty="dotted")
abline(v=mean(virusRes$log2FoldChange), col="red",lwd=2)

#--------------------------------------------------------------------------------------------


### BACTERIA

# Stress Random null distribution
logfoldrandoMean3 <- vector(mode="numeric")
for(i in 1:10000){
  logfoldrando3<-vector(mode='numeric')
  logfoldrando3 <- sample(x=resFULL$log2FoldChange, size=dim(bacteriaRes)[1], replace=T)
  logfoldrandoMean3[i]<-mean(logfoldrando3)
}

## Histogram - stress
hist(logfoldrandoMean3, 
     breaks=50,
     main="BACTERIA GENES",
     xlab="log(fold change)",
     ylim=c(0,800),
     xlim=c(-.3, .15),
     cex.axis=1,
     cex.lab=1.2,
     font.lab=2 
)

Interval <- quantile(x=logfoldrandoMean3, probs=c(0.025,0.975))
abline(v = Interval,col="black",lwd=2,lty="dotted")
abline(v=mean(bacteriaRes$log2FoldChange), col="red",lwd=2)

#--------------------------------------------------------------------------------------------
# set original par for multipanel plotting downstream
par(oPar)



##############################################################################################
### SET UP FUNCTION FOR DATA ANALYSIS:

##############################################################################################
# Function= pvalue -> takes two DFs and calculates summary stats and pval (z test)
# INPUT: nullDat=nullDF, obsDat = ResDF
# OUTPUT: xbar, mu, n, SD ,SE , z, pval
##############################################################################################
pvalue <- function(nullDat=logfoldrandoMean1, obsDat=stressRes){
  
  obsDat <- obsDat[,2]
  
  xbar <- mean(obsDat)
  mu <- mean(nullDat)
  n <- length(obsDat)
  SD <- sd(obsDat)
  SE <- SD/sqrt(n)
  
  z <- (xbar - mu)/(SE)
  
  pval <- 2*pnorm(-abs(z), lower.tail = TRUE)
  
  out <- list(xbar=xbar, mu=mu, n=n, SD=SD, SE=SE, z=z, pval=pval)
  
  return(out)
}

##############################################################################################
# END FUNCTION

### running test (pvalue) for each of 6 groups

#-------------------------------------------------------------------------------------------
immune <- pvalue(nullDat=logfoldrandoMean, obsDat=immuneRes)
#-------------------------------------------------------------------------------------------
stress <- pvalue(nullDat=logfoldrandoMean1, obsDat=stressRes)
#-------------------------------------------------------------------------------------------
virus <- pvalue(nullDat=logfoldrandoMean3, obsDat=virusRes)
#-------------------------------------------------------------------------------------------
bacteria <- pvalue(nullDat=logfoldrandoMean3, obsDat=bacteriaRes)
#-------------------------------------------------------------------------------------------

#############################################################################################
# creating data frame for figure: (genegroup, mean, SE etc.)
GeneGroup <- c("immune",
               "stress",
               "virus",
               "bacteria")


xBar <- c(immune$xbar,
          stress$xbar,
          virus$xbar,
          bacteria$xbar)


SE <- c(immune$SE,
        stress$SE,
        virus$SE,
        bacteria$SE)

pVal <- c(immune$pval,
          stress$pval,
          virus$pval,
          bacteria$pval)


DF3 <- data.frame(GeneGroup, xBar, SE, pVal)


###########################################################################################
# calculating confidence interval from null data:
meanCI <- mean(immune$mu,
               stress$mu,
               virus$mu,
               bacteria$mu)

# caclulating quantiles foe each group
IntervalIM <- quantile(x=logfoldrandoMean, probs=c(0.025,0.975))
IntervalST<- quantile(x=logfoldrandoMean1, probs=c(0.025,0.975))
IntervalVI <- quantile(x=logfoldrandoMean2, probs=c(0.025,0.975))
IntervalBA <- quantile(x=logfoldrandoMean3, probs=c(0.025,0.975))

lower <- mean(IntervalIM[1],IntervalST[1],IntervalVI[1],IntervalBA[1])
upper <- mean(IntervalIM[2],IntervalST[2],IntervalVI[2],IntervalBA[2])


###########################################################################################
# plotting using ggplot

#choosing color pallet
colors <- c(rep("dodgerblue4",4))

#Create a bar graph for with CI and SE bars
plot1 <- ggplot(DF3, aes(x=GeneGroup,
                         y=xBar,
                         fill=GeneGroup)) + 
  geom_bar(stat="identity",
           color = "black",
           position=position_dodge()) + labs(x="Gene Group", y = "log(Fold Change)") + geom_errorbar(aes(ymin=xBar-SE, ymax=xBar+SE), width=.2, position=position_dodge(.9)) 

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors) + coord_cartesian(ylim = c(-.3, .3)) + annotate("segment", x = .5,xend = 4.5 ,y = meanCI, yend = meanCI, col = "red") + annotate("segment", x = .5, xend = 4.5, y = upper, yend = upper, col = "red", linetype = 2, size = 0.3) + annotate("segment", x = .5, xend = 4.5, y = lower, yend = lower, col = "red", linetype = 2, size = 0.3) + annotate("segment", x = .5, xend = 4.5, y = 0, yend = 0, col = "black", size = 1) + annotate(geom = "text", x = 3, y = .2, label = "*",cex = 8) + theme(legend.position=c(3, 3))

#############################################################################################
# END SCRIPT
