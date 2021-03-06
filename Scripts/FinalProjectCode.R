###########################################################################################
# P. Alexander Burnham
# April 23, 2017
# Ecological Genimics
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

############################################################################################
# create fake data to test code:
HealthStat <- rep(c("Sick","Sick","Healthy","Sick","Healthy","Healthy","Sick","Sick","Healthy","Sick","Healthy","Healthy","Sick","Sick","Healthy","Sick","Healthy","Healthy","Sick","Healthy"),2)

ID <- rep(1:20, 2)

Location <- rep(c("I", "O", "I", "O", "O", "I", "O", "I", "O", "I"), 4)

FoldChange <- rnorm(40, mean = 0, sd = 2)

FoldChangeBack <- rnorm(40, mean = 0, sd = 1)

Gene <- c(rep("Gene1", 20), rep("Gene2", 20))

DF <- data.frame(ID, Gene, HealthStat, Location, FoldChange, FoldChangeBack)
###########################################################################################

# using ddply to get summary of fold chane by by gene and status:
DF1 <- ddply(DF, c("HealthStat", "Gene"), summarise, 
                   n = length(FoldChange),
                   mean = mean(FoldChange, na.rm=TRUE),
                   sd = sd(FoldChange, na.rm=TRUE),
                   se = sd / sqrt(n))

###########################################################################################
# calculating confidence interval from null data:
meanCI <- mean(DF$FoldChangeBack)
n <- length(DF$FoldChangeBack)
s <- sd(DF$FoldChangeBack)
MOE <- qnorm(0.975)*s/sqrt(n)
upper <- meanCI + MOE
lower <- meanCI - MOE

###########################################################################################
# plotting using ggplot

#choosing color pallet
colors <- c("slategray3", "dodgerblue4")

#Create a bar graph for with CI and SE bars
plot1 <- ggplot(DF1, aes(x=Gene,
                         y=mean, 
                         fill=HealthStat)) + 
  geom_bar(stat="identity",
           color = "black",
           position=position_dodge()) + labs(x="Gene", y = "log(Fold Change)") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Health Status:") + theme(legend.position=c(.2, .85)) + coord_cartesian(ylim = c(-2, 2)) + annotate("segment", x = .5,xend = 2.5,y = meanCI, yend = meanCI, col = "red") + annotate("segment", x = .5, xend = 2.5,y = upper, yend = upper, col = "red", linetype = 2, size = 0.3) + annotate("segment", x = .5, xend = 2.5, y = lower, yend = lower, col = "red", linetype = 2, size = 0.3)








###########################################################################################
# Lauren's null model code:

### Ecological Genomics Final Project
### April 24, 2017   


###########################################################################################
# subsettting data:

library("DESeq2")
library("ggplot2")

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


### reading in list of immune related genes (269)
immuneNames <- read.delim("RawData/immune.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(immuneNames) <- immuneNames$trinID
immuneRes<-resFULL[which(row.names(resFULL) %in% row.names(immuneNames)),]
dim(immuneRes) # 254 matched -  less than 321 because of low read counts


### Reading in stress related genes (276)
stressNames <- read.delim("RawData/Pisaster_stress.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(stressNames) <- stressNames$trinID
stressRes<-resFULL[which(row.names(resFULL) %in% row.names(stressNames)),]
dim(stressRes) 


### Reading in virus related genes (325)
virusNames <- read.delim("RawData/virus.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(virusNames) <- virusNames$trinID
virusRes<-resFULL[which(row.names(resFULL) %in% row.names(virusNames)),]
dim(virusRes) 


### Reading in bacteria related genes (1320)
bacteriaNames <- read.delim("RawData/bacteria.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(bacteriaNames) <- bacteriaNames$trinID
bacteriaRes<-resFULL[which(row.names(resFULL) %in% row.names(bacteriaNames)),]
dim(bacteriaRes)


### Reading in disease related genes (29)
diseaseNames <- read.delim("RawData/disease.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(diseaseNames) <- diseaseNames$trinID
diseaseRes<-resFULL[which(row.names(resFULL) %in% row.names(diseaseNames)),]
dim(diseaseRes)


### Reading in pathogen related genes (44)
pathogenNames <- read.delim("RawData/pathogen.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(pathogenNames) <- pathogenNames$trinID
pathogenRes<-resFULL[which(row.names(resFULL) %in% row.names(pathogenNames)),]
dim(pathogenRes)



###########################################################################################
# Random distribution of full data set 


### IMMUNE

# Immune Random null distribution
logfoldrandoMean <- vector(mode="numeric")
for(i in 1:10000){
  logfoldrando<-vector(mode='numeric')
  logfoldrando <- sample(x=resFULL$log2FoldChange, size=dim(immuneRes)[1], replace=T)
  logfoldrandoMean[i]<-mean(logfoldrando)
}

## Histogram immune
hist(logfoldrandoMean, breaks=50, main="Log2Fold Change Randomization Test")
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
hist(logfoldrandoMean1, breaks=50, main="Log2Fold Change Randomization Test: Stress")
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
hist(logfoldrandoMean2, breaks=50, main="Log2Fold Change Randomization Test: Stress")
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
hist(logfoldrandoMean3, breaks=50, main="Log2Fold Change Randomization Test: Stress")
Interval <- quantile(x=logfoldrandoMean3, probs=c(0.025,0.975))
abline(v = Interval,col="black",lwd=2,lty="dotted")
abline(v=mean(bacteriaRes$log2FoldChange), col="red",lwd=2)

#--------------------------------------------------------------------------------------------


### DISEASE

# Stress Random null distribution
logfoldrandoMean4 <- vector(mode="numeric")
for(i in 1:10000){
  logfoldrando4<-vector(mode='numeric')
  logfoldrando4 <- sample(x=resFULL$log2FoldChange, size=dim(diseaseRes)[1], replace=T)
  logfoldrandoMean4[i]<-mean(logfoldrando4)
}

## Histogram - stress
hist(logfoldrandoMean4, breaks=50, main="Log2Fold Change Randomization Test: Stress")
Interval <- quantile(x=logfoldrandoMean4, probs=c(0.025,0.975))
abline(v = Interval,col="black",lwd=2,lty="dotted")
abline(v=mean(diseaseRes$log2FoldChange), col="red",lwd=2)

#--------------------------------------------------------------------------------------------


### PATHOGEN

# Stress Random null distribution
logfoldrandoMean5 <- vector(mode="numeric")
for(i in 1:10000){
  logfoldrando5<-vector(mode='numeric')
  logfoldrando5 <- sample(x=resFULL$log2FoldChange, size=dim(pathogenRes)[1], replace=T)
  logfoldrandoMean5[i]<-mean(logfoldrando5)
}

## Histogram - stress
hist(logfoldrandoMean5, breaks=50, main="Log2Fold Change Randomization Test: Stress")
Interval <- quantile(x=logfoldrandoMean5, probs=c(0.025,0.975))
abline(v = Interval,col="black",lwd=2,lty="dotted")
abline(v=mean(pathogenRes$log2FoldChange), col="red",lwd=2)

#--------------------------------------------------------------------------------------------


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
disease <- pvalue(nullDat=logfoldrandoMean4, obsDat=diseaseRes)
#-------------------------------------------------------------------------------------------
pathogen <- pvalue(nullDat=logfoldrandoMean5, obsDat=pathogenRes)
#-------------------------------------------------------------------------------------------


GeneGroup <- c("immune",
                 "stress",
                 "virus",
                 "bacteria",
                 "disease",
                 "pathogen")

xBar <- c(immune$xbar,
                stress$xbar,
                virus$xbar,
                bacteria$xbar,
                disease$xbar,
                pathogen$xbar)

SE <- c(immune$SE,
          stress$SE,
          virus$SE,
          bacteria$SE,
          disease$SE,
          pathogen$SE)

pVal <- c(immune$pval,
        stress$pval,
        virus$pval,
        bacteria$pval,
        disease$pval,
        pathogen$pval)

DF3 <- data.frame(GeneGroup, xBar, SE, pVal)

# removeing disease and pathogen
DF3 <- DF3[-c(5,6),]

###########################################################################################
# calculating confidence interval from null data:
meanCI <- mean(immune$mu,
                  stress$mu,
                  virus$mu,
                  bacteria$mu)

Interval1 <- quantile(x=logfoldrandoMean, probs=c(0.025,0.975))
Interval2 <- quantile(x=logfoldrandoMean1, probs=c(0.025,0.975))
Interval3 <- quantile(x=logfoldrandoMean2, probs=c(0.025,0.975))
Interval4 <- quantile(x=logfoldrandoMean3, probs=c(0.025,0.975))
#Interval5 <- quantile(x=logfoldrandoMean4, probs=c(0.025,0.975))
#Interval6 <- quantile(x=logfoldrandoMean5, probs=c(0.025,0.975))

lower <- mean(Interval1[1],Interval2[1],Interval3[1],Interval4[1])
upper <- mean(Interval1[2],Interval2[2],Interval3[2],Interval4[2])


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













































########## Testing stuff


pvalrep<-vector(mode='numeric')
for(i in 1:1000){
  logfoldrando <- sample(x=resFULL$log2FoldChange, size=dim(upreg)[1], replace=T)
  pvalrep[i]<-t.test(x=logfoldrando, y=upreg$log2FoldChange)$p.value
}

hist(upreg$log2FoldChange)
summary(pvalrep)
hist(pvalrep, breaks=50)

### a priori group them 

hist(logfoldrando, breaks=30) 
length(logfoldrando)
length(immuneRes$log2FoldChange)
hist(immuneRes$log2FoldChange, breaks=20, col=rgb(1,0,0,0.25), add=T)
length(subset(pvalrep, pvalrep<0.05))
shapiro.test(immuneRes$log2FoldChange)
immuneRes$log2FoldChange

Interval975 <- quantile(x=logfoldrandoMean, probs=c(0.025,0.975))
abline(v = Interval975,col="black",lwd=2,lty="dotted")
#abline(v=mean(immuneRes$log2FoldChange),col="indianred",lwd=2)
abline(v=mean(upreg$log2FoldChange), col="red",lwd=2)


### Grouping immune related into up and down regulated

upreg<-immuneRes[which(immuneRes$log2FoldChange > 0.04681608),]
downreg<-immuneRes[which(immuneRes$log2FoldChange < -0.16077978),]
#Interval95 <- quantile(x=logfoldrandoMean, probs=c(0.05,0.95))
#abline(v = Interval95,col="black",lwd=2,lty="dashed")
dim(immuneRes)
dim(upreg)
dim(downreg)



#### See if any immune related genes match with the significantly expressed genes

summary(resFULL)
sigdiffexpress<-subset(resFULL, resFULL$padj<0.1)
dim(sigdiffexpress) ## 39

immuneSigDiff<-resFULL[which(row.names(sigdiffexpress) %in% row.names(immuneNames)),]
dim(immuneSigDiff)

immuneSigDiff<-immuneNames[which(row.names(immuneNames) %in% row.names(sigdiffexpress)),]
dim(immuneSigDiff)


## Loading in full annotated table 
geneNames <- read.delim("poch_uniprot_GO_nr.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
str(geneNames)

sigdiffexpress<-subset(resFULL, resFULL$padj<0.1)

row.names(geneNames)<-geneNames$trinID

immuneSigDiff<-geneNames[which(row.names(geneNames) %in% row.names(sigdiffexpress)),]

sigDiff<-as.data.frame(sigdiffexpress, row.names=row.names(sigdiffexpress))

annotatedSigDiff<-merge(immuneSigDiff, sigDiff, by= "trinID")
names(sigDiff)
names(annotatedSigDiff)

write.csv(annotatedSigDiff, file="AnnotatedSigDEGenes_Day06.csv")
















