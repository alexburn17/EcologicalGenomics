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

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Health Status:") + theme(legend.position=c(.2, .85)) + coord_cartesian(ylim = c(-2, 2)) + annotate("segment", 
           x = .5,
           xend = 2.5,
           y = meanCI, 
           yend = meanCI, 
           col = "red") + annotate("segment", 
                                   x = .5, xend = 2.5,
                                   y = upper, 
                                   yend = upper, 
                                   col = "red", 
                                   linetype = 2, 
                                   size = 0.3) + annotate("segment", 
                                                          x = .5, 
                                                          xend = 2.5, 
                                                          y = lower, 
                                                          yend = lower, 
                                                          col = "red", 
                                                          linetype = 2, 
                                                          size = 0.3)

###########################################################################################
# data analysis:

splitDF <- split(DF, DF$Gene)

model <- aov(data=DF, formula = FoldChange ~ HealthStat)
summary(model)

library(lme4)

nullMod <- lmer(data=splitDF$Gene1, formula = FoldChange ~ 1 + (1|ID) + (1|Location), REML = FALSE)
fullMod <- lmer(data=splitDF$Gene1, formula = FoldChange ~ HealthStat + (1|ID) + (1|Location), REML = FALSE)

anova(nullMod, fullMod)

###########################################################################################
# Lauren's null model code:

### Ecological Genomics Final Project
### April 24, 2017   

library("DESeq2")
library("ggplot2")


conds <- read.delim("RawData/cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
colData1 <- as.data.frame(conds)
colDataDay6<-subset(colData1, colData1$day=="day06")
colDataDay6
colData<-subset(colDataDay6, colDataDay6$indiv != "37")
colData # 7 sick # 11 healthy
# Day 15 H: 10, 24, 27, 31, 33, 34
# Day 15 S: 8, 9, 15, 19, 20, 23
## no 7, 22 (sick) in day 6
#sample log2foldchange of 279 genes with replacement - calculate average 1000 times to generate distribution; which proportion are sig on their own

colDataDay15<-subset(colData1, colData1$day=="day15")
colDataDay15

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

# TRINITY_DN45186_c2_g4_TRINITY_DN45186_c2_g4_i1_g.18786_m.18786 4.939978e-09; tumor suppressor


# reading in list of immune related genes (321)
immuneNames <- read.delim("Scripts/immune.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
row.names(immuneNames) <- immuneNames$trinID
immuneRes<-resFULL[which(row.names(resFULL) %in% row.names(immuneNames)),]
dim(immuneRes) # 254 -  less than 321 because low read counts
mean(immuneRes$log2FoldChange) # [1] 0.0199265
mean(immuneRes$padj[1:188]) # [1] 0.9695943


# Random null distribution
logfoldrandoMean <- vector(mode="numeric")
for(i in 1:10000){
  logfoldrando<-vector(mode='numeric')
  logfoldrando <- sample(x=resFULL$log2FoldChange, size=dim(immuneRes)[1], replace=T)
  logfoldrandoMean[i]<-mean(logfoldrando)
}

hist(logfoldrandoMean, breaks=50, main="Log2Fold Change Randomization Test")
hist()
Interval975 <- quantile(x=logfoldrandoMean, probs=c(0.025,0.975))
abline(v = Interval975,col="black",lwd=2,lty="dotted")
abline(v=mean(immuneRes$log2FoldChange),col="indianred",lwd=2)

#Interval95 <- quantile(x=logfoldrandoMean, probs=c(0.05,0.95))
#abline(v = Interval95,col="black",lwd=2,lty="dashed")

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

install.packages("topGO")
library(topGO)
library(ALL)
data(ALL)
data(geneList)


















