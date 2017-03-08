##################################################################################

# P. Alexander Burnham
# Eco Genomics 
# Homework 1
# 6, March 2017

##################################################################################
# Preliminaries:
# Clear memory of characters:
ls()
rm(list=ls())

# set working directory: (files are in /RawData)
setwd("~/EcologicalGenomics")

# plot my quality plots (qual, qualINT and qualSUB together in one pannel)
par(mfrow=c(1,3))

##################################################################################
# Loading packages, data and prepping the data

library("DESeq2")
library("ggplot2")

countsTable <- read.delim("RawData/countsdata_trim2.txt",  
                          stringsAsFactors=TRUE,
                          header = TRUE,
                          row.names=1)

countData <- as.matrix(countsTable)

head(countData)

conds <- read.delim("RawData/cols_data_trim.txt",
                    header=TRUE, 
                    stringsAsFactors=TRUE, 
                    row.names=1)

colData <- as.data.frame(conds)
head(colData)

# dimensions of the two data structures:
dim(countData)
dim(colData)

######################################################################################
# Subsetting colData by health status (H/S)

# colData
x<- split(colData, colData$location)
colDataInt <- x$int
colDataSub <- x$sub

# countData
# get names of rows for INT and SUB
i <- rownames(colDataInt)
s <- rownames(colDataSub)

# subsetting by INT and SUB
countDataInt <- countData[,i]
countDataSub <- countData[,s]
dim(countDataInt)
dim(countDataSub)

################################################################################
################# MODEL NUMBER 1: TEST EFFECT OF HEALTH CONTROLLING FOR LOCATION
################################################################################

# model looks for a main effect of health while controlling for the variance ascociated with location
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)

# check dimensions of data
dim(dds)
dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)

# For practice, let's work with fewer genes so that the models can run faster...
#dds <- dds[sample(nrow(dds), 1200), ]
#dim(dds)

#sets that "healthy is the reference
colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) 

# Run the model on the data
dds <- DESeq(dds) 

# summarise data
res <- results(dds)
res <- res[order(res$padj),]
head(res)
summary(res)



#--------------------------------------------------------------------------------
#--------------------------------------------DATA VISUALIZATION (FIGURES) MODEL 1
#--------------------------------------------------------------------------------

# Data quality assessment by sample clustering and visualization 
qual <- plotMA(res, main="DESeq2 Model 1 (controlling for location)", ylim=c(-2,2))

# PCA plot for health
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
plot <- plotPCA(vsd, intgroup=c("health"))
plot

# out of 12946 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 206, 1.6% 
# LFC < 0 (down)   : 58, 0.45% 
# outliers [1]     : 396, 3.1% 
# low counts [2]   : 7433, 57% 
# (mean count < 22)

################################################################################
######################## MODEL NUMBER 2: TEST EFFECT OF HEALTH WITHIN INTERTIDAL
################################################################################

# model looks for a main effect of health within intertidal population:
ddsINT <- DESeqDataSetFromMatrix(countData = countDataInt, colData = colDataInt, design = ~ health)

# check dimensions of data
dim(ddsINT)
ddsINT <- ddsINT[ rowSums(counts(ddsINT)) > 100, ]
dim(ddsINT)

# For practice, let's work with fewer genes so that the models can run faster...
#ddsINT <- ddsINT[sample(nrow(ddsINT), 1200), ]
#dim(ddsINT)

#sets that "healthy is the reference
colData(ddsINT)$health <- factor(colData(ddsINT)$health, levels=c("H","S")) 

# Run the model on the data
ddsINT <- DESeq(ddsINT) 

# summarise data
resINT <- results(ddsINT)
resINT <- resINT[order(resINT$padj),]
head(resINT)
summary(resINT)


#--------------------------------------------------------------------------------
#--------------------------------------------DATA VISUALIZATION (FIGURES) MODEL 2
#--------------------------------------------------------------------------------

# Data quality assessment by sample clustering and visualization 
qualINT <- plotMA(resINT, main="DESeq2 Model 2 (Intertidal Only)", ylim=c(-2,2))

# PCA plot for health
vsdINT <- varianceStabilizingTransformation(ddsINT, blind=FALSE)
plotINT <- plotPCA(vsdINT, intgroup=c("health"))
plotINT 

# out of 12396 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 204, 1.6% 
# LFC < 0 (down)   : 33, 0.27% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 9139, 74% 
# (mean count < 39)

################################################################################
########################## MODEL NUMBER 3: TEST EFFECT OF HEALTH WITHIN SUBTIDAL
################################################################################

# model looks for a main effect of health within subtidal population:
ddsSUB <- DESeqDataSetFromMatrix(countData = countDataSub, colData = colDataSub, design = ~ health)

# check dimensions of data
dim(ddsSUB)
ddsSUB <- ddsSUB[ rowSums(counts(ddsSUB)) > 100, ]
dim(ddsSUB)

# For practice, let's work with fewer genes so that the models can run faster...
#ddsSUB <- ddsSUB[sample(nrow(ddsSUB), 1200), ]
#dim(ddsSUB)

#sets that "healthy is the reference
colData(ddsSUB)$health <- factor(colData(ddsSUB)$health, levels=c("H","S")) 

# Run the model on the data
ddsSUB <- DESeq(ddsSUB) 

# summarise data
resSUB <- results(ddsSUB)
resSUB <- resSUB[order(resSUB$padj),]
head(resSUB)
summary(resSUB)


#--------------------------------------------------------------------------------
#--------------------------------------------DATA VISUALIZATION (FIGURES) MODEL 3
#--------------------------------------------------------------------------------

# Data quality assessment by sample clustering and visualization 
qualSUB <- plotMA(resSUB, main="DESeq2 Model 3 (Subtidal Only)", ylim=c(-2,2))

# PCA plot for health
vsdSUB <- varianceStabilizingTransformation(ddsSUB, blind=FALSE)
plotSUB <- plotPCA(vsdSUB, intgroup=c("health"))
plotSUB

# out of 12392 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 20, 0.16% 
# LFC < 0 (down)   : 113, 0.91% 
# outliers [1]     : 647, 5.2% 
# low counts [2]   : 4289, 35% 
# (mean count < 13)


################################################################################
########################################## FUNCTION TO PLOT GGPLOT FIGS TOGETHER
################################################################################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Plot figures together:
multiplot(plot, plotINT, plotSUB, cols=3)


