# set working directory: (files are in /RawData)
setwd("~/EcologicalGenomics")
library("DESeq2")

library("ggplot2")

countsTable <- read.delim("RawData/cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)

countData <- as.matrix(countsTable)
head(countData)
immunegenes <- read.csv("RawData/immunegenes.csv")
conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ health)
dim(dds)
dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)
colSums(counts(dds))
hist(colSums(counts(dds)), breaks = 80, xlim=c(0,max(colSums(counts(dds)))))
colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S"))

dds <- DESeq(dds, parallel = TRUE) 
norm.counts <- as.data.frame(counts(dds,normalized=TRUE))
dim(norm.counts)
head(row.names(norm.counts))

names <- rownames(norm.counts)
data <- cbind(names, norm.counts)
head(data)
dim(norm.counts)
dim(norm.counts1)
library(dplyr)
norm.counts1 <- add_rownames (norm.counts, "trinID")

newdata<- merge(immunegenes,norm.counts1, by="trinID")
dim(newdata)
dim(immunegenes)

newdata_counts<- as.matrix(newdata[,50:126])

heatmap(newdata_counts,hclustfun = hclust, scale = c("row"))
