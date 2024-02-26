#/ Loading libraries

library(ggplot2)
library(ggfortify) 
library(limma)
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(EnhancedVolcano)

#/ Setting working directory
setwd("")

#/ Reading data
# Read Quantile normalized and Batch corrected train data

train1.bc=readRDS("train1.batch-corrected.rds")
labels.train1= train1.bc[c(6859:6861),]

train1.bc.mat= as.matrix(train1.bc[-c(6859:6861),])
class(train1.bc.mat)="numeric"

train1.bc.log2=log2(train1.bc.mat)# NAN produced
train1.bc.log2[is.nan(train1.bc.log2)]=0 # Replace NAN with 0

#/ Conversion of class labels to factors
labels.train1.t=as.data.frame(t(labels.train1))
groups=labels.train1.t$Source

f.source= factor(groups,levels = c("SCZ","CNT"))

# Create design matrix

design = model.matrix(~ 0 + f.source) 
colnames(design)=c("SCZ","CNT")

# Limma
train1.log2.t= as.data.frame(t(train1.bc.log2))
fit = lmFit(train1.bc.log2,design)

#/ Comparison between the groups 

contrast.matrix = makeContrasts(SCZ-CNT,levels=design)
fit.con = contrasts.fit(fit,contrast.matrix)
fit.eb = eBayes(fit.con)

# Top-table for all genes

train1.tab = topTable(fit.eb, number=Inf,adjust.method="BH", confint=0.95)

topgenes = train1.tab[train1.tab[, "adj.P.Val"] < 0.05, ]

# Volcano plot using Enchancedvolcano function with P value

EnhancedVolcano(train1.tab,
                lab = rownames(train1.tab),
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = "-Log10 adj.P.Value",
                pCutoff = 5e-2,
                FCcutoff = 1,
                labSize = 4,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

# Save DE file
saveRDS(train1.tab, "train1.tab.rds")

# Note: Repeat this analysis for all the iterations of training and datasets.
