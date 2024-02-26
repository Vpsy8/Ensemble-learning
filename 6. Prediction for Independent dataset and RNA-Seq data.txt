#/ Loading libraries
library(plyr)
library(dplyr)
library(tibble)
library(tidyverse)
library(e1071)
library(janitor)
library(caret)
library(pamr)

# Best performing models, SVM-radial and PAM with 400 DEGs, were used for prediction analysis.

#/ Setting working directory
setwd("")

#/ Loading data
# 1. DEG Data
train1.DEG= readRDS("train1.tab.rds")

# 2. Quantile + Batch corrected data
train1.bc= as.data.frame(t(readRDS("train1.batch-corrected.rds")))
train1.bc = train1.bc %>% relocate(GSE_ID, Source, Batch)
labels.train1= as.factor(train1.bc[,c(2)])

testI1.bc= as.data.frame(t(readRDS("I1.GSE27383.rds")))
testI1.ph= read.csv("Pheno_27383.csv")
testI1.ph= testI1.ph[-c(1)]
colnames(testI1.ph)[1] = "GSE_ID"

testI1.bc= cbind(testI1.bc, testI1.ph)

testI1.bc = testI1.bc %>% relocate(GSE_ID, Source, Batch)
labels.testI1= as.factor(testI1.bc[,c(2)])

# Select the DEGs of interest w.r.t to the adj.P.value
abc= train1.DEG[order(train1.DEG$adj.P.Val, decreasing = F),]
top400DEG= abc[1:400,]

# Subset training and testing data for DEGs of interest

train1.top400= train1.bc[,rownames((top400DEG))]
train1.top400.mat= as.matrix(train1.top400)
class(train1.top400.mat)="numeric"

testI1.top400= testI1.bc[,rownames((top400DEG))]
testI1.top400.mat= as.matrix(testI1.top400)
class(testI1.top400.mat)="numeric"

############################ Prediction using SVM model ##########################

# Loading saved model

model.train1= readRDS("top400.SVM.train1.rds")

sink("top400.SVM.train1.txt")

#/ Testing the model on train set
p.train1 = predict(model.train1, train1.top400.mat) 

#/ Confusion Matrix for train set
train1.cm=confusionMatrix(p.train1, labels.train1, positive = “SCZ”)
train1.cm

#/ Testing the model on test set
p.testI1 = predict(model.train1, testI1.top400.mat)

#/ Confusion matrix for test set
testI1.cm=confusionMatrix(p.testI1, labels.testI1, positive = “SCZ”)
testI1.cm

sink()

#Save model and predictions

saveRDS(model.train1, "top400.SVM.train1.rds")

write.table(p.testI1, "top400.testI1.Pred.txt", sep="\t")

#/ Prediction probabilities using validation test data

pred1 = predict(model.train1, testI1.top400.mat, probability = TRUE)
ev_prob1 = attr(pred1, "probabilities")

#/ Save probablities

write.table(ev_prob1, "SVM.top400.ev.prob1.txt", sep="\t")

############################ Prediction using PAM model ##########################

#/ Setting working directory
setwd("")

# 3. Loading saved model
Model.train1= readRDS("PAM.top400.train1.rds")

# 4. Reading Delta Value

Delta = readRDS("D1.rds")

#/ Confusion matrix

sink("PAM.top400.train1.txt")

train1.cm=pamr.confusion(Model.train1, Delta, extra = F)

confusionMatrix(t(train1.cm), positive = “SCZ”)

#/ Prediction using test data
b= t(testI1.top400.mat)

p.testI1 = pamr.predict(Model.train1, b, Delta)

testI1.top400.table=table(true = labels.testI1, pred = p.testI1)

confusionMatrix(t(testI1.top400.table), positive = “SCZ”)

MCC.top400.testI1=mltools::mcc(TN=testI1.top400.table[1,1],
                      FP=testI1.top400.table[1,2],
                      TP=testI1.top400.table[2,2],
                      FN=testI1.top400.table[2,1])

MCC.top400.testI1

sink()

#/ Prediction probabilities using validation test data

ev_prob1 = pamr.predict(Model.train1, b, Delta, type = "posterior")

#/ Save test predictions
write.table(p.testI1, "PAM.top400.testI1.predictions.txt", sep="\t")

write.table(ev_prob1, "PAM.top400.ev.prob1.txt", sep="\t")

Note: Repeat the above steps for SVM and PAM models with the other train iterations and the DEGs of interest 

#############################################################################################
Cross-platform validation
#############################################################################################

# Prediction of RNA-Seq samples by SVM models 

setwd("")

train1.400DEGs.SVMradial=readRDS("train1.400DEGs.radial.rds")

setwd("")

train1.DEGs=readRDS("train1.tab.rds")

setwd("")

R1.RNA_seq=readRDS("R1.RNA_seq.rds")


setwd("")

RNA.seq.phdata=read.delim("RNA-Seq_phenodata_20_samples.txt", strip.white = T, row.names = "SampleID")

# Subset top DEGs from test data

abc=train1.DEGs[order(train1.DEGs$adj.P.Val, decreasing =F),]
top400.DEGs=abc[1:400,]
dim(top400.DEGs)

# Subset

R1.RNA_seq.df=R1.RNA_seq[rownames(top400.DEGs),]
R1.RNA_seq.final=as.matrix(t(R1.RNA_seq.df))# SVM needs columns as features

# subset test phenodata
R1.RNA.Seq.phData=RNA.seq.phdata[colnames(R1.RNA_seq.df),, drop=F]

# Factorize labels
R1.RNA_seq.flabels=factor(as.character(R1.RNA.Seq.phData$Source))

# Predict test data samples 
R1.400DEGs.predict=predict(train1.400DEGs.SVMradial, R1.RNA_seq.final)R1.400DEGs.cm=table(true = R1.RNA_seq.flabels, pred = R1.400DEGs.predict)
confusionMatrix(t(R1.400DEGs.cm), positive = “SCZ”)

MCC.top400.R1=mcc(TN=R1.400DEGs.cm[1,1],
                     FP=R1.400DEGs.cm[1,2],
                     TP=R1.400DEGs.cm[2,2],
                     FN=R1.400DEGs.cm[2,1])

MCC.top400.R1

##########################
# Save predictions
setwd("")

write.table(R1.400DEGs.predict, R1.400DEGs.predict.txt, sep="\t")

# Save data

setwd("")

sink("R1.RNA.Seq.20.samples.SVMradial.top400.DEGs.txt")

a= "R1.RNA.Seq.SVMradial.top400.DEGs"
a

confusionMatrix(R1.400DEGs.cm)

MCC.top400.R1

sink()

#####################################################################################

# Prediction of RNA-Seq samples by PAM models 

# Read PAM model
setwd("")

train1.top400.PAM=readRDS("top400.train1.model.rds")

# Read microarray DEGs

setwd("")

train1.DEGs=readRDS("train1.tab")

# Read RNA_seq batch corrected data

setwd("")

R1.RNA_Seq=readRDS("R1.RNA_seq.20.samples.rds")

# Read RNA-seq pheno data

setwd("")

RNA.seq.phdata=read.delim("RNA-Seq_phenodata_20_samples.txt", strip.white = T, row.names = "SampleID")

# Subset top DEGs from test data

abc=train1.DEGs[order(train1.DEGs$adj.P.Val, decreasing =F),]
top400.DEGs=abc[1:400,]
dim(top400.DEGs)

# Subset

R1.RNA_Seq.df=R1.RNA_Seq[rownames(top400.DEGs),]

# subset test phenodata
R1.RNA.Seq.phData=RNA.seq.phdata[colnames(R1.RNA_Seq.df),, drop=F]

# Factorize labels
R1.RNA_Seq.flabels=factor(as.character(R1.RNA.Seq.phData$Source))


# Predict test samples
Delta=0.1
top400.R1.predict=pamr.predict(train1.top400.PAM, R1.RNA_Seq.df, Delta)

top400.R1.predict.table=table(true = R1.RNA_Seq.flabels, pred = top400.R1.predict)

confusionMatrix(t(top400.R1.predict.table), positive = “SCZ”)

MCC.top400.R1=mcc(TN=top400.R1.predict.table[1,1],
                  FP=top400.R1.predict.table[1,2],
                  TP=top400.R1.predict.table[2,2],
                  FN=top400.R1.predict.table[2,1])

MCC.top400.R1

########################
## Save data

# Save predictions
setwd("")

write.table(top400.R1.predict, "R1.RNA.Seq.20.samples.PAM.top400DEGs.predictions.txt", sep="\t")

# Save data

setwd("")

sink("R1.RNA.Seq.20.samples.PAM.top400.DEGs.txt")

a= "R1.RNA.Seq.PAM.top400.DEGs"
a

confusionMatrix(top400.R1.predict.table)

MCC.top400.R1

sink()

##########
