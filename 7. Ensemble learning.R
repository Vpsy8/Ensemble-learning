#/ Loading libraries
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(caret)
library(mccr)
library(stats)
library(mlbench)
library(mltools)

#/ Setting working directory

setwd("")

################################# I1 #######################################

#/ Loading data

# 1. SVM test predictions

I1.top400.SVM=read.table("top400.testI1.Pred.txt", strip.white = T)

# 2. PAM test predictions

I1.top400.PAM=read.table("PAM.top400.testI1.predictions.txt", strip.white = T)

#/ Reading test data for sample IDs

test1.bc = readRDS("testI1.batch-corrected.rds")

#/ Subset phenodata

test1.bc.t= as.data.frame(t(test1.bc))
test1.bc.pheno= test1.bc.t[c(6859:6861)]

#/ Join predictions by ML models and true labels

I1.400.PAM=as.character(I1.top400.PAM$x)
I1.400.SVM=as.character(I1.top400.SVM$x)

I1.ensemble.df=as.data.frame(cbind( I1.400.PAM, I1.400.SVM))

#/ Apply boolean operator "AND"

I1.ensemble=I1.ensemble.df %>% mutate(Ensemble =
                                        case_when(I1.400.PAM == "SCZ" & I1.400.SVM=="SCZ" ~ "SCZ", 
                                                  I1.400.PAM == "SCZ" & I1.400.SVM=="CNT" ~ "CNT",
                                                  I1.400.PAM == "CNT" & I1.400.SVM=="SCZ" ~ "CNT",
                                                  I1.400.PAM == "CNT" & I1.400.SVM=="CNT" ~ "CNT"))

#/ Adding sample IDs to the ensemble data

I1.ensemble.phdata=cbind(test1.bc.pheno, I1.ensemble)

saveRDS(I1.ensemble.phdata, "I1.ensemble.rds")

#/ Predict ensemble results

I1.test.labels=as.factor(test1.bc.pheno$Source)
I1.ensemble.predicted=factor(as.character(I1.ensemble$Ensemble))

#/ Predict output

sink("Ensemble.top400.I1.txt")

I1.ensemble.cm=table(true = I1.test.labels, pred = I1.ensemble.predicted)
confusionMatrix(t(I1.ensemble.cm), positive = “SCZ”)

#/ MCC

I1.ensemble.MCC=mcc(TP=I1.ensemble.cm[1,1],
                    FN=I1.ensemble.cm[1,2],
                    TN=I1.ensemble.cm[2,2],
                    FP=I1.ensemble.cm[2,1])
I1.ensemble.MCC

sink()

# NOTE: Repeat the steps for all other iterations
# Similarly, ensemble learning for RNA-Seq samples can be performed
