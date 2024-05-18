#/ Loading libraries
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(preprocessCore)

#/ Setting working directory
setwd("")

#/ Reading data
Train1 = as.data.frame(readRDS("train1.rds"))Test1 = as.data.frame(readRDS("test1.rds"))

#/ Samples from each batch need to be quantile normalized separately

#/ Sub-setting first set of data based on GSE IDs

# Sub-setting train data

GSE18312.Train1.raw=Train1[Train1$GSE_ID=="GSE18312",]
GSE38481.Train1.raw=Train1[Train1$GSE_ID=="GSE38481",]
GSE38484.Train1.raw=Train1[Train1$GSE_ID=="GSE38484",]
GSE48072.Train1.raw=Train1[Train1$GSE_ID=="GSE48072",]
GSE54913.Train1.raw=Train1[Train1$GSE_ID=="GSE54913",]
Kum.Train1.raw=Train1[Train1$GSE_ID=="Kumarasinghe et al.",]

labels.Train1.GSE18312= GSE18312.Train1.raw[,c(1:3)]
labels.Train1.GSE38481= GSE38481.Train1.raw[,c(1:3)]
labels.Train1.GSE38484= GSE38484.Train1.raw[,c(1:3)]
labels.Train1.GSE48072= GSE48072.Train1.raw[,c(1:3)]
labels.Train1.GSE54913= GSE54913.Train1.raw[,c(1:3)]
labels.Train1.Kum= Kum.Train1.raw[,c(1:3)]

# Sub-setting Test data

GSE18312.Test1.raw=Test1[Test1$GSE_ID=="GSE18312",]
GSE38481.Test1.raw=Test1[Test1$GSE_ID=="GSE38481",]
GSE38484.Test1.raw=Test1[Test1$GSE_ID=="GSE38484",]
GSE48072.Test1.raw=Test1[Test1$GSE_ID=="GSE48072",]
GSE54913.Test1.raw=Test1[Test1$GSE_ID=="GSE54913",]
Kum.Test1.raw=Test1[Test1$GSE_ID=="Kumarasinghe et al.",]

labels.Test1.GSE18312= GSE18312.Test1.raw[,c(1:3)]
labels.Test1.GSE38481= GSE38481.Test1.raw[,c(1:3)]
labels.Test1.GSE38484= GSE38484.Test1.raw[,c(1:3)]
labels.Test1.GSE48072= GSE48072.Test1.raw[,c(1:3)]
labels.Test1.GSE54913= GSE54913.Test1.raw[,c(1:3)]
labels.Test1.Kum= Kum.Test1.raw[,c(1:3)]

############################ Quantile Normalization ############################

################################ Train 1 ########################################

# 1. GSE18312

# a. Train data
Train1.GSE18312.df=t(GSE18312.Train1.raw)
class(Train1.GSE18312.df)="numeric"
Train1.GSE18312.df= Train1.GSE18312.df[-c(1:3),]
Train1.GSE18312.df=log2(Train1.GSE18312.df)

# Check for rows with negative values
abc=apply(Train1.GSE18312.df, 1, function(row) any(row<0))
length(which(abc)) # 0 rows with negative values

Train1.GSE18312.df[Train1.GSE18312.df<0]=0

# Normalize
Train1.GSE18312.qnorm=normalize.quantiles(Train1.GSE18312.df, copy = T)
boxplot(Train1.GSE18312.qnorm)

# Save targets
Train1.GSE18312.targets= normalize.quantiles.determine.target(Train1.GSE18312.qnorm,target.length=NULL,subset=NULL)

# Add column and row names from original data
dimnames(Train1.GSE18312.qnorm)=dimnames(Train1.GSE18312.df)

# Add back the pheno-data

Train1.GSE18312.qnorm.ph= rbind(Train1.GSE18312.qnorm, t(labels.Train1.GSE18312))
#View(Train1.GSE18312.qnorm.ph)

# b. Test data

Test1.GSE18312.df=t(GSE18312.Test1.raw)
class(Test1.GSE18312.df)="numeric"
Test1.GSE18312.df= Test1.GSE18312.df[-c(1:3),]
Test1.GSE18312.df=log2(Test1.GSE18312.df)

# Check for rows with negative values
abc=apply(Test1.GSE18312.df, 1, function(row) any(row<0))
length(which(abc)) # 0 rows with negative values

Test1.GSE18312.df[Test1.GSE18312.df<0]=0

# Normalize
Test1.GSE18312.qnorm <- normalize.quantiles.use.target(Test1.GSE18312.df,Train1.GSE18312.targets,copy=TRUE,subset=NULL)

boxplot(Test1.GSE18312.qnorm)

# Add column and row names from original data
dimnames(Test1.GSE18312.qnorm)=dimnames(Test1.GSE18312.df)

# Add back the pheno-data

Test1.GSE18312.qnorm.ph= rbind(Test1.GSE18312.qnorm, t(labels.Test1.GSE18312))
#View(Test1.GSE18312.qnorm.ph)

# (Optional) Confirm if the Test data has been normalized with respect to the train data

#Alldata.GSE18312=cbind(Train1.GSE18312.qnorm, Test1.GSE18312.qnorm)
#boxplot(Alldata.GSE18312)

# 2.GSE38481
# a. Train data

Train1.GSE38481.df=t(GSE38481.Train1.raw)
class(Train1.GSE38481.df)="numeric"
Train1.GSE38481.df= Train1.GSE38481.df[-c(1:3),]
Train1.GSE38481.df=log2(Train1.GSE38481.df)

# Check for rows with negative values
abc=apply(Train1.GSE38481.df, 1, function(row) any(row<0))
length(which(abc)) # 400 rows with negative values

Train1.GSE38481.df[Train1.GSE38481.df<0]=0

# Normalize
Train1.GSE38481.qnorm=normalize.quantiles(Train1.GSE38481.df, copy = T)
boxplot(Train1.GSE38481.qnorm)

# Save targets
Train1.GSE38481.targets= normalize.quantiles.determine.target(Train1.GSE38481.qnorm,target.length=NULL,subset=NULL)

# Add column and row names from original data
dimnames(Train1.GSE38481.qnorm)=dimnames(Train1.GSE38481.df)

# Add back the pheno-data

Train1.GSE38481.qnorm.ph= rbind(Train1.GSE38481.qnorm, t(labels.Train1.GSE38481))
#View(Train1.GSE38481.qnorm.ph)

# b. Test data
Test1.GSE38481.df=t(GSE38481.Test1.raw)
class(Test1.GSE38481.df)="numeric"
Test1.GSE38481.df= Test1.GSE38481.df[-c(1:3),]
Test1.GSE38481.df=log2(Test1.GSE38481.df)

# Check for rows with negative values
abc=apply(Test1.GSE38481.df, 1, function(row) any(row<0))
length(which(abc)) # 101 rows with negative values

Test1.GSE38481.df[Test1.GSE38481.df<0]=0

# Normalize

Test1.GSE38481.qnorm <- normalize.quantiles.use.target(Test1.GSE38481.df,Train1.GSE38481.targets,copy=TRUE,subset=NULL)

boxplot(Test1.GSE38481.qnorm)

# Add column and row names from the original data
dimnames(Test1.GSE38481.qnorm)=dimnames(Test1.GSE38481.df)

# Add back the pheno-data

Test1.GSE38481.qnorm.ph= rbind(Test1.GSE38481.qnorm, t(labels.Test1.GSE38481))
#View(Test1.GSE38481.qnorm.ph)

# (Optional) Confirm if the Test data has been normalized with respect to the train data

#Alldata.GSE38481=cbind(Train1.GSE38481.qnorm, Test1.GSE38481.qnorm)
#boxplot(Alldata.GSE38481)

# 3. GSE38484
# a. Train data

Train1.GSE38484.df=t(GSE38484.Train1.raw)
class(Train1.GSE38484.df)="numeric"
Train1.GSE38484.df= Train1.GSE38484.df[-c(1:3),]
Train1.GSE38484.df=log2(Train1.GSE38484.df)

# Check for rows with negative values
abc=apply(Train1.GSE38484.df, 1, function(row) any(row<0))
length(which(abc)) # 1903 rows with negative values

Train1.GSE38484.df[Train1.GSE38484.df<0]=0

# Normalize
Train1.GSE38484.qnorm=normalize.quantiles(Train1.GSE38484.df, copy = T)
boxplot(Train1.GSE38484.qnorm)

# Save targets
Train1.GSE38484.targets= normalize.quantiles.determine.target(Train1.GSE38484.qnorm,target.length=NULL,subset=NULL)

# Add column and row names from the original data
dimnames(Train1.GSE38484.qnorm)=dimnames(Train1.GSE38484.df)

# Add back the pheno-data

Train1.GSE38484.qnorm.ph= rbind(Train1.GSE38484.qnorm, t(labels.Train1.GSE38484))
#View(Train1.GSE38484.qnorm.ph)

# b. Test data
Test1.GSE38484.df=t(GSE38484.Test1.raw)
class(Test1.GSE38484.df)="numeric"
Test1.GSE38484.df= Test1.GSE38484.df[-c(1:3),]
Test1.GSE38484.df=log2(Test1.GSE38484.df)

# Check for rows with negative values
abc=apply(Test1.GSE38484.df, 1, function(row) any(row<0))
length(which(abc)) # 1276 rows with negative values

Test1.GSE38484.df[Test1.GSE38484.df<0]=0

# Normalize

Test1.GSE38484.qnorm <- normalize.quantiles.use.target(Test1.GSE38484.df,Train1.GSE38484.targets,copy=TRUE,subset=NULL)

boxplot(Test1.GSE38484.qnorm)

# Add column and row names from original data
dimnames(Test1.GSE38484.qnorm)=dimnames(Test1.GSE38484.df)

# Add back the pheno-data

Test1.GSE38484.qnorm.ph= rbind(Test1.GSE38484.qnorm, t(labels.Test1.GSE38484))
#View(Test1.GSE38484.qnorm.ph)

# (Optional) Confirm if the Test data has been normalized with respect to the train data

#Alldata.GSE38484=cbind(Train1.GSE38484.qnorm, Test1.GSE38484.qnorm)
#boxplot(Alldata.GSE38484)

# 4. GSE48072
# a. Train data

Train1.GSE48072.df=t(GSE48072.Train1.raw)
class(Train1.GSE48072.df)="numeric"
Train1.GSE48072.df= Train1.GSE48072.df[-c(1:3),]
Train1.GSE48072.df=log2(Train1.GSE48072.df)

# Check for rows with negative values
abc=apply(Train1.GSE48072.df, 1, function(row) any(row<0))
length(which(abc)) # 1510 rows with negative values

Train1.GSE48072.df[Train1.GSE48072.df<0]=0

# Normalize
Train1.GSE48072.qnorm=normalize.quantiles(Train1.GSE48072.df, copy = T)
boxplot(Train1.GSE48072.qnorm)

# Save targets
Train1.GSE48072.targets= normalize.quantiles.determine.target(Train1.GSE48072.qnorm,target.length=NULL,subset=NULL)

# Add column and row names from the original data
dimnames(Train1.GSE48072.qnorm)=dimnames(Train1.GSE48072.df)

# Add back the pheno-data

Train1.GSE48072.qnorm.ph= rbind(Train1.GSE48072.qnorm, t(labels.Train1.GSE48072))
#View(Train1.GSE48072.qnorm.ph)

# b. Test data
Test1.GSE48072.df=t(GSE48072.Test1.raw)
class(Test1.GSE48072.df)="numeric"
Test1.GSE48072.df= Test1.GSE48072.df[-c(1:3),]
Test1.GSE48072.df=log2(Test1.GSE48072.df)

# Check for rows with negative values
abc=apply(Test1.GSE48072.df, 1, function(row) any(row<0))
length(which(abc)) # 866 rows with negative values

Test1.GSE48072.df[Test1.GSE48072.df<0]=0

# Normalize
Test1.GSE48072.qnorm <- normalize.quantiles.use.target(Test1.GSE48072.df,Train1.GSE48072.targets,copy=TRUE,subset=NULL)

boxplot(Test1.GSE48072.qnorm)

# Add column and row names from the original data
dimnames(Test1.GSE48072.qnorm)=dimnames(Test1.GSE48072.df)

# Add back the pheno-data

Test1.GSE48072.qnorm.ph= rbind(Test1.GSE48072.qnorm, t(labels.Test1.GSE48072))
#View(Test1.GSE48072.qnorm.ph)

# (Optional) Confirm if the Test data has been normalized with respect to the train data

#Alldata.GSE48072=cbind(Train1.GSE48072.qnorm, Test1.GSE48072.qnorm)
#boxplot(Alldata.GSE48072)

# 5. GSE54913
# a. Train data

Train1.GSE54913.df=t(GSE54913.Train1.raw)
class(Train1.GSE54913.df)="numeric"
Train1.GSE54913.df= Train1.GSE54913.df[-c(1:3),]
Train1.GSE54913.df=log2(Train1.GSE54913.df)

# Check for rows with negative values
abc=apply(Train1.GSE54913.df, 1, function(row) any(row<0))
length(which(abc)) # 170 rows with negative values

Train1.GSE54913.df[Train1.GSE54913.df<0]=0

# Normalize
Train1.GSE54913.qnorm=normalize.quantiles(Train1.GSE54913.df, copy = T)
boxplot(Train1.GSE54913.qnorm)

# Save targets
Train1.GSE54913.targets= normalize.quantiles.determine.target(Train1.GSE54913.qnorm,target.length=NULL,subset=NULL)

# Add column and row names from original data
dimnames(Train1.GSE54913.qnorm)=dimnames(Train1.GSE54913.df)

# Add back the pheno-data

Train1.GSE54913.qnorm.ph= rbind(Train1.GSE54913.qnorm, t(labels.Train1.GSE54913))
#View(Train1.GSE54913.qnorm.ph)

# b. Test data
Test1.GSE54913.df=t(GSE54913.Test1.raw)
class(Test1.GSE54913.df)="numeric"
Test1.GSE54913.df= Test1.GSE54913.df[-c(1:3),]
Test1.GSE54913.df=log2(Test1.GSE54913.df)

# Check for rows with negative values
abc=apply(Test1.GSE54913.df, 1, function(row) any(row<0))
length(which(abc)) # 4 rows with negative values

Test1.GSE54913.df[Test1.GSE54913.df<0]=0

# Normalize

Test1.GSE54913.qnorm <- normalize.quantiles.use.target(Test1.GSE54913.df,Train1.GSE54913.targets,copy=TRUE,subset=NULL)

boxplot(Test1.GSE54913.qnorm)

# Add column and row names from original data
dimnames(Test1.GSE54913.qnorm)=dimnames(Test1.GSE54913.df)

# Add back the pheno-data

Test1.GSE54913.qnorm.ph= rbind(Test1.GSE54913.qnorm, t(labels.Test1.GSE54913))
#View(Test1.GSE54913.qnorm.ph)

# (Optional) Confirm if the Test data has been normalized with respect to the train data

#Alldata.GSE54913=cbind(Train1.GSE54913.qnorm, Test1.GSE54913.qnorm)
#boxplot(Alldata.GSE54913)

# 6. KumaraSinghe et. al
# a. Train data

Train1.Kum.df=t(Kum.Train1.raw)
class(Train1.Kum.df)="numeric"
Train1.Kum.df= Train1.Kum.df[-c(1:3),]
Train1.Kum.df=log2(Train1.Kum.df)

# Check for rows with negative values
abc=apply(Train1.Kum.df, 1, function(row) any(row<0))
length(which(abc)) # 2521 rows with negative values

Train1.Kum.df[Train1.Kum.df<0]=0

# Normalize
Train1.Kum.qnorm=normalize.quantiles(Train1.Kum.df, copy = T)
boxplot(Train1.Kum.qnorm)

# Save targets
Train1.Kum.targets= normalize.quantiles.determine.target(Train1.Kum.qnorm,target.length=NULL,subset=NULL)

# Add column and row names from the original data
dimnames(Train1.Kum.qnorm)=dimnames(Train1.Kum.df)

# Add back the pheno-data

Train1.Kum.qnorm.ph= rbind(Train1.Kum.qnorm, t(labels.Train1.Kum))
#View(Train1.Kum.qnorm.ph)

# b. Test data
Test1.Kum.df=t(Kum.Test1.raw)
class(Test1.Kum.df)="numeric"
Test1.Kum.df= Test1.Kum.df[-c(1:3),]
Test1.Kum.df=log2(Test1.Kum.df)

# Check for rows with negative values
abc=apply(Test1.Kum.df, 1, function(row) any(row<0))
length(which(abc)) # 512 rows with negative values

Test1.Kum.df[Test1.Kum.df<0]=0

# Normalize
Test1.Kum.qnorm <- normalize.quantiles.use.target(Test1.Kum.df,Train1.Kum.targets,copy=TRUE,subset=NULL)

boxplot(Test1.Kum.qnorm)

# Add column and row names from the original data
dimnames(Test1.Kum.qnorm)=dimnames(Test1.Kum.df)

# Add back the pheno-data

Test1.Kum.qnorm.ph= rbind(Test1.Kum.qnorm, t(labels.Test1.Kum))
#View(Test1.Kum.qnorm.ph)

# (Optional) Confirm if the Test data has been normalized with respect to the train data

#Alldata.Kum=cbind(Train1.Kum.qnorm, Test1.Kum.qnorm)
#boxplot(Alldata.Kum)

#/ Join all q-normalize data-sets for train and Test 

Train1.qnorm=cbind(Train1.GSE18312.qnorm.ph, Train1.GSE38481.qnorm.ph, Train1.GSE38484.qnorm.ph, Train1.GSE48072.qnorm.ph,Train1.GSE54913.qnorm.ph, Train1.Kum.qnorm.ph)
#View(Train1.qnorm)
#colnames(Train1.qnorm)

Test1.qnorm=cbind(Test1.GSE18312.qnorm.ph, Test1.GSE38481.qnorm.ph, Test1.GSE38484.qnorm.ph, Test1.GSE48072.qnorm.ph, Test1.GSE54913.qnorm.ph, Test1.Kum.qnorm.ph)
#View(Test1.qnorm)
#colnames(Test1.qnorm)

# save q-normalize train and Test data

Train1.qnorm.df=as.data.frame(Train1.qnorm)
saveRDS(Train1.qnorm.df, "train1.qnorm.rds")

Test1.qnorm.df=as.data.frame(Test1.qnorm)
saveRDS(Test1.qnorm.df, "test1.qnorm.rds")

# Note: Repeat this analysis for all the iterations of training and testing datasets.

##########################################################
# Quantile normalization of Independent Dataset (GSE27383)
###########################################################
# Note: The independent dataset is pre-processed using Q-normalized and batch-corrected train data. Hence this step should be performed after batch correction of train data.

#/ Setting working directory

setwd("")

#/ Reading files 

# a. GSE27383

GSE27383.I1=readRDS("GSE27383.samples.final.1.rds")

# b. batch-corrected train data

train1.bc= as.data.frame(readRDS("train1.batch-corrected.rds"))

labels.train1= train1.bc[c(6859:6861),]
train1.bc.mod= train1.bc[-c(6859:6861),]
train1.bc.mat= as.matrix(train1.bc.mod)
class(train1.bc.mat)="numeric"

GSE27383.ordered= GSE27383.I1[rownames(train1.bc.mod),]

GSE27383.ordered[is.na(GSE27383.ordered)] <- 1

rownames(GSE27383.ordered)= rownames(train1.bc.mod)

GSE27383.ordered=as.matrix(GSE27383.ordered)
class(GSE27383.ordered)="numeric"

#/ Normalization

# Save targets
train1.bc.targets <- normalize.quantiles.determine.target(train1.bc.mat,
                                                           target.length=NULL,subset=NULL)

# Normalize quantile of test dataset 

GSE27383.Qnorm <- normalize.quantiles.use.target(GSE27383.ordered,train1.bc.targets,
                                                 copy=TRUE,subset=NULL)

# The above function does not retain dim names
dimnames(GSE27383.Qnorm)=dimnames(GSE27383.ordered)

#/ Saving files

saveRDS(GSE27383.Qnorm, "I1.GSE27383.Qnorm.rds")

# NOTE: Repeat for all iterations of the independent dataset


##########################################################
# Quantile normalization of RNA-Seq samples
##########################################################

# Note: The RNA-Seq samples are pre-processed using Q-normalized and batch-corrected train data. Hence, this step should be performed after batch correction of train data.


# Read microarray train data
setwd("")

train1.bcr.data=readRDS("train1.batch-corrected")

# Read RNA-Seq data (TPM values)
setwd("")

RNA.Seq.final=readRDS("RNA_Seq.20.samples.final.rds")

# pheno-info of microarray data used in train data

setwd("")

phdata=read.delim("Phenodata_Metafile_6 Batches.txt", strip.white=TRUE)

phdata.df=phdata[c(2,5,6,9)]# Keep relevant info
phdata.df=column_to_rownames(phdata.df,"Sample_ID") 

# Read RNA-seq pheno info
setwd("")

RNA.seq.phdata=read.delim("RNA-Seq_phenodata_20_samples.txt", strip.white = T, row.names = "SampleID")
RNA.seq.phdata.df=RNA.seq.phdata[c(1:3)]

## Quantile normalization

# Have row names in the same order for train and test data

RNA.Seq.final.order= RNA.Seq.final[rownames(train1.bcr.data),]
RNA.Seq.final.order=as.matrix(RNA.Seq.final.order)
class(RNA.Seq.final.order)="numeric"

# Save targets
train1.bcr.targets <- normalize.quantiles.determine.target(train1.bcr.data,
                                                           target.length=NULL,subset=NULL)

# Normalize quantile of RNA-Seq dataset (test)

RNA_seq.Qnorm <- normalize.quantiles.use.target(RNA.Seq.final.order,train1.bcr.targets,
                                               copy=TRUE,subset=NULL)
# The above function does not retain dim names
dimnames(RNA_seq.Qnorm)=dimnames(RNA.Seq.final.order)

# Save file

# Repeat the normalization of RNA-Seq samples using all the iterations of training data. 
