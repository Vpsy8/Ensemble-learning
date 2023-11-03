# Ensemble learning for higher diagnostic precision in schizophrenia.
Aim of the project: To achieve higher machine learning-based classification precision for psychiatric disorders using ensemble learning.
Data type: Microarray Gene Expression Datasets and RNA-Seq data
Packages used: Refer to scripts
About the scripts:
The gene expression data was processed using the scripts shared. The numbers before the title of the script indicate the sequence in which the analysis was performed. Below are the details of the analysis performed using the R scripts.
1. Pre-processing: The script includes the creation of a meta-file from the gene expression datasets followed by a random selection of training and testing datasets. The script also includes the processing of an independent dataset and RNA-Seq samples.
2. Quantile normalization: Training data was independently normalized and the quantile targets generated for each dataset/batch were used for normalizing samples of the same dataset/batch in testing data.
3. Batch correction: Different batches in train data were batch corrected using comBat. The test was batch-corrected with batch-corrected train data as a reference. 
4. Feature selection: Differential gene expression analysis as feature selection method. Feature genes identified from each train data were used to build machine learning models.
5. Model development and test data prediction: Models were built using different numbers of feature genes and varying kernels. Models were tested using test datasets. Models and test data predictions were saved for further analysis. 
6. Prediction for independent dataset and RNA-Seq data: Quantile normalized and batch-corrected samples from the independent dataset and RNA-Seq were predicted using microarray-based machine learning models.
7. Ensemble learning: Test data prediction of SVM and PAM models were ensemble using the Boolean operator “AND”.  Ensemble models were evaluated based on precision.

Contact:
Vipul : vipulwagh31@gmail.com
