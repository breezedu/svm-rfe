
##
##
#####################################
## 
## Feature selection with SVM-RFE
## 04/14/2019
## Jeff Du
## 
#####################################

#####################################
# install necessar packages
# library(devtools)
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("shiny")
# install_github("StatsWithR/statsr")
# BiocManager::install("sigFeature", version = "3.8")
########################################################

## this code will run RFE + cross-validation on all 784 genes, for any combination number of genes;
## 







#######################################################################
## 
## Use SVM-RFE from Caret package to 
## Selecte features (genes)
## since the input data is too small, only 43 samples
## it is relative difficult to run training and validation sets
## 
## 
#######################################################################


## install.packages("mlbench")
## install.packages("caret")


library("mlbench")
library("caret")
library("pROC")


count <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/test_miRNA/miRNA_OmicData.txt", row.names = 1, header = T, sep = "\t")

design <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/test_miRNA/miRNA_OmicData_Design.txt", row.names = 1, header = T, sep = "\t")

x <- as.data.frame( t( count ) )
y <- design$Progression.Logit

head(x)


#############################################################################
## as in the progression columns, there are two factors: yes and no
rfFuncs$summary <- twoClassSummary

## rfFuncs
## 
## control RFE set 
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 3,
                   verbose = TRUE,
                   saveDetails = TRUE,
                   returnResamp = "all")


## 
## training control set 
## 
## 1. merged progression column to the t(count) dataframe, so classProbs = true
## summary function: two class summary, based on progressin
## method, this could be cv (default) or svmRadial (for further model training)
## number = 10 folds
## 
trainctrl <- trainControl(classProbs= TRUE,
                          verboseIter = TRUE,
                          summaryFunction = twoClassSummary, 
                          method = "cv", 
                          number = 10 ,
                          returnResamp = "final", 
                          returnData = TRUE)


## 
## training control set with svm Radical
trainctrl_rf <- trainControl(classProbs= TRUE,
                             verboseIter = TRUE,
                             summaryFunction = twoClassSummary, 
                             method = "cv", 
                             number = 10 ,
                             returnResamp = "final", 
                             returnData = TRUE)


## 
## to repelicate the result, set random seed here; 
set.seed(123)

############################
## set tune grid
tunegrid <- expand.grid(.mtry=c(1:10))


##########################################################################
## this step might take ~60 mins, depends on the hardware of the computer

rfe_rf <- rfe(t(count), progress$Progressed, sizes=c(1:784),
              method="rf",
              rfeControl = ctrl, 
              metric = "ROC", 
              trControl = trainctrl,
              tuneGrid = tunegrid,
              preProc = c("center", "scale")
)


############################################################

subset <- c(1:784)

data.set <- as.data.frame( t(count) )
data.set$class <- design$Progression.Logit

rfe_rf <- rfe(data.set, data.set$class, 
              sizes=subset,
              method="svm",
              rfeControl = ctrl, 
              metric = "ROC", 
              trControl = trainctrl,
              tuneGrid = tunegrid,
              preProc = c("center", "scale")
              )


## check RFE results; 
rfe_rf 

## check optmized variables (genes )
rfe_rf$optVariables


####################################################################### 
## based on above results, the gene set with 13 genes showed optmized results
## 
selectedIndices <- rfe_rf$pred$Variables == rfe_rf$optsize

## require(pROC) 
rfe_rf$pred$yes


#####################################################
## brief ROC plot
ROC = plot.roc(rfe_rf$pred$obs[selectedIndices],
               rfe_rf$pred$yes[selectedIndices], 
               legacy.axes = TRUE
)

# check first 6 gene sets (gene #1-6) result; 
head( rfe_rf$results )

# check the last 6 group of gene setns (gene #779-784) result; 
tail( rfe_rf$results )


## 
## check optsize, to check the gene-set from the input dataset, 
## i.e. how many gene could there be, what genes are they;
## 

rfe_rf$optVariables
rfe_rf$optsize

## index the genes from the gene list;
selectedIndices <- rfe_rf$pred$Variables == rfe_rf$optsize

selectedIndices
## selectedIndices could be any digita from 1 to 784


#############################
## YOU can set any number of genes here
n_genes <- 
  
  ## or, just call the best geneset by assigning the rfe_rf$optsize to the n_genes;
  n_genes <- rfe_rf$optsize

selectedIndices <- rfe_rf$pred$Variables == n_genes

roc_obj <- roc(rfe_rf$pred$obs[selectedIndices], rfe_rf$pred$yes[selectedIndices], plot=TRUE, 
               legacy.axes=TRUE, percent=TRUE, 
               ci=TRUE,
               # of="se",
               main= paste0( n_genes," genes " ),
               xlab="False Positive Percentage", 
               ylab="True Postive Percentage", 
               col="darkblue", 
               lwd=4, 
               print.auc = T, 
               print.auc.x=45
)



legend("bottomright", 
       legend=c( paste( n_genes, "genes ROC-AUC 95% CI ", round(ci(roc_obj)[1], 2), "-", round(ci(roc_obj)[3], 2) ) ), 
       col=c("darkblue"), 
       lwd=4
)

rfe_rf$results[13, ]

auc(roc_obj)

ci(roc_obj)[1:3]


ci.coords(smooth(roc_obj), x=0.9, input = "sensitivity", ret=c("specificity", "ppv", "npv"))

ci.coords(roc_obj, x=0.9, input = "sensitivity", ret=c("specificity", "sensitivity", "ppv", "npv"))

coords(roc_obj, "local maximas", ret=c("threshold", "sens", "spec", "ppv", "npv"))

## get performance matrix specificity at 90% sensitivity, and PPV, NPV
t( coords( smooth( roc_obj ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))


coords(smooth(roc_obj), 0.5, ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                   "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity",
                                   "1-sensitivity", "npv", "ppv"))

#plot.roc.default(x = rfe_rf$pred$obs[selectedIndices], predictor = rfe_rf$pred$neg[selectedIndices],     legacy.axes = TRUE)
## Return more values:
all.values <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "npv", "ppv")

roc_results <- t(coords(roc_obj, "all", ret = all.values))

tail( roc_results)
head( roc_results)


##########################################
## training control with svmRadical



## 
## control RFE set 
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 3,
                   verbose = TRUE,
                   saveDetails = TRUE,
                   returnResamp = "all")


## 
## training control set 
trainctrl <- trainControl(classProbs= TRUE,
                          verboseIter = TRUE,
                          summaryFunction = twoClassSummary, 
                          method = "cv", 
                          number = 10 ,
                          returnResamp = "final", 
                          returnData = TRUE)


## 
## training control set with svm Radical
trainctrl_svm <- trainControl(classProbs= TRUE,
                              verboseIter = TRUE,
                              summaryFunction = twoClassSummary, 
                              method = "svmRadical", 
                              number = 10 ,
                              returnResamp = "final", 
                              returnData = TRUE)



rfe_rf_svm <- rfe(t(count), progression$Progression.Logit, sizes=c(1:784),
                  method="rf",
                  rfeControl = ctrl, 
                  metric = "ROC", 
                  trControl = trainctrl_svm,
                  tuneGrid = tunegrid,
                  preProc = c("center", "scale")
)


## check SVM-RFE results; 
rfe_rf_svm

## check optmized variables (genes )
rfe_rf_svm$optVariables


## based on above results, there are 13 genes showing 
selectedIndices <- rfe_rf_svm$pred$Variables == rfe_rf_svm$optsize

## require(pROC) 
rfe_rf_svm$pred$yes
ROC = plot.roc(rfe_rf_svm$pred$obs[selectedIndices],
               rfe_rf_svm$pred$yes[selectedIndices], 
               legacy.axes = TRUE
)

# check first 6 gene sets (gene #1-6) result; 
head( rfe_rf_svm$results )

# check the last 6 group of gene setns (gene #779-784) result; 
tail( rfe_rf_svm$results )


## plot ROC-AUC
selectedIndices <- rfe_rf_svm$pred$Variables == rfe_rf_svm$optsize
selectedIndices
## selectedIndices could be any digita from 1 to 784
rfe_rf_svm$optsize
n_genes <- 79

selectedIndices <- rfe_rf_svm$pred$Variables == n_genes

roc_obj <- roc(rfe_rf_svm$pred$obs[selectedIndices], rfe_rf_svm$pred$yes[selectedIndices], plot=TRUE, 
               legacy.axes=TRUE, percent=TRUE, 
               main= paste0( n_genes," genes " ),
               xlab="False Positive Percentage", 
               ylab="True Postive Percentage", 
               col="darkblue", 
               lwd=4, 
               print.auc = T, 
               print.auc.x=45
)





legend("bottomright", 
       legend=c( paste( n_genes, "genes SVM-RFE, ROC-AUC 95% CI ", round(ci(roc_obj)[1], 2), "-", round(ci(roc_obj)[3], 2) ) ), 
       col=c("darkblue"), 
       lwd=4
)

