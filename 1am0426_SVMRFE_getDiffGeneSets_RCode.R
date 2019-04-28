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



## on Lenove laptop
count <- read.table("D:/miRNA_Prj/miRNA_OmicData.txt", row.names = 1, header = T, sep = "\t")
progress <- read.table("D:/miRNA_Prj/miRNA_OmicData_Design.txt", row.names = 1, header = T, sep = "\t")

top.genes <- readLines("D:/miRNA_Prj/0414_cleanedInput/top_list.txt")
length(top.genes)
datamatrix <- count

count <- datamatrix[ row.names(datamatrix) %in% top.genes, ]

dim(count)

x <- t(count)
y <- progress$Progression.value

head(x)
head(y)


#############################################################################
## assuming (believe) this is the reference fitting set in Kossenkov's paper; 
library(e1071)



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


## 
## as in the progression columns, there are two factors: yes and no
## 
## as in the progression columns, there are two factors: yes and no
rfFuncs$summary <- twoClassSummary

## 
## control RFE set 
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 5,
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
                        ##  repeats = 3,
                          returnResamp = "final", 
                          returnData = TRUE)


## 
## training control set with svm Radical
trainctrl_knn <- trainControl(classProbs= TRUE,
                          verboseIter = TRUE,
                          summaryFunction = twoClassSummary, 
                          method = "knn", 
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

rfe_rf <- rfe(t(count), progress$Progression.Logit, sizes=c(1:49),
              method="rf",
              rfeControl = ctrl, 
              metric = "ROC", 
              trControl = trainctrl,
              tuneGrid = tunegrid,
              preProc = c("center", "scale")
              )


## check RFE results; 
rfe_rf 


##########################################################################
## this step might take ~60 mins, depends on the hardware of the computer
set.seed(123)
rfe_svm <- rfe(t(count), progress$Progression.Logit, sizes=c(1:784),
              method="svmRadial",
              rfeControl = ctrl, 
              metric = "ROC", 
              trControl = trainctrl,
              tuneGrid = tunegrid,
              preProc = c("center", "scale")
              )

plot_ROC(rfe_svm, " svm-RFE model") 

## check optmized variables (genes )
rfe_rf$optVariables
rfe_svm$optVariables



printOptSet(13, rfe_rf)

printOptSet(13, rfe_cv_svmR)

printOptSet <- function(set.size, rfe_model){
  
  ###############################################################################################
  ## 0426 update print out variable set
  ## 
  # set.size = 13 # suppoe you want set-size of 10
  
  rfe.vars <- rfe_model$variables 
  
  rfe.set <- rfe.vars[rfe.vars$Variables==set.size,  ] # selects variables of set-size (= 10 here)
  
  rfe.set
  
  #use aggregate to calculate mean ranking score (under column "Overall")
  lm.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)
  
  #order from highest to low, and select first 10:
  rfe.order <- order(lm.set[, c("x")], decreasing = TRUE)[1:set.size]

  print( rfe.set[rfe.order, ] )
  
  ###############################################################################################
  
  
}


###############################################################################################
## 0426 update print out variable set
## 
set.size = 13 # suppoe you want set-size of 10

rfe.vars <- rfe_rf$variables 

rfe.set <- rfe.vars[rfe.vars$Variables==set.size,  ] # selects variables of set-size (= 10 here)

rfe.set


#use aggregate to calculate mean ranking score (under column "Overall")
lm.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)

#order from highest to low, and select first 10:
rfe.order <- order(lm.set[, c("x")], decreasing = TRUE)[1:set.size]
rfe.set[rfe.order, ]
###############################################################################################

plot_ROC(rfe_rf, "RFE Random Forest")
plot_ROC(rfe_svm, " RFE SVM model")



####################################################################### 
## based on above results, the gene set with 13 genes showed optmized results
## 
rfe_rf$optsize
rfe_rf$bestSubset

selectedIndices <- rfe_rf$pred$Variables == rfe_rf$optsize
selectedIndices
## require(pROC) 
rfe_rf$pred$Yes


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
n_genes <- 13
  
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


################################################
## Updated code at 4/26/2019
## 




## on Lenove laptop
count <- read.table("D:/miRNA_Prj/miRNA_OmicData.txt", row.names = 1, header = T, sep = "\t")
progress <- read.table("D:/miRNA_Prj/miRNA_OmicData_Design.txt", row.names = 1, header = T, sep = "\t")

top.genes <- readLines("D:/miRNA_Prj/0414_cleanedInput/top_list.txt")
length(top.genes)
datamatrix <- count
count <- datamatrix[ row.names(datamatrix) %in% top.genes, ]

dim(count)

x <- t(count)
y <- progress$Progression.value

head(x)
head(y)

## install.packages("mlbench")
## install.packages("caret")


library("mlbench")
library("caret")
library("pROC")

## control RFE set 
#####################################################
# rfFuncs <- twoClassSummary

## 
## as in the progression columns, there are two factors: yes and no
rfFuncs$summary <- twoClassSummary


rfe_repeatedcv <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 3,
                   verbose = TRUE,
                   saveDetails = TRUE,
                   returnResamp = "all")

## control RFE set 
rfe_cv <- rfeControl(functions = rfFuncs,
                     method = "cv",
                     number = 10,
                     repeats = 3,
                     verbose = TRUE,
                     saveDetails = TRUE,
                     returnResamp = "all")



################################################################## 
## training control set 

## training control set with svm Radical
trainctrl_repCV <- trainControl(classProbs= TRUE,
                               verboseIter = TRUE,
                               summaryFunction = twoClassSummary, 
                               method = "repeatedcv", 
                               number = 10 ,
                               returnResamp = "final", 
                               returnData = TRUE)


################################################################## 
## training control set with svm Radical
trainctrl_svmR <- trainControl(classProbs= TRUE,
                              verboseIter = TRUE,
                              summaryFunction = twoClassSummary, 
                              method = "svmRadical", 
                              number = 10 ,
                              returnResamp = "final", 
                              returnData = TRUE)


trainctrl_svmLinear <- trainControl(classProbs= TRUE,
                              verboseIter = TRUE,
                              summaryFunction = twoClassSummary, 
                              method = "svmLinear", 
                              number = 10 ,
                              returnResamp = "final", 
                              returnData = TRUE)

trainctrl_mars <- trainControl( classProbs = TRUE,
                                verboseIter = TRUE,
                                summaryFunction = twoClassSummary,
                                method = "earth",
                                number = 10,
                                returnData = TRUE,
                                returnResamp = "final"
                                 )


trainctrl_rf <- trainControl( classProbs = TRUE,
                                verboseIter = TRUE,
                                summaryFunction = twoClassSummary,
                                method = "rf",
                                number = 10,
                                returnData = TRUE,
                                returnResamp = "final"
                              )

trainctrl_knn <- trainControl( classProbs = TRUE,
                              verboseIter = TRUE,
                              summaryFunction = twoClassSummary,
                              method = "knn",
                              number = 10,
                              returnData = TRUE,
                              returnResamp = "final"
                              )



############################
## set tune grid
tunegrid <- expand.grid(.mtry=c(1:10))


####################################################################
## run_RFE function
## pass a ref_ctr and a train_ctr, plus count data and class to the function
## return an model
## 
run_RFE <- function(data.count, class, rfe_ctr, train_ctr, method.assign){
  
  rfe_obj <- rfe(data.count, class, sizes=c(1:49),
                 summaryFunction = twoClassSummary,
                 
                     method = method.assign,
                     rfeControl = rfe_ctr, 
                     metric = "ROC", 
                     trControl = train_ctr,
                     tuneGrid = tunegrid,
                     preProc = c("center", "scale")
  )
  
  print(method.assign)
  
  return(rfe_obj)
  
} ### end run_RFE() function;


#####################################################################################################
### 
###    Plot_ROC function; 
###     Pass a rfe model and main string argument
###
#####################################################################################################
plot_ROC <- function(rfe_model, main.str, n_genes){
  
  ########################################
  ## one argument n_genes to this function; 
  ## or, just call the best geneset by assigning the rfe_rf$optsize to the n_genes;
  if( is.null( n_genes) )
    n_genes <- rfe_model$optsize
  
  selectedIndices <- rfe_model$pred$Variables == n_genes
  
  roc_obj <- roc(rfe_model$pred$obs[selectedIndices], rfe_model$pred$yes[selectedIndices], plot=TRUE, 
                 legacy.axes=TRUE, percent=TRUE, 
                 ci=TRUE,
                 # of="se",
                 main= paste0( n_genes, main.str ),
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
  
  print(rfe_model$optVariable)
  
} ### end of plot_ROC function
################################################################################################

class <- progress$Progression.Logit

count.data <- as.data.frame(  t(count) )
count.data$class <- progress$Progression.Logit
dim(count.data)


## 1 cv plus svmLinear
set.seed(100)
rfe_cv_svmL <- run_RFE(count.data[ ,1:49], count.data$class, rfe_cv, trainctrl_svmLinear, "svmLinear")

rfe_cv_svmL$optsize
#rfe_cv_svmL

plot_ROC(rfe_cv_svmL, " genes SVM_Linear and CV")


## 2 cv plus svmRadial
set.seed(100)
rfe_cv_svmR <- run_RFE(count.data[ ,1:49], class, rfe_cv, trainctrl_svmR, "svmRadial")

rfe_cv_svmR$optsize


plot_ROC(rfe_cv_svmR, " genes SVM_Radial and CV")

## 3 cv plus random forest
set.seed(100) 
rfe_cv_rf <- run_RFE(count.data[, 1:49], class, rfe_cv, trainctrl_rf, "rf")

rfe_cv_rf$optsize
rfe_cv_rf

plot_ROC(rfe_cv_rf, " genes Random Forest and CV")

## 4 cv Plus MARS
set.seed(100)
rfe_cv_MARS <- run_RFE(count.data[ , 1:49], class, rfe_cv, trainctrl_mars, "earth")

rfe_cv_MARS$optsize


plot_ROC(rfe_cv_MARS, " genes MARS and CV")

## 5 cv plus KNN
set.seed(100)
rfe_cv_knn <- run_RFE(count.data[, 1:49], class, rfe_cv, trainctrl_knn, "knn")

rfe_cv_knn$optsize
rfe_cv_knn

plot_ROC(rfe_cv_knn, " genes KNN and CV")




# 1
rfe_cv_svmL$optsize
rfe_cv_svmR$optsize
rfe_cv_rf$optsize
rfe_cv_MARS$optsize
rfe_cv_knn$optsize






plot_ROC(rfe_cv_svmL, " genes SVM_Linear and CV")

plot_ROC(rfe_cv_svmR, " genes SVM_Radial and CV")


plot_ROC(rfe_cv_rf, " genes Random Forest and CV")


plot_ROC(rfe_cv_MARS, " genes MARS and CV")


plot_ROC(rfe_cv_knn, " genes KNN and CV")












set.seed(100)
rfe_cv_svmR <- rfe(t(count), progress$Progression.Logit, sizes=c(1:49),
                  method="svmRadial",
                  rfeControl = ctrl_cv, 
                  metric = "ROC", 
                  trControl = trainctrl_svmR,
                  tuneGrid = tunegrid,
                  preProc = c("center", "scale")
                  )

rfe_cv_svmR$optsize



set.seed(100)
rfe_cv_svmL <- rfe(t(count), progress$Progression.Logit, sizes=c(1:49),
                   method="svmLinear",
                   rfeControl = rfe_cv, 
                   metric = "ROC", 
                   trControl = trainctrl_repCV,
                   tuneGrid = tunegrid,
                   preProc = c("center", "scale")
                  )



rfe_cv_svmL$optsize
plot_ROC(rfe_cv_svmL, " genes svmLinear and CV")



set.seed(100)
rfe_cv_knn <- rfe(t(count), progress$Progression.Logit, sizes=c(1:49),
                   method="knn",
                   rfeControl = rfe_cv, 
                   metric = "ROC", 
                   trControl = trainctrl_repCV,
                   tuneGrid = tunegrid,
                   preProc = c("center", "scale")
       )

rfe_cv_knn <- rfe(t(count), progress$Progression.Logit, sizes=c(1:49),
                  #method="knn",
                  rfeControl = rfe_cv, 
                  metric = "ROC", 
                  trControl = trainctrl_repCV,
                  tuneGrid = tunegrid,
                  preProc = c("center", "scale")
)

set.seed(123)
rfe_cv_knn <- rfe(t(count), progress$Progression.Logit, sizes=c(1:49),
                  method="knn",
                  rfeControl = rfe_cv, 
                  metric = "ROC", 
                #  trControl = trainctrl_repCV,
                  tuneGrid = tunegrid,
                  preProc = c("center", "scale")
)


set.seed(123)
rfe_cv_knn <- rfe(t(count), progress$Progression.Logit, sizes=c(1:49),
                  method="knn",
                  rfeControl = rfe_cv, 
                  metric = "ROC", 
                  trControl = trainctrl_repCV,
                  #tuneGrid = tunegrid,
                  tuneLength = 5,
                  preProc = c("center", "scale")
)



rfe_cv_knn$optsize
rfe_cv_knn

plot_ROC(rfe_cv_knn, " genes KNN and CV")



n_genes <- 13
## or, just call the best geneset by assigning the rfe_rf$optsize to the n_genes;
n_genes <- rfe_cv_svmR$optsize

selectedIndices <- rfe_rf$pred$Variables == n_genes

roc_obj <- roc(rfe_rf$pred$obs[selectedIndices], rfe_rf$pred$yes[selectedIndices], plot=TRUE, 
               legacy.axes=TRUE, percent=TRUE, 
               ci=TRUE,
               # of="se",
               main= paste0( n_genes," genes SVM_Radial and CV" ),
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

rfe_cv_svmR$results[n_genes, ] 
rfe_cv_svmR

## check SVM-RFE results; 
rfe_cv_svmR$optVariables

## check optmized variables (genes )
rfe_cv_svmR$variables$

## importance <- varImp(rfe_cv_svmR, scale=FALSE)
## predictors(rfe_cv_svmR)
  
## based on above results, there are 13 genes showing 
selectedIndices <- rfe_rcv_svm$pred$Variables == rfe_rcv_svm$optsize

## require(pROC) 
rfe_rcv_svm$pred$yes
ROC = plot.roc(rfe_rcv_svm$pred$obs[selectedIndices],
               rfe_rcv_svm$pred$yes[selectedIndices], 
               ci = TRUE,
               main = "SVM repeatedCV",
               legacy.axes = TRUE,
               print.auc = TRUE
              )

# check first 6 gene sets (gene #1-6) result; 
head( rfe_rf_svm$results )

# check the last 6 group of gene setns (gene #779-784) result; 
tail( rfe_rf_svm$results )



set.seed(100)

rfe_cv_svmL <- rfe(t(count), progress$Progression.Logit, sizes=c(1:49),
                   #method="knn",
                   rfeControl = ctrl_cv, 
                   metric = "ROC", 
                   trControl = trainctrl_svmLinear,
                   tuneGrid = tunegrid,
                   preProc = c("center", "scale")
                  )

rfe_cv_svmL$optsize
rfe_cv_svmL


## plot ROC-AUC
selectedIndices <- rfe_cv_svmL$pred$Variables == rfe_cv_svmL$optsize
selectedIndices
## selectedIndices could be any digita from 1 to 784
rfe_cv_svmL$optsize
n_genes <- rfe_cv_svmL$optsize

selectedIndices <- rfe_rf_svm$pred$Variables == n_genes

roc_obj <- roc(rfe_cv_svmL$pred$obs[selectedIndices], rfe_cv_svmL$pred$yes[selectedIndices], plot=TRUE, 
               legacy.axes=TRUE, percent=TRUE, 
               ci = TRUE,
               main= paste0( n_genes," genes svm_Linear and CV" ),
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









###############################################################################################
## 0426 update print out variable set
## 
set.size = 10 # suppoe you want set-size of 10

rfe.vars <- rfe_cv_svmR$variables 

rfe.set <- rfe.vars[rfe.vars$Variables==set.size,  ] # selects variables of set-size (= 10 here)

rfe.set


#use aggregate to calculate mean ranking score (under column "Overall")
lm.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)

#order from highest to low, and select first 10:
rfe.order <- order(lm.set[, c("x")], decreasing = TRUE)[1:set.size]
rfe.set[rfe.order, ]
###############################################################################################












######################################################################
## update
library(caret)
library(ggplot2)
library(mlbench)
library(pROC)

data(Sonar)

head(Sonar)

########
count <- read.table("D:/miRNA_Prj/0414_cleanedInput/count.txt", header = T, row.names = 1, sep = "\t")
progress <- read.table("D:/miRNA_Prj/0414_cleanedInput/CountData_Design.txt", header = T, row.names = 1, sep = "\t")

head( count )
count_t <- as.data.frame( t( count ) )
head( count_t )
dim(count_t)

summary(Sonar$Class)

count_t$Class <- progress$Progressed


ctrl <- trainControl(method="cv", 
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)

rfFit <- train(Class ~ ., data=count_t, 
               method="rf", 
               preProc=c("center", "scale"), 
               metric = "ROC", 
               trControl=ctrl)

# library(pROC)
# Select a parameter setting
selectedIndices <- rfFit$pred$mtry == 39

head(rfFit$pred)

selectedIndices

rfFit$pred$obs[selectedIndices]
rfFit$pred$Yes[selectedIndices]
# rfFit$pred$M
# Plot:

rfFit$results

plot.roc(rfFit$pred$obs[selectedIndices],
         rfFit$pred$Yes[selectedIndices],
         print.auc=TRUE, 
         print.auc.x=45, 
         legacy.axes = T)




roc(rfFit$pred$obs[selectedIndices], rfFit$pred$Yes[selectedIndices], plot=TRUE, 
    legacy.axes=TRUE, percent=TRUE, 
    xlab="False Positive Percentage", 
    ylab="True Postive Percentage", 
    col="darkblue", 
    lwd=4, 
    print.auc=TRUE, 
    print.auc.x=45
    )

rfFit$ 

plot.roc(rfFit$pred$obs[selectedIndices],
         rfFit$pred$No[selectedIndices],
         AUC=T,
         legacy.axes = T)




########

head(Sonar)
head(count2)


ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T,
                     savePredictions = T)

rfFit <- train(Class ~ ., data=Sonar, method="rf", preProc=c("center", "scale"), 
               trControl=ctrl)

# Select a parameter setting
selectedIndices <- rfFit$pred$mtry == 2



## plot ROC with AUC 
g <- ggplot(rfFit$pred[selectedIndices, ], aes(m=M, d=factor(obs, levels = c("R", "M")))) + 
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc()

g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))  



## Not run: 
data(BloodBrain)

x <- scale(bbbDescr[,-nearZeroVar(bbbDescr)])
x <- x[, -findCorrelation(cor(x), .8)]
x <- as.data.frame(x)

set.seed(1)
lmProfile <- rfe(x, logBBB,
                 sizes = c(2:25, 30, 35, 40, 45, 50, 55, 60, 65),
                 rfeControl = rfeControl(functions = lmFuncs,
                                         number = 200))
set.seed(1)
lmProfile2 <- rfe(x, logBBB,
                  sizes = c(2:25, 30, 35, 40, 45, 50, 55, 60, 65),
                  rfeControl = rfeControl(functions = lmFuncs,
                                          rerank = TRUE,
                                          number = 200))

xyplot(lmProfile$results$RMSE + lmProfile2$results$RMSE  ~
         lmProfile$results$Variables,
       type = c("g", "p", "l"),
       auto.key = TRUE)

rfProfile <- rfe(x, logBBB,
                 sizes = c(2, 5, 10, 20),
                 rfeControl = rfeControl(functions = rfFuncs))

bagProfile <- rfe(x, logBBB,
                  sizes = c(2, 5, 10, 20),
                  rfeControl = rfeControl(functions = treebagFuncs))

set.seed(1)
svmProfile <- rfe(x, logBBB,
                  sizes = c(2, 5, 10, 20),
                  rfeControl = rfeControl(functions = caretFuncs,
                                          number = 200),
                  ## pass options to train()
                  method = "svmRadial")

## classification

data(mdrr)
mdrrDescr <- mdrrDescr[,-nearZeroVar(mdrrDescr)]
mdrrDescr <- mdrrDescr[, -findCorrelation(cor(mdrrDescr), .8)]

set.seed(1)
inTrain <- createDataPartition(mdrrClass, p = .75, list = FALSE)[,1]

train <- mdrrDescr[ inTrain, ]
test  <- mdrrDescr[-inTrain, ]
trainClass <- mdrrClass[ inTrain]
testClass  <- mdrrClass[-inTrain]

set.seed(2)
ldaProfile <- rfe(train, trainClass,
                  sizes = c(1:10, 15, 30),
                  rfeControl = rfeControl(functions = ldaFuncs, method = "cv"))
plot(ldaProfile, type = c("o", "g"))

postResample(predict(ldaProfile, test), testClass)


## End(Not run)

#######################################
## Parallel Processing Example via multicore

## Not run: 
library(doMC)

## Note: if the underlying model also uses foreach, the
## number of cores specified above will double (along with
## the memory requirements)
registerDoMC(cores = 2)

set.seed(1)
lmProfile <- rfe(x, logBBB,
                 sizes = c(2:25, 30, 35, 40, 45, 50, 55, 60, 65),
                 rfeControl = rfeControl(functions = lmFuncs,
                                         number = 200))



## End(Not run)
