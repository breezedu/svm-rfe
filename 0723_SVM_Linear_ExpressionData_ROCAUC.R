


#############################################################################
## 
## Feature selection with SVM-RFE
## 04/14/2019 - 05/09/2019
## Jeff Du
##  
#############################################################################
## 
## use caret package to build classfication model
## use the DESeq2 generated significant list as pre-selected genes
## use pROC to plot ROC for training set
## 
#############################################################################
## 
#############################################################################
# install necessar packages
# library(devtools)
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("shiny")
# install_github("StatsWithR/statsr")
# BiocManager::install("sigFeature", version = "3.8")
# install.packages("mlbench")
# install.packages("caret")
# install.packages("ellipse")   ## for feature plot; 
#############################################################################

### load packages
### library(mlbench)
library(caret)
library(ggplot2)
# library(ellipse) ## one of the dependencies for pROC package
library(pROC)
library(doParallel)


set.seed(12345)

####################################################
## Data manipulation and pre-processing
## #1 deal with missing/NA data
## #2 remove constant rows
## #3 optimize training set and testing set
####################################################


########################################################################################################
## data of RNAseq for each tumor type from TCGA
## count <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/SVM_RFE_Prj/svm_rfe_prj/THCA_RNA_Set.RnaSeq_Transcript.Genes.txt",
##                    header = TRUE,
##                    row.names = 1,
##                    sep = "\t"
##                    )
########################################################################################################
## design <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/SVM_RFE_Prj/svm_rfe_prj/THCA_RNA_Set.RnaSeq_Transcript.Genes_Design.txt",
##                     header = TRUE,
##                     row.names = 1,
##                     sep = "\t"
##                     )
########################################################################################################


########################################################################################
## 
## mRNA data 
## 
## #1 input a 49 * 43 matrix, which is the count table of a mRNA dataset
## originally, there are 784 genes, after deseq2 One-Way-Test, only 49 of them showed
## slighly diffrencial expression;
## 
## #2 a design talbe, within the table there are several columns for survival plot
## and the classification ml model works on progression column;
## 
## #3 the goal is to trains a model to predict progrsssion Yes/No
## 
########################################################################################



########################################################
## Section 0: save ROC plot
## 
## function plot ROC cure and save as jpg file
## 
## pass an roc_object together with routine and file name
## Plot ROC curve, add AUC, 95% CI, save to target folder
## Also, pass the ML algorithm to the legend;

PlotROCSaveJPG <- function(roc.obj, path, file.name, ml.algo){
  
  jpeg(paste0(path, file.name, "_ROC_with_AUC.jpg"), 
       width = 8, height = 6, units = "in", res= 400)
  
  ## We can calculate the area under the curve...
  
  plot( 
    roc.obj,
    legacy.axes=TRUE, percent=TRUE, 
    main= ml.algo,
    xlab="False Positive Percentage", 
    ylab="True Postive Percentage", 
    col="darkblue", lwd=4, 
    print.auc=TRUE
  )  # end plotting
  
  # the 95% CI could be printed directly by print.auc=T; 
  #ci.95 <- paste0("95% CI = ", round(ci.auc(rocobj)[1], digits = 2), " - ", round(ci.auc(rocobj)[3], digits = 2) ) 
  
  legend("bottomright", 
         legend=c( paste("model:", ml.algo) ), 
         col=c("darkblue", "darkblue"), 
         lwd=4)
  
  dev.off()
  
} ## End PlotROCSaveJPG() function  
#####################################################################################################################



########################################################################################
## Section I 
## read in the count data and the clinical information table; 
## 
########################################################################################

count <- read.table("C:/Users/dug/OneDrive - QIAGEN GmbH/SVM_RFE_Prj/0507CleanedData/Count.txt", row.names = 1, header = T, sep = "\t")
design <- read.table("C:/Users/dug/OneDrive - QIAGEN GmbH/SVM_RFE_Prj/0507CleanedData/Count_Design.txt", row.names = 1, header = T, sep = "\t")
top.genes <- readLines("C:/Users/dug/OneDrive - QIAGEN GmbH/SVM_RFE_Prj/0507CleanedData/Deseq2_top45rows.txt")


#########
setwd("D:/miRNA_Prj/0507CleanedData/")

setwd("C:/Users/Jeff/OneDrive - QIAGEN GmbH/SVM_RFE_Prj/0507CleanedData/")

count <- read.table("Count.txt", row.names = 1, header = T, sep = "\t")
design <- read.table("Count_Design.txt", row.names = 1, header = T, sep = "\t")
top.genes <- readLines("Deseq2_top45rows.txt")


top.genes[1:10] 
length(top.genes) 


head(count)
colnames(count)


# rownames.ori <- rownames(clinical.data, do.NULL = T, prefix = "X")
# 
# rownames.new <- paste0('X', rownames.ori)
# row.names(clinical.data) <- rownames.new

count[1:5, 1:5]


############################################################
## remove rows with constant values (zeros)
## for mRNA data, remove rows with varience less than 0.0001
## 

datamatrix <- count

dim(datamatrix)

constantRow <- row.names(as.matrix(which(apply(datamatrix, MARGIN = 1, function(x) var(x) < 0.001) ) ) ) 

#exclude rows with constant value: 
count <- datamatrix[!row.names(datamatrix) %in% constantRow,] 

count <- datamatrix[ row.names(datamatrix) %in% top.genes, ]

print("These variables get constant values, thus got droped during PCA scale:") 
constantRow


############################################################
## within the input count table, each row represents a gene
## in the model, it will consider each column as one feature
## thus we need to transfer the per-row-a-gene into per-col-a gene
## use t() fucntion to directly transport the input table
my_data <- as.data.frame( t(count) )

dim(my_data)


##########################################
# > dim(my_data)
# [1] 41 45
# > dim(clinical.data)
# [1] 41 13
## [1] 43 49 ## 43 samples 49 genes

## for RnaSeq my_data$class <- design$sample.type

## for mRNA data: 

my_data$class <- design$Progressed
my_data$class <- ifelse(my_data$class==1, "yes", "no")

my_data[1:5, 1:5]



sort( row.names(my_data) )
sort( row.names(clinical.data) )


## in case there are multi factors in the class columns, remove them;
## subset the data, only keep primiary tumor and solid normal
my_data <- subset( my_data, class!='Metastatic') # my_data$class != "Metastatic", ]

## check factor levels left
summary(my_data$class) 

## drop the empty Metastatic level
my_data$class <- factor(my_data$class)




####################################################################################
## Section II 
## Create the training and test datasets
####################################################################################

set.seed(1234)

###################################################################
## Section II, Step 1: Get row numbers for the training data
## The limitation in this case is: dataset too small
## So, we split the original data into 75% training + 25% testing
###################################################################

trainRowNumbers <- createDataPartition(my_data$class, p=0.60, list=FALSE)

set.seed(321)
testRowNumbers <- createDataPartition(my_data$class, p=0.55, list=FALSE)
testData <- my_data[testRowNumbers, ]

###################################################################
## Section II, Step 2: Create the training  dataset
trainData <- my_data[trainRowNumbers,]

###################################################################
## Section II, Step 3: Create the test dataset
#testData <- my_data[-trainRowNumbers,]

dim(trainData)
dim(testData)

# Store X and Y for later use.
x = trainData[, 1:45]
y = trainData$class

x_test = testData[ , 1:45]
y_test = testData$class
dim(x)
dim(x_test)

x$class
# x$class should be NULL


###############################################################################################
## Section III, data, mamipulation, normalization and scalling
###############################################################################################


########################################################################
## 
## convert all the numeric variables to range between 0 and 1, 
## by setting method=range in preProcess(). 
## preProcess_range_model <- preProcess(trainData, method='range')
## trainData <- predict(preProcess_range_model, newdata = trainData)
## 
## Append the Y variable
## trainData$class <- y
########################################################################


########################################################################
## Section III, step 1 
## remove constant columns
## for top 500/1000 significant genes, this step could be ignored; 

x <- x[,apply(x, 2, var, na.rm=TRUE) > 0.0001]
dim(x)


########################################################################
## Section III, step 2
## convert all the numeric variables to range between 0 and 1, by setting method=range in preProcess(). 

trainData_range_model <- preProcess( x, method = c("center", "scale") )

trainData <- predict( trainData_range_model, newdata = x) 
trainData$class <- y

head(trainData)

########################################################################
## Section III, step 3
## convert all test data 

dim(x_test)

testData_range_model <- preProcess( x_test, method = c("center", "scale") )
testData <- predict( testData_range_model, newdata = x_test) 
testData$class <- y_test
# check head of the test data
head(testData)

## 
## check dim for both train and test data
dim( trainData )
dim( testData )


########################################################################
## Section IV briefly check the feature weight and density
## call featurePlot() function 
######################################################################## 
##
## 
## box plot, to see significantly differential expression genes 
featurePlot(x = trainData[, 1:15], 
            y = trainData$class, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free"))
)


## density plot, to visulize more important variables 
featurePlot(x = trainData[, 1:15], 
            y = trainData$class, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))



##########################################################################################################  
##### 
##### Section V: RFE 
##### The most important part of this code
##### 
##########################################################################################################  


##########################################################################################################  
##  RFE works in 3 broad steps:
##    
##  Step 1: Build a ML model on a training dataset and estimate the feature importances on the test dataset.
##          Here in this case, the features are genes; 
## 
##  Step 2: Keeping priority to the most important variables, 
##          iterate through by building models of given subset sizes, 
##          that is, subgroups of most important predictors determined from step 1. 
##          Ranking of the predictors is recalculated in each iteration.
##  
##  Step 3: The model performances are compared across different subset sizes 
##          to arrive at the optimal number and list of final predictors.
##  
##  Stop 4: It can be implemented using the rfe() function and you have the flexibility 
##          to control what algorithm rfe uses and how it cross validates by defining the rfeControl().
##########################################################################################################  


set.seed(1234)

options(warn=-1)

subsets <- c(1:45)

dim(trainData)

##########################################################################################################  
## Section V, Step 1
## build a control model on training dataset; estimate the feautre importances 

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",      ## cross-validation
                   number = 10,                ## fold of cross validation
                   repeats = 3,                ## repeats of cv
                   verbose = FALSE
)


#########################################
## Recursive Feature Elimination rfe() 

lmProfile <- rfe(x=trainData[, 1:45], 
                 y=trainData$class,
                 sizes = subsets,        ## model sizes (the number of most important genes to choose)
                 rfeControl = ctrl       ## use rfeControl output as reference: algorithm and cv to use;
)

## check lm_profile results
lmProfile 

## check the optimized variables (gens in this case)
lmProfile$optVariables

#####################################################
## tested by 4/15/2019 evening;
##################################################### 



################################################################################ 
#### Section V, step 2
#### 
#### train() the model and interpret the results
#### 
## Set the seed for reproducibility
set.seed(100)

################################################################################
## Train the model using MARS and predict on the training data itself.
## train() will cross validate the model 'earth' = Multivariate Adaptive Regression Splines (MARS); 
model_mars = train(class ~ ., 
                   data=trainData, 
                   method='earth'
                   #savePredictions = T
)

## check the results generated by train() with MARS algo:
model_mars
# model_mars$pred
# lmProfile$optsize

########################################
## fit the predicted model_mars
fitted <- predict(model_mars)

#######################################
## plot fit accuracies 
plot(model_mars, main="Model Accuracies with MARS") 

## 
## compute variable importance (genes more important to the model)

varimp_mars <- varImp(model_mars)
plot(varimp_mars, main="Variable Importance with MARS")
varimp_mars$importance 

## > varimp_mars$importance
## Overall
## PLAU                      100.00000
## CCL14                      81.73621
## IFITM1                     69.61859
## IL1RN                      42.42643
## ABCB1                       0.00000
## BST1                        0.00000
## C4B                         0.00000
## CCL17                       0.00000

##########################################################################################

##########################################################################################################  
#### Section V, Step 3  
#### Define the training control

fitControl <- trainControl(
  method = 'repeatedcv',           # k-fold cross validation
  number = 15,                     # number of folds
  repeats = 10,                    # number of repeats
  savePredictions = T,             # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary  #,  # results summary function
  # savePredictions = 'final'
) 


#################################################################################################
#### Section V, step 4
#### Tune hyper parameters by setting tuneLength

set.seed(12345)

#############################
#### Step 4.1 choose svm_linear method
model_svmLinear = train(class ~ ., 
                        data=trainData, 
                        method='svmLinear', 
                        tuneLength = 20, 
                        metric='ROC', 
                        trControl = fitControl
)

#############################
## briefly check the svmLinear results
model_svmLinear 
# model_svmLinear$pred$yes
# model_svmLinear$pred$no


varimp_svmLinear <- varImp(model_svmLinear)
varimp_svmLinear

varImp( model_svmLinear )$importance

plot(varimp_svmLinear, main="Variable Importance with svmLinear")

########################################################################################
## step 4.2
## plot roc for the training data
## We can calculate the area under the curve...
## Select a parameter setting
## selectedIndices <- model_mars2$pred

rocobj_svmlinear <- roc(model_svmLinear$pred$obs, model_svmLinear$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=5, 
                        print.auc=TRUE)


############################

###############################################################
## Step 4.3  Predict on testData and Compute the confusion matrix
## 

predicted2 <- predict(model_svmLinear, testData)

################################################################
## print out confusion matrix 

confusionMatrix(reference = testData$class, data = predicted2, mode='everything', positive='yes')

#######################################
# Confusion Matrix and Statistics
#
#             Reference
# Prediction  No Yes
#         No   7   1
#         Yes  0   2
# #####################################
# Accuracy : 0.9            
# 95% CI : (0.555, 0.9975)
# No Information Rate : 0.7            
# P-Value [Acc > NIR] : 0.1493  
# 
### save the ROC with svmRadiacl plot to local drive
path <- "C:/Users/dug/OneDrive - QIAGEN GmbH/SVM_RFE_Prj/0507CleanedData/Plots/"
path <- paste0(getwd(), "/60vs40Plots/" )
path

PlotROCSaveJPG(rocobj_svmlinear, path, "svmLinear_tune18_train60", "svmLinear" )


################################################################
## print out the whole prediction results
predicted2


######################################################################
#### Section V, step 5
#### Hyper Parameter Tuning using tuneGrid
#### Alternately, you can set the tuneGrid instead of tuneLength.

## Step 5.1: Define the tuneGrid
marsGrid <-  expand.grid(nprune = c(1:10), 
                         degree = c(1:3))

## Step 5.2: Tune hyper parameters by setting tuneGrid
set.seed(123)
model_marsG = train(class ~ ., 
                    data=trainData, 
                    method='earth', 
                    metric='ROC', 
                    tuneGrid = marsGrid, 
                    trControl = fitControl
)

## check tuned MARS model results
model_marsG

## Step 5.2 plot roc for tuned MARS model
rocobj_marsG<- roc(model_marsG$pred$obs, 
                   model_marsG$pred$yes, ci=TRUE,
                   plot=TRUE, 
                   legacy.axes=TRUE, percent=TRUE, 
                   main="MARS_with_Grid",
                   xlab="False Positive Percentage", 
                   ylab="True Postive Percentage", 
                   col="darkblue", lwd=4, 
                   print.auc=TRUE)

## Step 5.3: Predict on testData and Compute the confusion matrix 
predictedG <- predict(model_marsG, testData)
confusionMatrix(reference = testData$class, data = predictedG, mode='everything', positive='yes')

### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_marsG, path, "MARS_Grid_failed", "MARS_Grid" )


######################################################################
#### Section V, Step 6  
#### use random forest method
####
###################################################################### 
## the gbm algo does not work well on this 'small' training set;
## gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 5), 
##                        n.trees = (1:5)*5, 
##                        shrinkage = 0.1,
##                        n.minobsinnode = 15)

##################################################
### Step 6.1 train with RF
set.seed(100)
model_rf = train(class ~ ., 
                 data=trainData, 
                 method='rf', 
                 metric='ROC', 
                 #tuneGrid = gbmGrid, 
                 tuneLength = 1, 
                 trControl = fitControl
)

### Step 6.2 plot ROC 
rocobj_rf<- roc(model_rf$pred$obs, model_rf$pred$yes, 
                ci=TRUE,
                plot=TRUE, 
                legacy.axes=TRUE, percent=TRUE,
                main="Random Forest",
                xlab="False Positive Percentage", 
                ylab="True Postive Percentage", 
                col="darkblue", lwd=4, 
                print.auc=TRUE)



### Step 6.3 Predict on testData and Compute the confusion matrix
predicted_rf <- predict(model_rf, testData)
confusionMatrix(reference = testData$class, data = predicted_rf, mode='everything', positive='yes')

### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_rf, path, "RandomForest_tune1", "Random Forest" )



########################################################################################
#### Section V, step 7 Training SVM_Radial
#### based on the few training methods used, svm_Radial is the optimized one
############################################################

set.seed(100) 

### step 7.1 Train the model using SVMRadial
model_svmRadial = train(class ~ ., 
                        data=trainData, 
                        method='svmRadial', 
                        tuneLength=7,
                        trControl = fitControl
)

## check the prediction
## model_svmRadial$pred

########################################################
## step #7.2 create one roc object, with AUC and 95% CI

rocobj_svmRadial <- roc(model_svmRadial$pred$obs, 
                        model_svmRadial$pred$yes, 
                        ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE, 
                        main="svmRadial",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)

## plot ROC 

### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_svmRadial, path, "svmRadical_tune5", "svmRadical" )


###################################################################
## Step 7.3: Predict on testData and Compute the confusion matrix
predict_svmRadical <- predict(model_svmRadial, testData)
confusionMatrix(reference = testData$class, data = predict_svmRadical, mode='everything', positive='yes') 


#######################################################
#      Confusion Matrix and Statistics
#      
#                   Reference
#      Prediction   no yes
#             no    6   0
#             yes   1   3
#######################################################





#############################################
#### Section V, step 8 Training KNN

set.seed(100)

## step 8.1 Train the model using KNN
model_knn = train(class ~ ., 
                  data=trainData, 
                  method='knn', 
                  tuneLength=4, 
                  trControl = fitControl
)

## check the prediction
# model_knn$pred
model_knn$result

############################################################
## step #8.2 plot roc
rocobj_knn <- roc(model_knn$pred$obs, 
                  model_knn$pred$yes, 
                  ci=TRUE,
                  plot=TRUE, 
                  legacy.axes=TRUE, percent=TRUE, 
                  main="KNN",
                  xlab="False Positive Percentage", 
                  ylab="True Postive Percentage", 
                  col="darkblue", lwd=4, 
                  print.auc=TRUE
)


############################################################
## Step 8.3: Predict on testData and Compute the confusion matrix
predict_knn <- predict(model_knn, testData)
confusionMatrix(reference = testData$class, data = predict_knn, mode='everything', positive='yes') 

### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_knn, path, "KNN_tune4", "KNN" )




#############################################
#### Section V, step 9 Training adaBoost & xgbDART

set.seed(100)

##########################################3
## Unfortunately, adaboost failed 
## Train the model using adaboost
## 
## got an error: Error in { : task 1 failed - "object 'X.HLA.DPB1.' not found"
## hopythesis, the gene name HLA-DPB1 could not be parssed as table header, HLA_DPB1 works;

## change HLA-DPB1 column header into HLA_DPB1:
## trainData$HLA_DPB1 <- trainData$`HLA-DPB1`
## trainData$`HLA-DPB1` <- NULL

## 
model_ada = train(class ~ ., 
                  data=trainData, 
                  method='ada', 
                  tuneLength=2, 
                  trControl = fitControl
) 


############################################################
## step #9.2 plot roc for ada
## this step took about 90 mins

rocobj_ada <- roc(model_ada$pred$obs, 
                  model_ada$pred$yes, 
                  ci=TRUE,
                  plot=TRUE, 
                  legacy.axes=TRUE, percent=TRUE, 
                  main="ADA",
                  xlab="False Positive Percentage", 
                  ylab="True Postive Percentage", 
                  col="darkblue", lwd=4, 
                  print.auc=TRUE
)


############################################################
## Step #9.3: Predict on testData and Compute the confusion matrix
# testData$HLA_DPB1 <- testData$`HLA-DPB1`
# testData$`HLA-DPB1`<- NULL

predict_ada <- predict(model_ada, testData)
confusionMatrix(reference = testData$class, data = predict_ada, mode='everything', positive='yes') 


### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_ada, path, "ada_turne2", "ADA" )




#############################################
## Section V Step 10 train with sbgDART
## 
#############################################
#### Step #10. 1 Train the model using xgbDART 
#### this step might take 60+ mins for a small input data matrix
set.seed(100)

model_xgbDART = train(class ~ ., 
                      data=trainData, 
                      method='xgbDART', 
                      tuneLength=2, 
                      trControl = fitControl, 
                      verbose=F
)

model_xgbDART


############################################################
## step #9.2 plot roc
rocobj_xgbDART <- roc(model_xgbDART$pred$obs, 
                      model_xgbDART$pred$yes, 
                      ci=TRUE,
                      plot=TRUE, 
                      legacy.axes=TRUE, percent=TRUE, 
                      xlab="False Positive Percentage", 
                      ylab="True Postive Percentage", 
                      col="darkblue", lwd=4, 
                      print.auc=TRUE
)


############################################################
## Step #9.3: Predict on testData and Compute the confusion matrix
predict_xgbDART <- predict(model_xgbDART, testData)
confusionMatrix(reference = testData$class, data = predict_xgbDART, mode='everything', positive='yes') 


### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_xgbDART, path, "xgbDART", "xgbDART" )






#####################################################################################
#### Section VI: Compare models
#### call resample() function



# # Compare model performances using resample()

models_compare <- resamples(
  
  list(RF=model_rf, XGBDART=model_xgbDART,   
       SVM=model_svmRadial, SVML = model_svmLinear, knn = model_knn
  )
  
)

model_ada$resample        #150
model_rf$resample         #10
model_xgbDART$resample    #10
model_mars$resample       #25
model_marsG$resample      #150
model_svmLinear$resample  #10
model_svmRadial$resample  #10 
model_knn$resample        #10


# Summary of the models performances
summary(models_compare)

# > summary(models_compare)
# 
#      Call:
#        summary.resamples(object = models_compare)
#      
#      Models: RF, XGBDART, SVM, SVML, knn 
#      Number of resamples: 10 
#      
#      ROC 
#      Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#      RF      0.00 1.00000      1 0.8833333       1    1    0
#      XGBDART 0.75 1.00000      1 0.9750000       1    1    0
#      SVM     0.50 1.00000      1 0.9000000       1    1    0
#      SVML    0.00 0.62500      1 0.8000000       1    1    0
#      knn     0.50 0.90625      1 0.9125000       1    1    0
#      
#      Sens 
#      Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#      RF       0.5   1.000   1.00 0.9000000       1    1    0
#      XGBDART  0.5   0.625   1.00 0.8500000       1    1    0
#      SVM      0.5   0.500   0.75 0.7500000       1    1    0
#      SVML     0.5   0.750   1.00 0.8666667       1    1    0
#      knn      1.0   1.000   1.00 1.0000000       1    1    0
#      
#      Spec 
#      Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#      RF         0   0.000   1.00 0.60           1.00    1    0
#      XGBDART    0   0.625   1.00 0.75           1.00    1    0
#      SVM        0   1.000   1.00 0.85           1.00    1    0
#      SVML       0   0.000   0.75 0.55           1.00    1    0
#      knn        0   0.000   0.00 0.30           0.75    1    0
########################################################################      

## Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)





##################################################################################
####  Section VII
####  ensemble predictions from multiple models using caretEnsemble
#### 

library(caretEnsemble)

# Stacking Algorithms - Run multiple algos in one call.
trainControl <- trainControl(method="repeatedcv", 
                             number=25, 
                             repeats=10,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

algorithmList <- c('rf', 'knn', 'earth', 'xgbDART', 'svmRadial', 'svmLinear')

algorithmList <- c( 'rf', 'svmLinear')


################################################################################
## Run all algorithms in the list: 
set.seed(1234)

models <- caretList(class ~ ., 
                    data=trainData, 
                    trControl=trainControl, 
                    methodList=algorithmList) 


################################################################################
## check resample() results

results <- resamples(models)
summary(results)


#    > summary(results)
#    
#    Call:
#      summary.resamples(object = results)
#    
#    Models: rf, knn, earth, xgbDART, svmRadial, svmLinear 
#    Number of resamples: 30 
#    
#    Accuracy 
#    Min.   1st Qu.    Median      Mean 3rd Qu. Max. NA's
#    rf        0.3333333 0.6666667 0.7500000 0.7822222  1.0000    1    0
#    knn       0.3333333 0.6666667 0.6666667 0.7350000  0.7500    1    0
#    earth     0.0000000 0.6666667 0.6666667 0.7100000  1.0000    1    0
#    xgbDART   0.3333333 0.6666667 0.9000000 0.8211111  1.0000    1    0
#    svmRadial 0.3333333 0.6666667 0.7750000 0.8072222  1.0000    1    0
#    svmLinear 0.2500000 0.6666667 0.6666667 0.7355556  0.9375    1    0
#    
#    Kappa 
#              Min. 1st Qu.    Median      Mean 3rd Qu. Max. NA's
#    rf        -0.5     0.0 0.5000000 0.4181818   1.000    1    0
#    knn       -0.5     0.0 0.0000000 0.2681818   0.500    1    0
#    earth     -0.8     0.0 0.2000000 0.3515152   1.000    1    0
#    xgbDART   -0.5     0.0 0.7727273 0.5448485   1.000    1    0
#    svmRadial  0.0     0.4 0.5227273 0.6115152   1.000    1    0
#    svmLinear -0.5     0.0 0.2000000 0.3245455   0.875    1    0


################################################################################
## Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

## Save as multi_algo_Accuracy_Kappa_boxplot



################################################################################
## Plot multi ROCs in one plot
rocobj_models <- roc(models$rf$pred$obs, 
                     models$rf$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(models$svmRadial$pred$obs, 
                     models$svmRadial$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="green", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 44,
                     add = TRUE
)


rocobj_models <- roc(models$svmLinear$pred$obs, 
                     models$svmLinear$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 48,
                     add = TRUE
)

rocobj_models <- roc(models$xgbDART$pred$obs, 
                     models$xgbDART$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="black", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 52,
                     add = TRUE
)

rocobj_models <- roc(models$earth$pred$obs, 
                     models$earth$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="yellow", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 56,
                     add = TRUE
)


rocobj_models <- roc(models$knn$pred$obs, 
                     models$knn$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="pink", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 35.1,
                     add = TRUE
)


legend("bottomright", 
       legend=c( "rf", "svmRadial", "svmLinear", "xgbDART", "MARS", "knn" ), 
       col=c( "darkblue", "green", "red", "black", "yellow", "pink" ), 
       lwd=4
)


########################################################################
## Run all algorithms in the list: against validation dataset

set.seed(123)

valid.models <- caretList(class ~ ., 
                          data=testData, 
                          trControl=trainControl, 
                          #trControl = fitControl,
                          methodList=algorithmList) 



###################3
## Plot two ROC curves together, 
## compare ROC train and ROC validation, although this is a pretty odd comparison
## 


rocobj_models_svmLinear <- roc(models$svmLinear$pred$obs, 
                               models$svmLinear$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               main = "svm-Linear ROC-AUC Expression",
                               xlab=" 1 - Specificity", 
                               ylab="Sensitivity", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 35.1
                               #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$svmLinear$pred$obs, 
                                 valid.models$svmLinear$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 main = "svm-Linear ROC-AUC Expression",
                                 xlab="1 - Specifiity", 
                                 ylab="Sensitivity", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "svm Linear Training", "svm Linear Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)

## the variable importance would change a bit if we use different seeds; 
## for example, KLRC1 and IRF4 contricute almost equally to the model, 
## so the rank of those two genes might switch if we use another seed. 
plot( varImp(models$svmLinear) ) 


frame()

rocobj_models_svmRadial <- roc(models$svmRadial$pred$obs, 
                               models$svmRadial$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 35.1
                               #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$svmRadial$pred$obs, 
                                 valid.models$svmRadial$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "svm Radial Training", "svm Radial Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)





frame()

rocobj_models_svmLinear <- roc(models$rf$pred$obs, 
                               models$rf$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 35.1
                               #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$rf$pred$obs, 
                                 valid.models$rf$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "Random Forest Training", "Random Forest Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)





### KNN
frame()

rocobj_models_ <- roc(models$knn$pred$obs, 
                      models$knn$pred$yes, 
                      ci=TRUE,
                      plot=TRUE, 
                      legacy.axes=TRUE, percent=TRUE, 
                      xlab="False Positive Percentage", 
                      ylab="True Postive Percentage", 
                      col="darkblue", lwd=4, 
                      print.auc=TRUE,
                      print.auc.y = 35.1
                      #add = TRUE
)

rocobj_models_knn <- roc(valid.models$knn$pred$obs, 
                         valid.models$knn$pred$yes, 
                         ci=TRUE,
                         plot=TRUE, 
                         legacy.axes=TRUE, percent=TRUE, 
                         xlab="False Positive Percentage", 
                         ylab="True Postive Percentage", 
                         col="red", lwd=4, 
                         print.auc=TRUE,
                         print.auc.y = 25.1,
                         add = TRUE
)



legend("bottomright", 
       legend=c( "KNN Training", "KNN Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)



### xgbDART
frame()

rocobj_models_svmLinear <- roc(models$xgbDART$pred$obs, 
                               models$xgbDART$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 35.1
                               #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$xgbDART$pred$obs, 
                                 valid.models$xgbDART$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "xgbDART Training", "xgbDART Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)


## MARS
frame()

rocobj_models_mars <- roc(models$earth$pred$obs, 
                          models$earth$pred$yes, 
                          ci=TRUE,
                          plot=TRUE, 
                          legacy.axes=TRUE, percent=TRUE, 
                          xlab="False Positive Percentage", 
                          ylab="True Postive Percentage", 
                          col="darkblue", lwd=4, 
                          print.auc=TRUE,
                          print.auc.y = 35.1
                          #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$earth$pred$obs, 
                                 valid.models$earth$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "MARS Training", "MARS Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)



################################################################################
####  Section VIII 
####  combine the predictions of multiple models to form a final prediction


# Create the trainControl
set.seed(100)

stackControl <- trainControl(method="repeatedcv", 
                             number=20, 
                             repeats=5,
                             savePredictions=TRUE, 
                             classProbs=TRUE
)

# Ensemble the predictions of `models` to form a new combined prediction based on glm
stack.glm <- caretStack(models, method="glm", metric="Accuracy", trControl=stackControl)

print(stack.glm)

# Predict on testData 
stack_predicteds <- predict(stack.glm, newdata=testData) 
head(stack_predicteds) 

confusionMatrix(reference = testData$class, data = stack_predicteds, mode='everything', positive='Yes') 


########################################################################3
####
####         END OF THE MAIN CODE SECTION                           ####3
####
########################################################################3




