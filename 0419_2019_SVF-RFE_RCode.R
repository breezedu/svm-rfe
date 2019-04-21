## 
#################
## use caret package to plot ROC for training set and 10-fold validation set
## 

# install.packages("mlbench")
# install.packages("caret")
# install.packages("ellipse")   ## for feature plot; 
library(mlbench)
library(caret)
library(ggplot2)
library(ellipse)

set.seed(12345)

# Prepare data ------------------------------------------------------------



####################################################
## data of RNAseq for each tumor type from TCGA
## count <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/SVM_RFE_Prj/svm_rfe_prj/THCA_RNA_Set.RnaSeq_Transcript.Genes.txt",
##                    header = TRUE,
##                    row.names = 1,
##                    sep = "\t"
##                    )
##
## design <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/SVM_RFE_Prj/svm_rfe_prj/THCA_RNA_Set.RnaSeq_Transcript.Genes_Design.txt",
##                     header = TRUE,
##                     row.names = 1,
##                     sep = "\t"
##                     )
####################################################


########################################################################################
## 
## mRNA data
## 
########################################################################################
## read in the count data and the clinical information table; 
## 
count <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/test_miRNA/miRNA_OmicData.txt", row.names = 1, header = T, sep = "\t")

design <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/test_miRNA/miRNA_OmicData_Design.txt", row.names = 1, header = T, sep = "\t")

top.genes <- readLines("D:/WorkRecord/Companies/Qiagen_Sales/201904/test_miRNA/significant_genes_DESeq2.txt")

top.genes[1:10] 
length(top.genes) 

############################################################
## remove rows with constant values (zeros)
## for mRNA data, remove rows with varience less than 0.01
## 

datamatrix <- count

dim(datamatrix)

constantRow <- row.names(as.matrix(which(apply(datamatrix, MARGIN = 1, function(x) var(x) < 0.01) ) ) ) 

#exclude rows with constant value: 
count <- datamatrix[!row.names(datamatrix) %in% constantRow,] 

count <- datamatrix[ row.names(datamatrix) %in% top.genes, ]

print("These variables get constant values, thus got droped during PCA scale:") 
constantRow


my_data <- as.data.frame( t(count) )


dim(my_data)


## for RnaSeq my_data$class <- design$sample.type

## for mRNA data: 
my_data$class <- design$Progression.Logit

summary(my_data$class)

## subset the data, only keep primiary tumor and solid normal
my_data <- subset( my_data, class!='Metastatic') # my_data$class != "Metastatic", ]

## check factor levels left
summary(my_data$class) 

## drop the empty Metastatic level
my_data$class <- factor(my_data$class)








# Create the training and test datasets
set.seed(123)

# Step 1: Get row numbers for the training data
trainRowNumbers <- createDataPartition(my_data$class, p=0.8, list=FALSE)

# Step 2: Create the training  dataset
trainData <- my_data[trainRowNumbers,]

# Step 3: Create the test dataset
testData <- my_data[-trainRowNumbers,]

dim(trainData)
dim(testData)

# Store X and Y for later use.
x = trainData[, 1:42]
y = trainData$class

x_test = testData[ , 1:42]
y_test = testData$class
dim(x)

x$class

########################################################################
## convert all the numeric variables to range between 0 and 1, 
## by setting method=range in preProcess(). 
## preProcess_range_model <- preProcess(trainData, method='range')
## trainData <- predict(preProcess_range_model, newdata = trainData)
## 
## Append the Y variable
## trainData$class <- y
########################################################################


########################################################################
## remove constant columns
## for top 500/1000 significant genes, this step could be ignored; 
x <- x[,apply(x, 2, var, na.rm=TRUE) > 0.0001]
dim(x)


########################################################################
## convert all the numeric variables to range between 0 and 1, by setting method=range in preProcess(). 
trainData_range_model <- preProcess( x, method = c("center", "scale") )

trainData <- predict( trainData_range_model, newdata = x) 
trainData$class <- y

head(trainData)

########################################################################
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
## plot
## 
## box plot, to see significantly differential expression genes 
featurePlot(x = trainData[, 1:20], 
            y = trainData$class, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free"))
            )


## density plot, to visulize more important variables 
featurePlot(x = trainData[, 1:20], 
            y = trainData$class, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))



##########################################################################################################  
#  RFE works in 3 broad steps:
#    
#  Step 1: Build a ML model on a training dataset and estimate the feature importances on the test dataset.
#  
#  Step 2: Keeping priority to the most important variables, 
#   iterate through by building models of given subset sizes, 
#   that is, subgroups of most important predictors determined from step 1. 
#   Ranking of the predictors is recalculated in each iteration.
#  
#  Step 3: The model performances are compared across different subset sizes 
#   to arrive at the optimal number and list of final predictors.
#  
#  It can be implemented using the rfe() function and you have the flexibility 
#   to control what algorithm rfe uses and how it cross validates by defining the rfeControl().
##########################################################################################################  


set.seed(1234)

options(warn=-1)

subsets <- c(1:34)

dim(trainData)


## build a control model, on training dataset; estimate the feautre importances 
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",      ## cross-validation
                   repeats = 15,               ## fold of cv
                   verbose = FALSE
                   )


## Recursive Feature Elimination rfe() 
lmProfile <- rfe(x=trainData[, 1:42], y=trainData$class,
                 sizes = subsets,        ## model sizes (the number of most important genes to choose)
                 rfeControl = ctrl       ## use rfeControl output as reference: algorithm and cv to use;
                 )

## check lm_profile results
lmProfile 

## end by 4/15/2019 evening;
##################################################### 
 
 
 
 
#### 
#### train() the model and interpret the results
# Set the seed for reproducibility
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





fitted <- predict(model_mars)

## 
plot(model_mars, main="Model Accuracies with MARS") 
 
### 

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

 
## Define the training control
fitControl <- trainControl(
  method = 'repeatedcv',                   # k-fold cross validation
  number = 10,                     # number of folds
  repeats = 10,                     # number of repeats
  savePredictions = T,             # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary  #,  # results summary function
  # savePredictions = 'final'
) 





#################################################################################################
# Step 1: Tune hyper parameters by setting tuneLength
set.seed(123)

model_svmLinear = train(class ~ ., 
                    data=trainData, 
                    method='svmLinear', 
                    tuneLength = 9, 
                    metric='ROC', 
                    trControl = fitControl
                    )
model_svmLinear 

model_svmLinear$pred$yes

########################################################################################
## plot roc for the training data
## We can calculate the area under the curve...

library(pROC)

# Select a parameter setting
model_svmLinear
model_svmLinear$results

# selectedIndices <- model_mars2$pred


# Plot:
plot( 
  
    smooth( roc(model_svmLinear$pred$obs, model_svmLinear$pred$yes)
          )
    )


rocobj <- roc(model_svmLinear$pred$obs, model_svmLinear$pred$yes, ci=TRUE,
                      plot=TRUE, 
                      legacy.axes=TRUE, percent=TRUE, 
                      xlab="False Positive Percentage", 
                      ylab="True Postive Percentage", 
                      col="darkblue", lwd=4, 
                      print.auc=TRUE)
          

plot(
  
  smooth( rocobj ),
  #ci=TRUE,
  legacy.axes=TRUE, percent=TRUE, 
  xlab="False Positive Percentage", 
  ylab="True Postive Percentage", 
  col="darkblue", lwd=4, 
  print.auc=TRUE
  
)

#ci.95 <- paste0("95% CI = ", round(ci.auc(rocobj)[1], digits = 2), " - ", round(ci.auc(rocobj)[3], digits = 2) ) 

legend("topleft", 
      legend=c( paste("mRNA", "svmLinear optimzed") ), 
      col=c("darkblue", "darkblue"), 
      lwd=4)

# auc <- ci(rocobj)
# got this from previous plot
auc <- "95% CI: 60.4%-93.9%"

legend("bottomright", 
       legend=c( paste("AUC", auc) ), 
       col=c("darkblue", "darkblue"), 
       lwd=4)



# Step 2: Predict on testData and Compute the confusion matrix
predicted2 <- predict(model_svmLinear, testData)
dim(testData)
head(predicted2)
## 
confusionMatrix(reference = testData$class, data = predicted2, mode='everything', positive='yes')


predicted2

######################################################################
## Hyper Parameter Tuning using tuneGrid
## Alternately, you can set the tuneGrid instead of tuneLength.

# Step 1: Define the tuneGrid
marsGrid <-  expand.grid(nprune = c(2, 4, 6, 8, 10), 
                         degree = c(1, 2, 3))

# Step 2: Tune hyper parameters by setting tuneGrid
set.seed(123)
model_mars3 = train(class ~ ., 
                    data=trainData, 
                    method='earth', 
                    metric='ROC', 
                    tuneGrid = marsGrid, 
                    trControl = fitControl
                    )

model_mars3

rocobj_mars<- roc(model_mars3$pred$obs, model_mars3$pred$yes, ci=TRUE,
              plot=TRUE, 
              legacy.axes=TRUE, percent=TRUE, 
              xlab="False Positive Percentage", 
              ylab="True Postive Percentage", 
              col="darkblue", lwd=4, 
              print.auc=TRUE)

# Step 3: Predict on testData and Compute the confusion matrix
predicted3 <- predict(model_mars3, testData)
confusionMatrix(reference = testData$class, data = predicted3, mode='everything', positive='yes')



############################################################
## random forest method

######## 

# gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 5), 
#                        n.trees = (1:5)*5, 
#                        shrinkage = 0.1,
#                        n.minobsinnode = 15)

model_rf = train(class ~ ., 
                 data=trainData, 
                 method='rf', 
                 metric='ROC', 
                 #tuneGrid = gbmGrid, 
                 tuneLength = 5, 
                 trControl = fitControl
)


rocobj_rf<- roc(model_rf$pred$obs, model_rf$pred$yes, ci=TRUE,
                  plot=TRUE, 
                  legacy.axes=TRUE, percent=TRUE, 
                  xlab="False Positive Percentage", 
                  ylab="True Postive Percentage", 
                  col="darkblue", lwd=4, 
                  print.auc=TRUE)



# Step rf: Predict on testData and Compute the confusion matrix
predicted_rf <- predict(model_rf, testData)
confusionMatrix(reference = testData$class, data = predicted_rf, mode='everything', positive='yes')


#############################################
#### Training SVM
set.seed(125)

# Train the model using SVM
model_svmRadial = train(class ~ ., 
                        data=trainData, 
                        method='svmRadial', 
                        tuneLength=5, 
                        trControl = fitControl
                        )

model_svmRadial$pred

## plot roc
rocobj_svm <- roc(model_svmRadial$pred$obs, model_svmRadial$pred$yes, ci=TRUE,
               plot=TRUE, 
               legacy.axes=TRUE, percent=TRUE, 
               xlab="False Positive Percentage", 
               ylab="True Postive Percentage", 
               col="darkblue", lwd=4, 
               print.auc=TRUE)

plot(rocobj_svm,
     legacy.axes=TRUE, percent=TRUE, 
     xlab="False Positive Percentage", 
     ylab="True Postive Percentage", 
     col="darkblue", lwd=4, 
     print.auc=TRUE
     
     )

path <- "D:/WorkRecord/Companies/Qiagen_Sales/201904/SVM_RFE_Prj/"


### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_svm, path, "svmRadical_tune5", "svmRadical" )

# Step 3: Predict on testData and Compute the confusion matrix
predict_svm <- predict(model_svmRadial, testData)
confusionMatrix(reference = testData$class, data = predict_svm, mode='everything', positive='yes') 


#######################################################
#      Confusion Matrix and Statistics
#      
#      Reference
#      Prediction no yes
#      no   8   1
#      yes  0   3
#######################################################





























#############################################################################
### extra test and colde templates
### 
########################################################
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
}








# Cross Validation Definition ---------------------------------------------------

fitControl <-
  trainControl(
    method = "cv",
    number = 10,
    classProbs = T,
    savePredictions = T,
    summaryFunction = twoClassSummary
  )


# Training with Supporting Vector Machine ----------------------------------------------------------------


model <- train(
  Class ~ .,
  data = my_data,
  method = "rf",
  trControl = fitControl,
  metric = "ROC"
)


for_lift <- data.frame(Class = model$pred$obs, 
                       rf = model$pred$M, 
                       resample = model$pred$Resample
                       )

lift_df <-  data.frame()

for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(Class ~ rf, data = fold_df, class = "R")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}

lift_obj <- lift(Class ~ rf, data = for_lift, class = "R")


# Plot ROC ----------------------------------------------------------------

ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))

##################################################################

library('pROC')

model$results
model$pred$M
model$method

selectedIndices <- model$pred$mtry == 2

roc_train <- roc(model$pred$obs[selectedIndices],
                 model$pred$M[selectedIndices],
                 print.auc = T
                 )

# plot.roc( roc_train,  )


## We can calculate the area under the curve...
rocobj <- roc(model$pred$obs[selectedIndices], model$pred$M[selectedIndices], ci=TRUE,
              plot=TRUE, 
              legacy.axes=TRUE, percent=TRUE, 
              #  xlab="False Positive Percentage", 
              #  ylab="True Postive Percentage", 
              col="darkblue", lwd=4, 
              print.auc=TRUE
              )


