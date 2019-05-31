

##############################
## clinical model
## Jeff Du 05/13/2019
## 
#############################

####################################################
## 
## Section I
## Initial Setup, load package and data
## 
####################################################

library(caret)

getwd()

#########
setwd("D:/miRNA_Prj/0520_plot2RocAuc/")

clinical.data <- read.table("Clinical_Data.txt", 
                            header = T, row.names=1, sep="\t")


str(clinical.data)

dim(clinical.data)
clinical.data[1:5, 1:5]
clinical.data$Time2FollowUp <- NULL
clinical.data$Time2Progression <- NULL
clinical.data$AliveOrDead <- NULL

# clinical.data$locoregional <- NULL 

clinical.data$Progressed

summary(clinical.data$AliveOrDead)



####################################################
## 
## Sectopm II
## split the dataset into training and validation 
## 
#################################################### 

set.seed(123)

trainRowNum <- createDataPartition(clinical.data$Progressed, p=0.7, list = F)
testRowNum <- createDataPartition(clinical.data$Progressed, p=0.6, list = F)

train.data <- clinical.data[trainRowNum, ]
test.data <- clinical.data[-trainRowNum, ]
test.data <- clinical.data[testRowNum, ]
## save X and Y for later combination;
y <- train.data$Progressed
x <- train.data
x$Progressed <- NULL

dim(x)
dim(train.data)
length(y)

x_t <- test.data 
y_t <- test.data$Progressed
x_t$Progressed <- NULL 

##############################################
## 
## naive visulization of features in the training data
##
##############################################
dim(train.data)
dim(test.data)

library(skimr)
skim.train <- skim_to_wide(train.data)
skim.train[1:13, c(1:13)]

skim.test <- skim_to_wide(test.data)
skim.test[1:13, c(1:13)]
dim(test.data)
## good, there's no missing data



##############################################
## 
## create One-Hot Encoding variables 
## for the categories in the training data
##
##############################################

##
## One-Hot Encoding to create dummy variables by converting a categorical variable to as many binary variables
# dummie.model <- dummyVars(Progressed ~ ., data=train.data)
# 
# train.dataDum <- predict(dummie.model, newdata = train.data)
train.data <- data.frame(train.data)
test.data <- data.frame(test.data)

## cound not use train.data <- as.data.frame(train.dataDum)
## test.dataDum <- predict(dummie.model, newdata = test.data)

################################
## remove constant columns

#datamatrix <- train.data
#dim(datamatrix)

#datamatrix <- datamatrix[ , apply(datamatrix, 2, var, na.rm=T)!=0]
#dim(datamatrix)

#constantCol <- colnames(as.matrix( which(apply(datamatrix, MARGIN=2, function(x) var(x) < 0.001))))
#count <- datamatrix[ !colnames(datamatrix) %in% constantCol, ]
#dim(count)
#head(count)

#train.data <- datamatrix

## check the new data structure of the train.data;
str(train.data)
dim(train.data)

#str(test.dataDum)
#dim(test.dataDum)
#test.dataDum <- as.data.frame(test.dataDum)

preProcess.model <- preProcess( train.data, method = 'range' )

train.data <- predict(preProcess.model, newdata = train.data)
test.data <- predict(preProcess.model, newdata = test.data)

# preProcess.modelT <- preProcess( test.dataDum, method = 'range')
# test.dataPro <- predict(preProcess.modelT, newdata = test.dataDum)

## ad the progressed variable to the train.data
## train.data$Progress <- y

Y <- ifelse(y==1, "yes", "no")
Y
Y_T <- ifelse(y_t==1, "yes", "no")

train.data$class <- Y

test.data$class <- Y_T

dim(test.data)
dim(train.data)

num_col <- dim(train.data)[2]
num_col

t.col <- dim(test.data)[2]

## ZScore ? 
apply( train.data[, 1:num_col], 2, FUN=function(x){c('min' = min(x), 'max' = max(x) )})
apply(  test.data[,   1:t.col], 2, FUN=function(x){c('min' = min(x), 'max' = max(x) )})

head(train.data[, num_col-5:num_col])

train.data$class <- as.factor(train.data$class)
test.data$class <- as.factor((test.data$class))

featurePlot( x = train.data[ , 1:num_col-1],
             y = train.data$class,
             plot = 'box',
             strip = strip.custom(par.strip.text = list( cex = 0.7)),
             scales = list(x = list(relation="free"), 
                           y = list(relation="free")
             )
)


featurePlot(x = train.data[, 1:num_col-1], 
            y = train.data$class, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free"))
)



#######################################################################


set.seed(100)
options(warn=-1)

subsets <- c(1:num_col)

str(train.data)
train.data$class <- as.factor(train.data$class)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

lmProfile <- rfe(x=train.data[, 1:num_col], y=train.data$class,
                 sizes = subsets,
                 rfeControl = ctrl)




lmProfile 

varimp_lm <- varImp(lmProfile)

plot(varimp_lm, main="variable importance with lmProfile")

# colnames(train.data) 


set.seed(123)

fitControl <- trainControl(
  method = 'repeatedcv',           # k-fold cross validation
  number = 15,                     # number of folds
  repeats = 10,                    # number of repeats
  savePredictions = T,             # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary  #,  # results summary function
  # savePredictions = 'final'
) 

summary(train.data$Progressed)

train.data$Progressed <- NULL

set.seed(123)
#############################
#### Step 4.1 choose svm_linear method


clin.model_svmLinear = train(class ~ ., 
                        data=train.data, 
                        method='svmLinear', 
                        tuneLength = 5, 
                        metric='ROC', 
                        trControl = fitControl
                        )

#############################
## briefly check the svmLinear results
clin.model_svmLinear 

# model_svmLinear$pred$yes
# model_svmLinear$pred$no


varimp_svmLinear <- varImp(clin.model_svmLinear)
plot(varimp_svmLinear, main="Clinical Variable Importance with svmLinear")


########################################################################################
## step 4.2
## plot roc for the training data
## We can calculate the area under the curve...
## Select a parameter setting
## selectedIndices <- model_mars2$pred
library(pROC)



dim(train.data)
dim(test.data)

test.data$Progressed <- NULL
#test.data$Time2Progression <- NULL

### build testing model on validation set
clin.model_svmLinear.valid = train(class ~ ., 
                             data=test.data, 
                             method='svmLinear', 
                             tuneLength = 5, 
                             metric='ROC', 
                             trControl = fitControl
)


### plot two ROC-AUC in on figure

# roc from training set
rocobj_svmlinear <- roc(clin.model_svmLinear$pred$obs, 
                        clin.model_svmLinear$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear with Cinical Data",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)

# roc from validation set
rocobj_svmlinear.valid <- roc(clin.model_svmLinear.valid$pred$obs, 
                        clin.model_svmLinear.valid$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear with Cinical Data",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="red", lwd=4, 
                        print.auc=TRUE,
                        add = T)




################################################################################
## Plot multi ROCs in one plot
rocobj_clinical.svmLinear <- roc(clin.model_svmLinear$pred$obs, 
                               clin.model_svmLinear$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               main="svmLinear Clinical ROC",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 40
)

rocobj_clinical.svmLinear.valid <- roc(clin.model_svmLinear.valid$pred$obs, 
                                       clin.model_svmLinear.valid$pred$yes, 
                                     ci=TRUE,
                                     plot=TRUE, 
                                     legacy.axes=TRUE, percent=TRUE, 
                                     xlab="False Positive Percentage", 
                                     ylab="True Postive Percentage", 
                                     col="red", lwd=4, 
                                     print.auc=TRUE,
                                     print.auc.y = 30,
                                     add = TRUE
)




legend("bottomright", 
       legend=c( "svmLinear Training Set", "svmLinear Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)




#############################
#### Step 5.1 choose svm_Radial method
set.seed(100)
clin.model_svmRadial = train(class ~ ., 
                        data=train.data, 
                        method='svmRadial', 
                        #tuneLength = 4, 
                        metric='ROC', 
                        trControl = fitControl
                        )

set.seed(100)
clin.model_svmRadial.valid = train(class ~ ., 
                             data=test.data, 
                             method='svmRadial', 
                             # tuneLength = 4, 
                             metric='ROC', 
                             trControl = fitControl
                          )


#############################
## briefly check the svmLinear results
clin.model_svmRadial



varimp_svmLinear <- varImp(clin.model_svmRadial)
plot(varimp_svmLinear, main="Clinical Variable Importance with svmRadial")


########################################################################################
## step 4.2
## plot roc for the training data
## We can calculate the area under the curve...
## Select a parameter setting
## selectedIndices <- model_mars2$pred
# library(pROC)





################################################################################
## Plot multi ROCs in one plot
rocobj_clinical.svmRadial <- roc(clin.model_svmRadial$pred$obs, 
                                 clin.model_svmRadial$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 main="svmRadial Clinical ROC",
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="darkblue", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 40
)

rocobj_clinical.svmRadial.valid <- roc(clin.model_svmRadial.valid$pred$obs, 
                                       clin.model_svmRadial.valid$pred$yes, 
                                       ci=TRUE,
                                       plot=TRUE, 
                                       legacy.axes=TRUE, percent=TRUE, 
                                       xlab="False Positive Percentage", 
                                       ylab="True Postive Percentage", 
                                       col="red", lwd=4, 
                                       print.auc=TRUE,
                                       print.auc.y = 30,
                                       add = TRUE
)




legend("bottomright", 
       legend=c( "svmRadial Training Set", "svmRadial Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)




#########################################################################################
### clinical heatmap 
library(corrplot)

dim(train.data)
corrplot( cor(train.data[, 1:12]), method = "square", tl.cex = 0.8)

corrplot( cor(clinical.data), method = "square", tl.cex = 0.8)



predicted2 <- predict(clin.model_svmLinear, test.data) 

head(predicted2)

confusionMatrix(reference = test.data$class, data = predicted2, mode='everything', positive='yes')












#############################
#### Step 4.1 choose svm_linear method
set.seed(100)

clin.model_svmRadial = train(class ~ ., 
                        data=train.data, 
                        method='svmRadial', 
                        tuneLength = 8, 
                        metric='ROC', 
                        trControl = fitControl
)

#############################
## briefly check the svmLinear results
clin.model_svmRadial 
# model_svmLinear$pred$yes
# model_svmLinear$pred$no


varimp_svmRadial <- varImp(clin.model_svmRadial)
plot(varimp_svmRadial, main="Variable Importance with svmLinear")

rocobj_svmRadial <- roc(clin.model_svmRadial$pred$obs, 
                        clin.model_svmRadial$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmRadial with Cinical Data",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)

test.data3$Time2Progression <- NULL
model_svmRadial.valid = train(class ~ ., 
                        data=test.data, 
                        method='svmRadial', 
                        tuneLength = 8, 
                        metric='ROC', 
                        trControl = fitControl
)

#rocobj_svmRadial.valid
#test.data3$Time2Progression <- NULL

rocobj_svmRadial <- roc(model_svmRadial.valid$pred$obs, 
                        model_svmRadial.valid$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmRadial with Cinical Data",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)



predicted_svmRadial <- predict(model_svmRadial.valid, test.data) 

head(predicted_svmRadial)

confusionMatrix(reference = test.data$class, data = predicted_svmRadial, mode='everything', positive='yes')







##################################################################################
####  Section VII
####  ensemble predictions from multiple models using caretEnsemble
#### 

library(caretEnsemble)

fitControl <- trainControl(
  method = 'repeatedcv',           # k-fold cross validation
  number = 15,                     # number of folds
  repeats = 10,                    # number of repeats
  savePredictions = T,             # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary  #,  # results summary function
  # savePredictions = 'final'
) 



# Stacking Algorithms - Run multiple algos in one call.
trainControl <- trainControl(method="repeatedcv", 
                             number=25, 
                             repeats=20,
                             savePredictions=TRUE, 
                             classProbs=TRUE)


algorithmList <- c('rf', 'knn', 'earth', 'svmRadial', 'svmLinear')
algorithmList <- c('rf', 'knn', 'earth', 'xgbDART', 'svmRadial', 'svmLinear')


################################################################################
## Run all algorithms in the list: 
set.seed(100)


dim(train.data)
clin.models <- caretList( class ~ ., 
                          data=train.data, 
                          metric='ROC', 
                          trControl=fitControl, 
                          methodList=algorithmList
                        ) 




dim(test.data)
dim(train.data)

clin.models.valid <- caretList(class ~ ., 
                               data=test.data, 
                               metric='ROC',
                               trControl=fitControl, 
                               methodList=algorithmList
                               ) 


#################################

## Plot combination ROCs




################################################################################
## check resample() results

# models_compare <- resamples(
#  list(clin.models$rf, clin.models$knn, clin.models$svmRadial, clin.models$svmLinear)
#   )
# summary(models_compare)

clin.results <- resamples(clin.models)
summary(clin.results)


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
bwplot(clin.results, scales=scales)

## Save as multi_algo_Accuracy_Kappa_boxplot



################################################################################
## Plot TWO ROCs in one plot
rocobj_models <- roc(clin.models$rf$pred$obs, 
                     clin.models$rf$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main="Random Forest Clinical ROC",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(clin.models.valid$rf$pred$obs, 
                     clin.models.valid$rf$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 30,
                     add = TRUE
)




legend("bottomright", 
       legend=c( "Random Forest Training", "Random Forest Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)


### 


rocobj_models <- roc(clin.models$svmRadial$pred$obs, 
                     clin.models$svmRadial$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="green", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 44
                    )


################################################################################
## Plot multi ROCs in one plot
rocobj_models <- roc(clin.models$svmRadial$pred$obs, 
                     clin.models$svmRadial$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main="svmRadial Clinical ROC",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(clin.models.valid$svmRadial$pred$obs, 
                     clin.models.valid$svmRadial$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 30,
                     add = TRUE
)




legend("bottomright", 
       legend=c( "svmRadial Training Set", "svmRadial Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)


#######################


rocobj_models <- roc(clin.models$svmLinear$pred$obs, 
                     clin.models$svmLinear$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 48
                    )


## Plot multi ROCs in one plot
rocobj_models <- roc(clin.models$svmLinear$pred$obs, 
                     clin.models$svmLinear$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main="svmLinear Clinical ROC",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(clin.models.valid$svmLinear$pred$obs, 
                     clin.models.valid$svmLinear$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 30,
                     add = TRUE
)




legend("bottomright", 
       legend=c( "svmLinear Training Set", "svmLinear Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)

#######################################################


rocobj_models <- roc(clin.models$xgbDART$pred$obs, 
                     clin.models$xgbDART$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main="xgbDART Clinical ROC",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="black", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 52,
                     add = TRUE
)


## Plot multi ROCs in one plot
rocobj_models <- roc(clin.models$xgbDART$pred$obs, 
                     clin.models$xgbDART$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main="xgbDART Clinical ROC",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(clin.models.valid$xgbDART$pred$obs, 
                     clin.models.valid$xgbDART$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 30,
                     add = TRUE
)




legend("bottomright", 
       legend=c( "xgbDART Training Set", "xgbDART Validation Set" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)




####################################################
rocobj_models <- roc(clin.models$earth$pred$obs, 
                     clin.models$earth$pred$yes, 
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

## Plot multi ROCs in one plot
rocobj_models <- roc(clin.models$earth$pred$obs, 
                     clin.models$earth$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main="MARS Clinical ROC",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(clin.models.valid$earth$pred$obs, 
                     clin.models.valid$earth$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 30,
                     add = TRUE
)




legend("bottomright", 
       legend=c( "MARS Training Set", "MARS Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)



##############################################################
rocobj_models <- roc(clin.models$knn$pred$obs, 
                     clin.models$knn$pred$yes, 
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

###
## Plot multi ROCs in one plot
rocobj_models <- roc(clin.models$knn$pred$obs, 
                     clin.models$knn$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main="KNN Clinical ROC",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(clin.models.valid$knn$pred$obs, 
                     clin.models.valid$knn$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 30,
                     add = TRUE
)




legend("bottomright", 
       legend=c( "KNN Training Set", "KNN Validation Set" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)



############################################################################################
legend("bottomright", 
       legend=c( "RandomForest", "svmRadial", "svmLinear", "xgbDART", "MARS", "knn" ), 
       col=c( "darkblue", "green", "red", "black", "yellow", "pink" ), 
       lwd=4
)


### END
