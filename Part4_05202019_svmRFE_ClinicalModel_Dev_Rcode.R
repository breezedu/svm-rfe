

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
setwd("D:/miRNA_Prj/0507CleanedData/")

clinical.data <- read.table("Count_Design.txt", 
                            header = T, row.names=1, sep="\t")


str(clinical.data)

dim(clinical.data)
clinical.data[1:5, 1:5]

clinical.data$locoregional

clinical.data$Progressed

summary(clinical.data$TNM_Stage)

## remove results categories from the clinical table

clinical.data$progtype1_systemic_2locoregional <- NULL
clinical.data$Did_the_patient_develop_a_New_Tumor_Event <- NULL
clinical.data$Alive_or_Dead <- NULL
clinical.data$locoregional <- NULL
clinical.data$TTP <- NULL
clinical.data$New_Tumor_Event <- NULL
clinical.data$Time2FollowUp <- NULL
clinical.data$Alive_Dead <- NULL
clinical.data$TNM_Stage <- NULL
clinical.data$Record_ID <- NULL

head(clinical.data[ ,1:8])




####################################################
## 
## Sectopm II
## split the dataset into training and validation 
## 
#################################################### 

set.seed(123)

trainRowNum <- createDataPartition(clinical.data$Progressed, p=0.6, list = F)

train.data <- clinical.data[trainRowNum, ]
test.data <- clinical.data[-trainRowNum, ]

## save X and Y for later combination;
y <- train.data$Progressed
x <- train.data
x$Progressed <- NULL

dim(x)
dim(train.data)
length(y)

y_t <- test.data$Progressed

##############################################
## 
## naive visulization of features in the training data
##
##############################################

library(skimr)
skim.train <- skim_to_wide(train.data)
skim.train[1:20, c(1:5, 15:16)]

## good, there's no missing data



##############################################
## 
## create One-Hot Encoding variables 
## for the categories in the training data
##
##############################################

##
## One-Hot Encoding to create dummy variables by converting a categorical variable to as many binary variables
dummie.model <- dummyVars(Progressed ~ ., data=train.data)

train.dataDum <- predict(dummie.model, newdata = train.data)
train.data <- data.frame(train.dataDum)

## cound not use train.data <- as.data.frame(train.dataDum)

test.dataDum <- predict(dummie.model, newdata = test.data)

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

str(test.dataDum)
dim(test.dataDum)
test.dataDum <- as.data.frame(test.dataDum)

preProcess.model <- preProcess( train.data, method = 'range' )

train.data <- predict(preProcess.model, newdata = train.data)

# preProcess.modelT <- preProcess( test.dataDum, method = 'range')
# test.dataPro <- predict(preProcess.modelT, newdata = test.dataDum)

## ad the progressed variable to the train.data
## train.data$Progress <- y

Y <- ifelse(y==1, "yes", "no")
Y
Y_T <- ifelse(y_t==1, "yes", "no")

train.data$class <- Y
# test.dataPro$class <- Y_T

dim(train.data)
num_col <- dim(train.data)[2]

## ZScore ? 
apply( train.data[, 1:num_col], 2, FUN=function(x){c('min' = min(x), 'max' = max(x) )})

head(train.data[, 50:num_col])

train.data$class <- as.factor(train.data$class)

featurePlot( x = train.data[ , 1:20],
             y = train.data$class,
             plot = 'box',
             strip = strip.custom(par.strip.text = list( cex = 0.7)),
             scales = list(x = list(relation="free"), 
                           y = list(relation="free")
             )
)


featurePlot(x = train.data[, 1:20], 
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

rocobj_svmlinear <- roc(clin.model_svmLinear$pred$obs, 
                        clin.model_svmLinear$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear with Cinical Data",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)




#############################
#### Step 5.1 choose svm_Radial method
set.seed(123)
clin.model_svmRadial = train(class ~ ., 
                        data=train.data, 
                        method='svmRadial', 
                        #tuneLength = 8, 
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

rocobj_svmRadial <- roc(clin.model_svmRadial$pred$obs, 
                        clin.model_svmRadial$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmRadial with Cinical Data",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)



#########################################################################################
### clinical heatmap 
library(corrplot)

dim(train.data)
corrplot( cor(train.data[, 1:29]), method = "square", tl.cex = 0.5)



## there are several constant columns, remove them 
corr.data <- train.data

corr.data$histology_1_AD <- NULL
corr.data$Record_ID <- NULL

summary(corr.data$class)
corr.data$class <- ifelse(corr.data$class == "yes", "1", "0")

summary(corr.data$class)
corr.data$class <- as.numeric(corr.data$class)
corr.data$progression <- corr.data$class
corr.data$class <- NULL
dim(corr.data)
str(corr.data)

corrplot( cor( corr.data[, 1:28] ), method = "square", tl.cox = 0.15)

x <- train.data[, 1:53]
x$class <- y
x$progtype1_systemic_2locoregional <- NULL
x$TTP <- NULL
x$Did_the_patient_develop_a_New_Tumor_Event.Yes <- NULL

corrplot( cor(x), method = "square", tl.cex = 0.35)
colnames(x)


dim(train.data)
dim(test.data)
test.data2 <- predict(dummie.model, test.data)

test.data2 <- data.frame(test.data2)
str(test.data2)
str(train.dataDum)

# train.dataDum <- data.frame(train.dataDum)
# predict(preProcess.model, train.dataDum)
# test.data2$`Race.Black African American` <- NULL
# test.dataDum$`Race.Black African American` <- NULL
test.data3 <- predict(preProcess.model, test.data2)
test.data3$Race.Unknown.Refused <- rep( 0, nrow(test.data3))


colnames(test.data3)
colnames(train.data)


# str(train.data)
# str(test.data3)
test.data3$class <- Y_T
summary(test.data3$class)
test.data3$class <- as.factor(test.data3$class)
summary(test.data3$class)

predicted2 <- predict(model_svmLinear, test.data3) 

head(predicted2)

confusionMatrix(reference = test.data3$class, data = predicted2, mode='everything', positive='yes')












###########################################

dim(test.data3)
## ZScore ? 
apply( train.data[, 1:60], 2, FUN=function(x){c('min' = min(x), 'max' = max(x) )})

head(test.data3[, 50:60])

#############################
#### Step 4.1 choose svm_linear method
set.seed(100)

model_svmRadial = train(class ~ ., 
                        data=train.data, 
                        method='svmRadial', 
                        tuneLength = 8, 
                        metric='ROC', 
                        trControl = fitControl
)

#############################
## briefly check the svmLinear results
model_svmRadial 
# model_svmLinear$pred$yes
# model_svmLinear$pred$no


varimp_svmRadial <- varImp(model_svmRadial)
plot(varimp_svmRadial, main="Variable Importance with svmLinear")

rocobj_svmRadial <- roc(model_svmRadial$pred$obs, model_svmRadial$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmRadial with Cinical Data",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)


predicted_svmRadial <- predict(model_svmRadial, test.data3) 

head(predicted_svmRadial)

confusionMatrix(reference = test.data3$class, data = predicted_svmRadial, mode='everything', positive='yes')







##################################################################################
####  Section VII
####  ensemble predictions from multiple models using caretEnsemble
#### 

library(caretEnsemble)

# Stacking Algorithms - Run multiple algos in one call.
trainControl <- trainControl(method="repeatedcv", 
                             number=10, 
                             repeats=5,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

algorithmList <- c('rf', 'knn', 'earth', 'xgbDART', 'svmRadial', 'svmLinear')


################################################################################
## Run all algorithms in the list: 
set.seed(100)


dim(train.data)
clin.models <- caretList(class ~ ., 
                    data=train.data, 
                    trControl=trainControl, 
                    methodList=algorithmList) 


################################################################################
## check resample() results

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
## Plot multi ROCs in one plot
rocobj_models <- roc(clin.models$rf$pred$obs, 
                     clin.models$rf$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(clin.models$svmRadial$pred$obs, 
                     clin.models$svmRadial$pred$yes, 
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


rocobj_models <- roc(clin.models$svmLinear$pred$obs, 
                     clin.models$svmLinear$pred$yes, 
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

rocobj_models <- roc(clin.models$xgbDART$pred$obs, 
                     clin.models$xgbDART$pred$yes, 
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


legend("bottomright", 
       legend=c( "RandomForest", "svmRadial", "svmLinear", "xgbDART", "MARS", "knn" ), 
       col=c( "darkblue", "green", "red", "black", "yellow", "pink" ), 
       lwd=4
)





dim(test.data3)
dim(train.data)

clin.models <- caretList(class ~ ., 
                         data=test.data3, 
                         trControl=trainControl, 
                         methodList=algorithmList) 
