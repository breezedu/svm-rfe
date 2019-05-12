

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

## remove results categories from the clinical table

clinical.data$progtype1_systemic_2locoregional <- NULL
clinical.data$Did_the_patient_develop_a_New_Tumor_Event <- NULL
clinical.data$Alive_or_Dead <- NULL
clinical.data$locoregional <- NULL
clinical.data$TTP <- NULL
clinical.data$New_Tumor_Event <- NULL

head(clinical.data[ ,1:8])




####################################################
## 
## Sectopm II
## split the dataset into training and validation 
## 
#################################################### 

set.seed(123)

trainRowNum <- createDataPartition(clinical.data$Progressed, p=0.9, list = F)

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

lmProfile <- rfe(x=train.data[, 1:53], y=train.data$class,
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


#############################
#### Step 4.1 choose svm_linear method
model_svmLinear = train(class ~ ., 
                        data=train.data, 
                        method='svmLinear', 
                        tuneLength = 4, 
                        metric='ROC', 
                        trControl = fitControl
                        )

#############################
## briefly check the svmLinear results
model_svmLinear 
# model_svmLinear$pred$yes
# model_svmLinear$pred$no


varimp_svmLinear <- varImp(model_svmLinear)
plot(varimp_svmLinear, main="Variable Importance with svmLinear")


########################################################################################
## step 4.2
## plot roc for the training data
## We can calculate the area under the curve...
## Select a parameter setting
## selectedIndices <- model_mars2$pred
library(pROC)

rocobj_svmlinear <- roc(model_svmLinear$pred$obs, model_svmLinear$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear with Cinical Data",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)

library(corrplot)

corrplot( cor(train.data[, 1:53]), method = "square", tl.cex = 0.5)

## there are several constant columns, remove them
corr.data <- train.data
corr.data$histology_1_AD <- NULL
corr.data$TNM_Stage.T1b <- NULL
corr.data$TNM_Stage.T1bN0M0 <- NULL
dim(corr.data)

corrplot( cor( corr.data[, 1:54] ), method = "square", tl.cox = 0.2)

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

