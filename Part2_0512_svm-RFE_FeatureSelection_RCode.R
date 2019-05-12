#########################################################################
####
#####################################
## 
## Feature selection with SVM-RFE
## 04/14/2019
## Jeff Du
## 
#########################################################################
## 
#########################################################################
# install necessar packages
# library(devtools)
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("shiny")
# install_github("StatsWithR/statsr")
# BiocManager::install("sigFeature", version = "3.8")
#########################################################################




#########################################################################
## 
## Use SVM-RFE from Caret package to 
## Selecte features (genes)
## since the input data is too small, only 43 samples
## it is relative difficult to run training and validation sets
## 
## 
#########################################################################


## install.packages("mlbench")
## install.packages("caret")


library("mlbench")
library("caret")
library("pROC")


#########
setwd("D:/miRNA_Prj/0507CleanedData/")

count <- read.table("Count.txt", row.names = 1, header = T, sep = "\t")
design <- read.table("Count_Design.txt", row.names = 1, header = T, sep = "\t")
top.genes <- readLines("Deseq2_top45rows.txt")




#########################################################################
## 
## as in the progression columns, there are two factors: yes and no
rfFuncs$summary <- twoClassSummary

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
                          # repeats = 3, #repeats has no meaning for resampling;
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


rfe_rf <- rfe(t(count), progression$Progression.Logit, sizes=c(1:784),
              method="rf",
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


## 
## based on above results, there are 13 genes showing optmized results
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

progress

rfe_rf_svm <- train(t(count), progress$Progressed, sizes=c(1:784),
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

count_t$Class <- progress$Progressed


ctrl <- trainControl(method="cv", 
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)

rfFit <- train(Class ~ ., data=count_t, 
               method="svmRadial", 
               preProc=c("center", "scale"), 
               metric = "ROC", 
               trControl=ctrl)

# library(pROC)
# Select a parameter setting
selectedIndices <- rfFit$pred$mtry == 39

rfFit$resample

head(rfFit$pred)

selectedIndices

rfFit$pred$obs[selectedIndices]
rfFit$pred$Yes[selectedIndices]
# rfFit$pred$M
# Plot:

rfFit$results
rfFit$bestTune

gene_13 <- rfFit$pred$mtry == 13 
gene_13

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


## 04/15/2019
## plot roc for both training and the validation set

rfe_rf_svm

for_lift <- data.frame(Class = rfe_rf_svm$pred$obs, rf = rfe_rf_svm$pred$yes, resample = rfe_rf_svm$pred$Resample)

for_lift_f10r3 <- for_lift[for_lift$resample == "Fold10.Rep3", ]
for_lift <- for_lift_f10r3

lift_df <-  data.frame()

dim(for_lift)
head( for_lift )

for_lift$resample


for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(Class ~ rf, data = fold_df, class = "yes")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}


## only plot roc for resample: Fold10.Rep3 
## 

lift_obj_data <- lift(Class ~ rf, data = fold_df, class = "yes")$data

lift_obj <- lift(Class ~ rf, data = for_lift, class = "yes")


# Plot ROC ----------------------------------------------------------------



dim(lift_df)

roc()

ggplot(lift_df) +
  geom_line(aes(1 - Sp, Sn, color = fold)) +
  scale_color_discrete(guide = guide_legend(title = "Fold"))


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



#################################
#################################

install.packages("caret")
library(caret)
T2DData<-read.table("T2DRecords.csv",header=T,sep=",")
head(T2DData)
T2DData$T2D<-as.factor(T2DData$T2D)
levels(T2DData$T2D)<-c('No','Yes')
set.seed(330)# Use any positive integer to "shuffle" the data records
partitionRule<-createDataPartition(T2DData$T2D, p = 0.7, list = FALSE)
trainingSet<-T2DData[partitionRule,]
testingSet<-T2DData[-partitionRule,]

splitRule<-trainControl(method="repeatedcv",repeats=3,number=10,classProbs = TRUE, summaryFunction = twoClassSummary)
gbmModel<-train(T2D ~ .,data = trainingSet, trControl=splitRule, method = "gbm",preProc = c("center", "scale"),metric = "ROC")

gbmTest<- predict(gbmModel, newdata = testingSet)
confusionMatrix(data=gbmTest,testingSet$T2D)



glmModel<-train(T2D ~ .,data = trainingSet, trControl=splitRule, method = "glm",preProc = c("center", "scale"),metric = "ROC")
glmTest<- predict(glmModel, newdata = testingSet)
confusionMatrix(data=glmTest,testingSet$T2D)

resamps <- resamples(list(GBM = gbmModel,GLM = glmModel))
summary(resamps)
trellis.par.set(caretTheme())
dotplot(resamps, metric = "ROC")



###################################

