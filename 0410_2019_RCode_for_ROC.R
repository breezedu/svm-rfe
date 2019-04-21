###################################################
##  mRNA count data analysis with DE & ROC/AUC  
##  Jeff Du
##  04/10/2019
########################################
## install required packages
## install.packages('pROC')
## install.packages('randomForest')
## 
## load pROC and randomForest


library(pROC)
library(randomForest)

#######################################################################################
## Now we will decide if a mRNA gene is 'correlated' to progression or not.
## NOTE: This method for classifying a sample as progression Yes/No
## the two category was based on the progression column in the second row from the counts table.
## 

########################################################################################
## read in the count data and the clinical information table; 
## on HP laptop
count <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/test_miRNA/miRNA_OmicData.txt", row.names = 1, header = T, sep = "\t")

progression <- read.table("D:/WorkRecord/Companies/Qiagen_Sales/201904/test_miRNA/miRNA_OmicData_Design.txt", row.names = 1, header = T, sep = "\t")

## on Lenove laptop
count <- read.table("D:/miRNA_Prj/miRNA_OmicData.txt", row.names = 1, header = T, sep = "\t")
progression <- read.table("D:/miRNA_Prj/miRNA_OmicData_Design.txt", row.names = 1, header = T, sep = "\t")

## 
## check the data frames
dim(count)
# head(count)

head(progression)
dim(progression)
count[ "FCGR2B", ]


gene.name <- "HLA-DRB4"
## gene.name <- "IL8" 
## gene.name <- "CCL20"

plotRocCalAuc(gene.name)

top.genes <- c( "HLA-DRB4", "IL8", "CCL20", ""  )


##########################################################
### get the top 43 genes from DESeq2 analysis
### set the p-value cutoff at 0.05
### saved them into a txt file
### 
### 
top.genes <- readLines("D:/WorkRecord/Companies/Qiagen_Sales/201904/test_miRNA/significant_genes_DESeq2.txt")

length(top.genes)

for(gene.name in top.genes){ 
  
  plotRocCalAuc(gene.name)
  
  }


save_plot_path <- "D:/miRNA_Prj/plots/"


##########################################################################################
###
###    The function to plot AOC curves 
###
##########################################################################################
plotRocCalAuc<-function(gene.name, save_plot_path) {
  
  ######################################################################################
  ## each row.name from the count table represents one mRNA gene
  ## pass count[row.name=gene.name, all_col] from count table to the logistic regression
  ## 
  gene1 <- as.numeric( count[ gene.name, ] )
  
  
  progress <- as.numeric(  progression$Progression.value )
  
  # length( gene1 )
  # length( progress)
  # str(progress)
  #
  plot(x=gene1, y=progress)
  
  ## fit a logistic regression to the data...
  glm.fit=glm( progress ~ gene1, family = binomial )
  
  ## fit line, to visulize the model and the fit line
  lines(gene1, glm.fit$fitted.values)
  
  ## plot a simple ROC curve; 
  roc(progress, glm.fit$fitted.values, plot=T)

  #######################################
  ##
  ## draw ROC and AUC using pROC
  ##
  #######################################
  ## 
  ## NOTE: By default, the graphs come out looking flat
  ## So, it is good to configure R so that it plots the whole ROC as a square.
  ##
  par(pty = "s") ## pty sets the aspect ratio of the plot region. Two options:
  ##                "s" - creates a square plotting region
  ##                "m" - (the default) creates a maximal plotting region 
  roc(progress, glm.fit$fitted.values, plot=TRUE)
  
  ## 
  ## To use 1-specificity or the False Positive Rate on the x-axis, have to set "legacy.axes = T".
  roc(progress, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE)
  
  ## Rename X and Y axes...
  roc(progress, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
      xlab="False Positive Percentage", 
      ylab="True Postive Percentage")
  
  ## Modify the color/width of the ROC line...
  roc(progress, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
      xlab="False Positive Percentage", 
      ylab="True Postive Percentage", 
      col="darkblue", lwd=4)
  
  ## If we want to find out the optimal threshold we can store the
  ## data used to make the ROC graph in a variable...
  roc.obj <- roc(progress, glm.fit$fitted.values, legacy.axes=TRUE, ci=T)
  str(roc.obj)
  
  ## and then extract just the information that we want from that variable.
  roc.df <- data.frame(
    tpp=roc.obj$sensitivities*100, ## tpp = true positive percentage
    fpp=(1 - roc.obj$specificities)*100, ## fpp = false positive precentage
    thresholds=roc.obj$thresholds)
  
  head(roc.df) ## head() will show us the values for the upper right-hand corner
  ## of the ROC graph, when the threshold is so low 
  ## (negative infinity) that every single sample is called "obese".
  ## Thus TPP = 100% and FPP = 100%
  
  tail(roc.df) ## tail() will show us the values for the lower left-hand corner
  ## of the ROC graph, when the threshold is so high (infinity) 
  ## that every single sample is called "not obese". 
  ## Thus, TPP = 0% and FPP = 0%
  
  
  
  ## now let's look at the thresholds between TPP 60% and 80%
  roc.df[roc.df$tpp > 60 & roc.df$tpp < 80,]
  
  print(gene.name)


    
  jpeg(paste0(save_plot_path, gene.name, "_ROC_with_AUC.jpg"), 
       width = 8, height = 6, units = "in", res= 400)
  
  ## We can calculate the area under the curve...
  rocobj <- roc(progress, glm.fit$fitted.values, ci=TRUE,
                plot=TRUE, 
                legacy.axes=TRUE, percent=TRUE, 
                xlab="False Positive Percentage", 
                ylab="True Postive Percentage", 
                col="darkblue", lwd=4, 
                print.auc=TRUE)
  
  
  #ci.95 <- paste0("95% CI = ", round(ci.auc(rocobj)[1], digits = 2), " - ", round(ci.auc(rocobj)[3], digits = 2) ) 
  
  legend("bottomright", 
         legend=c( paste("miRNA", gene.name) ), 
         col=c("darkblue", "darkblue"), 
         lwd=4)
  
    
  dev.off()
  
  
  rocobj$ci
  
  
  ### plot and save ROC with partial AUC 
  jpeg(paste0(save_plot_path, gene.name, "_ROC_PartialAUC.jpg"), 
       width = 8, height = 6, units = "in", res= 400)
  ## ...and the partial area under the curve.
  roc(progress, glm.fit$fitted.values, plot=TRUE, 
      legacy.axes=TRUE, percent=TRUE, 
      xlab="False Positive Percentage", 
      ylab="True Postive Percentage", 
      col="#377eb8", 
      lwd=4, 
      print.auc=TRUE, print.auc.x=45, 
      partial.auc=c(100, 90), 
      auc.polygon = TRUE, 
      auc.polygon.col = "#377eb822")
  
  legend("bottomright", legend=c( paste("miRNA", gene.name) ), col=c("#377eb8"), lwd=4)
  dev.off()
  
  #########################################################
  ##
  ## Next, we can fit the gene data with a random forest...
  ##
  #######################################
  rf.model <- randomForest(factor(progress) ~ gene1)
  
  ## ROC for random forest
  roc(progress, rf.model$votes[,1], plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
      xlab="False Positive Percentage", 
      ylab="True Postive Percentage", 
      col="#4daf4a", lwd=4, print.auc=TRUE)
  
  ##########################################################
  ##
  ## compare logistic regression and random forest ROC graphs..
  ##
  #######################################
  
  jpeg(paste0(save_plot_path, gene.name, "_ROC_plot_Logistic_Forest.jpg"), 
       width = 8, height = 6, units = "in", res= 400)
  
  # jpeg(file="save_plot2.png")
  
  roc(progress, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, 
      xlab="False Positive Percentage", 
      ylab="True Postive Percentage", 
      col="red", lwd=4, print.auc=TRUE)
  
  # 1. Open jpeg file
  plot.roc(progress, rf.model$votes[,1], percent=TRUE, col="blue", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
  legend("bottomright", legend=c("Logisitic Regression", "Random Forest" ), col=c("red", "blue"), lwd=4)
  legend("topleft", legend=c( paste("gene: ", gene.name) ), col=c("black"), lwd=4)
  dev.off()
  
  
  
} ## end function plotROC with Cal AUC; 




