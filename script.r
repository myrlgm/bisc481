######################################
# 24/10/2016
# Assignment #3
# BISC 481
######################################

## Install packages
# Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()
# DNAshapeR
biocLite("DNAshapeR")
# for writeXStringSet
biocLite("Biostrings")
# Caret
install.packages("caret")

## Initialization
library(DNAshapeR)
library(caret)
library(ROCR)
library(Biostrings)

## Encode Mad feature vectors
MadFeaturesMer <- encodeSeqShape(
  "gcPBM/Mad.txt.fa",
  getShape("gcPBM/Mad.txt.fa"),
  c("1-mer")
)

MadFeaturesMerShape <- encodeSeqShape(
  "gcPBM/Mad.txt.fa",
  getShape("gcPBM/Mad.txt.fa"),
  c("1-mer", "1-shape")
)

## Encode Max feature vectors
MaxFeaturesMer <- encodeSeqShape(
  "gcPBM/Max.txt.fa",
  getShape("gcPBM/Max.txt.fa"),
  c("1-mer")
)

MaxFeaturesMerShape <- encodeSeqShape(
  "gcPBM/Max.txt.fa",
  getShape("gcPBM/Max.txt.fa"),
  c("1-mer", "1-shape")
)

## Encode Myc feature vectors
MycFeaturesMer <- encodeSeqShape(
  "gcPBM/Myc.txt.fa",
  getShape("gcPBM/Myc.txt.fa"),
  c("1-mer")
)

MycFeaturesMerShape <- encodeSeqShape(
  "gcPBM/Myc.txt.fa",
  getShape("gcPBM/Myc.txt.fa"),
  c("1-mer", "1-shape")
)

# train control parameters
tc = trainControl (method = "cv", number = 10,
                   savePredictions = TRUE , allowParallel = TRUE, classProbs = TRUE)

# Mad linear fits
MadFitMer = train(affinity~.,
            data = 
              df <- data.frame(
                affinity=read.table("gcPBM/Mad.txt")$V2,
                MadFeaturesMer),
            trControl = tc,
            method = "lm")
summary(MadFitMer)

MadFitMerShape = train(affinity~.,
                  data = 
                    df <- data.frame(
                      affinity=read.table("gcPBM/Mad.txt")$V2,
                      MadFeaturesMerShape),
                  trControl = tc,
                  method = "lm")
summary(MadFitMerShape)

# Max linear fits
MaxFitMer = train(affinity~.,
                  data = 
                    df <- data.frame(
                      affinity=read.table("gcPBM/Max.txt")$V2,
                      MaxFeaturesMer),
                  trControl = tc,
                  method = "lm")
summary(MaxFitMer)

MaxFitMerShape = train(affinity~.,
                  data = 
                    df <- data.frame(
                      affinity=read.table("gcPBM/Max.txt")$V2,
                      MaxFeaturesMerShape),
                  trControl = tc,
                  method = "lm")
summary(MaxFitMerShape)

# Myc linear fits
MycFitMer = train(affinity~.,
                  data = 
                    df <- data.frame(
                      affinity=read.table("gcPBM/Myc.txt")$V2,
                      MycFeaturesMer),
                  trControl = tc,
                  method = "lm")
summary(MycFitMer)

MycFitMerShape = train(affinity~.,
                  data = 
                    df <- data.frame(
                      affinity=read.table("gcPBM/Myc.txt")$V2,
                      MycFeaturesMerShape),
                  trControl = tc,
                  method = "lm")
summary(MycFitMerShape)

# the R^2 values were plotted in Excel.

# bound sequences and their plots
bound <- getShape("CTCF/bound_500.fa")
plotShape(bound$MGW)
plotShape(bound$ProT)
plotShape(bound$Roll)
plotShape(bound$HelT)

# unbound sequences and their plots
unbound <- getShape("CTCF/unbound_500.fa")
plotShape(unbound$MGW)
plotShape(unbound$ProT)
plotShape(unbound$Roll)
plotShape(unbound$HelT)

# get smaller datasets for this step because it takes longer
boundTxt <- data.frame(
  seq = paste(readDNAStringSet("CTCF/bound_30.fa")),
  isBound = "Y"
  )

unboundTxt <- data.frame(
  seq = paste(readDNAStringSet("CTCF/unbound_30.fa")),
  isBound = "N"
  )

# merge two datasets
writeXStringSet(
  c(readDNAStringSet("CTCF/bound_30.fa"), readDNAStringSet("CTCF/unbound_30.fa")),
  "CTCF/ctcf.fa"
)

boundUnboundData <- rbind(boundTxt, unboundTxt)

# both bound and unbound
bothShapes <- getShape("CTCF/ctcf.fa")

# encode feature vectors
merFeatureVector <- encodeSeqShape("CTCF/ctcf.fa", bothShapes, c("1-mer"))
merFrame <- data.frame(isBound = boundUnboundData$isBound, merFeatureVector)

merShapeFeatureVector <- encodeSeqShape("CTCF/ctcf.fa", bothShapes, c("1-mer", "1-shape"))
merShapeFrame <- data.frame(isBound = boundUnboundData$isBound, merShapeFeatureVector)

# perform prediction
merClassifier <- train(isBound ~ ., data = merFrame, trControl = tc,
               method = "glm", family = binomial, metric ="ROC")
merShapeClassifier <- train(isBound ~ ., data = merShapeFrame, trControl = tc,
                       method = "glm", family = binomial, metric ="ROC")

## Plot 1-mer AUROC
merPrediction <- prediction( merClassifier$pred$Y, merClassifier$pred$obs )
merPerformance <- performance( merPrediction, "tpr", "fpr" )
plot(merPerformance)

## Plot 1-mer+shape AUROC
merShapePrediction <- prediction( merShapeClassifier$pred$Y, merShapeClassifier$pred$obs )
merShapePerformance <- performance( merShapePrediction, "tpr", "fpr" )
plot(merShapePerformance)

## Caluculate 1-mer AUROC
merAuc <- performance(merPrediction, "auc")
merAuc <- unlist(slot(merAuc, "y.values"))

## Caluculate 1-mer+shape AUROC
merShapeAuc <- performance(merShapePrediction, "auc")
merShapeAuc <- unlist(slot(merShapeAuc, "y.values"))

