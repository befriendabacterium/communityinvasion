# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")

# LOAD PACKAGES --------------------------------------------------

#install packages
#install.packages('ape')
#install.packages('ade4')
#install.packages('caret')
#install.packages('ggsci')
#install.packages('picante')
#install.packages('plyr')
#install.packages('RColorBrewer')
#install.packages('randomForest')
#install.packages('scales')
#install.packages('vegan')
#install.packages('lavaan')
#install.packages('semPlot')

#load libraries
library(ape)
library(ade4)
library(caret)
library(ggsci)
library(picante)
library(plotrix) #for std.error function
library(plyr)
library(randomForest)
library(RColorBrewer)
library(scales)
library(vegan)
library(lavaan)
library(semPlot)
library(tibble)

# READ IN INVASION DATA ---------------------------------------------------

source('src/X_data_prepper.R')

# RESULTS 3: Run full random forests with functional groups, diversity metrics and phenotypic variables -----------------------------------------------------------------

##sbw25
#24h
set.seed(1234)
sbw25.24h.ranf_full<-randomForest(cbind(composition,divandgrowth[,-4]),invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
saveRDS(sbw25.24h.ranf_full,"outputs/randomforest/full/sbw25.24h.ranf_full.RDS")
sbw25.24h.ranf_full<-readRDS("outputs/randomforest/full/sbw25.24h.ranf_full.RDS")
plot(sbw25.24h.ranf_full)
sbw25.24h.ranf_full
varImp(sbw25.24h.ranf_full)
varImpPlot(sbw25.24h.ranf_full, type=1)
#saveRDS(sbw25.24h.ranf_full,"outputs/randomforest/full/sbw25.24h.ranf_full.RDS")

#96h
set.seed(1234)
sbw25.96h.ranf_full<-randomForest(cbind(composition,divandgrowth[,-4]),invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
saveRDS(sbw25.96h.ranf_full,"outputs/randomforest/full/sbw25.96h.ranf_full.RDS")
sbw25.96h.ranf_full<-readRDS("outputs/randomforest/full/sbw25.96h.ranf_full.RDS")
plot(sbw25.96h.ranf_full)
sbw25.96h.ranf_full
varImp(sbw25.96h.ranf_full)
varImpPlot(sbw25.96h.ranf_full, type=1)

#7d
set.seed(1234)
sbw25.7d.ranf_full<-randomForest(cbind(composition,divandgrowth[,-4]),invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
saveRDS(sbw25.7d.ranf_full,"outputs/randomforest/full/sbw25.7d.ranf_full.RDS")
sbw25.7d.ranf_full<-readRDS("outputs/randomforest/full/sbw25.7d.ranf_full.RDS")
plot(sbw25.7d.ranf_full)
sbw25.7d.ranf_full
varImp(sbw25.7d.ranf_full)
varImpPlot(sbw25.7d.ranf_full, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
kt2440.24h.ranf_full<-randomForest(cbind(composition,divandgrowth[,-3]),invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
saveRDS(kt2440.24h.ranf_full,"outputs/randomforest/full/kt2440.24h.ranf_full.RDS")
kt2440.24h.ranf_full<-readRDS("outputs/randomforest/full/kt2440.24h.ranf_full.RDS")
plot(kt2440.24h.ranf_full)
kt2440.24h.ranf_full
varImp(kt2440.24h.ranf_full)
varImpPlot(kt2440.24h.ranf_full, type=1)
#saveRDS(kt2440.24h.ranf_full,"outputs/randomforest/full/kt2440.24h.ranf_full.RDS")

#96h
set.seed(1234)
kt2440.96h.ranf_full<-randomForest(cbind(composition,divandgrowth[,-3]),invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
saveRDS(kt2440.96h.ranf_full,"outputs/randomforest/full/kt2440.96h.ranf_full.RDS")
kt2440.96h.ranf_full<-readRDS("outputs/randomforest/full/kt2440.96h.ranf_full.RDS")
plot(kt2440.96h.ranf_full)
kt2440.96h.ranf_full
varImp(kt2440.96h.ranf_full)
varImpPlot(kt2440.96h.ranf_full, type=1)

#7d
set.seed(1234)
kt2440.7d.ranf_full<-randomForest(cbind(composition,divandgrowth[,-3]),invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
saveRDS(kt2440.7d.ranf_full,"outputs/randomforest/full/kt2440.7d.ranf_full.RDS")
kt2440.7d.ranf_full<-readRDS("outputs/randomforest/full/kt2440.7d.ranf_full.RDS")
plot(kt2440.7d.ranf_full)
kt2440.7d.ranf_full
varImp(kt2440.7d.ranf_full)
varImpPlot(kt2440.7d.ranf_full, type=1)

# COMPARE MODELS ----------------------------------------------------------

#bind all random forests into a list
all_randomforests<-list(sbw25.24h.ranf_full,sbw25.96h.ranf_full,sbw25.7d.ranf_full,
                        kt2440.24h.ranf_full,kt2440.96h.ranf_full,kt2440.7d.ranf_full)

#extract variance explained
extracted_varexps<-as.numeric(sapply(strsplit(capture.output(all_randomforests), "% Var explained: "), "[", 2))
extracted_varexps<-extracted_varexps[!is.na(extracted_varexps)]
extracted_varexps_df<-cbind(
                            invader=c('sbw25','sbw25','sbw25','kt2440','kt2440','kt2440'),
                            timepoint=c(24,96,168,24,96,168),
                            varexp=extracted_varexps
                            )

#extract variable importance
varimps<-cbind(
  varImp(sbw25.24h.ranf_full, type=1),
  varImp(sbw25.96h.ranf_full, type=1),
  varImp(sbw25.7d.ranf_full, type=1),
  NA,
  varImp(kt2440.24h.ranf_full, type=1),
  varImp(kt2440.96h.ranf_full, type=1),
  varImp(kt2440.7d.ranf_full, type=1)
)

rownames(varimps)[which(rownames(varimps)=='sbw25_phydist')]<-'phydist'

write.csv(varimps,'outputs/randomforest/full/varimps.csv')
write.csv(extracted_varexps_df,'outputs/randomforest/full/extracted_varexps_df.csv') 
saveRDS(all_randomforests,'outputs/randomforest/full/all_randomforests_list.RDS')

