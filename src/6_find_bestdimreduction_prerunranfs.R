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

#read in invasion data
invexp.data<-read.csv('inputs/3_ready/invasion/invasiondata_cfu_assaymeans_matched.csv', row.names = 1)
invexp.data<-log10(invexp.data+1)

# RESULTS 3: Random forest at OTU level -----------------------------------------------------------------

#read in relevant OTU table
composition<-read.csv('inputs/3_ready/composition/composition_otu_matched.csv', row.names = 1)
#add total seqs column
composition<-cbind(composition,Total=rowSums(composition))
#add total seqs column composition<-cbind(composition,Total=rowSums(composition)) #log transform
composition<-log10(composition+1)


##sbw25
#24h
set.seed(1234)
#sbw25.24h.ranf_OTU<-randomForest(composition,invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.24h.ranf_OTU,"outputs/randomforest/OTU/sbw25.24h.ranf_OTU.RDS")
sbw25.24h.ranf_OTU<-readRDS("outputs/randomforest/OTU/sbw25.24h.ranf_OTU.RDS")
plot(sbw25.24h.ranf_OTU)
sbw25.24h.ranf_OTU
varImp(sbw25.24h.ranf_OTU)
varImpPlot(sbw25.24h.ranf_OTU, type=1)
#saveRDS(sbw25.24h.ranf_OTU,"outputs/randomforest/OTU/sbw25.24h.ranf_OTU.RDS")

#96h
set.seed(1234)
#sbw25.96h.ranf_OTU<-randomForest(composition,invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.96h.ranf_OTU,"outputs/randomforest/OTU/sbw25.96h.ranf_OTU.RDS")
sbw25.96h.ranf_OTU<-readRDS("outputs/randomforest/OTU/sbw25.96h.ranf_OTU.RDS")
plot(sbw25.96h.ranf_OTU)
sbw25.96h.ranf_OTU
varImp(sbw25.96h.ranf_OTU)
varImpPlot(sbw25.96h.ranf_OTU, type=1)

#7d
set.seed(1234)
#sbw25.7d.ranf_OTU<-randomForest(composition,invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.7d.ranf_OTU,"outputs/randomforest/OTU/sbw25.7d.ranf_OTU.RDS")
sbw25.7d.ranf_OTU<-readRDS("outputs/randomforest/OTU/sbw25.7d.ranf_OTU.RDS")
plot(sbw25.7d.ranf_OTU)
sbw25.7d.ranf_OTU
varImp(sbw25.7d.ranf_OTU)
varImpPlot(sbw25.7d.ranf_OTU, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
#kt2440.24h.ranf_OTU<-randomForest(composition,invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.24h.ranf_OTU,"outputs/randomforest/OTU/kt2440.24h.ranf_OTU.RDS")
kt2440.24h.ranf_OTU<-readRDS("outputs/randomforest/OTU/kt2440.24h.ranf_OTU.RDS")
plot(kt2440.24h.ranf_OTU)
kt2440.24h.ranf_OTU
varImp(kt2440.24h.ranf_OTU)
varImpPlot(kt2440.24h.ranf_OTU, type=1)
#saveRDS(kt2440.24h.ranf_OTU,"outputs/randomforest/OTU/kt2440.24h.ranf_OTU.RDS")

#96h
set.seed(1234)
#kt2440.96h.ranf_OTU<-randomForest(composition,invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.96h.ranf_OTU,"outputs/randomforest/OTU/kt2440.96h.ranf_OTU.RDS")
kt2440.96h.ranf_OTU<-readRDS("outputs/randomforest/OTU/kt2440.96h.ranf_OTU.RDS")
plot(kt2440.96h.ranf_OTU)
kt2440.96h.ranf_OTU
varImp(kt2440.96h.ranf_OTU)
varImpPlot(kt2440.96h.ranf_OTU, type=1)

#7d
set.seed(1234)
#kt2440.7d.ranf_OTU<-randomForest(composition,invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.7d.ranf_OTU,"outputs/randomforest/OTU/kt2440.7d.ranf_OTU.RDS")
kt2440.7d.ranf_OTU<-readRDS("outputs/randomforest/OTU/kt2440.7d.ranf_OTU.RDS")
plot(kt2440.7d.ranf_OTU)
kt2440.7d.ranf_OTU
varImp(kt2440.7d.ranf_OTU)
varImpPlot(kt2440.7d.ranf_OTU, type=1)

# RESULTS 3: Random forest at genus level -----------------------------------------------------------------

#read in relevant genus table
composition<-read.csv('inputs/3_ready/composition/composition_genus_matched.csv', row.names = 1)
#add total seqs column 
composition<-cbind(composition,Total=rowSums(composition)) 
#log transform
composition<-log10(composition+1)

##sbw25
#24h
set.seed(1234)
#sbw25.24h.ranf_genus<-randomForest(composition,invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.24h.ranf_genus,"outputs/randomforest/genus/sbw25.24h.ranf_genus.RDS")
sbw25.24h.ranf_genus<-readRDS("outputs/randomforest/genus/sbw25.24h.ranf_genus.RDS")
plot(sbw25.24h.ranf_genus)
sbw25.24h.ranf_genus
varImp(sbw25.24h.ranf_genus)
varImpPlot(sbw25.24h.ranf_genus, type=1)
#saveRDS(sbw25.24h.ranf_genus,"outputs/randomforest/genus/sbw25.24h.ranf_genus.RDS")

#96h
set.seed(1234)
#sbw25.96h.ranf_genus<-randomForest(composition,invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.96h.ranf_genus,"outputs/randomforest/genus/sbw25.96h.ranf_genus.RDS")
sbw25.96h.ranf_genus<-readRDS("outputs/randomforest/genus/sbw25.96h.ranf_genus.RDS")
plot(sbw25.96h.ranf_genus)
sbw25.96h.ranf_genus
varImp(sbw25.96h.ranf_genus)
varImpPlot(sbw25.96h.ranf_genus, type=1)

#7d
set.seed(1234)
#sbw25.7d.ranf_genus<-randomForest(composition,invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.7d.ranf_genus,"outputs/randomforest/genus/sbw25.7d.ranf_genus.RDS")
sbw25.7d.ranf_genus<-readRDS("outputs/randomforest/genus/sbw25.7d.ranf_genus.RDS")
plot(sbw25.7d.ranf_genus)
sbw25.7d.ranf_genus
varImp(sbw25.7d.ranf_genus)
varImpPlot(sbw25.7d.ranf_genus, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
#kt2440.24h.ranf_genus<-randomForest(composition,invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.24h.ranf_genus,"outputs/randomforest/genus/kt2440.24h.ranf_genus.RDS")
kt2440.24h.ranf_genus<-readRDS("outputs/randomforest/genus/kt2440.24h.ranf_genus.RDS")
plot(kt2440.24h.ranf_genus)
kt2440.24h.ranf_genus
varImp(kt2440.24h.ranf_genus)
varImpPlot(kt2440.24h.ranf_genus, type=1)
#saveRDS(kt2440.24h.ranf_genus,"outputs/randomforest/genus/kt2440.24h.ranf_genus.RDS")

#96h
set.seed(1234)
#kt2440.96h.ranf_genus<-randomForest(composition,invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.96h.ranf_genus,"outputs/randomforest/genus/kt2440.96h.ranf_genus.RDS")
kt2440.96h.ranf_genus<-readRDS("outputs/randomforest/genus/kt2440.96h.ranf_genus.RDS")
plot(kt2440.96h.ranf_genus)
kt2440.96h.ranf_genus
varImp(kt2440.96h.ranf_genus)
varImpPlot(kt2440.96h.ranf_genus, type=1)

#7d
set.seed(1234)
#kt2440.7d.ranf_genus<-randomForest(composition,invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.7d.ranf_genus,"outputs/randomforest/genus/kt2440.7d.ranf_genus.RDS")
kt2440.7d.ranf_genus<-readRDS("outputs/randomforest/genus/kt2440.7d.ranf_genus.RDS")
plot(kt2440.7d.ranf_genus)
kt2440.7d.ranf_genus
varImp(kt2440.7d.ranf_genus)
varImpPlot(kt2440.7d.ranf_genus, type=1)

# RESULTS 3: Random forest at family level -----------------------------------------------------------------

#read in relevant family table
composition<-read.csv('inputs/3_ready/composition/composition_family_matched.csv', row.names = 1)
#add total seqs column 
composition<-cbind(composition,Total=rowSums(composition)) 
#log transform
composition<-log10(composition+1)

##sbw25
#24h
set.seed(1234)
#sbw25.24h.ranf_family<-randomForest(composition,invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.24h.ranf_family,"outputs/randomforest/family/sbw25.24h.ranf_family.RDS")
sbw25.24h.ranf_family<-readRDS("outputs/randomforest/family/sbw25.24h.ranf_family.RDS")
plot(sbw25.24h.ranf_family)
sbw25.24h.ranf_family
varImp(sbw25.24h.ranf_family)
varImpPlot(sbw25.24h.ranf_family, type=1)
#saveRDS(sbw25.24h.ranf_family,"outputs/randomforest/family/sbw25.24h.ranf_family.RDS")

#96h
set.seed(1234)
#sbw25.96h.ranf_family<-randomForest(composition,invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.96h.ranf_family,"outputs/randomforest/family/sbw25.96h.ranf_family.RDS")
sbw25.96h.ranf_family<-readRDS("outputs/randomforest/family/sbw25.96h.ranf_family.RDS")
plot(sbw25.96h.ranf_family)
sbw25.96h.ranf_family
varImp(sbw25.96h.ranf_family)
varImpPlot(sbw25.96h.ranf_family, type=1)

#7d
set.seed(1234)
#sbw25.7d.ranf_family<-randomForest(composition,invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.7d.ranf_family,"outputs/randomforest/family/sbw25.7d.ranf_family.RDS")
sbw25.7d.ranf_family<-readRDS("outputs/randomforest/family/sbw25.7d.ranf_family.RDS")
plot(sbw25.7d.ranf_family)
sbw25.7d.ranf_family
varImp(sbw25.7d.ranf_family)
varImpPlot(sbw25.7d.ranf_family, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
#kt2440.24h.ranf_family<-randomForest(composition,invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.24h.ranf_family,"outputs/randomforest/family/kt2440.24h.ranf_family.RDS")
kt2440.24h.ranf_family<-readRDS("outputs/randomforest/family/kt2440.24h.ranf_family.RDS")
plot(kt2440.24h.ranf_family)
kt2440.24h.ranf_family
varImp(kt2440.24h.ranf_family)
varImpPlot(kt2440.24h.ranf_family, type=1)
#saveRDS(kt2440.24h.ranf_family,"outputs/randomforest/family/kt2440.24h.ranf_family.RDS")

#96h
set.seed(1234)
#kt2440.96h.ranf_family<-randomForest(composition,invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.96h.ranf_family,"outputs/randomforest/family/kt2440.96h.ranf_family.RDS")
kt2440.96h.ranf_family<-readRDS("outputs/randomforest/family/kt2440.96h.ranf_family.RDS")
plot(kt2440.96h.ranf_family)
kt2440.96h.ranf_family
varImp(kt2440.96h.ranf_family)
varImpPlot(kt2440.96h.ranf_family, type=1)

#7d
set.seed(1234)
#kt2440.7d.ranf_family<-randomForest(composition,invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.7d.ranf_family,"outputs/randomforest/family/kt2440.7d.ranf_family.RDS")
kt2440.7d.ranf_family<-readRDS("outputs/randomforest/family/kt2440.7d.ranf_family.RDS")
plot(kt2440.7d.ranf_family)
kt2440.7d.ranf_family
varImp(kt2440.7d.ranf_family)
varImpPlot(kt2440.7d.ranf_family, type=1)

# RESULTS 3: Random forest at order level -----------------------------------------------------------------

#read in relevant order table
composition<-read.csv('inputs/3_ready/composition/composition_order_matched.csv', row.names = 1)
#add total seqs column 
composition<-cbind(composition,Total=rowSums(composition)) 
#log transform
composition<-log10(composition+1)

##sbw25
#24h
set.seed(1234)
#sbw25.24h.ranf_order<-randomForest(composition,invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.24h.ranf_order,"outputs/randomforest/order/sbw25.24h.ranf_order.RDS")
sbw25.24h.ranf_order<-readRDS("outputs/randomforest/order/sbw25.24h.ranf_order.RDS")
plot(sbw25.24h.ranf_order)
sbw25.24h.ranf_order
varImp(sbw25.24h.ranf_order)
varImpPlot(sbw25.24h.ranf_order, type=1)
#saveRDS(sbw25.24h.ranf_order,"outputs/randomforest/order/sbw25.24h.ranf_order.RDS")

#96h
set.seed(1234)
#sbw25.96h.ranf_order<-randomForest(composition,invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.96h.ranf_order,"outputs/randomforest/order/sbw25.96h.ranf_order.RDS")
sbw25.96h.ranf_order<-readRDS("outputs/randomforest/order/sbw25.96h.ranf_order.RDS")
plot(sbw25.96h.ranf_order)
sbw25.96h.ranf_order
varImp(sbw25.96h.ranf_order)
varImpPlot(sbw25.96h.ranf_order, type=1)

#7d
set.seed(1234)
#sbw25.7d.ranf_order<-randomForest(composition,invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.7d.ranf_order,"outputs/randomforest/order/sbw25.7d.ranf_order.RDS")
sbw25.7d.ranf_order<-readRDS("outputs/randomforest/order/sbw25.7d.ranf_order.RDS")
plot(sbw25.7d.ranf_order)
sbw25.7d.ranf_order
varImp(sbw25.7d.ranf_order)
varImpPlot(sbw25.7d.ranf_order, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
#kt2440.24h.ranf_order<-randomForest(composition,invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.24h.ranf_order,"outputs/randomforest/order/kt2440.24h.ranf_order.RDS")
kt2440.24h.ranf_order<-readRDS("outputs/randomforest/order/kt2440.24h.ranf_order.RDS")
plot(kt2440.24h.ranf_order)
kt2440.24h.ranf_order
varImp(kt2440.24h.ranf_order)
varImpPlot(kt2440.24h.ranf_order, type=1)
#saveRDS(kt2440.24h.ranf_order,"outputs/randomforest/order/kt2440.24h.ranf_order.RDS")

#96h
set.seed(1234)
#kt2440.96h.ranf_order<-randomForest(composition,invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.96h.ranf_order,"outputs/randomforest/order/kt2440.96h.ranf_order.RDS")
kt2440.96h.ranf_order<-readRDS("outputs/randomforest/order/kt2440.96h.ranf_order.RDS")
plot(kt2440.96h.ranf_order)
kt2440.96h.ranf_order
varImp(kt2440.96h.ranf_order)
varImpPlot(kt2440.96h.ranf_order, type=1)

#7d
set.seed(1234)
#kt2440.7d.ranf_order<-randomForest(composition,invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.7d.ranf_order,"outputs/randomforest/order/kt2440.7d.ranf_order.RDS")
kt2440.7d.ranf_order<-readRDS("outputs/randomforest/order/kt2440.7d.ranf_order.RDS")
plot(kt2440.7d.ranf_order)
kt2440.7d.ranf_order
varImp(kt2440.7d.ranf_order)
varImpPlot(kt2440.7d.ranf_order, type=1)

# RESULTS 3: Random forest at class level -----------------------------------------------------------------

#read in relevant class table
composition<-read.csv('inputs/3_ready/composition/composition_class_matched.csv', row.names = 1)
#add total seqs column 
composition<-cbind(composition,Total=rowSums(composition)) 
#log transform
composition<-log10(composition+1)

##sbw25
#24h
set.seed(1234)
#sbw25.24h.ranf_class<-randomForest(composition,invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.24h.ranf_class,"outputs/randomforest/class/sbw25.24h.ranf_class.RDS")
sbw25.24h.ranf_class<-readRDS("outputs/randomforest/class/sbw25.24h.ranf_class.RDS")
plot(sbw25.24h.ranf_class)
sbw25.24h.ranf_class
varImp(sbw25.24h.ranf_class)
varImpPlot(sbw25.24h.ranf_class, type=1)
#saveRDS(sbw25.24h.ranf_class,"outputs/randomforest/class/sbw25.24h.ranf_class.RDS")

#96h
set.seed(1234)
#sbw25.96h.ranf_class<-randomForest(composition,invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.96h.ranf_class,"outputs/randomforest/class/sbw25.96h.ranf_class.RDS")
sbw25.96h.ranf_class<-readRDS("outputs/randomforest/class/sbw25.96h.ranf_class.RDS")
plot(sbw25.96h.ranf_class)
sbw25.96h.ranf_class
varImp(sbw25.96h.ranf_class)
varImpPlot(sbw25.96h.ranf_class, type=1)

#7d
set.seed(1234)
#sbw25.7d.ranf_class<-randomForest(composition,invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.7d.ranf_class,"outputs/randomforest/class/sbw25.7d.ranf_class.RDS")
sbw25.7d.ranf_class<-readRDS("outputs/randomforest/class/sbw25.7d.ranf_class.RDS")
plot(sbw25.7d.ranf_class)
sbw25.7d.ranf_class
varImp(sbw25.7d.ranf_class)
varImpPlot(sbw25.7d.ranf_class, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
#kt2440.24h.ranf_class<-randomForest(composition,invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.24h.ranf_class,"outputs/randomforest/class/kt2440.24h.ranf_class.RDS")
kt2440.24h.ranf_class<-readRDS("outputs/randomforest/class/kt2440.24h.ranf_class.RDS")
plot(kt2440.24h.ranf_class)
kt2440.24h.ranf_class
varImp(kt2440.24h.ranf_class)
varImpPlot(kt2440.24h.ranf_class, type=1)
#saveRDS(kt2440.24h.ranf_class,"outputs/randomforest/class/kt2440.24h.ranf_class.RDS")

#96h
set.seed(1234)
#kt2440.96h.ranf_class<-randomForest(composition,invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.96h.ranf_class,"outputs/randomforest/class/kt2440.96h.ranf_class.RDS")
kt2440.96h.ranf_class<-readRDS("outputs/randomforest/class/kt2440.96h.ranf_class.RDS")
plot(kt2440.96h.ranf_class)
kt2440.96h.ranf_class
varImp(kt2440.96h.ranf_class)
varImpPlot(kt2440.96h.ranf_class, type=1)

#7d
set.seed(1234)
#kt2440.7d.ranf_class<-randomForest(composition,invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.7d.ranf_class,"outputs/randomforest/class/kt2440.7d.ranf_class.RDS")
kt2440.7d.ranf_class<-readRDS("outputs/randomforest/class/kt2440.7d.ranf_class.RDS")
plot(kt2440.7d.ranf_class)
kt2440.7d.ranf_class
varImp(kt2440.7d.ranf_class)
varImpPlot(kt2440.7d.ranf_class, type=1)

# RESULTS 3: Random forest at phylum level -----------------------------------------------------------------

#read in relevant phylum table
composition<-read.csv('inputs/3_ready/composition/composition_phylum_matched.csv', row.names = 1)
#add total seqs column
composition<-cbind(composition,Total=rowSums(composition))
#log transform
composition<-log10(composition+1)

##sbw25
#24h
set.seed(1234)
#sbw25.24h.ranf_phylum<-randomForest(composition,invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.24h.ranf_phylum,"outputs/randomforest/phylum/sbw25.24h.ranf_phylum.RDS")
sbw25.24h.ranf_phylum<-readRDS("outputs/randomforest/phylum/sbw25.24h.ranf_phylum.RDS")
plot(sbw25.24h.ranf_phylum)
sbw25.24h.ranf_phylum
varImp(sbw25.24h.ranf_phylum)
varImpPlot(sbw25.24h.ranf_phylum, type=1)
#saveRDS(sbw25.24h.ranf_phylum,"outputs/randomforest/phylum/sbw25.24h.ranf_phylum.RDS")

#96h
set.seed(1234)
#sbw25.96h.ranf_phylum<-randomForest(composition,invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.96h.ranf_phylum,"outputs/randomforest/phylum/sbw25.96h.ranf_phylum.RDS")
sbw25.96h.ranf_phylum<-readRDS("outputs/randomforest/phylum/sbw25.96h.ranf_phylum.RDS")
plot(sbw25.96h.ranf_phylum)
sbw25.96h.ranf_phylum
varImp(sbw25.96h.ranf_phylum)
varImpPlot(sbw25.96h.ranf_phylum, type=1)

#7d
set.seed(1234)
#sbw25.7d.ranf_phylum<-randomForest(composition,invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.7d.ranf_phylum,"outputs/randomforest/phylum/sbw25.7d.ranf_phylum.RDS")
sbw25.7d.ranf_phylum<-readRDS("outputs/randomforest/phylum/sbw25.7d.ranf_phylum.RDS")
plot(sbw25.7d.ranf_phylum)
sbw25.7d.ranf_phylum
varImp(sbw25.7d.ranf_phylum)
varImpPlot(sbw25.7d.ranf_phylum, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
#kt2440.24h.ranf_phylum<-randomForest(composition,invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.24h.ranf_phylum,"outputs/randomforest/phylum/kt2440.24h.ranf_phylum.RDS")
kt2440.24h.ranf_phylum<-readRDS("outputs/randomforest/phylum/kt2440.24h.ranf_phylum.RDS")
plot(kt2440.24h.ranf_phylum)
kt2440.24h.ranf_phylum
varImp(kt2440.24h.ranf_phylum)
varImpPlot(kt2440.24h.ranf_phylum, type=1)
#saveRDS(kt2440.24h.ranf_phylum,"outputs/randomforest/phylum/kt2440.24h.ranf_phylum.RDS")

#96h
set.seed(1234)
#kt2440.96h.ranf_phylum<-randomForest(composition,invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.96h.ranf_phylum,"outputs/randomforest/phylum/kt2440.96h.ranf_phylum.RDS")
kt2440.96h.ranf_phylum<-readRDS("outputs/randomforest/phylum/kt2440.96h.ranf_phylum.RDS")
plot(kt2440.96h.ranf_phylum)
kt2440.96h.ranf_phylum
varImp(kt2440.96h.ranf_phylum)
varImpPlot(kt2440.96h.ranf_phylum, type=1)

#7d
set.seed(1234)
#kt2440.7d.ranf_phylum<-randomForest(composition,invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.7d.ranf_phylum,"outputs/randomforest/phylum/kt2440.7d.ranf_phylum.RDS")
kt2440.7d.ranf_phylum<-readRDS("outputs/randomforest/phylum/kt2440.7d.ranf_phylum.RDS")
plot(kt2440.7d.ranf_phylum)
kt2440.7d.ranf_phylum
varImp(kt2440.7d.ranf_phylum)
varImpPlot(kt2440.7d.ranf_phylum, type=1)

# RESULTS 3: Random forest at funcgroups level -----------------------------------------------------------------

#read in relevant phylum table
composition<-read.table('inputs/3_ready/composition/funcgroups/SamplePropsOFtaxaClus_Time0_NL_Average_StopStep-91_AbundMean.dat', row.names = 1, header=T)

##sbw25
#24h
set.seed(1234)
#sbw25.24h.ranf_funcgroups<-randomForest(composition,invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.24h.ranf_funcgroups,"outputs/randomforest/funcgroups/sbw25.24h.ranf_funcgroups.RDS")
sbw25.24h.ranf_funcgroups<-readRDS("outputs/randomforest/funcgroups/sbw25.24h.ranf_funcgroups.RDS")
plot(sbw25.24h.ranf_funcgroups)
sbw25.24h.ranf_funcgroups
varImp(sbw25.24h.ranf_funcgroups)
varImpPlot(sbw25.24h.ranf_funcgroups, type=1)
#saveRDS(sbw25.24h.ranf_funcgroups,"outputs/randomforest/funcgroups/sbw25.24h.ranf_funcgroups.RDS")

#96h
set.seed(1234)
#sbw25.96h.ranf_funcgroups<-randomForest(composition,invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.96h.ranf_funcgroups,"outputs/randomforest/funcgroups/sbw25.96h.ranf_funcgroups.RDS")
sbw25.96h.ranf_funcgroups<-readRDS("outputs/randomforest/funcgroups/sbw25.96h.ranf_funcgroups.RDS")
plot(sbw25.96h.ranf_funcgroups)
sbw25.96h.ranf_funcgroups
varImp(sbw25.96h.ranf_funcgroups)
varImpPlot(sbw25.96h.ranf_funcgroups, type=1)

#7d
set.seed(1234)
#sbw25.7d.ranf_funcgroups<-randomForest(composition,invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.7d.ranf_funcgroups,"outputs/randomforest/funcgroups/sbw25.7d.ranf_funcgroups.RDS")
sbw25.7d.ranf_funcgroups<-readRDS("outputs/randomforest/funcgroups/sbw25.7d.ranf_funcgroups.RDS")
plot(sbw25.7d.ranf_funcgroups)
sbw25.7d.ranf_funcgroups
varImp(sbw25.7d.ranf_funcgroups)
varImpPlot(sbw25.7d.ranf_funcgroups, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
#kt2440.24h.ranf_funcgroups<-randomForest(composition,invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.24h.ranf_funcgroups,"outputs/randomforest/funcgroups/kt2440.24h.ranf_funcgroups.RDS")
kt2440.24h.ranf_funcgroups<-readRDS("outputs/randomforest/funcgroups/kt2440.24h.ranf_funcgroups.RDS")
plot(kt2440.24h.ranf_funcgroups)
kt2440.24h.ranf_funcgroups
varImp(kt2440.24h.ranf_funcgroups)
varImpPlot(kt2440.24h.ranf_funcgroups, type=1)
#saveRDS(kt2440.24h.ranf_funcgroups,"outputs/randomforest/funcgroups/kt2440.24h.ranf_funcgroups.RDS")

#96h
set.seed(1234)
#kt2440.96h.ranf_funcgroups<-randomForest(composition,invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.96h.ranf_funcgroups,"outputs/randomforest/funcgroups/kt2440.96h.ranf_funcgroups.RDS")
kt2440.96h.ranf_funcgroups<-readRDS("outputs/randomforest/funcgroups/kt2440.96h.ranf_funcgroups.RDS")
plot(kt2440.96h.ranf_funcgroups)
kt2440.96h.ranf_funcgroups
varImp(kt2440.96h.ranf_funcgroups)
varImpPlot(kt2440.96h.ranf_funcgroups, type=1)

#7d
set.seed(1234)
#kt2440.7d.ranf_funcgroups<-randomForest(composition,invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.7d.ranf_funcgroups,"outputs/randomforest/funcgroups/kt2440.7d.ranf_funcgroups.RDS")
kt2440.7d.ranf_funcgroups<-readRDS("outputs/randomforest/funcgroups/kt2440.7d.ranf_funcgroups.RDS")
plot(kt2440.7d.ranf_funcgroups)
kt2440.7d.ranf_funcgroups
varImp(kt2440.7d.ranf_funcgroups)
varImpPlot(kt2440.7d.ranf_funcgroups, type=1)

# RESULTS 3: Random forest at PCoA level -----------------------------------------------------------------

#read in relevant pcoa
pcoa<-readRDS('inputs/3_ready/composition/PCoA/communities_pcoa.RDS')
composition<-pcoa$li[,1:5]

##sbw25
#24h
set.seed(1234)
#sbw25.24h.ranf_pcoa<-randomForest(composition,invexp.data$sbw25.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.24h.ranf_pcoa,"outputs/randomforest/pcoa/sbw25.24h.ranf_pcoa.RDS")
sbw25.24h.ranf_pcoa<-readRDS("outputs/randomforest/pcoa/sbw25.24h.ranf_pcoa.RDS")
plot(sbw25.24h.ranf_pcoa)
sbw25.24h.ranf_pcoa
varImp(sbw25.24h.ranf_pcoa)
varImpPlot(sbw25.24h.ranf_pcoa, type=1)
#saveRDS(sbw25.24h.ranf_pcoa,"outputs/randomforest/pcoa/sbw25.24h.ranf_pcoa.RDS")

#96h
set.seed(1234)
#sbw25.96h.ranf_pcoa<-randomForest(composition,invexp.data$sbw25.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.96h.ranf_pcoa,"outputs/randomforest/pcoa/sbw25.96h.ranf_pcoa.RDS")
sbw25.96h.ranf_pcoa<-readRDS("outputs/randomforest/pcoa/sbw25.96h.ranf_pcoa.RDS")
plot(sbw25.96h.ranf_pcoa)
sbw25.96h.ranf_pcoa
varImp(sbw25.96h.ranf_pcoa)
varImpPlot(sbw25.96h.ranf_pcoa, type=1)

#7d
set.seed(1234)
#sbw25.7d.ranf_pcoa<-randomForest(composition,invexp.data$sbw25.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(sbw25.7d.ranf_pcoa,"outputs/randomforest/pcoa/sbw25.7d.ranf_pcoa.RDS")
sbw25.7d.ranf_pcoa<-readRDS("outputs/randomforest/pcoa/sbw25.7d.ranf_pcoa.RDS")
plot(sbw25.7d.ranf_pcoa)
sbw25.7d.ranf_pcoa
varImp(sbw25.7d.ranf_pcoa)
varImpPlot(sbw25.7d.ranf_pcoa, type=1)

##kt2440

dev.off()

#24h
set.seed(1234)
#kt2440.24h.ranf_pcoa<-randomForest(composition,invexp.data$kt2440.cfu.24h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.24h.ranf_pcoa,"outputs/randomforest/pcoa/kt2440.24h.ranf_pcoa.RDS")
kt2440.24h.ranf_pcoa<-readRDS("outputs/randomforest/pcoa/kt2440.24h.ranf_pcoa.RDS")
plot(kt2440.24h.ranf_pcoa)
kt2440.24h.ranf_pcoa
varImp(kt2440.24h.ranf_pcoa)
varImpPlot(kt2440.24h.ranf_pcoa, type=1)
#saveRDS(kt2440.24h.ranf_pcoa,"outputs/randomforest/pcoa/kt2440.24h.ranf_pcoa.RDS")

#96h
set.seed(1234)
#kt2440.96h.ranf_pcoa<-randomForest(composition,invexp.data$kt2440.cfu.96h,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.96h.ranf_pcoa,"outputs/randomforest/pcoa/kt2440.96h.ranf_pcoa.RDS")
kt2440.96h.ranf_pcoa<-readRDS("outputs/randomforest/pcoa/kt2440.96h.ranf_pcoa.RDS")
plot(kt2440.96h.ranf_pcoa)
kt2440.96h.ranf_pcoa
varImp(kt2440.96h.ranf_pcoa)
varImpPlot(kt2440.96h.ranf_pcoa, type=1)

#7d
set.seed(1234)
#kt2440.7d.ranf_pcoa<-randomForest(composition,invexp.data$kt2440.cfu.7d,importance=T, proximity = T, ntree=1000)
#saveRDS(kt2440.7d.ranf_pcoa,"outputs/randomforest/pcoa/kt2440.7d.ranf_pcoa.RDS")
kt2440.7d.ranf_pcoa<-readRDS("outputs/randomforest/pcoa/kt2440.7d.ranf_pcoa.RDS")
plot(kt2440.7d.ranf_pcoa)
kt2440.7d.ranf_pcoa
varImp(kt2440.7d.ranf_pcoa)
varImpPlot(kt2440.7d.ranf_pcoa, type=1)

# COMPARE MODELS ----------------------------------------------------------

sbw25_24h<-list(sbw25.24h.ranf_OTU,sbw25.24h.ranf_genus,sbw25.24h.ranf_family,sbw25.24h.ranf_order,sbw25.24h.ranf_class,sbw25.24h.ranf_phylum,sbw25.24h.ranf_funcgroups, sbw25.24h.ranf_pcoa)
sbw25_96h<-list(sbw25.96h.ranf_OTU,sbw25.96h.ranf_genus,sbw25.96h.ranf_family,sbw25.96h.ranf_order,sbw25.96h.ranf_class,sbw25.96h.ranf_phylum,sbw25.96h.ranf_funcgroups, sbw25.96h.ranf_pcoa)
sbw25_7d<-list(sbw25.7d.ranf_OTU,sbw25.7d.ranf_genus,sbw25.7d.ranf_family,sbw25.7d.ranf_order,sbw25.7d.ranf_class,sbw25.7d.ranf_phylum,sbw25.7d.ranf_funcgroups, sbw25.7d.ranf_pcoa)
kt2440_24h<-list(kt2440.24h.ranf_OTU,kt2440.24h.ranf_genus,kt2440.24h.ranf_family,kt2440.24h.ranf_order,kt2440.24h.ranf_class,kt2440.24h.ranf_phylum,kt2440.24h.ranf_funcgroups, sbw25.7d.ranf_pcoa)
kt2440_96h<-list(kt2440.96h.ranf_OTU,kt2440.96h.ranf_genus,kt2440.96h.ranf_family,kt2440.96h.ranf_order,kt2440.96h.ranf_class,kt2440.96h.ranf_phylum,kt2440.96h.ranf_funcgroups, kt2440.24h.ranf_pcoa)
kt2440_7d<-list(kt2440.7d.ranf_OTU,kt2440.7d.ranf_genus,kt2440.7d.ranf_family,kt2440.7d.ranf_order,kt2440.7d.ranf_class,kt2440.7d.ranf_phylum,kt2440.7d.ranf_funcgroups, kt2440.96h.ranf_pcoa)
all_randomforests<-list(sbw25_24h,sbw25_96h,sbw25_7d,kt2440_24h,kt2440_96h,kt2440_7d)

dimred_types<-c('OTU','genus','family','order','class','phylum', 'funcgroups', 'pcoa')

extracted_varexps_all<-c()

for (tp in 1:6){
  extracted_varexps<-as.numeric(sapply(strsplit(capture.output(all_randomforests[[tp]]), "% Var explained: "), "[", 2))
  extracted_varexps<-extracted_varexps[!is.na(extracted_varexps)]
  extracted_varexps_all<-rbind(extracted_varexps_all, extracted_varexps)
}

colnames(extracted_varexps_all)<-dimred_types
dimred_performancerank<-order(colMeans(extracted_varexps_all), decreasing=T)
extracted_varexps_all<-extracted_varexps_all[,dimred_performancerank]

plot(NA,xlim=c(1,8),ylim=c(0,50), type='n',xaxt='n', ylab='Variance explained (%)',xlab='Taxonomic rank')
axis(1,1:8,dimred_types[dimred_performancerank])

for (tp in 1:6){
  points(extracted_varexps_all[tp,]~jitter(1:8))
}

invader_names<-rep(c('sbw25','kt2440'),each=3)
timepoint_names<-c(24,96,168,24,96,168)
extracted_varexps_df<-as.data.frame(cbind(invader_names,timepoint_names,extracted_varexps_all))
rownames(extracted_varexps_df)<-NULL

write.csv(extracted_varexps_df,'outputs/randomforest/dimredcomparison_extractedvarexps.csv') 
saveRDS(all_randomforests,"outputs/randomforest/dimredcomparison_all_randomforests_list.RDS")
