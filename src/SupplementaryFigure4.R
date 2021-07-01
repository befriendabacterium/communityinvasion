# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

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

# LOAD & PREP  DATA ---------------------------------------------------------------

#read in relevant pcoa
obs.pcoa<-readRDS('inputs/3_ready/composition/composition_pcoa.RDS')

# SupFig5 -----------------------------------------------------------------

tiff('outputs/figures/SupplementaryFigure4.tiff',  res=300, units='in', width=9, height=9)

plot(obs.pcoa$li[,2]~obs.pcoa$li[,1], 
     xlim=c(-60,60),ylim=c(-60,60),
     xlab='PC1', ylab='PC2', pch=16)

dev.off()
