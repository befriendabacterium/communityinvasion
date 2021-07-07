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
#install.packages('DescTools')
#install.packages('shape')

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
library(DescTools)
library(shape)

# LOAD DATA ---------------------------------------------------------------

sem_modelcomparison<-read.csv('outputs/sem/modelcomparison.csv')
sem_coeffs_all<-read.csv('outputs/sem/all_model_coeffs.csv')

#round the numeric cols to 2 dp
sem_modelcomparison[,sapply(sem_modelcomparison, is.numeric)]<-round(sem_modelcomparison[,sapply(sem_modelcomparison, is.numeric)],2)
sem_coeffs_all[,sapply(sem_coeffs_all, is.numeric)]<-round(sem_coeffs_all[,sapply(sem_coeffs_all, is.numeric)],2)

# FIGURE 3 PLOT -----------------------------------------------------------

tiff('outputs/figures/Figure4.tiff', res=300, units='in', width=12, height=4)
layout(matrix(c(1,2,3), nrow=1, byrow=T))
par(mar=c(0,0,0,0))
plot.new()

comp_col<-"blue"
prod_col<-"red"

# NO MEDIATION ------------------------------------------------------------

sem_coeffs_nomediation<-sem_coeffs_all[sem_coeffs_all$Model=='No mediation',]
#plot title
text(0.225,0.95, 'No mediation', cex=2.5)

#PRODUCTIVITY~COMPOSITION
#arrow
shape::Arrows(0.1,0.5,0.7,0.84, lwd=5, col=comp_col, arr.type = 'triangle', arr.length=0, arr.width = 0)
#add arrowhead
shape::Arrows(0.1,0.5,0.7,0.84, lwd=5, arr.col=comp_col,arr.type = 'triangle', arr.length=0.65, arr.width = 0.65, lcol='transparent')
#add coefficient
text(0.5,0.82, sem_coeffs_nomediation$Std.coefficient[sem_coeffs_nomediation$Path=='B'], cex=1.5, col=comp_col)

#INVASION~COMPOSITION
shape::Arrows(0.1,0.5,0.69,0.12, lwd=5, col=comp_col,arr.type = 'triangle', arr.length=0, arr.width = 0, lty='dashed')
#add arrowhead
shape::Arrows(0.1,0.5,0.69,0.12, lwd=5, arr.col=comp_col,arr.type = 'triangle', arr.length=0.65, arr.width = 0.65, lcol='transparent')
#add coefficient
text(0.5,0.18, sem_coeffs_nomediation$Std.coefficient[sem_coeffs_nomediation$Path=='C'], cex=1.5, col=comp_col)

#COMPOSITION MODULE
points(0.1,0.5, cex=10, pch=16, col=comp_col)
text(0.1,0.5, 'C', cex=2, col='white')

#PRODUCTIVITY MODULE
points(0.8,0.9, cex=10, pch=19, col=prod_col)
text(0.8,0.9, 'P', cex=2)

#invasion module
points(0.8,0.1, cex=10, pch=19, col='grey')
text(0.8,0.1, 'I', cex=2)

#MODEL SUMMARY
text(0,0.3,paste('Df = ',sem_modelcomparison$df[sem_modelcomparison$model=='No mediation']), cex=1.25, pos=4)
text(0,0.25,paste('AIC = ',sem_modelcomparison$AIC[sem_modelcomparison$model=='No mediation']), cex=1.25, pos=4)
text(0,0.20,paste('Delta AIC = ',sem_modelcomparison$Delta_AIC[sem_modelcomparison$model=='No mediation']), cex=1.25, pos=4)
text(0,0.15,substitute(paste("R"^2," (Invasion) = ",rsq), list(rsq=format(sem_modelcomparison$rsq_inv[sem_modelcomparison$model=='No mediation'],nsmall=2))), cex=1.25, pos=4)
text(0,0.1,paste('CFI = ',sem_modelcomparison$cfi[sem_modelcomparison$model=='No mediation']), cex=1.25, pos=4)

# PARTIAL MEDIATION ------------------------------------------------------------

sem_coeffs_partial<-sem_coeffs_all[sem_coeffs_all$Model=='Partial mediation',]

#colour mixing
comp_prod_mix<-MixColor(comp_col, prod_col, 0.5)
prod_via_comp<-sem_coeffs_partial$Std.coefficient[sem_coeffs_partial$Path=='D = A * B']

plot.new()

#plot title
text(0.275,0.95, 'Partial mediation', cex=2.5)

#PRODUCTIVITY~COMPOSITION
#arrow
shape::Arrows(0.1,0.5,0.7,0.84, lwd=5, col=comp_col, arr.type = 'triangle', arr.length=0, arr.width = 0)
#add arrowhead
shape::Arrows(0.1,0.5,0.7,0.84, lwd=5, arr.col=comp_col,arr.type = 'triangle', arr.length=0.65, arr.width = 0.65, lcol='transparent')
#add coefficient
text(0.5,0.82, sem_coeffs_partial$Std.coefficient[sem_coeffs_partial$Path=='B'], cex=1.5, col=comp_col)

#INVASION~COMPOSITION
shape::Arrows(0.1,0.5,0.69,0.12, lwd=5, col=comp_col,arr.type = 'triangle', arr.length=0, arr.width = 0, lty='dashed')
#add arrowhead
shape::Arrows(0.1,0.5,0.69,0.12, lwd=5, arr.col=comp_col,arr.type = 'triangle', arr.length=0.65, arr.width = 0.65, lcol='transparent')
#add coefficient
text(0.5,0.18, sem_coeffs_partial$Std.coefficient[sem_coeffs_partial$Path=='C'], cex=1.5, col=comp_col)

#COMPOSITION MODULE
points(0.1,0.5, cex=10, pch=16, col=comp_col)
text(0.1,0.5, 'C', cex=2, col='white')

#INVASION~PRODUCTIVITY|COMPOSITION
#arrow
shape::Arrows(0.76,0.9,0.76,0.21, lwd=5, lty='dashed', col=comp_prod_mix, arr.type='triangle', arr.width = 0, arr.length = 0)
#arrowhead
shape::Arrows(0.76,0.9,0.76,0.21, lwd=5, lty='dashed', arr.col = comp_prod_mix, lcol = "transparent", arr.type = 'triangle', arr.length=0.6, arr.width = 0.6)
#add coefficient
text(0.685,0.65, sem_coeffs_partial$Std.coefficient[sem_coeffs_partial$Path=='D = A * B'], cex=1.5, col=comp_prod_mix)

#INVASION~PRODUCTIVITY ONLY
shape::Arrows(0.84,0.9,0.84,0.21, lwd=5, lty='dashed', col=prod_col, arr.type='triangle', arr.width = 0, arr.length = 0)
#arrowhead
shape::Arrows(0.84,0.9,0.84,0.21, lwd=5, lty='dashed', arr.col = prod_col, arr.type='triangle', arr.length=0.6, arr.width = 0.6, lcol = "transparent")
#add coefficient
text(0.925,0.65, sem_coeffs_partial$Std.coefficient[sem_coeffs_partial$Path=='E = A - D'], cex=1.5, col=prod_col)

#Add blank arrow to create white bg for the coefficient
shape::Arrows(0.75,0.375,0.85,0.375, lwd=25, arr.width=0, arr.length=0, col='white')
#add coefficient
text(0.795,0.375, sem_coeffs_partial$Std.coefficient[sem_coeffs_partial$Path=='A'], cex=1.5, col='black')

#PRODUCTIVITY MODULE
points(0.8,0.9, cex=10, pch=19, col=prod_col)
text(0.8,0.9, 'P', cex=2)

#INVASION MODULE
points(0.8,0.1, cex=10, pch=19, col='grey')
text(0.8,0.1, 'I', cex=2)

#MODEL SUMMARY
text(0,0.3,paste('Df = ',sem_modelcomparison$df[sem_modelcomparison$model=='Partial']), cex=1.25, pos=4)
text(0,0.25,paste('AIC = ',sem_modelcomparison$AIC[sem_modelcomparison$model=='Partial']), cex=1.25, pos=4)
text(0,0.20,paste('Delta AIC = ',sem_modelcomparison$Delta_AIC[sem_modelcomparison$model=='Partial']), cex=1.25, pos=4)
text(0,0.15,substitute(paste("R"^2," (Invasion) = ",rsq), list(rsq=format(sem_modelcomparison$rsq_inv[sem_modelcomparison$model=='Partial'],nsmall=2))), cex=1.25, pos=4)
text(0,0.1,paste('CFI = ',sem_modelcomparison$cfi[sem_modelcomparison$model=='Partial']), cex=1.25, pos=4)

# COMPLETE MEDIATION ------------------------------------------------------------

sem_coeffs_complete<-sem_coeffs_all[sem_coeffs_all$Model=='Complete mediation',]

plot.new()

#plot title
text(0.315,0.95, 'Complete mediation', cex=2.5)

#PRODUCTIVITY~COMPOSITION
#arrow
shape::Arrows(0.1,0.5,0.7,0.84, lwd=5, col=comp_col, arr.type = 'triangle', arr.length=0, arr.width = 0)
#add arrowhead
shape::Arrows(0.1,0.5,0.7,0.84, lwd=5, arr.col=comp_col,arr.type = 'triangle', arr.length=0.65, arr.width = 0.65, lcol='transparent')
#add coefficient
text(0.5,0.82, sem_coeffs_complete$Std.coefficient[sem_coeffs_partial$Path=='B'], cex=1.5, col=comp_col)

#COMPOSITION MODULE
points(0.1,0.5, cex=10, pch=16, col=comp_col)
text(0.1,0.5, 'C', cex=2, col='white')

#INVASION~PRODUCTIVITY|COMPOSITION
#arrow
shape::Arrows(0.76,0.9,0.76,0.21, lwd=5, lty='dashed', col=comp_prod_mix, arr.type='triangle', arr.width = 0, arr.length = 0)
#arrowhead
shape::Arrows(0.76,0.9,0.76,0.21, lwd=5, lty='dashed', arr.col = comp_prod_mix, lcol = "transparent", arr.type = 'triangle', arr.length=0.6, arr.width = 0.6)
#add coefficient
text(0.685,0.65, sem_coeffs_complete$Std.coefficient[sem_coeffs_complete$Path=='D = A * B'], cex=1.5, col=comp_prod_mix)

#INVASION~PRODUCTIVITY ONLY
shape::Arrows(0.84,0.9,0.84,0.21, lwd=5, lty='dashed', col=prod_col, arr.type='triangle', arr.width = 0, arr.length = 0)
#arrowhead
shape::Arrows(0.84,0.9,0.84,0.21, lwd=5, lty='dashed', arr.col = prod_col, arr.type='triangle', arr.length=0.6, arr.width = 0.6, lcol = "transparent")
#add coefficient
text(0.925,0.65, sem_coeffs_complete$Std.coefficient[sem_coeffs_complete$Path=='E = A - D'], cex=1.5, col=prod_col)

#Add blank arrow to create white bg for the coefficient
shape::Arrows(0.75,0.375,0.85,0.375, lwd=25, arr.width=0, arr.length=0, col='white')
#add coefficient
text(0.795,0.375,sem_coeffs_complete$Std.coefficient[sem_coeffs_complete$Path=='A'], cex=1.5, col='black')

#PRODUCTIVITY MODULE
points(0.8,0.9, cex=10, pch=19, col=prod_col)
text(0.8,0.9, 'P', cex=2)

#INVASION MODULE
points(0.8,0.1, cex=10, pch=19, col='grey')
text(0.8,0.1, 'I', cex=2)

#MODEL SUMMARY
text(0,0.3,paste('Df = ',sem_modelcomparison$df[sem_modelcomparison$model=='Complete']), cex=1.25, pos=4)
text(0,0.25,paste('AIC = ',sem_modelcomparison$AIC[sem_modelcomparison$model=='Complete']), cex=1.25, pos=4)
text(0,0.20,paste('Delta AIC = ',sem_modelcomparison$Delta_AIC[sem_modelcomparison$model=='Complete']), cex=1.25, pos=4)
text(0,0.15,substitute(paste("R"^2," (Invasion) = ",rsq), list(rsq=format(sem_modelcomparison$rsq_inv[sem_modelcomparison$model=='Complete'],nsmall=2))), cex=1.25, pos=4)
text(0,0.1,paste('CFI = ',sem_modelcomparison$cfi[sem_modelcomparison$model=='Complete']), cex=1.25, pos=4)

dev.off()
