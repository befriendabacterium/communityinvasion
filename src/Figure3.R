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
library(tibble)

# LOAD DATA ---------------------------------------------------------------
source('src/X_data_prepper.R')

#colour palette for functional group 4
heatcols<-RColorBrewer::brewer.pal(11,'RdBu')
heatcols<-rev(grDevices::colorRampPalette(heatcols)(678))
palette_heat<-heatcols[as.factor(composition$C18)]

# FIGURE 3 PANEL ----------------------------------------------------------

grDevices::tiff('outputs/figures/Figure3.tiff', res=300, units='in', width=12.5, height=8.5)

layout(matrix(c(1,1,2,2,2,3,3,3,4,4,4,
                1,1,5,5,5,6,6,6,7,7,7), nrow=2, byrow=T))
par(mar=c(7,0,2,3))

# Y AXIS LABEL ------------------------------------------------------------

plot.new()
text(0.5,0.55,expression(paste(italic("P. fluorescens"))), cex=2)
text(0.5,0.50," invader survival ", cex=2)
text(0.5,0.45,"(cells/ml at 24 h)", cex=2)

# Plot C18  ---------------------------------------------------------------

log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))
log.ticks.x<-sapply(seq(0,6,1), function(i) as.expression(bquote(10^ .(i))))
plot(invexp.data$sbw25.cfu.24h~composition$C18,
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),
     xlab='',ylab='',
     cex.axis=1.5,
     xaxt='n',yaxt='n')
axis(1,seq(0,6,1),log.ticks.x, cex.axis=1.4, padj=0.3)
axis(2,seq(3,7,1),log.ticks.y[4:8], las=2, cex.axis=1.5)
mtext('Starting abundance of C18', side=1, line=3.5, cex=0.8)
mtext('(number of sequences in community inoculum)', side=1, line=5.5, cex=0.8)

model<-lm(invexp.data$sbw25.cfu.24h~composition$C18)
text(4.5,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)
summary(model)

#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-composition$C18
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

# PLOT CPM7 ---------------------------------------------------------------

log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))
log.ticks.x<-sapply(seq(3,6,1), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$sbw25.cfu.24h~divandgrowth$CPM7, 
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),xlim=c(3,6.25),
     xlab='',ylab='',
     xaxt='n',yaxt='n')
axis(1,seq(3,6,1),log.ticks.x, cex.axis=1.5, padj=0.3)

mtext('Cell density of the resident community at 7 days', side=1, line=3.5, cex=0.8)
mtext('(cells/ml)', side=1, line=5.5, cex=0.8)


model<-lm(invexp.data$sbw25.cfu.24h~divandgrowth$CPM7)
text(5.6,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)

summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$CPM7
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

# PLOT CO2 at 7 days ---------------------------------------------------------------
log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$sbw25.cfu.24h~divandgrowth$mgCO2.7, 
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),xlim=c(0,1.1),
     xlab='',
     ylab='',
     xaxt='n', yaxt='n')
axis(1,seq(0,1,0.2),seq(0,1,0.2),cex.axis=1.5, padj=0.3)

mtext("Respiration of the resident community at 7 days", side=1, line=3.5, cex=0.8)
mtext(expression("(mg CO"[2]*"/ml)"), side=1, line=5.5, cex=0.8)

model<-lm(invexp.data$sbw25.cfu.24h~divandgrowth$mgCO2.7)
text(0.88,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)

summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$mgCO2.7
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

par(mar=c(8,0,1,3))

# Plot Simpson's diversity ---------------------------------------------------------------

log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$sbw25.cfu.24h~divandgrowth$simp, 
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),xlim=c(0,1),
     xlab='',ylab='',
     xaxt='n',yaxt='n')
axis(1,seq(0,1,0.2),seq(0,1,0.2),cex.axis=1.5, padj=0.3)
axis(2,seq(3,7,1),log.ticks.y[4:8], las=2, cex.axis=1.5)
mtext('Simpson\'s Diversity Index (OTU-level)', side=1, line=3.5, cex=0.8)
mtext('in resident community inoculum',side=1,line=5.5, cex=0.8)
model<-lm(invexp.data$sbw25.cfu.24h~divandgrowth$simp)
text(0.8,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)

summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$simp
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

# Plot Rao's diversity ---------------------------------------------------------------

log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$sbw25.cfu.24h~divandgrowth$rao,
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),xlim=c(0,0.2),
     xlab='', ylab=expression(paste(italic("P. fluorescens"), " invasion at 24 hours (cells/ml)")),
     yaxt='n')
#axis(1,seq(0,0.06,0.02),seq(0,0.06,0.02),cex.axis=1.5, padj=0.3)
#axis(2,seq(3,7,1),log.ticks.y[4:8], las=2)
mtext('Rao\'s Phylogenetic Diversity (OTU-level)', side=1, line=3.5, cex=0.8)
mtext('(mean estimated number of 16S substitutions per 16S)', side=1, line=5.5, cex=0.8)

model<-lm(invexp.data$sbw25.cfu.24h~divandgrowth$rao)
text(0.175,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)
summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$rao
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

# PLOT Phylogenetic distance ---------------------------------------------------------------

log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$sbw25.cfu.24h~divandgrowth$sbw25_phydist, 
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),
     xlab='',
     ylab='',
     xaxt='n',yaxt='n')
axis(1,log10(c(0.125,0.25,0.5,1,2)),c(0.0125,0.25,0.5,1,2), cex.axis=1.4, padj=0.3)
#axis(2,seq(3,7,1),log.ticks.y[4:8], las=2)
mtext('Community OTUs-invader phylogenetic distance', side=1, line=3.5, cex=0.8)
mtext('(weighted mean estimated number of 16S substitutions)', side=1, line=5.5, cex=0.8)

model<-lm(invexp.data$sbw25.cfu.24h~divandgrowth$sbw25_phydist)
text(0.12,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)
summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$sbw25_phydist
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

dev.off()