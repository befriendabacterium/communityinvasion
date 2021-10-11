# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

source('src/X_data_prepper.R')

OTUtable<-read.csv('inputs/3_ready/composition/composition_otu_matched.csv', row.names=1)

# SUPP FIG 1 --------------------------------------------------------------

grDevices::tiff('outputs/figures/Figure2_Supplement1.tiff', res=300, units='in', width=11.69, height=8.27)

OTUtable_rel<-OTUtable/rowSums(OTUtable)*100
OTUmeans<-as.matrix(sort(apply(OTUtable_rel,2,mean),decreasing = T))
OTUse<-as.matrix(sort(apply(OTUtable_rel,2,plotrix::std.error),decreasing = T))

sub<-1:581
spcols<-c(rep('red',20),rep('grey',561))

plot(OTUmeans[sub]~sub, cex=1,
     col=spcols, pch=20,
     ylim=c(0,20), xlim=c(min(sub),max(sub)), 
     xlab='Rank', ylab='Abundance (%)')
lines(OTUmeans[1:20]~c(1:20), col='red')
lines(OTUmeans[21:581]~c(21:581), col='grey')

arrows(sub, as.numeric(OTUmeans)[sub]+as.numeric(OTUse)[sub], sub, as.numeric(OTUmeans)[sub]-as.numeric(OTUse)[sub], angle=90, length = 0, col=spcols)
text(75, seq(1,20,length.out = 20), rev(rownames(OTUmeans)[1:20]), cex=0.9, adj=0)
arrows(1:20+5, as.numeric(OTUmeans)[1:20], #from
       73, rev(seq(1,20,length.out = 20)), #to
       angle=90, length = 0, col=spcols)

dev.off()

