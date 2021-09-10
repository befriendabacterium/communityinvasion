# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

# LOAD & PREP  DATA ---------------------------------------------------------------

#read in relevant pcoa
obs.pcoa<-readRDS('inputs/3_ready/composition/composition_pcoa.RDS')

# SupFig5 -----------------------------------------------------------------

grDevices::tiff('outputs/figures/SupplementaryFigure5.tiff',  res=300, units='in', width=9, height=9)

plot(obs.pcoa$li[,2]~obs.pcoa$li[,1], 
     xlim=c(-60,60),ylim=c(-60,60),
     xlab='PC1', ylab='PC2', pch=16)

dev.off()
