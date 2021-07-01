# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")

# LOAD PACKAGES ------------------------------------------------------------

library(randomForest)

# READ IN DATA ------------------------------------------------------------

extracted_varexps_df<-read.csv('outputs/randomforest/extracted_varexps_df.csv') 
all_randomforests<-readRDS("outputs/randomforest/all_randomforests_list.RDS")

# PLOT DATA ---------------------------------------------------------------
  
dimred_types<-c('OTU','genus','family','order','class','phylum', 'funcgroups', 'pcoa')
nvars<-c(581,360,177,90,42,18,20,5)


extracted_varexps_all<-c()

for (tp in 1:6){
  extracted_varexps<-as.numeric(sapply(strsplit(capture.output(all_randomforests[[tp]]), "% Var explained: "), "[", 2))
  extracted_varexps<-extracted_varexps[!is.na(extracted_varexps)]
  extracted_varexps_all<-rbind(extracted_varexps_all, extracted_varexps)
}

colnames(extracted_varexps_all)<-dimred_types

plot(extracted_varexps_all[tp,]~nvars, type='n', ylim=c(0,50))
for (tp in 1:6){
  points(extracted_varexps_all[tp,]~nvars, ylim=c(0,50))
  #lines(extracted_varexps_all[tp,]~nvars)
}

varexp_mean<-colMeans(extracted_varexps_all)
varexp_se<-plotrix::std.error(extracted_varexps_all)

# PLOT --------------------------------------------------------------------

tiff('outputs/figures/SupplementaryFigure1.tiff', res=300, units='in', width=11.7, height=8.3)

colours<-RColorBrewer::brewer.pal(8, 'Set2')

layout(matrix(c(1,1,1,1,1,1,2,2), nrow=1, byrow=T))
par(mar=c(10,10,5,0))

plot(varexp_mean~nvars, ylim=c(0,32), 
     col=colours, pch=20, cex=3, 
     cex.lab=2, cex.axis=1.5,
     xlab='',ylab='')

mtext('Number of dimensions of composition/explanatory variables permuted', side=1, line=5, cex=1.5) 
mtext('Mean variance explained across random forests (%)', side=2, line=5, cex=1.5) 

arrows(nvars,varexp_mean+varexp_se,
       nvars,varexp_mean-varexp_se,
       length=0, col=colours)

par(mar=c(0,0,0,0))
plot.new()
legend.labs<-c('Principal coordinates', 'Phylum', 'Functional groups', 'Class', 'Order', 'Family', 'Genus', 'OTU (no reduction)')
legend(0,0.75,legend.labs,fill=colours[order(nvars)], cex=2, bty='n', title='Dimensionality \n reduction')

dev.off()

