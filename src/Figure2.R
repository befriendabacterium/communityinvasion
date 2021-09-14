# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

# READ IN DATA ------------------------------------------------------------

library(randomForest)

#read in randomforest data
extracted_varexps_df<-read.csv('outputs/randomforest/full/extracted_varexps_df.csv') 
all_randomforests_list<-readRDS('outputs/randomforest/full/all_randomforests_list.RDS')
varimps<-read.csv('outputs/randomforest/full/varimps.csv', row.names = 1)

# PLOT --------------------------------------------------------------------

tiff('outputs/figures/Figure2.tiff', res=300, units='in', width=8.3, height=11.7)

layout(matrix(c(1,1,1,1,1,1,
                2,2,2,2,2,2,
                2,2,2,2,2,2,
                2,2,2,2,2,2), nrow=4, byrow=T))
par(mar=c(1,14.5,2,19))

barcols<-rep(c('darkgrey','lightgrey','orange'),3)


# FIGURE 1A ---------------------------------------------------------------

bars<-c(extracted_varexps_df$varexp[extracted_varexps_df$invader=='sbw25'],extracted_varexps_df$varexp[extracted_varexps_df$invader=='kt2440'])
bars<-cbind(bars[1:3],
            bars[4:6])

varplot<-barplot(bars,
                 xlab='Random forest model type',
                 ylab='',
                 ylim=c(0,100),
                 xaxt='n',
                 beside=T)

mtext('Variance explained (%)', cex=1.25, side=2, line=3)

text(3,90,expression(paste(italic("P. fluorescens"))), cex=2)
text(7,90,expression(paste(italic("P. putida"))), cex=2)

#axis(1,varplot,c(rep(types,3),rep(NA,3),rep(types,3)), tick=F, line=-1, cex.axis=0.8)
newmids<-varplot[2,]
axis(1,seq(1.5,7.5,1),c("24h","96h","168h",NA,"24h","96h","168h"), tick=F, line=-0.5, cex=0.5)

# FIGURE 1B ---------------------------------------------------------------

par(mar=c(2,15.5,1,20))

diversitymetrics<-c('simp','rao','phydist')
growthvariables<-c('CPM7','CPM14','mgCO2.7','mgCO2.14','ATP7','ATP14','mG7','mG14','mN7','mN14','mP7','mP14','mX7','mX14')
orderabundances<-rownames(varimps)[!rownames(varimps)%in%c(diversitymetrics,growthvariables,'Total')]

#DIVERSITY
varimps_diversitymetrics<-varimps[rownames(varimps)%in%diversitymetrics,]
rownames(varimps_diversitymetrics)<-c('OTU diversity (Simpson\'s)','Phylogenetic diversity (Rao\'s)', 'Phylogenetic distance from invader')

#COMPOSITION
varimps_orderabundances<-varimps[rownames(varimps)%in%orderabundances,]
#replace C label (Cluster) with FG label (Functional group)
rownames(varimps_orderabundances)<-gsub('C',rownames(varimps_orderabundances),replacement='FG')

#only do top 10 because of space limitations
#varimps_orderabundances_top10<-varimps_orderabundances[order(rowSums(varimps_orderabundances, na.rm=T), decreasing = T)[1:10],]

#TOTAL ABUNDANCE
varimps_totalabundance<-varimps[rownames(varimps)%in%'Total',]

#GROWTH
varimps_growthvariables<-varimps[rownames(varimps)%in%growthvariables,]
rownames(varimps_growthvariables)<-c('Respiration at 7 days','Respiration at 14 days','Cell density at 7 days','Cell density at 14 days','ATP yield at 7 days','ATP yield at 14 days','Cellulose degradation at 7 days','Cellulose degradation at 14 days','Chitin degradation at 7 days','Chitin degradation at 14 days','Phosphate cycling at 7 days','Phosphate cycling at 14 days','Xylose degradation at 7 days','Xylose degradation at 14 days')

varimps_organised<-rbind(varimps_diversitymetrics,varimps_orderabundances,varimps_totalabundance,NA,varimps_growthvariables)
matrix<-t(as.matrix(varimps_organised))
colnames(matrix)[24]<-NA
image(matrix[,ncol(matrix):1], rev=F, xaxt='n',yaxt='n', frame=F)
axis(2,seq(1,0,length.out=ncol(matrix)),colnames(matrix), tick=F, las=2,line=0, cex.axis=1)
AHMbook::image_scale(matrix,hcl.colors(12, "YlOrRd", rev = TRUE), x=1.3, digits=0, labels='range', size=c(0.1,0.05), cex.legend = 1.5)
mtext('IncMSE (%)', side=4,line=5, las=2, padj=-15, cex=1.5)

dev.off()