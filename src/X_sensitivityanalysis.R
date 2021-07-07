# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")

source('src/X_data_prepper.R')

invexp.belowdetection<-read.csv('inputs/3_ready/invasion/invasiondata_belowdetectionlimit_matched.csv', row.names = 1)
detectionlims<-readRDS('inputs/3_ready/invasion/detectionlims.RDS')
detectionlims<-rep(detectionlims,each=3)
names(detectionlims)<-c('sbw','sbw','sbw','kt','kt','kt')

# RANDOMLY SHUFFLE BELOW DETECTION LIMIT VALUES---------------------------------------------------------

dir.create('inputs/3_ready/invasion/belowdetection_shuffled/')

set.seed(1234)
#for 999 iterations
for (i in 1:999){
  print(i)
  invexp.data_temp<-invexp.data
  #go across each of 6 columns
  for (c in 1:ncol(invexp.data_temp)){
    
    #randomly sample between 0 and the detection limit for number of NAs/missing samples
    tail<-sample(invexp.data_temp[,c],size=sum(invexp.belowdetection[,c]), replace=F)
    
    #replace each NA/missing sample with one of the randomly generated numbers
    invexp.data_temp[,c][invexp.belowdetection[,c]]<-tail
    #write input data to csv
    write.csv(invexp.data_temp, paste('inputs/3_ready/invasion/belowdetection_shuffled/invasiondata_cfu_assaymeans_matched_shuffled_',i,'.csv',sep=''))
  }
  print(c)
}

# DO SOME TEST PLOTS OF 999TH ITERATION TO SEE HOW ACTUAL AND INFERRED COMPARE

breaky=seq(0,7,by=0.1)

layout(matrix(1:6, nrow=2, byrow=F))
for (c in 1:6){
  testhist1<-hist(invexp.data[,c], breaks=breaky, plot=F)
  testhist2<-hist(invexp.data_temp[,c], breaks=breaky, plot=F)
  histcols<-testhist2$mids<log10(detectionlims[c]+1)
  histcols[histcols==T]<-'red'
  histcols[histcols==F]<-'grey'
  plot(testhist1, col='grey',xlim=c(0,7), ylim=c(0,300), main=colnames(invexp.data_temp)[c])
  plot(testhist2, col=rep(histcols),xlim=c(0,7), ylim=c(0,300),main=colnames(invexp.data_temp)[c])
}

# VIEW DISTRIBUTIONS -------------------------------------------------------------------------

#before NA removal
hist(invexp.data)
     
#change below detection lims to NAs
invexp.data[as.matrix(invexp.belowdetection)]<-NA

#after NA removal
hist(invexp.data)

#for first experiment (sbw25 at 24h), plot predicted densities over histogram
hist(invexp.data[,1], xlim=c(0,7), probability=T)
#compute density estimates for the histogram of the invasion data
densities<-density(invexp.data[,1], na.rm = T, adjust = 1, n=1024)
lines(densities)

# RUN SENSITIVITY ANALYSIS ------------------------------------------------

#load dummy environment variable to tell data prepper this is sensitivity analysis
sensitivityanalysis<-T

filespath<-'inputs/3_ready/invasion/belowdetection_shuffled/'
files<-list.files(filespath)
filenos<-as.numeric(gsub(".*?([0-9]+).*", "\\1", files))
files<-files[order(filenos)]

aic.ranks_all<-c()

for (f in 1:length(files)){
  
  print(files[f])
  invexp.data_temp<-read.csv(paste(filespath,files[f], sep=''), row.names = 1)
  source('src/9_launch_SEM_multiple_simpsonsdiversity.R')
  aic.rank_temp<-cbind(sourcefile=files[f],aic.rank)
  aic.ranks_all<-rbind(aic.ranks_all,aic.rank_temp)
}

write.csv(aic.ranks_all,"outputs/sem/aicranks_sensitivityanalysis.csv")
aggregs<-aggregate(aic.ranks_all,list(aic.ranks_all$Modnames), mean, na.rm=T)
options(scipen=1000)
aggregs


#check what proportion of model comparisons supported Partial vs Complete
aic.ranks_all<-read.csv('outputs/sem/aicranks_sensitivityanalysis.csv')
#total 
tapply(aic.ranks_all$Delta_AICc==0,aic.ranks_all$Modnames,sum)/999
#mean delta AIC
tapply(aic.ranks_all$Delta_AICc==0,aic.ranks_all$Modnames,mean)
#mean AIC
tapply(aic.ranks_all$AICc,aic.ranks_all$Modnames,mean)
