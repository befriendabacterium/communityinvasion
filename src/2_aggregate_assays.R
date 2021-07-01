# INTRO -------------------------------------------------------------------

#This code takes the raw data and re-aggregates it, checking whether it matches to the old aggregated data

# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

# READ IN RAW DATA ------------------------------------------------------------

#raw invasion data
inv<-read.csv('inputs/1_raw/invasion/invasiondata_luminescence_4assays.csv')
#raw and cleaned growth data
gr<-read.csv('inputs/2_wrangled/growth/growthdata_cleaned.csv')

# AGGREGATE LAB ASSAYS (Growth and invasion) ---------------------------------------

#calculate mean of growth assays across 4 assays done for each community
gr_agg<-aggregate(gr,list(gr$Community), mean, na.rm=T)
#remove junk columns
gr_agg<-gr_agg[,-c(2:4)]
colnames(gr_agg)[1]<-'Community'
write.csv(gr_agg, "inputs/2_wrangled/growth/growthdata_cleaned_assaymeans.csv", row.names = F)

#calculate mean of invasion assays across 4 assays done for each community
inv_agg<-aggregate(inv,list(inv$Community), mean, na.rm=T)
#remove junk columns
inv_agg<-inv_agg[,-c(2:7)]
#rename Community column
colnames(inv_agg)[1]<-'Community'
write.csv(inv_agg,"inputs/2_wrangled/invasion/invasiondata_luminescence_assaymeans.csv")
