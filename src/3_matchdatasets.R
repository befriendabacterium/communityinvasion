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

#OTU table
composition_otu_cleaned<-read.csv('inputs/2_wrangled/composition/composition_otu_cleaned.csv')
#raw and cleaned invasion data
gr<-read.csv('inputs/2_wrangled/growth/growthdata_cleaned_assaymeans.csv')
inv<-read.csv('inputs/2_wrangled/invasion/invasiondata_luminescence_assaymeans.csv')
  
# AGGREGATE INVASION ASSAYS ---------------------------------------

# FIND THE COMMUNITIES THAT MATCH ACROSS DATASETS AND SUBSET DATA --------

#find the communities that match across the 3 comp analyses,1 growth assays, and 1 invasion assays
matchingcomms<-Reduce(intersect, list(composition_otu_cleaned$Community, gr$Community, inv$Community)) 
#subset the 5 dataframes to those that do match 
composition_otu_matched<-composition_otu_cleaned[composition_otu_cleaned$Community%in%matchingcomms,]

gr<-gr[gr$Community%in%matchingcomms,]
inv<-inv[inv$Community%in%matchingcomms,]

write.csv(composition_otu_matched, 'inputs/3_ready/composition/composition_otu_matched.csv', row.names = F)
write.csv(gr,'inputs/3_ready/growth/growthdata_assaymeans_matched.csv', row.names = F)
write.csv(inv,'inputs/3_ready/invasion/invasiondata_luminescence_assaymeans_matched.csv', row.names = F)

