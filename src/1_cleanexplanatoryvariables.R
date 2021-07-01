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

# COMPOSITION -------------------------------------------------------------

#read in raw OTU table
comp<-read.csv('inputs/1_raw/composition/Functional_Redundancy_Treehole_genera.csv', row.names = 1)

#remove OTUs with less than 100 reads across all samples
comp<-comp[,-which(colSums(comp)<100)]
#remove OTUs occurring in less than 10 samples
comp<-comp[,-which(colSums(comp>0)<10)]
#remove communities with less than 10,000 sequences
comp<-comp[which(rowSums(comp[,-1])>10000),]

#add community column name
comp<-cbind(Community=rownames(comp),comp)

#write the cleaned data to a CSV
write.csv(comp, 'inputs/2_wrangled/composition/composition_otu_cleaned.csv', row.names = F)

# GROWTH ---------------------------------------------------------

#read in raw growth assay data (productivity, metabolic and enzyme measurements)
growthmeasures<-read.csv('inputs/1_raw/growth/20151016_Functions_remainder.csv')

#ditch ATP at 1 and 4 days because don't have these for other growth measures (redundant)
growthmeasures<-growthmeasures[,-which(colnames(growthmeasures)%in%c("ATP1","ATP4"))]

#ditch the respiration per cell columns (redundant)
growthmeasures<-growthmeasures[,-which(colnames(growthmeasures)%in%c("pgRPC.7","pgRPC.14"))]

#write the cleaned data to a CSV
write.csv(growthmeasures,'inputs/2_wrangled/growth/growthdata_cleaned.csv', row.names = F)

