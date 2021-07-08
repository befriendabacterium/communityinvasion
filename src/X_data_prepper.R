# READ IN DATA -----------------------------------------------------

if(exists("sensitivityanalysis")==T){
  #imputed data
  print('data below detection limit inferred')
  invexp.data<-invexp.data_temp
} else {
  #original data
  print('original CFU data')
  invexp.data<-read.csv('inputs/3_ready/invasion/invasiondata_cfu_assaymeans_matched.csv', row.names = 1)
  invexp.data<-log10(invexp.data+1)
}

#diversity metrics
diversitymetrics<-read.csv('inputs/3_ready/composition/diversitymetrics.csv', row.names = 1)

#OTU abundances
#ALL OTUS
composition<-read.csv('inputs/3_ready/composition/composition_funcgroups_matched.csv', header=T, row.names=1)
composition$Total<-log10(composition$Total)
#match composition to taxa
composition<-composition[match(rownames(invexp.data),rownames(composition),),]
#write.csv(composition,'inputs/3_ready/composition/funcgroups/meanabundance_of_functionalgroups_matched.csv', row.names = F)

#growth metrics
growthmeasures<-read.csv('inputs/3_ready/growth/growthdata_assaymeans_matched.csv', row.names = 1)

# LOG AND COMBINE RELEVANT DIVERSITY METRICS AND GROWTH MEASURES INTO THREE DATAFRAMES (GENERAL, SBW25-SPECIFIC, KT2440-SPECIFIC) ---------------------------------------------------------

#identify column indexes NOT to log
dontlog<-which(colnames(diversitymetrics)%in%c('community','simp','rao'))
diversitymetrics[,-dontlog]<-log10(diversitymetrics[,-dontlog])
dontlog<-which(colnames(growthmeasures)%in%c('Community','mgCO2.7','mgCO2.14'))
growthmeasures[,-dontlog]<-log10(growthmeasures[,-dontlog])

#combine dataframes appropriately
divandgrowth<-cbind(diversitymetrics,growthmeasures)

# MAKE DATAFRAME FOR SEM WITH SELECTED VARIABLES (doesn't contain phydist) --------------------------------------------------

df.all<-cbind(invexp.data,composition,divandgrowth)
write.csv(divandgrowth,'inputs/3_ready/divandgrowth_data.csv')

#plot(invexp.data$sbw25.cfu.24h~sumabundance_of_functionalgroups$C4)
#model<-lm(invexp.data$sbw25.cfu.24h~sumabundance_of_functionalgroups$C4)
#summary(model)
#abline(model)
