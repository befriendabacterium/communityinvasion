rm(list=ls())
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
this.dir<-getwd()

source('src/9_launch_SEM_multiple_simpsonsdiversity.R')

# MAKE MODEL FIT TABLE ----------------------------------------------------

modnames=c("No mediation","Partial","Complete")
fittable<-c()

for (m in 1:length(model.fit.list)){
  fit_temp<-c(model=modnames[m],lavaan::fitMeasures(model.fit.list[[m]], c('df','cfi')))
  summ<-lavaan::summary(model.fit.list[[m]], rsq=T)$PE
  #extract r squared
  rsq_temp<-summ[nrow(summ),]$est
  fit_temp<-c(fit_temp,rsq_inv=rsq_temp)
  fittable<-rbind(fittable,fit_temp)

}

fittable<-as.data.frame(fittable, row.names = F)
fittable

aiccols<-which(colnames(aic.rank)%in%c("AICc","Delta_AICc"))
aic.rank<-as.data.frame(aic.rank)
fittable<-cbind(fittable,aic.rank[,aiccols])
fittable

write.csv(fittable, 'outputs/sem/modelcomparison.csv')

# GET MODEL COEFFICIENTS ------------------------------------------

pathdescs_needed<-c("Invasion ~ Composition","Productivity ~ Composition","Invasion ~ Productivity")

#NO MEDIATION MODEL
nomediationmodel_ests<-lavaan::parameterEstimates(model.fit.list[[1]], standardized = T)

path_descs_all<-paste(nomediationmodel_ests$lhs,nomediationmodel_ests$op,nomediationmodel_ests$rhs)  

#get paths interested in
estimates_table<-nomediationmodel_ests[path_descs_all%in%pathdescs_needed,]
estimates_table<-estimates_table[c(1,3,2),]
nomediationmodel_stdcoeffsall<-estimates_table$std.all

nomediationmodel_stdcoeffsall_df<-
  cbind("No mediation",
        c('A','B','C'),
        c('Invasion ~ Productivity','Productivity ~ Composition', 'Invasion ~ Composition'),
        nomediationmodel_stdcoeffsall
  )

nomediationmodel_stdcoeffsall_df

colnames(nomediationmodel_stdcoeffsall_df)<-c('Model','Path','Interpretation','Std.coefficient')
nomediationmodel_stdcoeffsall_df<-as.data.frame(nomediationmodel_stdcoeffsall_df)
nomediationmodel_stdcoeffsall_df

#PARTIAL MEDIATION MODEL
partialmodel_ests<-lavaan::parameterEstimates(model.fit.list[[2]], standardized = T)

path_descs_all<-paste(partialmodel_ests$lhs,partialmodel_ests$op,partialmodel_ests$rhs)  
pathdescs_needed<-c("Invasion ~ Composition","Productivity ~ Composition","Invasion ~ Productivity")

#get paths interested in
estimates_table<-partialmodel_ests[path_descs_all%in%pathdescs_needed,]
estimates_table<-estimates_table[c(1,3,2),]
#composition-mediated effects of productivity
D<-estimates_table$std.all[1]*estimates_table$std.all[2]
#productivity-only effects
E<-estimates_table$std.all[1]-D
partialmodel_stdcoeffsall<-c(estimates_table$std.all,D,E)
partialmodel_stdcoeffsall_df<-
  cbind("Partial mediation",
          c('A','B','C','D = A * B','E = A - D'),
        c('Invasion ~ Productivity','Productivity ~ Composition','Invasion ~ Composition','Invasion ~ Productivity | Composition','Productivity-only effect on invasion'),
        partialmodel_stdcoeffsall
  )

colnames(partialmodel_stdcoeffsall_df)<-c('Model','Path','Interpretation','Std.coefficient')
partialmodel_stdcoeffsall_df<-as.data.frame(partialmodel_stdcoeffsall_df)
partialmodel_stdcoeffsall_df

#COMPLETE MEDIAION MODEL
completemodel_ests<-lavaan::parameterEstimates(model.fit.list[[3]], standardized = T)

path_descs_all<-paste(completemodel_ests$lhs,completemodel_ests$op,completemodel_ests$rhs)  
pathdescs_needed<-c("Invasion ~ Composition","Productivity ~ Composition","Invasion ~ Productivity")

#get paths interested in
estimates_table<-completemodel_ests[path_descs_all%in%pathdescs_needed,]
estimates_table[3,]<-NA
#composition-mediated effects of productivity
D<-estimates_table$std.all[1]*estimates_table$std.all[2]
#productivity-only effects
E<-estimates_table$std.all[1]-D

completemodel_stdcoeffsall<-c(estimates_table$std.all,D,E)
completemodel_stdcoeffsall_df<-
  cbind("Complete mediation",
        c('A','B','C','D = A * B','E = A - D'),
        c('Invasion ~ Productivity','Productivity ~ Composition','Invasion ~ Composition','Invasion ~ Productivity | Composition','Productivity-only effect on invasion'),
        completemodel_stdcoeffsall
  )

colnames(completemodel_stdcoeffsall_df)<-c('Model','Path','Interpretation','Std.coefficient')
completemodel_stdcoeffsall_df<-as.data.frame(completemodel_stdcoeffsall_df)
completemodel_stdcoeffsall_df



all_model_coeffs<-
  rbind(
  nomediationmodel_stdcoeffsall_df,
  partialmodel_stdcoeffsall_df,
  completemodel_stdcoeffsall_df)

all_model_coeffs

write.csv(all_model_coeffs,'outputs/sem/all_model_coeffs.csv')
