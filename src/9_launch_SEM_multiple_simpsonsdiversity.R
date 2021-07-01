#################################################
#      SEM models with Lavaan                   #
#################################################
# author = Alberto Pascual-García
# web = apascualgarcia.github.io
# date = June 24th, 2020 (ETH-Zürich)
# description = This script is a wrapper of SEM_Lavaan_NetGroups_2pub.R and
#     it permits to run it several times with different models. However, if you
#     are interested in comparing several models, since the anova
#     function below to analyse different SEM has a rigid sintax, the script cannot be used without recoding the call of
#     this function. 
########## 
#rm(list=ls())
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
this.dir<-getwd()

#-- Load packages needed
#library(gplots) 
#library(MASS)
#require(lavaan)
#library(rstudioapi)    

#### START EDITING
# Models and files --------------------------------------------------------
# --- Provide the files for your models in a list
# Basic model with diversity 
models=c("CompProdInvasion_NoMediation_Diversity.lav",
         "CompProdInvasion_PartialMediation_Diversity.lav",
         "CompProdInvasion_CompleteMediation_Diversity.lav")

label_out="MediationTest" # a label to identify your results

# Optimized model
#models=c("CompProdInvasion_NoMediation_mod6.2.10.lav",
#         "CompProdInvasion_PartialMediation_mod6.2.10.lav",
#         "CompProdInvasion_CompleteMediation_mod6.2.10.lav")

# .... The output files saved with the model are created according to "selectModel" variable)

# ... Comparison of models
comparison=1 # if you want to compare models (=1), =0 otherwise. Note that this part of the code should be recoded (see below)

###### STOP EDITING

model.fit.list=list()
k=0
for(selectModel in models){
  k=k+1
  cat("** Running model ",selectModel,"\n")
  if(k == 1){
    source("src/SEM_Lavaan_process_input.R")
  }
    source("src/SEM_Lavaan_NetGroups_2pub.R")
    model.fit.list[[k]]=fit
}

if(comparison == 1){
  # Here you should choose the models you want to compare
  # ..... For nested models
  # 1="Independent",2="NoMediation",3="Partial",4="Complete"
  # df(Partial) < df(Complete) = df(Independent) < df(NoMed)
  anova.test=lavaan::lavTestLRT(model.fit.list[[2]],model.fit.list[[3]]) # Partial and complete
  
  
  # ..... Rank by parsimony (AIC), we expect a deltaAIC > 2 per degree of freedom lost
  library(AICcmodavg) 
  aic.rank=AICcmodavg::aictab(model.fit.list,
                  modnames=c("NoMediation","Partial","Complete"))
  
}

orderneeded<-match(c('NoMediation','Partial','Complete'),aic.rank$Modnames)
aic.rank<-aic.rank[orderneeded,]
aic.rank

setwd(pathOut)
fileOut=paste("MultipleTests_results_",label_out,".txt",sep="")
sink(fileOut)
cat("******** \n")
cat("Anova tests \n")
cat("1=NoMediation,2=Partial,3=Complete \n")
print(anova.test)
cat(" \n")
cat("******** \n")
cat("Rank by parsimony (AIC), we expect a deltaAIC > 2 per degree of freedom lost \n")
print(aic.rank)
sink()

#move back to starting directory
setwd("../..")
getwd()
