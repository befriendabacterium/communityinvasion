# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)

#if rstudioapi package not installed, install. Needed for setting working directory
if('rstudioapi'%in%rownames(installed.packages())==F){
  install.packages('rstudioapi')
}

#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")

# RUN CODE ----------------------------------------------------------------

#analysis
source('src/0_acquirepackages.R')
source('src/0_downloaddata.R')
source('src/1_cleanexplanatoryvariables.R')
source('src/2_aggregate_assays.R')
source('src/3_matchdatasets.R')
source('src/4_convertinvasionlumtocfu.R')
source('src/5_calculatedimensionalityreductions.R')
source('src/6_find_bestdimreduction_rerunranfs.R')
source('src/7_calculatediversitymetrics.R')
source('src/8_randomforests.R')
source('src/9_launch_SEM_multiple_simpsonsdiversity.R')
source('src/10_coeffs_extractor.R')

#figures
source('src/Figure1.R')
source('src/Figure2.R')
source('src/Figure3.R')
source('src/Figure4.R')

#supplementary figures
source('src/Figure2_Supplement1.R')
source('src/Figure2_Supplement2.R')
source('src/Figure2_Supplement3.R')
source('src/Figure3_Supplement1.R')
source('src/Figure3_Supplement2.R')
source('src/Figure3_Supplement3.R')
source('src/Figure4_Supplements1to3.R')

#sensitivity analysis - takes a long time hence hashed out and at end
source('src/X_sensitivityanalysis.R')                         