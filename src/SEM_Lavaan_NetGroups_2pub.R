#################################################
#      SEM models with Lavaan                   #
#################################################
# author = Alberto Pascual-García
# web = apascualgarcia.github.io
# date = June 24th, 2020 (ETH-Zürich)
# description = This script takes the invasion data [1], a processed matrix with information
#   regarding the abundances of functional groups found in [2] for the experiments described in [3] 
#   and the productivities of the communities, and it explores with structural equation models the relation between
#   composition (encoded in the functional groups), productivity and invasion.
# edit = The user must define the SEM models in lavaan format in external files located in the folder lavaanModels.
#   Then, the variable selectModel should point to that file. The input data should be located in the respective directories.
#   The paths to the file relative to the repository are provided, should be edited otherwise. Note that the function to
#   retrieve the files works only with rstudio, should be provided otherwise.
# usage =If only a single model is run, the model should be launched from the script: launch_SEM_single.R
#   If multiple models are run, the model should be launched from the script: launch_SEM_multiple.R
# results = The script will run each model and store it in a .mod file, and summary of results of the fit and the MI indexes
#   in another file (.fit). It will then create a representation for the plain model (excluding variances), the structural model
#   and the fulll model. If several models are analysed with the "multiple" version, an analysis of the models will
#   be provided. 
########## 

# Compute model -------------------

myModel <- readLines(paste(pathModel,selectModel,sep=""))
myModel

fit <- lavaan::sem(myModel, data=df.all, check.gradient = FALSE) 
summary(fit, standardized=T,rsq=T, fit.measures = TRUE)
mod.fit=lavaan::modindices(fit,sort. = TRUE)
mod.fit

#,do.fit=FALSE) # Use it when there are problems fitting the model
# inspectCov=as.data.frame(fitted(fit)$cov) # Then use this function to inspect the problem, for instance
# variables very close to zero in the covariance are likely to be problematic, you can
# fix starting values as well to control this.
#,test = "bootstrap",)

# ...... save model
#setwd(pathOut) # fix the output path
outModel=paste(pathOut,selectModel,".mod",sep="") # We save the model
save(myModel,file=outModel)

outFit=paste(pathOut,selectModel,".fit",sep="") # and the fit
sink(outFit) # Saving the summary is tricky, first I divert the output to the file
options(max.print=1000) # I skip the limitation in lines shown
summary(fit,standardized=T,rsq=T, fit.measures = TRUE) # I print the file
print(mod.fit)#[mod.fit$mi>3,]) # Also add modifications proposed
sink() # And close the pipe

# ..... Plot model
# include path to the script
#setwd(pathCode)
source(paste(pathCode,"/SEM_PathPlot.R",sep = ""))
