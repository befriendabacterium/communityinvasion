# This is auxiliary code to load data for the script 
# SEM_Lavaan_NetGroups_2pub.R

# Read Files -------------------------------------------------------------
# --- Paths are automatically set relative to the root of the repo
pathModel=paste(this.dir,"/inputs/4_semmodelspecifications/",sep="")
pathCode=paste(this.dir,"/src/",sep="")
pathOut=paste(this.dir,"/outputs/sem/",sep="")

#-- Read a single object

source('src/X_data_prepper.R')

# --- Modify some column names to facilitate model specification

idx=which(colnames(df.all)=="community") # make life easier changing the colnames
colnames(df.all)[idx]="Community"
idx=which(colnames(df.all)=="kt2440.cfu.24h") # make life easier changing the colnames
colnames(df.all)[idx]="Putida24h"
idx=which(colnames(df.all)=="kt2440.cfu.96h") # make life easier changing the colnames
colnames(df.all)[idx]="Putida96h"
idx=which(colnames(df.all)=="kt2440.cfu.7d") # make life easier changing the colnames
colnames(df.all)[idx]="Putida168h"
idx=which(colnames(df.all)=="sbw25.cfu.24h") # make life easier changing the colnames
colnames(df.all)[idx]="Fluoresc24h"
idx=which(colnames(df.all)=="sbw25.cfu.96h") # make life easier changing the colnames
colnames(df.all)[idx]="Fluoresc96h"
idx=which(colnames(df.all)=="sbw25.cfu.7d") # make life easier changing the colnames
colnames(df.all)[idx]="Fluoresc168h"

# Rescale data ---------------------------

Nvars=dim(df.all)[2]
all.vars=matrix(0,nrow=1,ncol=Nvars)
sapply(df.all,class)
numeric=which(sapply(df.all,class)=="numeric")
checkvars<-sapply(df.all[,numeric],var, na.rm=T) # check the variances
checkvars
select.cols=colnames(df.all)[numeric]

# In ill.scaled covariance matrix, the ratio of the largest to smallest variance 
# is greater than say, 100.0 (see e.g. book Rex Kline, pp 81)
min(checkvars[1:6])/min(checkvars[-c(1:6)])
mean(checkvars[1:6])/min(checkvars[-c(1:6)])
max(checkvars[1:6])/min(checkvars[-c(1:6)])

# The invasion variables variances are only ~40x larger (average 17x, min 2.5x) larger than the functional groups, so no need to rescale further
