# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")

# EMPTY EXISTING DATA/REVERT TO GITHUB REPO STATE -----------------------------------------------------

#empty the inputs and outputs directories of any data, keeping all READMEs (that's what the grep does) and jpgs
file.remove(
  grep(list.files(c('inputs','outputs'), recursive = T, full.names = T), pattern='.md|.jpg', invert=T, value=T)
  )

# DOWNLOAD NEW DATA -------------------------------------------------------

#change this to '2_preanalysis' to run the code from Step 8, to '3_end' to download the end result of running the code
whichpoint<-'1_start'

my_project <- osfr::osf_ls_files(osfr::osf_retrieve_node("hc57w"))
data_folder <- my_project[which(my_project$name==whichpoint),]

#download inputs and outputs folders
osfr::osf_download(osfr::osf_ls_files(data_folder),getwd(), recurse = T, conflicts='overwrite')
