####################################
# SEM_PAthPlot.R
####################################
# This script contains only the semPlot function, its use
# is just to store different plot configurations, but it should
# be considered the final lines of SEM_Lavaan_latent.R (i.e.
# it should be run that script first). I separated it to store
# different options of the plot depending on the model.
# NOTE: Some more variables may be needed if you run it from
# SEM_Lavaan_NetGroups.R, e.g. nodelabels
# Zürich, May 28th, 2018.
# A.P-G.
######################################

#N.B. THIS SCRIPT AS IT STANDS IS CURRENTLY ONLY EXECUTABLE VIA THE SEM LAUNCHER SCRIPTS IN COMMUNITYINVASION/SRC
#IF YOU WANT TO RUN IT MANUALLY, UNHASH THE FOLLOWING LINES TO SET THE OUTPUT DIRECTORY AND THE MODEL/FIGURE YOU WANT TO MAKE

#pathOut<-"~/GitHub/communityinvasion/sem_analysis/output_data/"
#selectModel<-"CompProdInvasion_PartialMediation_Diversity.lav"

#select the name of the model you want to plot
#Figure 6 (best model, default) = "CompProdInvasion_PartialMediation_Diversity.lav"
#Figure 7 = "CompProdInvasion_CompleteMediation_Diversity.lav"
#Figure 8 = "CompProdInvasion_NoMediation_Diversity.lav"


#setwd(pathOut)
# --- Models:
# Energy2Resp2Yield.lav, with no groups
# EnergyAndResp2Yield.lav
# Plot options:
# base = structural and indicator variables, no residuals
# struct = structural variables only, no residuals
# full = structural and indicator variables, residuals included
plotVec=c("struct","base","full")

for(plotType in plotVec){
  print(paste("Creating plot for options: ",plotType,sep=""))
  if(plotType == "base"){
    structOpt=FALSE
    residOpt=FALSE
    widthOpt=25
    heightOpt=30
  }else if(plotType == "struct"){
    structOpt=TRUE
    residOpt=FALSE
    widthOpt=7
    heightOpt=7
  }else if(plotType == "full"){
    structOpt=FALSE
    residOpt=TRUE
    widthOpt=25
    heightOpt=30
  }
  
  fileSEM=paste(pathOut,"Plot","SEMpath",selectModel,plotType,".pdf",sep="")
  grDevices::pdf(file=fileSEM,width = widthOpt,height = heightOpt)
  #newLabels=c("mgCO2.7"="CO2","CPM7"="Cells","ATP7"="ATP","mG7"="G","mX7"="X","mN7"="N",
  #            "mP7"="P","Respiration"="Resp","Yield"="Yld","Uptake"="Upt")
  
  
  semPlot::semPaths(fit, 
           #layoutSplit = FALSE, # Logical that can be used to split computing of layout between structural and measurment models. This is very useful in more complicated models where the structural part is best shown by using a spring layout.
           #measurementLayout = "tree3", # with layoutSplit=TRUE
           what="est", # What should the edges indicate in the path
           #nodeLabels=newLabels,
           intercepts = FALSE,
           whatLabels = "est", 
           residuals = residOpt,
           structural = structOpt,
           #combineGroups = FALSE,
           #curvature=1.5,
           #centerlLevels=FALSE, # For tree2, should each level be centered?
           #panelGroups = FALSE,
           label.cex=1,
           edge.label.cex = 1.5,
           rotation=4, # exogenous at the top (default=1, at the bottom)
           #optimizeLatRes = TRUE,
           #sizeLat=7,
           #sizeLat2=5
           #springLevels = TRUE,
           #reorder = FALSE,
           #title=TRUE,
           #title.color=colorCodes,
           #maximum=0.1,
           #filetype="pdf",
           #filename=PlotFitName,
           #filename=c(gr1,gr2,gr3,gr4,gr5,gr6), # This doesn't work, I saved them manually
           layout="tree")#tree")#tree3")# "tree3")#"circle2") #)
  
  
  dev.off()
  
}

# # --- Model Milestone 3, separated by groups and manually positioning each variable.
# layout.matrix.coor=matrix(nrow=6,ncol=2) 
# # I obtain the ides looking at the values assigned in the list ggraph.groups (APG. I do not understand how)
# # 1=cells,2=resp,3=energy,4=G,5=X,6=N,7=P (This was for FinalModel_Milestone3, i.e. with latent)
# # 1=cells,2=resp,3=energy,4=X,5=G,6=P,7=N (With manifest it changes, so I reconfigure the matrix below like this)
# newLabels=c("Cells","CO2","ATP","X","G","P","N")
# # The area is defined between -1 and 1 for x and y
# #coor.x=c(-0.5,0.5,0.12,0.25,-0.75,-0.25,0.75)
# #coor.y=c(0.35,0.85,-0.15,-0.75,-1,-0.75,-1) # individual groups
# #coor.y=c(-0.35,-0.125,-0.5,-0.75,-0.9,-0.75,-0.9) # combined groups
# # This is the new reconfiguration for the manifest model
# coor.x=c(-0.5,0.5,0.12,-0.4,0.75,0.4,-0.75)
# coor.y=c(0.35,0.85,-0.15,-1,-0.65,-1,-0.65) # individual groups
# layout.matrix.coor=cbind(coor.x,coor.y)
# 
# #fileSEM="Plot.FinalModel_Manifest-Milestone3.lav.PanelByGroups.pdf"
# fileSEM=paste("Plot","SEMpath",selectModel,"PanelByGroups","pdf",sep=".")
# grDevices::pdf(file=fileSEM,width = 15,height = 6)
# semPlot::semPaths(fit, 
#          #layoutSplit = FALSE, # Logical that can be used to split computing of layout between structural and measurment models. This is very useful in more complicated models where the structural part is best shown by using a spring layout.
#          #measurementLayout = "tree3", # with layoutSplit=TRUE
#          nodeLabels=newLabels,
#          what="std.all", # What should the edges indicate in the path
#          intercepts = FALSE,
#          whatLabels = "std.all", 
#          residuals = FALSE,
#          #structural = FALSE,
#          combineGroups = FALSE,
#          #curvature=1.5,
#          #centerlLevels=FALSE, # For tree2, should each level be centered?
#          panelGroups = TRUE,
#          label.cex=1.7,
#          edge.label.cex = 2.5,
#          #optimizeLatRes = TRUE,
#          sizeMan=12,
#          sizeMan2=8,
#          #springLevels = TRUE,
#          reorder = FALSE,
#          #title=TRUE,
#          #title.color=colorCodes,
#          #maximum=0.1,
#          #filetype="pdf",
#          #filename=PlotFitName,
#          #filename=c(gr1,gr2,gr3,gr4,gr5,gr6), # This doesn't work, I saved them manually
#          layout=layout.matrix.coor)# "tree3")#"circle2") #)
#          #layout="tree2")
# dev.off()
# 
# 
# fileSEM="Plot.FinalModel_Manifest-Milestone1.lav.pdf"
# grDevices::pdf(file=fileSEM,width = 6,height = 12)
# semPlot::semPaths(fit, 
#          #layoutSplit = FALSE, # Logical that can be used to split computing of layout between structural and measurment models. This is very useful in more complicated models where the structural part is best shown by using a spring layout.
#          #measurementLayout = "tree3", # with layoutSplit=TRUE
#          nodeLabels=newLabels,
#          what="std.all", # What should the edges indicate in the path
#          intercepts = FALSE,
#          whatLabels = "std.all", 
#          residuals = FALSE,
#          #structural = FALSE,
#          exoCov=FALSE,
#          combineGroups = FALSE,
#          #curvature=1.5,
#          #centerlLevels=FALSE, # For tree2, should each level be centered?
#          panelGroups = TRUE,
#          label.cex=1.2,
#          edge.label.cex = 1.5,
#          #optimizeLatRes = TRUE,
#          sizeMan=12,
#          sizeMan2=8,
#          #springLevels = TRUE,
#          reorder = FALSE,
#          #title=TRUE,
#          #title.color=colorCodes,
#          #maximum=0.1,
#          #filetype="pdf",
#          #filename=PlotFitName,
#          #filename=c(gr1,gr2,gr3,gr4,gr5,gr6), # This doesn't work, I saved them manually
#          layout=layout.matrix.coor)# "tree3")#"circle2") #)
# #layout="tree2")
# dev.off()
# 
# 
