# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

# LOAD DATA --------------------------------------------------------------------

composition<-read.csv('inputs/3_ready/composition/composition_otu_matched.csv', row.names = 1)
taxonomy_df<-read.csv("inputs/external/taxonomy/taxonomyNCBI.csv")
taxonomy_df

# LOAD FUNCTIONS ----------------------------------------------------------

# Function to compute the Shannon-Jensen distance 
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {    
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

# CALCULATE ABUNDANCES AT EACH PHYLOGENETIC LEVEL -------------------------

#transpose OTU table to allow aggregation
composition_transposed<-t(composition)

#genus
composition_genus<-t(aggregate(composition_transposed, list(taxonomy_df$genus),sum))
#set first row as col names
colnames(composition_genus)<-composition_genus[1,]
#remove first row and coerce to df
composition_genus<-composition_genus[-1,]
#make a Community id column
composition_genus<-cbind(Community=rownames(composition_genus),composition_genus)
write.csv(composition_genus,"inputs/3_ready/composition/composition_genus_matched.csv", row.names=F)

#family
composition_family<-t(aggregate(composition_transposed, list(taxonomy_df$family),sum))
#set first row as col names
colnames(composition_family)<-composition_family[1,]
#remove first row
composition_family<-composition_family[-1,]
#make a Community id column
composition_family<-cbind(Community=rownames(composition_family),composition_family)
write.csv(composition_family,"inputs/3_ready/composition/composition_family_matched.csv", row.names=F)

#order
composition_order<-t(aggregate(composition_transposed, list(taxonomy_df$order),sum))
#set first row as col names
colnames(composition_order)<-composition_order[1,]
#remove first row
composition_order<-composition_order[-1,]
#make a Community id column
composition_order<-cbind(Community=rownames(composition_order),composition_order)
write.csv(composition_order,"inputs/3_ready/composition/composition_order_matched.csv", row.names=F)

#class
composition_class<-t(aggregate(composition_transposed, list(taxonomy_df$class),sum))
#set first row as col names
colnames(composition_class)<-composition_class[1,]
#remove first row
composition_class<-composition_class[-1,]
#make a Community id column
composition_class<-cbind(Community=rownames(composition_class),composition_class)
write.csv(composition_class,"inputs/3_ready/composition/composition_class_matched.csv", row.names=F)

#phylum
composition_phylum<-t(aggregate(composition_transposed, list(taxonomy_df$phylum),sum))
#set first row as col names
colnames(composition_phylum)<-composition_phylum[1,]
#remove first row
composition_phylum<-composition_phylum[-1,]
#make a Community id column
composition_phylum<-cbind(Community=rownames(composition_phylum),composition_phylum)
write.csv(composition_phylum,"inputs/3_ready/composition/composition_phylum_matched.csv", row.names=F)

# FETCH FUNCTIONLINK CALCULATIONS AND MOVE TO 3_READY FOLDER ------------------------------------------------------

funcgroups<-read.table('inputs/external/functionInk/SamplePropsOFtaxaClus_Time0_NL_Average_StopStep-91_AbundMean.dat', header=T)
colnames(funcgroups)[colnames(funcgroups)=='Sample']<-'Community'
funcgroups<-funcgroups[match(rownames(composition),funcgroups$Community),]

write.csv(funcgroups,'inputs/3_ready/composition/composition_funcgroups_matched.csv', row.names = F)

# CALCULATE PCOA ----------------------------------------------------------

#calculate jensen shannon distance
jensenshannon_dist<-dist.JSD(t(composition))
#calculate the pcoa
obs.pcoa<-ade4::dudi.pco(jensenshannon_dist, scannf=F, nf=17)

saveRDS(obs.pcoa, 'inputs/3_ready/composition/composition_pcoa.RDS')

