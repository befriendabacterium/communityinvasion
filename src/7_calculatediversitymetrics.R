# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

# LOAD PACKAGES --------------------------------------------------

#install packages
#install.packages('ape')
#install.packages('ade4')
#install.packages('caret')
#install.packages('ggsci')
#install.packages('picante')
#install.packages('plyr')
#install.packages('RColorBrewer')
#install.packages('randomForest')
#install.packages('scales')
#install.packages("vegan")

#load libraries
library(ape)
library(ade4)
library(caret)
library(ggsci)
library(picante)
library(plyr)
library(randomForest)
library(RColorBrewer)
library(scales)
library(vegan)

# LOAD DATA ---------------------------------------------------------------

#load composition
composition<-read.csv('inputs/3_ready/composition/composition_otu_matched.csv', row.names = 1, header=T)

#load tree of life
tree<-ape::read.tree('inputs/external/geneious/tree101_rooted_invdistNAMES.newick')
#select first tree (best fit)
tree<-tree[[1]]
#root tree in H.salinarum
tree<-root(tree, outgroup = 583, resolve.root = T) #roots tree and resolves root cos must have two branches off it
#checks rooted
is.rooted(tree)
#load taxonomy table for OTUs
taxonomy<-read.csv('inputs/external/taxonomy/sp.names.final.csv', stringsAsFactors = T)

#check taxonomy species names are in otu table, in same order
count(taxonomy$species==colnames(composition))
#check taxonomy OTU names are in tree table, in any order
count(taxonomy$otu%in%tree$tip.label)
#rename OTUtable cols with OTU numbers so match those in tree (for sake of analysis)
colnames(composition)<-taxonomy$otu

#load invasion data
inv<-read.csv('inputs/3_ready/invasion/invasiondata_cfu_assaymeans_matched.csv')

# CALCULATE SIMPSON'S DIVERSITY -------------------------------------------

simp<-vegan::diversity(composition,index='simpson')
plot(log10(inv$sbw25.cfu.24h+1)~simp)
model<-lm(log10(inv$sbw25.cfu.24h+1)~simp)
summary(model)
abline(model)

# PHYLOGENETIC DIVERSITY AT THE ORDER LEVEL (where there is coverage in tree of life) ------------------------------------------------

#calculate phylogenetic distance matrix from tree of life
phydist<-cophenetic(tree)
#make dataframe
phydist.df<-as.data.frame(phydist)

#CALCULATE PHYLOGENETIC DIVERSITY
phylo_simp<-SYNCSA::rao.diversity(composition, phylodist = phydist, standardize = T)
phylo_simp$PhyRao

#CALCULATE PHYLOGENETIC DISTANCE BETWEEN INVADER AND SPECIES IN COMMUNITY

#make distance matrix of invaders to comm members only
invdists_nonum<-as.data.frame(cbind(phydist.df$SBW25_trimmed, phydist.df$KT2440_trimmed))
#set colnames to invader names
colnames(invdists_nonum)<-c('SBW25','KT2440')
#set rownames to same as phylogenetic distance matrix
rownames(invdists_nonum)<-colnames(phydist.df)
#remove two invader and root distances
invdists_nonum<-invdists_nonum[!rownames(invdists_nonum)%in%c('KT2440_trimmed','SBW25_trimmed','H.salinarum'),] 
#reorder invdists_nonum to order of OTUcols
invdists_nonum<-invdists_nonum[match(colnames(composition),rownames(invdists_nonum)),]
#rename invdists_nonum with species names of OTUS
rownames(invdists_nonum)<-taxonomy$species
#check if p.fluorescrens and putida are in the phylogenetic distance vectors for each invader
invdists_nonum[which(rownames(invdists_nonum)=='pseudomonas.fluorescens'),]
invdists_nonum[which(rownames(invdists_nonum)=='pseudomonas.putida'),]

#calculate overall phylogenetic distance of community by getting mean branch length between invader and community members (weighted by OTU abundance) per community
community.sbw25.wdists_nonum<-c()
for (r in 1:nrow(composition)){
  #For each community, calculated mean phylogenetic distance between community and SBW25, weighted by OTU abundance
  temp.rowsum<-weighted.mean(invdists_nonum$SBW25, composition[r,])
  community.sbw25.wdists_nonum<-c(community.sbw25.wdists_nonum,temp.rowsum)
}

#calculate overall phylogenetic distance of community by getting mean branch length between invader and community members (weighted by OTU abundance) per community
community.kt2440.wdists_nonum<-c()
for (r in 1:nrow(composition)){
  #For each community, calculated mean phylogenetic distance between community and SBW25, weighted by OTU abundance
  temp.rowsum<-weighted.mean(invdists_nonum$KT2440, composition[r,])
  community.kt2440.wdists_nonum<-c(community.kt2440.wdists_nonum,temp.rowsum)
}

plot(log10(inv$sbw25.cfu.24h+1)~log10(community.sbw25.wdists_nonum))
model<-lm(log10(inv$sbw25.cfu.24h+1)~log10(community.sbw25.wdists_nonum))
summary(model)
abline(model)

plot(log10(inv$sbw25.cfu.24h+1)~log10(community.kt2440.wdists_nonum))
model<-lm(log10(inv$sbw25.cfu.24h+1)~log10(community.kt2440.wdists_nonum))
summary(model)
abline(model)

# WRITE DIVERSITY METRICS -------------------------------------------

diversitymetrics<-
  
  cbind(community = rownames(composition),
        simp = simp,
        rao = phylo_simp$PhyRao,
        sbw25_phydist = community.sbw25.wdists_nonum,
        kt2440_phydist  = community.kt2440.wdists_nonum
  )

write.csv(diversitymetrics, 'inputs/3_ready/composition/diversitymetrics.csv', row.names = F)
