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
library(plotrix) #for std.error function
library(plyr)
library(randomForest)
library(RColorBrewer)
library(scales)
library(vegan)

# LOAD FUNCTIONS ----------------------------------------------------------

#standard error function

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))

addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}

# LOAD DATA ---------------------------------------------------------------

#load data of invaders growing in monoculture
invmono.data<-read.csv("inputs/1_raw/invasion/invadercalibrationgrowthcurve_96h_8reps.csv")
#load all invasion data and metadata
invexp.data<-read.csv("inputs/3_ready/invasion/invasiondata_cfu_assaymeans_matched.csv", row.names = 1)
detectionthresholds<-readRDS('inputs/1_raw/invasion/detectionlims.RDS')

# DATA PREP ---------------------------------------------------------------

#ESTIMATE INOCULUM MEANS AND SES
#estimate the inoculum means by getting densities (per ml) per invader at 96 hours, dividing by 25 (because only 40ul inoculum of a 96 hr culture added per ml) and getting mean
inocula_mean<-aggregate(invmono.data$cfu.96/25, by=list(invmono.data$invader), mean)
#append standard error to this dataframe
inocula_mean[,3]<-aggregate(invmono.data$cfu.96/25, by=list(invmono.data$invader), std.error)[,2]
#drop first column
inocula_mean<-inocula_mean[,-1]
#invaders as row names
rownames(inocula_mean)<-c('KT2440','SBW25')
#mean/se as col names
colnames(inocula_mean)<-c('Mean','SE')
#convert to scientific notation for purposes of putting in text
sapply(inocula_mean[,-1], formatC, format = "e", digits = 2)

#CALCULATE MONOCULTURE MEANS AND SES
mono_mean<-aggregate(invmono.data[,c('cfu.24','cfu.96','cfu.7d')], by=list(invmono.data$invader), mean)
mono_mean<-mono_mean[,-1]
rownames(mono_mean)<-c('KT2440','SBW25')
mono_se<-aggregate(invmono.data[,c('cfu.24','cfu.96','cfu.7d')], by=list(invmono.data$invader), std.error)
#convert to scientific notation for purposes of putting in text (standard error stays same format)
as.data.frame(sapply(mono_mean[,-1], formatC, format = "e", digits = 2))
mono_se


#CALCULATE COMMUNITY MEANS AND SES
comm_mean<-rbind(colMeans(invexp.data[,4:6]+1),colMeans(invexp.data[,1:3]+1))
colnames(comm_mean)<-c(24,96,168)
rownames(comm_mean)<-c("KT2440", "SBW25")
formatC(comm_mean, format = "e", digits = 2)
comm_se<-rbind(apply(invexp.data[,4:6]+1,2,std.error),apply(invexp.data[,1:3]+1,2,std.error))
rownames(comm_se)<-c("KT2440", "SBW25")
comm_se

#DETECTION THRESHOLDS
sbw.detectionthresh<-detectionthresholds[1]
formatC(sbw.detectionthresh, format = "e", digits = 2)
kt.detectionthresh<-detectionthresholds[2]
formatC(kt.detectionthresh, format = "e", digits = 2)
#count numbers below detection threshold
apply(invexp.data[,1:3]>sbw.detectionthresh,2,plyr::count)
apply(invexp.data[,4:6]<kt.detectionthresh,2,plyr::count)

# PLOT OPTIONS ---------------------------------------------------------------

#dev.off()
tiff('outputs/figures/Figure1.tiff', res=300, units='in', width=12, height=10)

#turn off plot options
#layout
layout(matrix(c(1,1,1,1,
                1,1,1,1,
                1,1,1,1,
                2,2,3,3,
                2,2,3,3,
                2,2,3,3,
                4,4,4,4),
                7, 4, byrow = T))
par(mar=c(0, 6, 2, 0))

# FIGURE 1A: SCHEMATIC ---------------------------------------------------------------

#load images
schematic<-imager::load.image('inputs/external/images/schematic.jpg')
plot(schematic, axes=F)
mtext('A', side=2, las=2, padj=-3, line=1, cex=4)

# FIGURE 1B: Relationship between community composition, overall growth performance and invasion resistance (FIGURE 3, LEFT) -----------------------------------------------------------------------

par(mar=c(0, 7, 2, 1))

#make timepoints vector
tps<-c(24,96,168)
#make y axis tick marks (log10) for invasion success (cells/ml)
log.ticks.y<-sapply(c(1:8), function(i) as.expression(bquote(10^ .(i))))
#log all the inocula, monoculture and community means
inocula_mean<-log10(inocula_mean+1)
mono_mean<-log10(mono_mean+1)
comm_mean<-log10(comm_mean+1)

#plot
plot(as.numeric(mono_mean[2,])~tps, xlim=c(0,180), ylim=c(1,9), yaxt='n', xaxt='n',
     xlab='',ylab='',
     pch=22, cex=2,bg='white', col="black", type='n')

mtext('B', side=2, las=2, padj=-4, line=1, cex=4)

text(40,8.75,expression(italic('P. fluorescens')), cex=3)
#x axis tick marks
axis(1, at=c(0,24,96,168),c(0,24,96,168),las=1, cex.axis=1.5, padj=1)
#y axis tick marks
axis(2, at=c(1,2,3,4,5,6,7,8),log.ticks.y,las=2, cex.axis=1.3)
mtext('Invader survival (cells/ml)', side=2, line=4, cex=1.5)

#add monoculture points, jittered
points(log10(c(invmono.data[invmono.data$invader=="SBW25",c('cfu.24','cfu.96','cfu.7d')][,1],
               invmono.data[invmono.data$invader=="SBW25",c('cfu.24','cfu.96','cfu.7d')][,2],
               invmono.data[invmono.data$invader=="SBW25",c('cfu.24','cfu.96','cfu.7d')][,3]))~jitter(c(rep(24,8),rep(96,8),rep(168,8)),0.5), 
       xlim=c(0,170), ylim=c(3,7), pch=18, cex=1, col=alpha("black",0.2))

points(as.numeric(mono_mean[2,])~tps,
       pch=23, cex=2, bg='white', col="black")

#draw lines betwen 24,96 and 168 hours for monoculture
lines(as.numeric(mono_mean[2,])~tps, xlim=c(0,170), col='black', lwd=1.5)

#add community data, coloured by partition
points(log10(c(invexp.data[,1],invexp.data[,2],invexp.data[,3])+1)~jitter(c(rep(24,nrow(invexp.data)),rep(96,nrow(invexp.data)),rep(168,nrow(invexp.data))),0.5), xlim=c(0,170), ylim=c(3,7),
       pch=16, cex=1, col=alpha("black",0.2))
#add community timepoint means as black points
points(as.numeric(comm_mean[2,])~tps,pch=21, cex=2, bg='white',col='black')
#add lines between community timepoint means
lines(comm_mean[2,]~tps, lwd=1.5, lty=1)

#add points and dotted line between inoculum mean and monoculture at 96hr
points(c(inocula_mean[2,1],mono_mean[2,1])~c(0,24), pch=4, cex=1,col=c("black",NA))
lines(c(inocula_mean[2,1],mono_mean[2,1])~c(0,24),col="black", lwd=1.5,lty=2)
#add points and dotted line between inoculum mean and community at 96hr
points(c(inocula_mean[2,1],comm_mean[2,1])~c(0,24), pch=4, cex=2,col=c("black",NA))
lines(c(inocula_mean[2,1],comm_mean[2,1])~c(0,24),col="black", lwd=1.5, lty=2)
#add line at detection threshold
abline(h=log10(sbw.detectionthresh), lwd=2, col=alpha("black",0.2))
#label numbers below detection threshold
belowdet_n<-apply(invexp.data[,1:3]<sbw.detectionthresh,2,sum)
text(96,2,"Numbers below detection limit (n = 678):",col='black', cex=1.25)
text(24,1.5,belowdet_n[1], col='black', cex=1.25)
text(96,1.5,belowdet_n[2], col='black', cex=1.25)
text(168,1.5,belowdet_n[3], col='black', cex=1.25)

# FIGURE 1C: Relationship between community composition, overall growth performance and invasion resistance (FIGURE 3, RIGHT) --------

par(mar=c(0, 2, 2, 6))

#plot
plot(as.numeric(mono_mean[2,])~tps, xlim=c(0,180), ylim=c(1,9), yaxt='n', xaxt='n',
     xlab='',ylab='',
     pch=22, cex=2,bg='white', col="black", type='n')
text(27,8.75,expression(italic('P. putida')), cex=3)
#x axis tick marks
axis(1, at=c(0,24,96,168),c(0,24,96,168),las=1, cex.axis=1.5, padj=1)

#y axis tick marks
axis(2, at=c(1,2,3,4,5,6,7,8),NA,las=2, cex.axis=1.3)

#add monoculture points, jittered
points(log10(c(invmono.data[invmono.data$invader=="KT2440",c('cfu.24','cfu.96','cfu.7d')][,1],
               invmono.data[invmono.data$invader=="KT2440",c('cfu.24','cfu.96','cfu.7d')][,2],
               invmono.data[invmono.data$invader=="KT2440",c('cfu.24','cfu.96','cfu.7d')][,3]))~jitter(c(rep(24,8),rep(96,8),rep(168,8)),0.5), 
       xlim=c(0,170), ylim=c(3,7), pch=18, cex=1, col=alpha("black",0.2))

points(as.numeric(mono_mean[2,])~tps,
       pch=23, cex=2, bg='white', col="black")

#draw lines betwen 24,96 and 168 hours for monoculture
lines(as.numeric(mono_mean[2,])~tps, xlim=c(0,170), col='black', lwd=1.5)

#add community data, coloured by partition
points(log10(c(invexp.data[,4:6][,1]+1,invexp.data[,4:6][,2]+1,invexp.data[,4:6][,3]+1))~jitter(c(rep(24,nrow(invexp.data)),rep(96,nrow(invexp.data)),rep(168,nrow(invexp.data))),0.5), xlim=c(0,170), ylim=c(3,7),
       pch=16, cex=2, col=alpha("black",0.2))

#add community timepoint means as black points
points(as.numeric(comm_mean[1,])~tps,pch=21, cex=2, bg='white',col='black')
#add lines between community timepoint means
lines(comm_mean[1,]~tps, lwd=1.5, lty=1)

#add points and dotted line between inoculum mean and monoculture at 96hr
points(c(inocula_mean[1,1],mono_mean[1,1])~c(0,24), pch=4, cex=1,col=c("black",NA))
lines(c(inocula_mean[1,1],mono_mean[1,1])~c(0,24), pch=4, cex=1,col=c("black",NA), lty=2)
#add points and dotted line between inoculum mean and community at 96hr
points(c(inocula_mean[1,1],comm_mean[1,1])~c(0,24), pch=4, cex=1,col=c("black",NA))
lines(c(inocula_mean[1,1],comm_mean[1,1])~c(0,24), pch=4, cex=1,col=c("black",NA), lty=2)
#add line at detection threshold
abline(h=log10(kt.detectionthresh), lwd=2, col=alpha("black",0.2))

#label numbers below detection threshold
belowdet_n<-apply(invexp.data[,1:3]<kt.detectionthresh,2,sum)
text(96,2,"Numbers below detection limit (n = 678):",col='black', cex=1.25)
text(24,1.5,belowdet_n[1], col='black', cex=1.25)
text(96,1.5,belowdet_n[2], col='black', cex=1.25)
text(168,1.5,belowdet_n[3], col='black', cex=1.25)

# X AXIS LABEL ------------------------------------------------------------

plot.new()
#x axis label
text(0.5,0.5,"Time since invasion (hours)", cex=2.5)

rm(inocula_mean)
rm(mono_mean)
rm(mono_se)
rm(comm_mean)
rm(comm_se)

dev.off()