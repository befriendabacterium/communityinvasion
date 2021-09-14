# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

# READ IN INVASION DATA ---------------------------------------------------

#read in calibration assay data
btiptg<-read.csv('inputs/1_raw/invasion/invadercalibrationgrowthcurve_96h_8reps.csv')
#read in raw invasion data
inv<-read.csv('inputs/3_ready/invasion/invasiondata_luminescence_assaymeans_matched.csv', row.names = 1)

# SBW25 calibration -------------------------------------------------------

btiptg.sbw25<-btiptg[btiptg$invader=='SBW25',]
tps<-c(12,14,16,18,20,22,24,96,168)

# SBW25 calibration -------------------------------------------------------

#sbw25
#cfu
log10.ticks.y<-sapply(c(4:8), function(i) as.expression(bquote(10^ .(i))))
plot(log10(colMeans(btiptg.sbw25[,7:15]))~tps, type='n', ylab='Growth (cfu/ml)', xlab='Time (hours)', yaxt='n', ylim=c(4,8))
axis(2, 4:8, log10.ticks.y, las=2)
points(log10(colMeans(btiptg.sbw25[,7:15]))~tps, pch=20)
lines(log10(colMeans(btiptg.sbw25[,7:15]))~tps)

#lum
log10.ticks.y<-sapply(c(2:5), function(i) as.expression(bquote(10^ .(i))))
plot(log10(colMeans(btiptg.sbw25[,16:24]))~tps, type='n', ylab='Luminescence (lumens)', xlab='Time (hours)', yaxt='n',ylim=c(2,5))
axis(2, 2:5, log10.ticks.y, las=2)
points(log10(colMeans(btiptg.sbw25[,16:24]))~tps, pch=20)
lines(log10(colMeans(btiptg.sbw25[,16:24]))~tps)

#loop to make vector of 12h-16h growth data (active phase when lum is linearly related to cfu)
all.sbw25.cfus<-c()
for (i in 1:8){
  temp<-as.numeric(btiptg.sbw25[i,7:9])
  all.sbw25.cfus<-c(all.sbw25.cfus, temp)
}

#loop to make vector of 12h-16h lux data (active phase when lum is linearly related to cfu)
all.sbw25.lums<-c()
for (j in 1:8){
  temp<-as.numeric(btiptg.sbw25[j,16:18])
  all.sbw25.lums<-c(all.sbw25.lums, temp)
}

all.sbw25.cfus[21]<-NA #remove contaminated sample
all.sbw25.lums[21]<-NA #remove contaminated sample

#plot log10 cfu vs log10 lum
log10.ticks.y<-sapply(c(4:7), function(i) as.expression(bquote(10^ .(i))))
log10.ticks.x<-sapply(c(1:4), function(i) as.expression(bquote(10^ .(i))))
plot(log10(all.sbw25.cfus)~log10(all.sbw25.lums), cex=0.2, xlim=c(1,4), ylim=c(4,7), yaxt='n', xaxt='n', xlab='Luminescence (lumens)', ylab='Number of cells (cfu/mL)')
axis(2, 4:7, log10.ticks.y, las=2)
axis(1,1:4, log10.ticks.x,las=1)
#linear model to fit cfu vs log10 lum
fit.sbw25<-lm(log10(all.sbw25.cfus)~log10(all.sbw25.lums))
summary(fit.sbw25)
#plot linear model line
abline(fit.sbw25)

# KT2440 CALIBRATION ------------------------------------------------------

btiptg.kt2440<-btiptg[btiptg$invader=='KT2440',]

#sbw25
#cfu
log10.ticks.y<-sapply(c(4:8), function(i) as.expression(bquote(10^ .(i))))
plot(log10(colMeans(btiptg.kt2440[,7:15]))~tps, ylab='Growth (cfu/ml)', xlab='Time (hours)', yaxt='n', ylim=c(4,8))
axis(2, 4:8, log10.ticks.y, las=2)
points(log10(colMeans(btiptg.kt2440[,7:15]))~tps, pch=20)
lines(log10(colMeans(btiptg.kt2440[,7:15]))~tps)
#lum
log10.ticks.y<-sapply(c(2:5), function(i) as.expression(bquote(10^ .(i))))
plot(log10(colMeans(btiptg.kt2440[,16:24]))~tps, type='n', ylab='Luminescence (lumens)', xlab='Time (hours)', yaxt='n',ylim=c(2,5))
axis(2, 2:5, log10.ticks.y, las=2)
points(log10(colMeans(btiptg.kt2440[,16:24]))~tps,pch=20)
lines(log10(colMeans(btiptg.kt2440[,16:24]))~tps)

#loop to make vector of 12h-16h growth data (active phase when lum is linearly related to cfu)
all.kt2440.cfus<-c()
for (i in 1:8){
  temp<-as.numeric(btiptg.kt2440[i,7:9])
  all.kt2440.cfus<-c(all.kt2440.cfus, temp)
}

#loop to make vector of 12h-16h lum data (active phase when lum is linearly related to cfu)
all.kt2440.lums<-c()
for (j in 1:8){
  temp<-as.numeric(btiptg.kt2440[j,16:18])
  all.kt2440.lums<-c(all.kt2440.lums, temp)
}
#plot linear model line
all.kt2440.cfus[13]<-NA #remove contaminated samples
all.kt2440.lums[13]<-NA #remove contaminated samples


#plot log10 cfu vs log10 lum
log10.ticks.y<-sapply(c(4:7), function(i) as.expression(bquote(10^ .(i))))
plot(log10(all.kt2440.cfus)~log10(all.kt2440.lums), cex=0.2,  xlim=c(1,4), ylim=c(4,7), yaxt='n', xaxt='n', xlab='Luminescence (lumens)', ylab='Number of cells (cfu/mL)')
axis(2, 4:7, log10.ticks.y, las=2)
axis(1,1:4, log10.ticks.x,las=1)
#linear model to fit cfu vs log10 lum
fit.kt2440<-lm(log10(all.kt2440.cfus)~log10(all.kt2440.lums))
summary(fit.kt2440)
#plot linear model line
abline(fit.kt2440)

# GRAPH SHOWING BOTH SBW AND KT CELLS-LUM RELATIONSHIP --------------------

plot.new
par(mar=c(5,7,1,2))

plot(log10(all.sbw25.cfus)~log10(all.sbw25.lums), cex=1, yaxt='n',xaxt='n', ylab='',ylim=c(4,7), xlim=c(1,4), xlab='Luminescence (lumens)')
axis(2, 4:7, log10.ticks.y, las=2)
axis(1,1:4, log10.ticks.x,las=1)
points(log10(all.kt2440.cfus)~log10(all.kt2440.lums), cex=1, ylim=c(0,17), xlim=c(0,9), xlab='log10 luminescence (lux)', ylab='Number of cells (cfu/mL)', col='red')
abline(fit.sbw25)
abline(fit.kt2440,col='red')
title(ylab='Number of cells (cfu/mL)', line=4)

# CONVERT LUM TO CFU FOR ACTUAL DATA------------------------------------------------------

#make a dataframe of which communities are below detection limit
belowdetectionlimit_df<-inv<=12
colSums(belowdetectionlimit_df)
belowdetectionlimit_df[,1]<-inv$Community

#convert sbw25s
summary(fit.sbw25)
sbw25.cfu.24h<-10^(fit.sbw25$coefficients[2]*log10(inv[,2])+fit.sbw25$coefficients[1])
sbw25.cfu.96h<-10^(fit.sbw25$coefficients[2]*log10(inv[,3])+fit.sbw25$coefficients[1])
sbw25.cfu.7d<-10^(fit.sbw25$coefficients[2]*log10(inv[,4])+fit.sbw25$coefficients[1])
detection.limit.sbw25<-10^(fit.sbw25$coefficients[2]*log10(12)+fit.sbw25$coefficients[1]) #background lux is 12 lumens, 
detection.limit.sbw25

#convert kt2440
summary(fit.kt2440)
kt2440.cfu.24h<-10^(fit.kt2440$coefficients[2]*log10(inv[,5])+fit.kt2440$coefficients[1])
kt2440.cfu.96h<-10^(fit.kt2440$coefficients[2]*log10(inv[,6])+fit.kt2440$coefficients[1])
kt2440.cfu.7d<-10^(fit.kt2440$coefficients[2]*log10(inv[,7])+fit.kt2440$coefficients[1])
detection.limit.kt2440<-10^(fit.kt2440$coefficients[2]*log10(12)+fit.kt2440$coefficients[1]) #background lux is 12 lumens
detection.limit.kt2440 #dl per 1ml

detectionlims<-c(detection.limit.sbw25,detection.limit.kt2440)
saveRDS(detectionlims,'inputs/3_ready/invasion/detectionlims.RDS')

#make dataframe of converted lum-cfu values
inv.cfu.cols<-as.data.frame(cbind(Community=inv$Community,
                                  sbw25.cfu.24h,sbw25.cfu.96h,sbw25.cfu.7d,
                                  kt2440.cfu.24h,kt2440.cfu.96h,kt2440.cfu.7d))

colnames(belowdetectionlimit_df)<-colnames(inv.cfu.cols)
write.csv(inv.cfu.cols, 'inputs/3_ready/invasion/invasiondata_cfu_assaymeans_matched.csv', row.names = F)
write.csv(belowdetectionlimit_df,'inputs/3_ready/invasion/invasiondata_belowdetectionlimit_matched.csv', row.names = F)

# END ---------------------------------------------------------------------
