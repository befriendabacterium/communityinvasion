# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")
#check working directory
getwd()

# LOAD & PREP DATA ---------------------------------------------------------------

source('src/X_data_prepper.R')
#colour palette for functional group 4

heatcols<-RColorBrewer::brewer.pal(11,'RdBu')
heatcols<-rev(grDevices::colorRampPalette(heatcols)(680))
palette_heat<-heatcols[as.factor(composition$C4)]

# SUPPLEMENTARY FIGURE 2 PANEL ----------------------------------------------------------

# Y AXIS LABEL -----------------------------------------------------------

tiff('outputs/figures/Figure3_Supplement3.tiff', res=300, units='in', width=11, height=11)


layout(matrix(c(1,2,2,3,3,
                1,2,2,3,3,
                1,2,2,3,3,
                1,2,2,3,3,
                1,2,2,3,3,
                1,2,2,3,3,
                1,2,2,3,3,
                1,2,2,3,3,
                1,2,2,3,3,
                1,2,2,3,3,
                1,4,4,5,5,
                1,4,4,5,5,
                1,4,4,5,5,
                1,4,4,5,5,
                1,4,4,5,5,
                1,4,4,5,5,
                1,4,4,5,5,
                1,4,4,5,5,
                1,4,4,5,5,
                1,4,4,5,5,
                6,6,6,6,6), nrow=21, byrow=T))
par(mar=c(5,0,3,3))

plot.new()
text(0.5,0.55,expression(paste(italic("P. putida"))), cex=1.75)
text(0.5,0.50," invader survival ", cex=1.75)
text(0.5,0.45,"(cells/ml at 24 h)", cex=1.75)

# PLOT Cellulose degradation at 7 days ---------------------------------------------------------------
log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$kt2440.cfu.24h~divandgrowth$mG7, 
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),xlim=c(0.5,3.1),
     xlab='',
     ylab='',
     yaxt='n',
     cex.lab=2)
axis(2,seq(3,7,1),log.ticks.y[4:8], las=2, cex.axis = 2)
mtext(expression(paste('Cellulose degradation at 7 days')),side=1,line=4)
mtext(expression(paste('(mg/ml ', beta,'-1-4-glucosidase)')), side=1, line=6)

model<-lm(invexp.data$kt2440.cfu.24h~divandgrowth$mG7)
text(2.75,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)

summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$mG7
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = scales::alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

# PLOT Chitin degradation at 7 days ---------------------------------------------------------------
log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$kt2440.cfu.24h~divandgrowth$mN7,
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),xlim=c(0.5,1.8),
     xlab='',
     ylab='',
     yaxt='n',
     cex.lab=1)

mtext(expression(paste('Chitinase degradation at 7 days')),side=1,line=4)
mtext(expression(paste('(mg/ml ','N-acetyl glucosamine)')), side=1, line=6)

model<-lm(invexp.data$kt2440.cfu.24h~divandgrowth$mN7)
text(1.6,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)

summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$mN7
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = scales::alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

# PLOT Phosphate degradation at 7 days ---------------------------------------------------------------
log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$kt2440.cfu.24h~divandgrowth$mP7, 
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),xlim=c(1.3,3.5),
     xlab='',
     ylab='',
     yaxt='n',
     
     cex.lab=1)
mtext(expression(paste('Phosphate degradation at 7 days')),side=1,line=4)
mtext(expression(paste('(mg/ml ','phosphotase)')), side=1, line=6)


axis(2,seq(3,7,1),log.ticks.y[4:8], las=2, cex.axis=2)

model<-lm(invexp.data$kt2440.cfu.24h~divandgrowth$mP7)
text(3.2,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)

summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$mP7
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = scales::alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

# PLOT Xylose degradation at 7 days ---------------------------------------------------------------
log.ticks.y<-sapply(c(0:10), function(i) as.expression(bquote(10^ .(i))))

plot(invexp.data$kt2440.cfu.24h~divandgrowth$mX7, 
     cex=1.5, pch=21, bg=palette_heat, 
     ylim=c(3,7),xlim=c(0.5,1.8),
     xlab=expression(paste('')),
     ylab='',
     yaxt='n',
     cex.lab=1)

mtext(expression(paste('Xylose degradation at 7 days')),side=1,line=4)
mtext(expression(paste('(mg/ml ','xylosidase)')), side=1, line=6)


model<-lm(invexp.data$kt2440.cfu.24h~divandgrowth$mX7)
text(1.6,6.9,substitute(paste("R"^2," = ",rsq), list(rsq=round(summary(model)$r.squared,2))), cex=2.5)
summary(model)
#Prediction intervals
pred.int =  predict(model,interval="prediction")

#Confidence intervals
conf.int =  predict(model,interval="confidence")
fitted.values = conf.int[,1]
conf.lower = conf.int[,2]
conf.upper = conf.int[,3]

x<-divandgrowth$mX7
ix <- sort(x,index.return=T)$ix
lines(x,fitted.values,col='black', lwd=2)
polygon(c(rev(x[ix]), x[ix]), c(rev(conf.int[ ix,3]), conf.int[ ix,2]), col = scales::alpha('grey',0.5), border = NA)
#polygon(c(rev(x[ix]), x[ix]), c(rev(pred.int[ ix,3]), pred.int[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

dev.off()
