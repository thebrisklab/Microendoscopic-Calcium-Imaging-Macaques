###############################
# Author: Benjamin Risk
# 
# Create plots of Z-normalized Jaccard and distance
# NOTE: This file needs to be run after JaccardIndexAnalysisAllSessions.R
###############################
library(readxl)
library(dplyr)
library(corrplot)
library(ggplotify)
library(devEMF)

#---------------->
# helper functions:
mat.na = function(mymat,p.value=FALSE) {
  if(p.value) {
    mymat[is.na(mymat)]=1
    mymat[abs(mymat)==Inf]=1
    mymat
  } else {
    mymat[is.na(mymat)]=0
    mymat[abs(mymat)==Inf]=0
    mymat
  }
}

# this function subsets to cells with spikes
# (e.g., greater than 0 during the spontaneous sesssion)
# and replaces data with color limits of 
# plot for compatibility with corrplot
mat.thresh = function(mymat,col.lim=NULL,threshold=NULL) {
  diag(mymat)=NA
  temp = apply(mymat,2,function(x) all(is.na(x)))
  cells = colnames(mymat)[!temp]
  newmat = mymat[cells,cells]
  if(!is.null(col.lim)) {
    newmat[newmat>=col.lim[2]]=col.lim[2]
    newmat[newmat<=col.lim[1]]=col.lim[1]
  }
  if(!is.null(threshold)) {
    newmat[abs(newmat)<threshold]=0
  }
  newmat
}

#<-----------------------

# Load Jaccard matrices:
# Z Jaccard is NA if there were 0 spikes in one of the cells
# Use z-normalized Jaccard

pdf('Figures/MonkeyQ_SMA_JaccardDistance.pdf')
par(mfrow=c(1,2))

load('./Results/Jac_Time_Series_Quartz_L_2022_09_15.RData')
corrplot(mat.thresh(p.spon$z.jaccard,col.lim=c(-5,5),threshold = NULL), method = 'circle', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey Q - SMA',cl.length=6,tl.cex=0.5,order='hclust',col.lim=c(-5,5))

mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<0,na.rm=TRUE)
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>=0,na.rm=TRUE)

jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.spon, span=1)
#model0 = mgcv::gam(z.jaccard~s(distance),data=jaccard.data.spon)

plot(z.jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
# lines(predict(model0), x=jaccard.data.spon$distance, col="dodgerblue",lwd=2)
dev.off()

pdf('Figures/MonkeyU_SMA_JaccardDistance.pdf')
par(mfrow=c(1,2))

load('./Results/Jac_Time_Series_Ulrik_L_2023_03_21.RData')
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<0,na.rm=TRUE)
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>=0,na.rm=TRUE)


corrplot(mat.thresh(p.spon$z.jaccard,col.lim=c(-5,5),threshold = NULL), method = 'circle', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey U - SMA',cl.length=6,tl.cex=0.5,order='hclust',col.lim=c(-5,5))

jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.spon, span=1)
plot(z.jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()

pdf('Figures/MonkeyF_SMA_JaccardDistance.pdf')
par(mfrow=c(1,2))

load("Results/Jac_Fuji_L_S_2024_02_23_09_04_09.RData")
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<0,na.rm=TRUE)
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>=0,na.rm=TRUE)


corrplot(mat.thresh(p.spon$z.jaccard,col.lim=c(-5,5),threshold = NULL), method = 'circle', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey F - SMA',cl.length=6,tl.cex=0.5,order='hclust',col.lim=c(-5,5))

jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.spon, span=1)
plot(z.jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)

dev.off()

model.fuji = lm(z.jaccard~distance, data = jaccard.data.spon)
plot(model.fuji)

####### Create plot for M1:


#emf('Figures/M1_JaccardDistance.emf')
#par(mfrow=c(2,2))

pdf('Figures/MonkeyU_M1_JaccardDistance.pdf')
par(mfrow=c(1,2))

load("Results/Jac_Time_Series_Ulrik_R_2023_03_01.RData")
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<0,na.rm=TRUE)
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>=0,na.rm=TRUE)

corrplot(mat.thresh(p.spon$z.jaccard,col.lim=c(-5,5),threshold = NULL), method = 'circle', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey U - M1',cl.length=6,tl.cex=0.5,order='hclust',col.lim=c(-5,5))

jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.spon, span=1)
plot(z.jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()

pdf('Figures/MonkeyV_M1_JaccardDistance.pdf')
par(mfrow=c(1,2))
load("Results/Jac_Time_Series_Vader_R_2023_03_16.RData")
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<0,na.rm=TRUE)
mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>=0,na.rm=TRUE)

corrplot(mat.thresh(p.spon$z.jaccard,col.lim=c(-5,5),threshold = NULL), method = 'circle', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey V - M1',cl.length=6,tl.cex=0.5,order='hclust',col.lim=c(-5,5))

jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.spon, span=1)
plot(z.jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()


pdf('Figures/SMA_unnormalizedJaccardDistance.pdf')
par(mfrow=c(3,2))

load('./Results/Jac_Time_Series_Quartz_L_2022_09_15.RData')

jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(jaccard ~ distance, data=jaccard.data.spon, span=1)
#model0 = mgcv::gam(z.jaccard~s(distance),data=jaccard.data.spon)
plot(jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
# lines(predict(model0), x=jaccard.data.spon$distance, col="dodgerblue",lwd=2)

load('./Results/Jac_Time_Series_Ulrik_L_2023_03_21.RData')
jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(jaccard ~ distance, data=jaccard.data.spon, span=1)
plot(jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
load("Results/Jac_Fuji_L_S_2024_02_23_09_04_09.RData")
jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(jaccard ~ distance, data=jaccard.data.spon, span=1)
plot(jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()

pdf('Figures/M1_unnormalized_JaccardDistance.pdf')
par(mfrow=c(2,2))
load("Results/Jac_Time_Series_Ulrik_R_2023_03_01.RData")
jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(jaccard ~ distance, data=jaccard.data.spon, span=1)
plot(jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
load("Results/Jac_Time_Series_Vader_R_2023_03_16.RData")
jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(jaccard ~ distance, data=jaccard.data.spon, span=1)
plot(jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()


##############
#############
######################
# plots for Ulrik:
# These plots appear in the iScience manuscript:
load("Results/Jac_Time_Series_Ulrik_R_2023_03_16.RData")

pdf('Figures/MonkeyU_3_16_zJaccard_M1_Spon_v2.pdf')
par(mfrow=c(1,2))

corrplot(abs(mat.thresh(p.spon$z.jaccard,col.lim=c(-3,3),threshold = NULL)), method = 'color', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey U - M1 - Spontaneous',cl.length=6,tl.cex=0.5,col.lim=c(0,3),tl.col='black')


jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.spon, span=3)
plot(z.jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19,xlim=c(0,600))
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()

summary(lm(z.jaccard~distance,data=jaccard.data.spon))


pdf('Figures/MonkeyU_3_16_zJaccard_M1_TS_v2.pdf')
par(mfrow=c(1,2))
corrplot(abs(mat.thresh(p.ts$z.jaccard,col.lim=c(-3,3),threshold = NULL)), method = 'color', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey U - M1 - Reaching Task',cl.length=6,tl.cex=0.5,col.lim=c(0,3),na.label=' ',tl.col='black')
jaccard.data.ts = jaccard.data.ts[order(jaccard.data.ts$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.ts, span=3)
plot(z.jaccard~distance,data = jaccard.data.ts, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19,xlim=c(0,600))
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()
summary(lm(z.jaccard~distance,data=jaccard.data.ts))



load("Results/Jac_Time_Series_Ulrik_L_2023_03_14.RData")

# additional edits after submission: 
# use abs(z-jaccard) in loess plots:
pdf('Figures/MonkeyU_3_14_zJaccard_SMA_Spon_v2.pdf')
par(mfrow=c(1,2))
corrplot(abs(mat.thresh(p.spon$z.jaccard,col.lim=c(-3,3),threshold = NULL)), method = 'color', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey U - SMA - Spontaneous',cl.length=6,tl.cex=0.5,col.lim=c(0,3),tl.col='black')
jaccard.data.spon = jaccard.data.spon[order(jaccard.data.spon$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.spon, span=3)
plot(z.jaccard~distance,data = jaccard.data.spon, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()
summary(lm(z.jaccard~distance,data=jaccard.data.spon))


#emf('Figures/MonkeyU_3_14_zJaccard_SMA_TS.emf')
pdf('Figures/MonkeyU_3_14_zJaccard_SMA_TS_v2.pdf')
par(mfrow=c(1,2))
corrplot(abs(mat.thresh(p.ts$z.jaccard,col.lim=c(-3,3),threshold = NULL)), method = 'color', type = 'lower', diag=FALSE,is.corr=FALSE,mar=c(2,1,2,2),main='Monkey U - SMA - Reaching Task',cl.length=6,tl.cex=0.5,col.lim=c(0,3),na.label=' ',tl.col='black')
jaccard.data.ts = jaccard.data.ts[order(jaccard.data.ts$distance),]
loess.smooth = loess(z.jaccard ~ distance, data=jaccard.data.ts, span=3)
plot(z.jaccard~distance,data = jaccard.data.ts, type="p", ylab="Z Jaccard",xlab=expression(paste('Distance (',mu,'m)')), pch=19)
lines(predict(loess.smooth), x=loess.smooth$x, col="dodgerblue",lwd=2)
dev.off()
summary(lm(z.jaccard~distance,data=jaccard.data.ts))
