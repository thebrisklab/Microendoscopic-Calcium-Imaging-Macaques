########################
# Author: Benjamin Risk
# Updated 2/3/25
# Calculate the proportion of cells
# that are synchronized using the Jaccard
# index. This code runs the function to calculate
# the z-Jaccard on each dataset. It calls jaccard.wrapper,
# which is contained in permfun.R.
# This code then creates boxplots of the proportion
# synchronized in spontaneous and task conditions 
# in M1 and SMA across all sessions. 

library(corrplot)
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(devEMF)
source('permfun.R')

# read in metadata sheet:
mdata = read_excel('Metadata_all_monkeys.xlsx')


##################
# Estimate proportion of cells synchronized
# grab file stub names for SMA:
monkey_files = unique(mdata$fname[mdata$Side%in%c('L (SMA)','SMA')])

# remove files already analyzed:
monkey_files = monkey_files[!monkey_files%in%c('Time_Series_Quartz_L_2022_09_15','Fuji_L_S_2024_02_23_09_04_09','Time_Series_Ulrik_L_2023_03_21')]

# remove fuji file with only 1 cell:
monkey_files = monkey_files[monkey_files!='Fuji_L_S_2024_03_19_10_21_14']

# set to 0 if  you have the results already saved:
# set to 1 if you want to run the full analysis (it takes awhile):
if (0) {
for (temp_file in monkey_files) {
  a = proc.time()
  message(paste0("Starting ",temp_file,'\n'))
  jaccard.wrapper(temp_file,nperm=10)
  proc.time()-a
}
}

# nrow(spikedata[spikedata$session=='Spontaneous',])
# nrow(spikedata[spikedata$session=='Task',])
# only 24, did we use all these data in the plots? 

#############################
# grab file stub names for M1:
monkey_files = unique(mdata$fname[mdata$Side%in%c('R (M1)')])
# remove files already run:
monkey_files = monkey_files[!monkey_files%in%c('Time_Series_Ulrik_R_2023_03_01','Time_Series_Vader_R_2023_03_16')]

# read in time series data to find the rest versus task period:
# set to 1 if running complete analysis
# set to 0 if you have saved results from previously compiling code
if (0) {
for (temp_file in monkey_files) {
  a = proc.time()
  jaccard.wrapper(temp_file,nperm=1000)
  proc.time()-a
}
}

#########################################
########################################
## Collect data from all sessions
# monkey_files = unique(mdata$fname)
templist = dir('./Results')
templist2 = substr(templist,5,stop=nchar(templist))
monkey_files = substr(templist2,1,nchar(templist2)-6)

allsessions_summary=NULL  
threshold = 1.96
## FUJI, we only have spontaneous:
filestub = "Fuji_L_S_2024_02_13_10_32_32"
monkey='F'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Fuji_L_S_2024_02_20_09_11_32"
monkey='F'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Fuji_L_S_2024_02_23_09_04_09"
monkey='F'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Fuji_L_S_2024_02_27_10_16_24"
monkey='F'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Fuji_L_S_2024_02_29_10_45_22"
monkey='F'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Fuji_L_S_2024_03_05_09_36_11"
monkey='F'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Fuji_L_S_2024_03_07_09_17_10"
monkey='F'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Fuji_L_S_2024_03_12_10_24_36"
monkey='F'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)


#################################
## Quartz SMA:
filestub = "Time_Series_Quartz_L_2022_09_15"
monkey='Q'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Quartz_L_2022_09_22"
monkey='Q'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Quartz_L_2022_10_03"
monkey='Q'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Quartz_L_2022_10_28"
monkey='Q'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Quartz_L_2022_12_01"
monkey='Q'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Quartz_L_2022_12_02"
monkey='Q'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Quartz_L_2022_12_06"
monkey='Q'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

###################
# Ulrik SMA:
filestub = "Time_Series_Ulrik_L_2023_02_27"
monkey='U'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Ulrik_L_TimeSeries_2023_03_07"
monkey='U'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)


filestub = "Time_Series_Ulrik_L_2023_03_14"
monkey='U'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Ulrik_L_2023_03_21"
monkey='U'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Ulrik_L_2023_03_28"
monkey='U'
region='SMA'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)



###################
## M1:

# Ulrik:
filestub = "Time_Series_Ulrik_R_2023_02_24"
monkey='U'
region='M1'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Ulrik_R_2023_03_01"
monkey='U'
region='M1'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Ulrik_R_2023_03_16"
monkey='U'
region='M1'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Ulrik_R_2023_03_23"
monkey='U'
region='M1'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

###############
# Vader
filestub = "Time_Series_Vader_R_2023_03_09"
monkey='V'
region='M1'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)

filestub = "Time_Series_Vader_R_2023_03_16"
monkey='V'
region='M1'
load(paste0('Results/Jac_',filestub,'.RData'))
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Spontaneous','ncell'=dim(p.spon$z.jaccard)[1],'prop_sync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Spontaneous'))
allsessions_summary=rbind(allsessions_summary,temp)
temp=data.frame('filestub'=filestub,'monkey'=monkey,'region'=region,'condition'='Reaching Task','ncell'=dim(p.ts$z.jaccard)[1],'prop_sync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]>threshold,na.rm=TRUE),'prop_nsync'=mean(p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)]<(-1*threshold),na.rm=TRUE),'nbin'=sum(spikedata$session=='Task'))
allsessions_summary=rbind(allsessions_summary,temp)
allsessions_summary$monkey = factor(allsessions_summary$monkey)
############




###########
##############
#################
# Plot results
#emf('Figures/AllSessions_ProportionSynchronized.emf')

# change order of SMA and M1:
allsessions_summary$region=relevel(factor(allsessions_summary$region),ref='SMA')
# 2/3/25: This version includes both positive and negative synchrony:
allsessions_summary$condition <- factor(allsessions_summary$condition, levels = rev(sort(unique(allsessions_summary$condition))))

pdf('Figures/AllSessions_ProportionSynchronized.pdf')
allsessions_summary$prop_pnsync = allsessions_summary$prop_sync+allsessions_summary$prop_nsync
allsessions_summary%>%ggplot(aes(x=condition,y=prop_pnsync)) + geom_boxplot(outlier.shape = NA) +geom_jitter(aes(colour=monkey),alpha=0.8,width=0.025,height=0)+geom_line(aes(group = filestub,colour=monkey), alpha = 0.5) + facet_wrap(~region) + xlab('') + theme(text = element_text(size=15))+ylab('Proportion Synchronized')
dev.off()


# 2/3/25: Positive and negative synchrony:
temp.data = allsessions_summary%>%select(c('filestub','monkey','region','condition','prop_pnsync'))
allsession_wide =  temp.data%>%pivot_wider(names_from = condition, values_from = prop_pnsync)
m1 = allsession_wide%>%filter(region=='M1')
wilcox.test(x=m1$Spontaneous,y=m1$`Reaching Task`,paired=TRUE)

sma = allsession_wide%>%filter(region=='SMA')
wilcox.test(x=sma$Spontaneous,y=sma$`Reaching Task`,paired=TRUE)

# We defined the proportion of cells that were synchronized as the fraction of cell pairs with Z-Jaccard greater than or less than 1.96. For the sessions in which we had both spontaneous and touchscreen conditions, we compared whether there was a change in the proportion of cells synchronized using a Wilcoxon signed rank test.

# The proportion of cells that were synchronized tended to increase during the touchscreen, while the proportion tended to decrease in SMA (Figure XX), although the trends were not significant (p=0.22 and 0.07 in Wilcoxon signed rank, respectively). 

########
# Added 1/8/25 in response to reviewer:
# Summarize number of bins across sessions:

summary(allsessions_summary$nbin[allsessions_summary$condition=='Spontaneous'])

summary(allsessions_summary$nbin[allsessions_summary$condition=='Touchscreen'])
