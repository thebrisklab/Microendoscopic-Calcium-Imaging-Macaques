trace_readdat = function(infile,reject=TRUE) {
  tracedat = read.csv(infile)
  acc_rej = tracedat[1,]
  col_keep = acc_rej%in%c(" accepted","Time(s)/Cell Status")
  if(reject) {
    tracedat = tracedat[-1,col_keep]
  } else {
    tracedat = tracedat[-1,]
  }
  tracedat = data.frame(sapply(tracedat,as.numeric))
  tracedat
}

trace_readdatplus = function(infile,reject=TRUE) {
    tracedat = trace_readdat(infile,reject=reject)
    end_session = length(tracedat$X)
    end_spont = which(tracedat$X[2:end_session]-tracedat$X[1:(end_session-1)]>1)
    if (length(end_spont)>1) {
      warning("There are more than one gaps in the time series -- using min and max of the time breaks")
      end_spont1 = min(end_spont)
      begin_task = max(end_spont)+1
      tracedat = tracedat[c(1:end_spont1,begin_task:end_session),]
      tracedat$session='Task'
      tracedat$session[1:end_spont1]='Spontaneous'
    } else { 
      tracedat$session='Task'
      tracedat$session[1:end_spont]='Spontaneous'
    }
    maxspon = max(tracedat$X[tracedat$session=='Spontaneous'])
    mintask = min(tracedat$X[tracedat$session=='Task'])
    tracedat$newX = tracedat$X
    tracedat$newX[tracedat$session=='Task']=(tracedat$newX[tracedat$session=='Task']-
      (mintask-maxspon))
    tracedat
}


corrmat_ts_spon_function = function(infile,reject=TRUE) {
  tracedat = trace_readdat(infile,reject=reject)
  end_session = length(tracedat$X)
  end_spont = which(tracedat$X[2:end_session]-tracedat$X[1:(end_session-1)]>1)
  if (length(end_spont)>1) {
    error("Error: there are more than one gaps in the time series -- can't detect spontaneous versus task break")
  }
  trace_spon = tracedat[1:end_spont,]
  trace_ts = tracedat[(end_spont+1):end_session,]
  corrmat_spon = cor(trace_spon[,-1],use='pairwise.complete.obs')
  corrmat_ts = cor(trace_ts[,-1],use='pairwise.complete.obs')
  return(list('corrmat_spon'=corrmat_spon,'corrmat_ts'=corrmat_ts))
}


corrmat_ts_spon_perm_function = function(infile,reject=TRUE,nperm=100,ncores=8) {
  seed.list = as.list(sample(1:100000,nperm,replace=FALSE))
  library(parallel)
  tracedat = trace_readdat(infile,reject=reject)
  end_session = length(tracedat$X)
  end_spont = which(tracedat$X[2:end_session]-tracedat$X[1:(end_session-1)]>1)
  if (length(end_spont)>1) {
    error("Error: there are more than one gaps in the time series -- can't detect spontaneous versus task break")
  }
  trace_spon = tracedat[1:end_spont,]
  trace_ts = tracedat[(end_spont+1):end_session,]
  corrmat_spon = cor(trace_spon[,-1],use='pairwise.complete.obs')
  corrmat_ts = cor(trace_ts[,-1],use='pairwise.complete.obs')
  diff_mat = corrmat_spon-corrmat_ts
  z_diff = atanh(corrmat_spon)-atanh(corrmat_ts)
  z_diff_vector = z_diff[lower.tri(z_diff)]
  n_cell = ncol(trace_spon)-1;
  mat.perm.corr=NULL
  perm.fun.temp = function(x) {
    set.seed(x)
    circ_shift = sample(1:end_session,size=1)
    new_trace = tracedat[c((circ_shift+1):end_session,1:circ_shift),]
    trace_spon = new_trace[1:end_spont,]
    trace_ts = new_trace[(end_spont+1):end_session,]
    corrmat_p_spon = cor(trace_spon[,-1],use='pairwise.complete.obs')
    corrmat_p_ts = cor(trace_ts[,-1],use='pairwise.complete.obs')
    z_perm_diff = atanh(corrmat_p_spon)-atanh(corrmat_p_ts)
    z_perm_diff[lower.tri(z_perm_diff)]
  }
  list.perm.corr = mclapply(X=seed.list,FUN=perm.fun.temp,mc.cores=ncores)
  # list.perm.corr = lapply(X=seed.list,FUN=perm.fun.temp)
  mat.perm.corr = matrix(unlist(list.perm.corr),nrow=choose(n_cell,2))
  # calculate edge-wise (not corrected for fwer) p-values:
  p.values.vec=apply(X = abs(mat.perm.corr)>abs(z_diff_vector),1,FUN=mean)
  p.values.mat = matrix(0,n_cell,n_cell)
  p.values.mat[lower.tri(p.values.mat)]=p.values.vec
  p.values.mat = p.values.mat+t(p.values.mat)
  diag(p.values.mat)=1

  p.values.fdr = p.adjust(p.values.vec,method='fdr')
  p.values.mat.fdr = matrix(0,n_cell,n_cell)
  p.values.mat.fdr[lower.tri(p.values.mat.fdr)]=p.values.fdr
  p.values.mat.fdr = p.values.mat.fdr+t(p.values.mat.fdr)
  diag(p.values.mat.fdr)=1
  diff_mat.fdr.05 = diff_mat
  diff_mat.fdr.05[p.values.mat.fdr>0.05]=0
  diff_mat.fdr.05[lower.tri(diff_mat.fdr.05)]=diff_mat[lower.tri(diff_mat)]
  
  diff_mat.fdr.20 = diff_mat
  diff_mat.fdr.20[p.values.mat.fdr>0.20]=0
  diff_mat.fdr.20[lower.tri(diff_mat.fdr.20)]=diff_mat[lower.tri(diff_mat)]
  
  
  
  # calculate fwer-corrected p-values:
  # get vector of max statistics for all permutations: 
  max.vec = apply(abs(mat.perm.corr),2,max)
  fwer.p.values = sapply(z_diff_vector,function(x) mean(max.vec>abs(x)))
  fwer.p.values.mat = matrix(NA,n_cell,n_cell)
  fwer.p.values.mat[lower.tri(fwer.p.values.mat)]=fwer.p.values
  # fwer.p.values.mat[upper.tri(fwer.p.values.mat)]=diff_mat[upper.tri(diff_mat)]
  return(list('corrmat_spon'=corrmat_spon,'corrmat_ts'=corrmat_ts,'diff_mat'=diff_mat,'p.values.mat'=p.values.mat,'fwer.p.values.mat'=fwer.p.values.mat,'diff_mat.fdr.05'=diff_mat.fdr.05,'diff_mat.fdr.20'=diff_mat.fdr.20))
}


#######
createCorrPlots = function(out, fileend, reorder=FALSE) {
  require(corrplot)
  if (reorder==TRUE) {
    ave_corr = (out$corrmat_spon+out$corrmat_ts)/2
    temp_corrplot = corrplot(ave_corr,order='hclust')
    # extract ordering:
    temp_data = data.frame('cell'=rownames(temp_corrplot$corr),order=c(1:nrow(temp_corrplot$corr)))
    temp_data = temp_data[order(temp_data$cell),]
    temp_data$orig_order = c(1:nrow(temp_data))           
    temp_data = temp_data[order(temp_data$order),]
    
    corr_spon = out$corrmat_spon[temp_data$orig_order,temp_data$orig_order]
    corr_ts = out$corrmat_ts[temp_data$orig_order,temp_data$orig_order]
    corr_diff = out$diff_mat[temp_data$orig_order,temp_data$orig_order]
    temp = out$diff_mat.fdr.20
    temp[lower.tri(temp)]=0
    temp = temp+t(temp)
    temp = temp[temp_data$orig_order,temp_data$orig_order]
    temp[lower.tri(temp)] = corr_diff[lower.tri(corr_diff)]
    corr_diff_upper = temp
  } else {
    corr_spon = out$corrmat_spon
    corr_ts = out$corrmat_ts
    corr_diff = out$diff_mat
    corr_diff_upper = out$diff_mat.fdr.20
  }
  
  pdf(file=paste0('./Figures/CorrelationPlot_',fileend,'_Spontaneous.pdf'))
    corrplot(corr_spon)
  dev.off()
  
  pdf(file=paste0('./Figures/CorrelationPlot_',fileend,'_Task.pdf'))
  corrplot(corr_ts)
  dev.off()

  pdf(file=paste0('./Figures/CorrelationPlot_',fileend,'_Diff.pdf'))
  corrplot(corr_diff)
  dev.off()
  
  pdf(file=paste0('./Figures/CorrelationPlot_',fileend,'_Diff_upper_fdr20.pdf'))
  corrplot(corr_diff_upper)
  dev.off()
  
  
  png(file=paste0('Figures/CorrelationPlot_',fileend,'_four_plots.png'),height=1600,width=1600,pointsize=30)
  par(mfrow=c(2,2))
  corrplot(corr_spon,main='Spontaneous')
  corrplot(corr_ts,main='Touchscreen')
  corrplot(corr_diff,main='Difference')
  corrplot(corr_diff_upper,main='Significance (draft)')
  dev.off()
  
}


spike_rate_perm = function(infile,nperm,binsize=0.3,statistic=c('difference','tstatistic','cv','skewness')) {
 
  #--->code included in spike_time_to_bins
  spike_file = paste0('Data/',infile,'_Spikes.csv')
  spikes = read.csv(spike_file)
  names(spikes)[1] = 'Time'
  spikes$Cell.Name = trimws(spikes$Cell.Name)
  
  denoised_file = paste0('Data/',infile,'_Denoised.csv')
  denoised =  trace_readdatplus(denoised_file)
  
  midtime = mean(c(max(denoised$X[denoised$session=='Spontaneous']),min(denoised$X[denoised$session=='Task'])))
  spontime = max(denoised$X[denoised$session=='Spontaneous'])
  tasktime = max(denoised$X)-min(denoised$X[denoised$session=='Task'])
  spikes$session = NULL
  spikes$session[spikes$Time<midtime]='Spontaneous'
  spikes$session[spikes$Time>midtime]='Task'
  
  message('Total number of spikes in spontaneous versus task:')
  print(table(spikes$session,useNA='always'))
  
  maxtime = max(denoised$X)
  # remove the time without recordings to create continuous time series:
  deadtime = maxtime - (spontime+tasktime)
  
  spikes$Time[spikes$session=='Task']= spikes$Time[spikes$session=='Task'] - deadtime
  
  newmaxtime = max(denoised$X)-deadtime
  timegrid = seq(0,newmaxtime,by=binsize)
  ntime = length(timegrid)
  cellname = names(table(spikes$Cell.Name))
  ncell = length(table(spikes$Cell.Name))
  # drop the bin that has partial time at end of task:
  spikes_timegrid = data.frame('time'=timegrid[1:(ntime-1)])
  spikes_timegrid$session='Spontaneous'
  spikes_timegrid$session[spikes_timegrid$time>spontime]='Task'
  #<---- the above code is also included in spike_times_to_bins
# cellname
# ntime
# timegrid
# spikes
# cell.spikes
  
  count=1
  for (cell in cellname) {
    spikecounts = numeric(ntime-1)
    count = count+1
    cell.spikes = spikes[spikes$Cell.Name==cell,]
    for (t in 2:ntime) {
      binstarttime = timegrid[t-1]
      binendtime = timegrid[t]
      spikecounts[t-1] = sum(cell.spikes$Time>=binstarttime & cell.spikes$Time<binendtime)
    }
    spikes_timegrid[[count]]=spikecounts 
  }
  names(spikes_timegrid)[2:(ncell+1)]=cellname
  spikes_timegrid$session='Spontaneous'
  spikes_timegrid$session[spikes_timegrid$time>spontime]='Task'
  # identify time point that has partial spontaneous and task:
  which.point = ceiling(spontime/binsize)
  # should be equivalent:
  # which.point = which(spontime>timegrid[1:(ntime-1)] & spontime<timegrid[2:ntime])
  # delete this mixed time bin:
  # check it's the correct time bin:
  # timegrid[ceiling(spontime/binsize)]
  spikes_timegrid = spikes_timegrid[-which.point,]
  
  message('Number of spikes in each cell:')
  print(table(spikes$Cell.Name))
  
  observed_tstat = numeric(ncell)
 
  
   # naive t-statistics for each cell:
  count=1
  spon_rate=NULL
  task_rate=NULL
  ttest_pvalue=NULL
  for (cid in 2:(ncell+1)) { 
    model.temp = t.test(spikes_timegrid[,cid]~spikes_timegrid$session)
    observed_tstat[count] = model.temp$statistic
    spon_rate[count] = model.temp$estimate[1]
    task_rate[count] = model.temp$estimate[2]
    ttest_pvalue[count] = model.temp$p.value
    count = count+1
  }
  
  # construct test statistic from circularly shifted data:
  # extract spike data only for the function:
  spikedata = spikes_timegrid[,-c(1,ncol(spikes_timegrid))]
  if(nrow(spikedata)<nperm) message(paste0('Maximum possible permutations is ',nrow(spikedata),'; enumerating all possible shuffles'))
  pvalue = apply(spikedata,2,perm.stat.test,session=spikes_timegrid$session,nperm=nperm,statistic=statistic)
  
  # TO DO 3/6/24: Needs to be updated to correspond to the different statistics
  results.df = data.frame('fname'=infile,'spontaneous'=spon_rate/binsize,'touchscreen'=task_rate/binsize,'difference'=task_rate/binsize-spon_rate/binsize,ttest_pvalue,pvalue)
  
  results.df$`IDPS cell ID`=names(table(spikes$Cell.Name))
  results.df
}

# function that converts data from the spike times to a regular grid
# with each data point containing the number of spikes during that time interval
spike_times_to_bins = function(infile,binsize=0.2) {
    spike_file = paste0('Data/',infile,'_Spikes.csv')
    spikes = read.csv(spike_file)
    names(spikes)[1] = 'Time'
    spikes$Cell.Name = trimws(spikes$Cell.Name)
    
    denoised_file = paste0('Data/',infile,'_Denoised.csv')
    denoised =  trace_readdatplus(denoised_file)
    
    midtime = mean(c(max(denoised$X[denoised$session=='Spontaneous']),min(denoised$X[denoised$session=='Task'])))
    spontime = max(denoised$X[denoised$session=='Spontaneous'])
    tasktime = max(denoised$X)-min(denoised$X[denoised$session=='Task'])
    spikes$session = NULL
    spikes$session[spikes$Time<midtime]='Spontaneous'
    spikes$session[spikes$Time>midtime]='Task'
    
    message('Total number of spikes in spontaneous versus task:')
    print(table(spikes$session,useNA='always'))
    
    maxtime = max(denoised$X)
    # remove the time without recordings to create continuous time series:
    deadtime = maxtime - (spontime+tasktime)
    
    spikes$Time[spikes$session=='Task']= spikes$Time[spikes$session=='Task'] - deadtime
    
    newmaxtime = max(denoised$X)-deadtime
    timegrid = seq(0,newmaxtime,by=binsize)
    ntime = length(timegrid)
    cellname = names(table(spikes$Cell.Name))
    ncell = length(table(spikes$Cell.Name))
    # drop the bin that has partial time at end of task:
    spikes_timegrid = data.frame('time'=timegrid[1:(ntime-1)])
    spikes_timegrid$session='Spontaneous'
    spikes_timegrid$session[spikes_timegrid$time>spontime]='Task'

    count=1
    for (cell in cellname) {
      spikecounts = numeric(ntime-1)
      count = count+1
      cell.spikes = spikes[spikes$Cell.Name==cell,]
      for (t in 2:ntime) {
        binstarttime = timegrid[t-1]
        binendtime = timegrid[t]
        spikecounts[t-1] = sum(cell.spikes$Time>=binstarttime & cell.spikes$Time<binendtime)
      }
      spikes_timegrid[[count]]=spikecounts 
    }
    names(spikes_timegrid)[2:(ncell+1)]=cellname
    spikes_timegrid$session='Spontaneous'
    spikes_timegrid$session[spikes_timegrid$time>spontime]='Task'
    # identify time point that has partial spontaneous and task:
    which.point = ceiling(spontime/binsize)
    # should be equivalent:
    # which.point = which(spontime>timegrid[1:(ntime-1)] & spontime<timegrid[2:ntime])
    # delete this mixed time bin:
    # check it's the correct time bin:
    # timegrid[ceiling(spontime/binsize)]
    spikes_timegrid = spikes_timegrid[-which.point,]
    
    message('Number of bins with spikes in each cell:')
    print(table(spikes$Cell.Name))
    
    # identify index where switch from spontaneous to task:
    last_spon_index = which(spikes_timegrid$session[1:(ntime-1)]!=spikes_timegrid$session[2:ntime])
    count = ncol(spikes_timegrid)
    ntime = nrow(spikes_timegrid)
    for (cell in cellname) {
      fwd_smooth = rep(NA,ntime)
      count = count+1
      spikecounts = spikes_timegrid[,cell]
      for (t in 1:(ntime-4)) {
        fwd_smooth[t] = (sum(spikecounts[t:(t+4)])!=0)
      }
      fwd_smooth[(last_spon_index-3):last_spon_index]=NA 
      spikes_timegrid[[count]]=fwd_smooth
      names(spikes_timegrid)[[count]]=paste0('fwd5_',cell)
    }
    spikes_timegrid
}
###############################################

spike_rate_analysis = function(monkey,side,binsize,nperm,seed,statistic=c('difference','tstatistic','cv','skewness'),fdr=0.20,firstappearance=FALSE) {
  set.seed(seed)
  temp_files = unique(mdata$fname[mdata$Monkey==monkey & mdata$Side==side])
  temp_all = NULL  
  for (mfile in temp_files) {
    message(paste0('Processing ',mfile))
    temp_session = spike_rate_perm(infile=mfile,nperm=nperm,binsize = binsize,statistic=statistic)
    temp_all = rbind(temp_all,temp_session)
    message('******')
  }
  # save(temp_all,file=paste0(monkey,'_',side,'_spikes_all_bin',binsize,'s_',seed,'.RData'))
  
  colsSelect = c('fname','Global cell ID','IDPS cell ID','spontaneous','touchscreen','difference','ttest_pvalue','pvalue','X_centroid','Y_centroid')
  
  # read in metadata sheet:
  mdata = read_excel('Metadata_all_monkeys.xlsx')
  temp_all2 = merge(temp_all,mdata)
  
  # subset to the first cell appearance:
  temp_first = temp_all2%>%filter(`Cell occurr. No.`==1)%>%select(all_of(colsSelect))
  temp_first$pvalue.fdr = p.adjust(temp_first$pvalue,method = 'BH')
  temp_first$ttest_pvalue.fdr = p.adjust(temp_first$ttest_pvalue,method='BH')

  # Identify cells that were significantly higher in TS, Spontaneous, or no difference
  # according to the permutation test:
  temp_first$State='No Difference'
  temp_first$State[temp_first$pvalue.fdr<fdr & temp_first$difference>0]='Touchscreen'
  temp_first$State[temp_first$pvalue.fdr<fdr & temp_first$difference<0]='Spontaneous'
  
  # Identify cells that were significantly higher in TS, Spontaneous, or no difference
  # according to the t-test:
  temp_first$State_ttest='No Difference'
  temp_first$State_ttest[temp_first$ttest_pvalue.fdr<fdr & temp_first$difference>0]='Touchscreen'
  temp_first$State_ttest[temp_first$ttest_pvalue.fdr<fdr & temp_first$difference<0]='Spontaneous'
  
  message('Significant cells using circular shift test:')
  print(table(temp_first$State))
  message('**************\n')
  
  message('Significant cells using t-test:')
  print(table(temp_first$State_ttest))
  message('**************\n')
  
  spikes_long = temp_first%>%pivot_longer(names_to='Session',values_to='Rate',cols = c('spontaneous','touchscreen'))
  spikes_long2 = spikes_long%>%rename(id=`Global cell ID`)%>%select(id,Rate,Session,State,State_ttest)
  

  spikes_long2$State = factor(spikes_long2$State,levels=c('No Difference','Spontaneous','Touchscreen'))
  spikes_long2$State_ttest = factor(spikes_long2$State,levels=c('No Difference','Spontaneous','Touchscreen'))
  
  p1=ggplot(spikes_long2,aes(x = Session,y = Rate,group = id,color=State)) +  geom_line(size = 1,alpha = 0.25) +  geom_point(size = 1,alpha=0.5)+scale_color_manual(values = c("gray80","#FC4E07", "blue"))+ggtitle(paste0(monkey,'_',side,' - circshift test'))+ylab('Spikes/sec')+theme_light()
  
  p2=ggplot(spikes_long2,aes(x = Session,y = Rate,group = id,color=State_ttest)) +  geom_line(size = 1,alpha = 0.25) +  geom_point(size = 1,alpha=0.25)+scale_color_manual(values = c("gray80","#FC4E07", "blue"))+ggtitle(paste0(monkey,'_',side,' - t test'))+ylab('Spikes/sec')+theme_light()
  
  a= grid.arrange(p1,p2,nrow=1)
  
  ggsave(paste0('Figures/',monkey,'_',side,'_cell_first_appearance_',statistic,'.png'),
         a,
         width=8,
         height=3,
         dpi=1200)
  
  message('Compare whether there is an overall difference in firing rate across cells:')
  print(wilcox.test(x = temp_first$spontaneous, y = temp_first$touchscreen, paired=TRUE))
  message('**********\n')

  return(list('spikes_first'=temp_first,'spikes_first_long'=spikes_long))
}


################################################
#' Calculate the pvalue for a cell for spontaneous versus touchscreen
#'
#' @description 
#' This function calculates the shuffled-based permutation test. 
#' 
#' @details
#' This function uses permute::shuffleSet to generate the permuted data. It is often that case that the exact permutation statistic can be calculated from calcium spike data because the number of possible permutations is equal to the length of the binned dataset. 
#' @param x vector in which each element corresponds to the number of spikes in a time interval. The vector concatenates the spontaneous and touchscreen.
#' @param session vector of the session labels. Statistic is calculated for each session. The names of the session are either "Spontaneous" or "Task"
#' @returns  p-value (scalar)
perm.stat.test <- function(x, session, nperm,statistic=c('difference','tstatistic','cv')) {
  ## t-statistic function
  # i is vector of shuffled indices (length equal to length of data)
  tStat <- function(i, tdata, grp) {
    grp <- grp[i] #shuffle the data
    t.test(tdata[grp == "Spontaneous"],tdata[grp == "Task"])$statistic
  }
  
  ## mean difference function
  meanDif <- function(i, tdata, grp) {
    grp <- grp[i]
    mean(tdata[grp == "Spontaneous"]) - mean(tdata[grp == "Task"])
  }
  
  ## cv function
  cv <- function(i, tdata, grp) {
    grp <- grp[i]
    grp1 = tdata[grp=='Spontaneous']
    grp2 = tdata[grp=='Task']
    cv.stat = sd(grp1)/mean(grp1) - sd(grp2)/mean(grp2)
    if (is.na(cv.stat)) cv.stat=0
    cv.stat
    }
    
  ## skewness function
  skewness <- function(i, tdata, grp) {
    grp <- grp[i]
    grp1 = tdata[grp=='Spontaneous']
    grp2 = tdata[grp=='Task']
    skewness.stat = mean(((grp1 - mean(grp1))/sd(grp1))^3) - mean(((grp2 - mean(grp2))/sd(grp2))^3)
    if (is.na(skewness.stat)) skewness.stat=0
    skewness.stat
  }
  
  ## check x and session are of same length
  stopifnot(all.equal(length(x), length(session)))
  ## number of observations
  N <- nobs(x)
  ## generate the required set of permutations
  CTRL <- how(within = Within(type = "series"))
  pset <- shuffleSet(N, nset = nperm, control = CTRL, quietly=TRUE)
  
  ## iterate over the set of permutations applying 
  ## desired statistic
  if (statistic == 'difference') {
    funstat = meanDif
  } else if (statistic == 'tstatistic') {
    funstat = tStat 
  } else if (statistic == 'cv') {
    funstat = cv
  } else if (statistic =='skewness') {
    funstat = skewness
  }
  D <- apply(pset, 1, funstat, tdata = x, grp = session)
  ## add on the observed statistic
  D <- c(funstat(seq_len(N), tdata = x, session), D)
  ## compute & return the p-value
  mean(abs(D) >= abs(D[1])) # how many >= to the observed statistic?
}
###########################################
#' Calculate the Jaccard Index
#'
#' @description 
#' This function calculates the Jaccard similarity index from a dataset output from the file spike_times_to_bins
#' 
#' @details
#' <to be added>
#' @param spike_timegrid dataset output from the file spike_times_to_bins.
#' @param session character indicate either "spontaneous" or "task"
#' @returns a symmetric matrix of dimension ncell x ncell with jaccard similarity index between cells.

jaccard_calc = function(logicmat,return='short') {
  if (!return%in%c('short','long')) stop('set "return" to "short" or "long"')
  ntime = nrow(logicmat)
  ncell = ncol(logicmat)
  nevents = apply(logicmat,2,sum)
  mat_nevents = matrix(nevents,ncol=ncell,nrow=ncell)
  # calculate co-occurring events:
  mat_cooccur = t(logicmat)%*%logicmat
  rownames(mat_cooccur) = substr(rownames(mat_cooccur),6,nchar(rownames(mat_cooccur)))
  colnames(mat_cooccur) = rownames(mat_cooccur)
  jaccard = mat_cooccur/(mat_nevents + t(mat_nevents) - mat_cooccur)
  prop_row = mat_cooccur/mat_nevents
  prop_col = mat_cooccur/t(mat_nevents)
  # 17 april 2024
  if (return=='long') return(list('jaccard'=jaccard,'prop_row'=prop_row,'prop_col'=prop_col))
  else {
    jaccard
  }
}

#####################################

jaccard_calc_x1x2 = function(x1,x2) {
  ntime = length(x1)
  logicmat = cbind(x1,x2)
  nevents = apply(logicmat,2,sum)
  mat_nevents = matrix(nevents,ncol=2,nrow=2)
  # calculate co-occurring events:
  mat_cooccur = t(logicmat)%*%logicmat
  rownames(mat_cooccur) = substr(rownames(mat_cooccur),6,nchar(rownames(mat_cooccur)))
  colnames(mat_cooccur) = rownames(mat_cooccur)
  ((mat_cooccur)/(mat_nevents + t(mat_nevents) - mat_cooccur))[2,1]
}

################################################
#' Calculate the pvalue for a cell for a similarity measure
#'
#' @description 
#' This function calculates the shuffled-based permutation test for the significance
#' of the similarity (e.g., Jaccard) between two cells. 
#' 
#' @details
#' This function uses permute::shuffleSet to generate the permuted data. It is often that case that the exact permutation statistic can be calculated from calcium spike data because the number of possible permutations is equal to the length of the binned dataset. 
#' @param x1 vector for cell1 in which each element corresponds to whether there was a spike in the time interval
#' @param x2 vector for cell2 
#' @returns  p-value (scalar)
perm.similarity.test.x1x2 <- function(x1, x2, nperm, statistic='jaccard') {
  require(permute)
  ## calculate jaccard for two cells:
  jaccard <- function(i, x1, x2) {
    x2 <- x2[i]
    stat = jaccard_calc_x1x2(x1,x2)
   # if (is.na(stat)) stat=0
   stat
  }
  
  ## number of observations
  N <- length(x1)
  ## generate the required set of permutations
  CTRL <- how(within = Within(type = "series"))
  pset <- shuffleSet(N, nset = nperm, control = CTRL, quietly=TRUE)
  
  ## iterate over the set of permutations applying 
  ## desired statistic
  if (statistic == 'jaccard') {
    funstat = jaccard
  } else {
    message("Only Jaccard implemented at this time")
  }
  
  D <- apply(pset, 1, funstat, x1 = x1, x2 = x2)
  obs.D = funstat(seq_len(N), x1=x1, x2=x2)
  sd.perm = sd(D)
  mean.perm = mean(D)
  norm.jacc =  obs.D/mean.perm # normalization used in Parker
  z.jacc = (obs.D - mean.perm)/sd.perm
  z.perm = (D - mean.perm)/sd.perm
  ## add on the observed statistic
  # D <- c(obs.D, D)
  # B Risk:
  # use z-statistic instead of D in order to be able to 
  # detect Jaccard values that are significantly LOWER
  # note: need to append original statistic to get full set of permutations
  z.perm = c(z.jacc,z.perm)
  
  ## compute & return the p-value
  #p.value = mean(abs(D) >= abs(obs.D)) # how many >= to the observed statistic?
  p.value = mean(abs(z.perm) >= abs(z.jacc)) # how many >= to the observed statistic?
  return(list('jacc'=obs.D, 'p.value'=p.value, 'z.jacc'=z.jacc, 'norm.jacc'=norm.jacc))
}

# Main function:
# add normalized Jaccard, like Parker
perm.similarity.test.mat <- function(spike_logical, nperm, statistic='jaccard') {
  ncell = ncol(spike_logical)
  pmat = matrix(0,ncell,ncell)
  jacc = pmat
  norm.jacc = pmat
  z.jacc = pmat
  for (i in 1:(ncell-1)) {
    for (j in (i+1):ncell) {
      out = perm.similarity.test.x1x2(x1=spike_logical[,i], x2=spike_logical[,j], nperm=nperm,statistic=statistic)
      pmat[i,j] = out$p.value
      jacc[i,j] = out$jacc
      norm.jacc[i,j]= out$norm.jacc
      z.jacc[i,j]=out$z.jacc
    }
  }
  rownames(pmat) = substr(colnames(spike_logical),6,nchar(colnames(spike_logical)))
  colnames(pmat) = rownames(pmat)
  pvec = pmat[upper.tri(pmat)]
  fdr.pvec = p.adjust(pvec,method='fdr')
  fdr.pmat = matrix(0,ncell,ncell)
  fdr.pmat[upper.tri(pmat)] = fdr.pvec
  pmat = pmat + t(pmat)
  fdr.pmat = fdr.pmat+t(fdr.pmat)
  diag(pmat) = 1
  diag(fdr.pmat)=1
  colnames(fdr.pmat) = colnames(pmat)
  rownames(fdr.pmat) = rownames(pmat)
  jacc = jacc+t(jacc)
  norm.jacc = norm.jacc+t(norm.jacc)
  z.jacc = z.jacc+t(z.jacc)
  colnames(jacc) = colnames(pmat)
  rownames(jacc)=rownames(pmat)
  colnames(norm.jacc) = colnames(pmat)
  rownames(norm.jacc) = rownames(pmat)
  colnames(z.jacc) = colnames(pmat)
  rownames(z.jacc) = rownames(pmat)
  return(list('jaccard'=jacc,'norm.jaccard'=norm.jacc,'z.jaccard'=z.jacc,'unadjusted.p'=pmat,'fdr.p'=fdr.pmat))
}




###################################
####################

# grab traces for comparison:
jaccard.wrapper = function(filestub,nperm=1000) {
  temp_file = filestub
  denoised = trace_readdatplus(paste0('Data/',temp_file,'_Denoised.csv'))
  
  # create spike data:
  spikedata = spike_times_to_bins(infile=temp_file,binsize = 0.2)
  # find the index where logical data start:
  myindex = which(names(spikedata)=='session')+1
  # remove the NAs; this corresponds to last four observations
  # when using 5 bin forward smoothing
  spike_logical = spikedata[!is.na(spikedata[[myindex]]),]
  # subset to either spontaneous or task:
  spike_spon = as.matrix(spike_logical[spike_logical$session=='Spontaneous',myindex:ncol(spike_logical)])
  spike_ts = as.matrix(spike_logical[spike_logical$session=='Task',myindex:ncol(spike_logical)])
  
  # calculate jaccard:
  spon = jaccard_calc(spike_spon)
  ts = jaccard_calc(spike_ts)
  spon = spon - diag(diag(spon))
  ts = ts - diag(diag(ts))
  spon[is.na(spon)]=0
  ts[is.na(ts)]=0
  
  p.spon = perm.similarity.test.mat(spike_spon,nperm=nperm)
  p.ts = perm.similarity.test.mat(spike_ts,nperm=nperm)
  
  # Examine relationship between firing rates and distance:
  # 1 pixel = 3.2 micrometer um when plotting the Jaccard analysis with the distance between cells.
  centroids = mdata%>%filter(fname==temp_file & `IDPS cell ID`%in%colnames(spon))%>%select(c(`IDPS cell ID`,X_centroid,Y_centroid))
  # check ordering matches:
  
  if (all(centroids$`IDPS cell ID`==colnames(spon))) {
  distmat = as.matrix(dist(centroids[,2:3],upper=TRUE))
  # convert to micrometers:
  distmat = distmat*3.2
  colnames(distmat) = colnames(spon)
  rownames(distmat) = colnames(spon)
  
  cor.spon = cor(denoised[denoised$session=='Spontaneous',names(denoised)%in%colnames(spon)],use = 'complete.obs')
  diag(cor.spon)=0
  cor.ts = cor(denoised[denoised$session=='Task',names(denoised)%in%colnames(ts)],use = 'complete.obs')
  diag(cor.ts)=0
  
  # approach using z-jaccard:
  centroids = mdata%>%filter(fname==temp_file & `IDPS cell ID`%in%colnames(p.spon$jaccard))%>%select(c(`IDPS cell ID`,X_centroid,Y_centroid,`S Rate`))
  
  # check ordering matches:
  all(centroids$`IDPS cell ID`==colnames(p.spon$jaccard))
  
  distmat = as.matrix(dist(centroids[,2:3],upper=TRUE))
  # convert to micrometers:
  distmat = distmat*3.2
  colnames(distmat) = colnames(spon)
  rownames(distmat) = colnames(spon)
  # Z Jaccard is NA if there were 0 spikes in one of the cells
  # Use z-normalized Jaccard
  jaccard.data.spon = data.frame(z.jaccard=p.spon$z.jaccard[lower.tri(p.spon$z.jaccard)],jaccard=spon[lower.tri(spon)],correlation=cor.spon[lower.tri(cor.spon)],"distance"=distmat[lower.tri(distmat)])
  
  jaccard.data.ts = data.frame(z.jaccard=p.ts$z.jaccard[lower.tri(p.ts$z.jaccard)],jaccard=ts[lower.tri(ts)],correlation=cor.ts[lower.tri(cor.ts)],distance=distmat[lower.tri(distmat)])
  
  save(spon,ts,p.spon,p.ts,denoised,spikedata,jaccard.data.spon,jaccard.data.ts,file=paste0('Results/Jac_',temp_file,'.RData'))
  }
  else{
    warning('Issue with mismatched cells in datasets and Metadata_all_monkeys.xlsx -- no distance analysis')
    save(spon,ts,p.spon,p.ts,denoised,spikedata,file=paste0('Results/Jac_',temp_file,'.RData'))
    }
}


