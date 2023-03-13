cv_fun <- function(k,index,datalm , datalong,datasurv,w) {
  # source("fun.R")
  
  require(data.table)
  require(plyr)
  require(survival)
  require(splines)
  require(doSNOW)
  require(foreach)
  require(JMbayes2)
  library(tidyverse)
  library(magrittr)
  
  
  
  
  # Define landmark times
  tLM <- seq(365, 365 * 5, by = 90)
  
  datasurv_test <- datasurv[datasurv$id%in%index[[k]],]
  datalong_test <- datalong[datalong$id%in%index[[k]],]
  datalm_test <- datalm[datalm$id%in%index[[k]],]
  
  datasurv_train <- datasurv[!(datasurv$id%in%index[[k]]),]
  datalong_train <- datalong[!(datalong$id%in%index[[k]]),]
  datalm_train <- datalm[datalm$id%in%index[[k]],]
  
  datasurv_train$age_stand <- scale(datasurv_train$age)
  datasurv_test$age_stand <- scale(datasurv_test$age)
  
  
  datalm_train$age_stand <- scale(datalm_train$age)
  datalm_train$k.dtrend_stand <- scale(datalm_train$k.dtrend * 100)
  datalm_train$k.trend_stand <- scale(datalm_train$k.trend)
  datalm_train$k_stand <- scale(datalm_train$k)
  datalm_train$hypo <- ifelse(datalm_train$k<3.5,1,0)
  datalm_train$hyper <- ifelse(datalm_train$k>5,1,0)
  
  datalm_test$age_stand <- scale(datalm_test$age)
  datalm_test$k.dtrend_stand <- scale(datalm_test$k.dtrend * 100)
  datalm_test$k.trend_stand <- scale(datalm_test$k.trend)
  datalm_test$k_stand <- scale(datalm_test$k)
  datalm_test$hypo <- ifelse(datalm_test$k<3.5,1,0)
  datalm_test$hyper <- ifelse(datalm_test$k>5,1,0)
  
  datasurv_train %<>%arrange(id)
  datalong_train %<>%arrange(id,time)
  
  print("test")
  ########
  # LOCF #
  ########
  
  
  ## Simple (ipl) locf
  
  LMsupercox_locf <-
    coxph(
      Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3 +
        tv(k,ord = 4,knots = c(2,3,4,6,7)) + strata(time) + cluster(id),
      data = datalm_train
    )
  
  
  ########
  # LOCF 2#
  ########
  
  
  ## Simple (ipl) locf
  
  LMsupercox_locf2 <-
    coxph(
      Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3 +
        hypo+hyper + strata(time) + cluster(id),
      data = datalm_train
    )
  
  #########
  # MIXED #
  #########
  
  
  ## Simple (ipl) l
  LMsupercox_mixed <-
    coxph(
      Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34 + lvef_classi +  somma_comor_nc3 +
        tv(k.trend,ord = 4,knots = c(2,3,4,6,7)) + strata(time) + cluster(id),
      data = datalm_train,
      method = "breslow"
    )
  
  
  ####################
  # MIXED & WAVELETS #
  ####################
  
  
  ## Simple (ipl)
  LMsupercox_mixed_wave <-
    coxph(
      Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34   + lvef_classi +  somma_comor_nc3 +
        tv(k.trend,ord = 4,knots = c(2,3,4,6,7))+
      
        b14_d2 +
        b30_d2 +
        b90_d2 +
        b180_d2 +
        over180_d2 +
        strata(time) + cluster(id),
      data = datalm_train,
      method = "breslow"
    )
  
  
  ###############
  # JOINT MODEL #
  ###############
  
  fp <- lme(k ~ gender + age_stand + irc + nyha_34 + tv(time, ord=4,knots = c(0, 218,
                                                                        752, 1447,4040)),
            data = potassio,
            random = ~tv(k,ord = 4,knots = c(2,3,4,6,7)) | id,
            control = lmeControl(opt = 'optim'))
  
  CoxFit <-
    coxph(
      Surv(fup_days, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3,
      data = index_visit,
      x = TRUE
    )
  
  
  
  fitJoint <-jm(CoxFit,
                 fp,
                 time_var = "time",
                 functional_forms = ~ value(tv(k,ord = 4,knots = c(2,3,4,6,7))),
                 control = list(n_iter = 50000,
                                n_burnin = 10000))
  print("ok models")
  
  
  
  
  
  
  # Obtain performance measures
  
  
  # null model
  p_null <- llply(
    c(365,725,1085),
    perfnull,
    data = datalm_test,
    w = w,
    .parallel=F
    # .paropts = list(.options.snow=opts)
  )
  
  print("ok null")
  
  # LOCF
  p_locf <- lapply(
    c(365,725,1085),
    perfLM,
    object = LMsupercox_locf,
    data = datalm_test,
    w = w
  )
  
  
  print("ok locf")
  
  
  # LOCF 2
  p_locf2 <- lapply(
    c(365,725,1085),
    perfLM,
    object = LMsupercox_locf2,
    data = datalm_test,
    w = w
  )
  
  
  print("ok locf2")
  
  #MIXED
  p_mixed <- lapply(
    c(365,725,1085),
    perfLM,
    object = LMsupercox_mixed,
    data = datalm_test,
    w = w
  )
  
  
  
  print("ok mixed")
  
  # MIXED-WAVELETS
  p_mixed_wave <- lapply(
    c(365,725,1085),
    perfLM,
    object = LMsupercox_mixed_wave,
    data = datalm_test,
    w =w
  )
  
  
  
  print("ok wave")
  
  # JOINT
  #debug(perfJM)
  p_joint <-
    lapply(
      c(365,725,1085),
      perfJM,
      object = fitJoint,
      data = datalong_test,
      dataevent=datasurv_test,
      w = w
    )
  
  print("ok joint")
  
  tranf_perf <- function(obj){
    extract<- function(landmark) lapply(obj, function(y) y[[landmark]])
    
    res<- sapply(1:9, function(my.year)
      do.call(rbind, extract(my.year)))   
    
    return(res)
    
  }
  
  p_null_long <- tranf_perf(p_null)
  p_locf_long <- tranf_perf(p_locf)
  p_locf2_long <- tranf_perf(p_locf2)
  p_mixed_long <- tranf_perf(p_mixed)
  p_mixed_wave_long <- tranf_perf(p_mixed_wave)
  p_joint_long <- tranf_perf(p_joint)
  
  
  
  
  data_perf <- data.frame(tLM=rep(c(365,725,1085) / 365, 6),
                          bs=c(p_null_long[[1]],p_locf_long[[1]],p_locf2_long[[1]],p_mixed_long[[1]],p_mixed_wave_long[[1]],p_joint_long[[1]]),
                          se.bs=c(p_null_long[[2]],p_locf_long[[2]],p_locf2_long[[2]],p_mixed_long[[2]],p_mixed_wave_long[[2]],p_joint_long[[2]]),
                          auc=c(NA,NA,NA,p_locf_long[[4]],p_locf2_long[[4]],p_mixed_long[[4]],p_mixed_wave_long[[4]],p_joint_long[[4]]),
                          se.auc=c(NA,NA,NA,p_locf_long[[5]],p_locf2_long[[5]],p_mixed_long[[5]],p_mixed_wave_long[[5]],p_joint_long[[5]]),
                          model=rep(0:5,each=3),
                          k=rep(k,3*6))
  
  
  # data_perf <- data.frame(tLM=rep(c(365,725,1085) / 365, 5),
  #                         bs=c(p_null[[1]],p_locf_long[[1]],p_locf2_long[[1]],p_mixed_long[[1]],p_mixed_wave_long[[1]]),
  #                         se.bs=c(p_null[[2]],p_locf_long[[2]],p_locf2_long[[2]],p_mixed_long[[2]],p_mixed_wave_long[[2]]),
  #                         auc=c(p_null[[4]],p_locf_long[[4]],p_locf2_long[[4]],p_mixed_long[[4]],p_mixed_wave_long[[4]]),
  #                         se.auc=c(p_null[[5]],p_locf_long[[5]],p_locf2_long[[5]],p_mixed_long[[5]],p_mixed_wave_long[[5]]),
  #                         model=rep(0:4,each=3),
  #                         k=rep(k,3*5))
  
  
  
  
  return(data_perf)
  
}


perfJM <- function(Tstart,w,data,object,dataevent){
  require(timeROC)
  require(JMbayes) 
  
  n <- length(unique(data$id))
  datasel <- data[data$time<=Tstart & data$fup_days>Tstart ,]
  Thoriz <- Tstart+w-0.5
  datasel$dht_w <- datasel$dht
  datasel$fup_days_w <- datasel$fup_days
  ids <- unique(datasel$id)
  
  datasel$fup_days <- Tstart
  datasel$dht <- 0
  datasel <- as.data.frame(datasel)
  dataeventsel <- dataevent %>% filter(id%in%ids)
  
  pred <- predict(object,
                  process = "event",
                  type_pred = "response",
                  type="subject_specific",
                  newdata =datasel ,
                  times = Thoriz,
                  #FtTimes = last_m,
                  return_newdata =  F)
  
  surv <- as.numeric(1-pred$pred)
  surv <- surv[seq(2,length(surv),by=2)]
  
  datasel %<>% distinct(id,.keep_all = T)  
  
  print(nrow(dataeventsel))
  print(length(surv))
 
  roc <- timeROC(dataeventsel$fup_days,
                 delta=dataeventsel$dht,
                 marker=1-surv,
                 cause=1,
                 weighting="marginal",
                 times=c(Thoriz),
                 iid=TRUE)
  
  
  
  bs <- BS(timepoints=c(Thoriz),
           times=dataeventsel$fup_days,
           status=dataeventsel$dht,
           pred=as.matrix(1-surv),
           cause=1)
  
  
  est.bs <- bs$BS # BS estimate
  se.bs <- bs$sd  # BS s.e. estimate
  Mat.iid.BS <-  c(rep(0,n-roc$n),as.vector(bs$iid))/(roc$n/n)  # BS iid decomposition
  est.auc <- roc$AUC[2] # AUC estimate
  sd.auc <- roc$inference$vect_sd_1[2]  # AUC s.e. estimate
  Mat.iid.AUC <- c( rep(0,n-roc$n),roc$inference$mat_iid_rep_1[,2]) # AUC iid decomposition   
  matrixstat <- roc$Stats[2,] # proportions of cases, controls, censored subjects within the prediction window and survivors at s+t
  CumInc <- roc$CumulativeIncidence[2]  # marginal cumulative incidence  = P(T<s+t,eta=1|T>s)
  Surv <- roc$survProb[2] # marginal survival probability = P(T>s+t|T>s)
  return(list(est.bs,se.bs,Mat.iid.BS,est.auc,sd.auc,Mat.iid.AUC,matrixstat,CumInc,Surv))
  
}



# Expected Brier score estimator iid representation 

# {{{ Inpus ;
# pred : prediction, e.g. P(T<t,cause=1|history) A MATRIX 
# timepoints : vector of time points for which we aim to compute the iid representations
# times : the vector of observed time-to-event
# status : 1=uncensored, 0=censored
# cause : cause for which we aim to compute the expected Brier score estimator
# }}}

# {{{ Outputs :
# matrix with iid representation for all time points
# }}}


# {{{ Functions 
# main function
BS <- function(timepoints,times,status,pred,cause=1,compute.iid=TRUE){ 
  
  require(pec)
  start_computation_time <- Sys.time()
  # define useful objects
  n <- length(times)
  n_times <- length(timepoints)
  timepoints <- timepoints[order(timepoints)]
  times_names <- paste("t=",timepoints,sep="")
  # output initialisation 
  BS <- rep(NA,n_times)
  CumInci <- rep(NA,n_times)
  surv <- rep(NA,n_times)
  Stats <- matrix(NA,nrow=n_times,ncol=4)
  hit1_all <- matrix(NA,nrow=n,ncol=n_times)
  hit2_all <- matrix(NA,nrow=n,ncol=n_times)
  epsilon_i <- matrix(NA,nrow=n,ncol=n_times)
  #adds name to outputs
  names(BS) <- times_names
  names(CumInci) <- times_names
  names(surv) <- times_names
  colnames(Stats) <- c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats) <- times_names
  colnames(epsilon_i) <- times_names
  colnames(hit1_all) <-  times_names
  colnames(hit2_all)  <- times_names 
  # we need to order to use the ipcw() function of the pec package
  #browser()
  order_T <- order(times)
  times <-  times[order_T]
  delta  <-  status[order_T]
  pred <-  pred[order_T,,drop=FALSE]
  #compute KM weights
  weights <- pec::ipcw(Surv(failure_time,status)~1,
                       data=data.frame(failure_time=times,status=as.numeric(delta!=0)),
                       method="marginal",times=timepoints,subjectTimes=times,subjectTimesLag=1)
  Mat_data <- cbind(times,delta,as.numeric(delta==0))
  colnames(Mat_data) <- c("T","delta","indic_Cens")
  # computate weights of cases
  Weights_cases_all <- 1/(weights$IPCW.subjectTimes*n)
  # compute KM censoring estimator iid representation
  if (compute.iid){ MatInt0TcidhatMCksurEff <- Compute.iid.KM(times,delta!=0)}
  # loop on all time points
  for(t in 1:n_times){
    Cases <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]==cause)
    Controls_1 <- (Mat_data[,"T"]> timepoints[t] )
    Controls_2 <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    # compute weights
    Weights_controls_1 <- rep(1/(weights$IPCW.times[t]*n),times=n)
    Weights_cases <- Weights_cases_all
    Weights_controls_2 <- Weights_cases_all
    Weights_cases[!Cases] <- 0
    Weights_controls_1[!Controls_1] <- 0
    Weights_controls_2[!Controls_2] <- 0   
    #compute outputs
    CumInci[t] <- c(sum(Weights_cases))
    surv[t] <- c(sum(Weights_controls_1))
    Stats[t,] <- c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2)) 
    hit1_all[,t] <- (Weights_controls_1*((pred[,t])^2))*n
    hit2_all[,t] <- (Weights_cases*((1-pred[,t])^2) + Weights_controls_2*((pred[,t])^2))*n
    BS[t] <- (sum(hit1_all[,t]) +sum(hit2_all[,t]))/n
    if (compute.iid){
      # compute 
      Int0tdMCsurEffARisk <- MatInt0TcidhatMCksurEff[max(which(Mat_data[,"T"]<=timepoints[t])),]   
      #browser()
      epsilon_i[,t] <- hit1_all[,t]+hit2_all[,t]-BS[t] + mean(hit1_all[,t])*Int0tdMCsurEffARisk +  colMeans(MatInt0TcidhatMCksurEff*hit2_all[,t])
    }
  } 
  #compute mean and sd of iid representation
  sd_all <- rep(NA,n_times)
  mean_all <- rep(NA,n_times)
  if (compute.iid){sd_all <- apply(epsilon_i,2,sd)/sqrt(n)
  mean_all <- apply(epsilon_i,2,mean)}
  #compute a table to print 
  print.tab <- cbind(Stats,BS,sd_all,mean_all)
  colnames(print.tab) <- c(colnames(Stats),"BS","sd","mean_iid")
  #compute the computation time
  stop_computation_time <- Sys.time()
  computation_time=difftime(stop_computation_time,start_computation_time,units="secs")
  
  out <- list(BS=BS,iid=epsilon_i,sd=sd_all,res=(hit1_all+hit2_all),
              CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,print.tab=print.tab,timepoints=timepoints,
              computation_time=difftime(stop_computation_time,start_computation_time,units="secs"))
  class(out) <- "ipcwEBS"
  out 
}
# print function
print.ipcwEBS <- function(x,No.lines=5,...){
  tab_ou_print <- round(cbind(x$Stats,x$BS*100,x$sd*100),2)
  colnames(tab_ou_print) <- c("Cases","survivors","Other events","Censored","EBS(*100)","se(*100)")
  l <- length(x$timepoints)
  if(l<=No.lines){  print(tab_ou_print) }
  else{print(tab_ou_print[unique(round(quantile(1:length(x$timepoints),probs=seq(0,1,length.out=No.lines)),0)),])}
  cat("\n")
  cat("Total computation time :",round(x$computation_time,2)," secs. \n")
}
# }}}

# function to compute iid decomposition of KM estimator of the censoring survival distribution

# {{{ Input
#times : observed time
#status : 1 if non censored, 0 if censored
# }}}

# {{{ Output
#iid.mat : matrix of the iid representation of KM for all time in vector times
# }}}

# {{{ Function
Compute.iid.KM <- function(times,status){
  #browser()
  times <- times[order(times)]
  status <- status[order(times)] 
  n <- length(times)
  mat.data<-cbind(times,as.numeric(status==0))
  colnames(mat.data)<-c("T","indic.Cens")
  # compute the empirical survival function corresponding to the counting process 1(\tilde{eta}=0, \tilde{T}<=t)
  hatSdeltaCensTc<-1-cumsum(mat.data[,c("indic.Cens")])/n  
  # Build the matrix required for computing  dM_C(u) for all time u (all observed times \tilde{T}_i)
  temp1 <- cbind(mat.data[,c("T","indic.Cens")],1-(1:n)/n,hatSdeltaCensTc)
  temp1 <- rbind(c(0,0,1,1),temp1) # Add the first row corresponding to time t=0
  colnames(temp1)<-c("T","indic.Cens","hatSTc","hatSdeltaCensTc")
  # compute hazard function of the censoring
  lambdaC<-(temp1[-1,"indic.Cens"])/(n:1)  
  # Add the column of the hazard function of the censoring (equal to 0 at time t=0)
  temp1<-cbind(temp1,c(0,lambdaC))
  colnames(temp1)[ncol(temp1)]<-"lambdaC"
  # Cumulative hazard of censoring
  LambdaC<-cumsum(lambdaC)         
  # Add the column of the cumulative hazard function of the censoring (equal to 0 at time t=0)
  temp1 <- cbind(temp1,c(0,LambdaC))
  colnames(temp1)[ncol(temp1)]<-"LambdaC"
  temp2<-temp1[-1,]
  # compute  martingale of censoring \hat{M}_{C_i}(u) for all time u (all observed times \tilde{T}_i) using previous matrix
  # We obtain a matrix. Each column contains the vector of M_{C_i}(\tilde{T}_j) for  all j.
  hatMC<-matrix(NA,n,n)
  for (i in 1:n){
    hatMC[,i] <-temp2[i,2]*as.numeric(temp2[i,1]<=temp2[,"T"])- c(temp2[0:i,"LambdaC"], rep(temp2[i,6],(n-i)))
  }  
  # In order to draw martingale paths
  #matplot(mat.data[,"T"],hatMC,type="l")
  #lines(mat.data[,"T"],rowMeans(hatMC),lwd=5)  
  # Compute d \hat{M}_{C_i} (u) for all time u (all observed times \tilde{T}_i)
  dhatMC<-rbind(hatMC[1,],hatMC[-1,]-hatMC[-nrow(hatMC),])
  # Compute d \hat{M}_{C_i} (u)/(S_{\tilde{T}}(u)) for all time u (all observed times \tilde{T}_i)
  # We need this for integrals in the martingale representation of the Kaplan-Meier estimator of the censoring survival function
  # function to divide d \hat{M}_{C_i} (u) by (S_{\tilde{T}}(u))
  MulhatSTc<-function(v){
    n <- length(v)
    v/c(1,1-(1:(n-1))/n)      # c(1,1-(1:(n-1))/n) is the at risk probability (S_{\tilde{T}}(u))
  }
  # apply the function for each column (corresponding to the
  # vector M_{C_i}(u)  for all time u (all observed times \tilde{T}_i), 
  # time \tilde{T}_i corresponds to the i-th row of the matrix)
  dhatMCdivST<-apply(dhatMC,2,MulhatSTc)
  # Compute \int_0^{\tilde{T}_j} d{ \hat{M}_{C_l} (u) } / (S_{\tilde{T}}(u)) for each subject l, we compute for all time \tilde{T}_j.
  # l=column, j=row
  MatInt0TcidhatMCksurEff<-apply(dhatMCdivST,2,cumsum)  # (Remark : on of the row corresponds to the previous step...) 
  colnames(MatInt0TcidhatMCksurEff)<-paste("M_{C_",1:length(times),"}",sep="")
  rownames(MatInt0TcidhatMCksurEff)<-times  
  return(MatInt0TcidhatMCksurEff)  
}




perfLM <- function(Tstart,w,data,object){
  require(timeROC)
  require(survival)
  require(splines)
  
  n <- length(unique(data$id))
  datasel <- data[data$time==Tstart,]
  Thoriz <- Tstart+w-0.5
  year <- as.character(unique(datasel$years))
  year2 <- setdiff(c("second", "third" , "fourth" ,"fifth" ),year)
  Tstart2 <- sample(unique(data$Tstart[data$years==year2[1]]),1)
  Tstart3 <- sample(unique(data$Tstart[data$years==year2[2]]),1)
  Tstart4 <- sample(unique(data$Tstart[data$years==year2[3]]),1)
  datasel2 <- data[data$Tstart%in%c(Tstart,Tstart2,Tstart3,Tstart4),]
  sfit <- survfit(object, newdata = datasel2)
  pred <- summary(sfit, times = Thoriz,extend=T)
  surv <- pred$surv[which(datasel2$time==Tstart)]
  datasel$pred <- pred$surv[which(datasel2$time==Tstart)]
  
  
  
  roc <- timeROC(T=datasel$fup_days_w,
                 delta=datasel$dht,
                 marker=1-surv,
                 cause=1,
                 weighting="marginal",
                 times=c(Thoriz),
                 iid=TRUE)
  
  bs <- BS(timepoints=c(Thoriz),
           times=datasel$fup_days_w,
           status=datasel$dht,
           pred=as.matrix(1-surv),
           cause=1)
  
  
  est.bs <- bs$BS # BS estimate
  se.bs <- bs$sd  # BS s.e. estimate
  Mat.iid.BS <-  c(rep(0,n-roc$n),as.vector(bs$iid))/(roc$n/n)  # BS iid decomposition
  est.auc <- roc$AUC[2] # AUC estimate
  sd.auc <- roc$inference$vect_sd_1[2]  # AUC s.e. estimate
  Mat.iid.AUC <- c( rep(0,n-roc$n),roc$inference$mat_iid_rep_1[,2]) # AUC iid decomposition 
  matrixstat <- roc$Stats[2,] # proportions of cases, controls, censored subjects within the prediction window and survivors at s+t
  CumInc <- roc$CumulativeIncidence[2]  # marginal cumulative incidence  = P(T<s+t,eta=1|T>s)
  Surv <- roc$survProb[2] # marginal survival probability = P(T>s+t|T>s)
  return(list(est.bs,se.bs,Mat.iid.BS,est.auc,sd.auc,Mat.iid.AUC,matrixstat,CumInc,Surv))
}



perfnull <- function(Tstart,w,data){
  require(timeROC)
  require(survival)
  require(splines)
  
  datasel <- data[data$time==Tstart,]
  Thoriz <- Tstart+w-0.5
  nLM <- n <- nrow(datasel)
  KMw <- summary(survfit(Surv(fup_days, dht) ~ 1, data=datasel), times=Thoriz)$surv
  pseudovals <- as.numeric(pseudosurv(datasel$fup_days_w, datasel$dht, Thoriz)$pseudo)
  datasel$null <- (nLM*KMw - pseudovals) / (nLM-1)

  roc <- timeROC(T=datasel$fup_days,
                 delta=datasel$dht,
                 marker=1-datasel$null,
                 cause=1,
                 weighting="marginal",
                 times=c(Thoriz),
                 iid=TRUE)
  
  bs <- BS(timepoints=c(Thoriz),
           times=datasel$fup_days,
           status=datasel$dht,
           pred=as.matrix(1-datasel$null),
           cause=1)
  
  
  est.bs <- bs$BS # BS estimate
  
  se.bs <- bs$sd  # BS s.e. estimate
  Mat.iid.BS <-  c(rep(0,n-roc$n),as.vector(bs$iid))/(roc$n/n)  # BS iid decomposition
  est.auc <- roc$AUC[2] # AUC estimate
  sd.auc <- roc$inference$vect_sd_1[2]  # AUC s.e. estimate
  Mat.iid.AUC <- c( rep(0,n-roc$n),roc$inference$mat_iid_rep_1[,2]) # AUC iid decomposition 
  matrixstat <- roc$Stats[2,] # proportions of cases, controls, censored subjects within the prediction window and survivors at s+t
  CumInc <- roc$CumulativeIncidence[2]  # marginal cumulative incidence  = P(T<s+t,eta=1|T>s)
  Surv <- roc$survProb[2] # marginal survival probability = P(T>s+t|T>s)
  return(list(est.bs,se.bs,Mat.iid.BS,est.auc,sd.auc,Mat.iid.AUC,matrixstat,CumInc,Surv))
}


PercRedPredErr <- function(Score1, Score0) (Score0 - Score1) / Score0
cloglog <- function(x) log(-log(x))
