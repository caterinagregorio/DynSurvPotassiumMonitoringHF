load("data.RData")
load("landmark_data.RData")


source("funcv.R")
source("fun.R")
library(doSNOW)
library(foreach)
library(doParallel)
library(nlme)
library(JMbayes2)
library(survival)
library(pseudo)
library(magrittr)
library(tidyverse)

# Create K=10 random folds
K <- 10
IDs <- unique(LMdata$id)
n <- length(IDs)
set.seed(123)
splits <- split(IDs, sample(rep(seq_len(K), length.out = n)))
tLM <- c(365,725,1085)

index_visit$age_stand <- scale(index_visit$age)[,1]
LMdata$age_stand <- scale(LMdata$age)[,1]
LMdata$iperk <- ifelse(LMdata$k>5,1,0)
LMdata$ipok <- ifelse(LMdata$k<3.5,1,0)

LMdata$k.dtrend_stand <- scale(LMdata$k.dtrend*100)
LMdata$k.trend_stand <- scale(LMdata$k.trend)
LMdata$years <- cut(LMdata$Tstart,breaks = c(365,365*2,365*3,365*4,365*5),labels = c("second","third","fourth","fifth"),include.lowest = T)

potassio %<>%left_join(index_visit[,c(1,58)]) 

pred <- list(K)
for (k in 1:10){
  gc()
  pred[[k]] <- cv_fun(k, 
                       splits,
                       LMdata,
                       potassio,
                       index_visit,
                       w=180)
  print(k)
  print(pred[[k]])
}


##############################################################################################

perf_data <- do.call("rbind",pred)
perf_summ<- perf_data%>% 
  dplyr::group_by(tLM,model) %>% 
  dplyr::summarize(aucm=mean(auc),
                   se.auc=sd(auc),
                   bsm=mean(bs),
                   se.bs=sd(bs))


save(perf_data,perf_summ,pred,file="Results cv.RData")


fred <- Vectorize(function(x,nullm){(nullm-x)/nullm*100},"x")

perf_summ$tLM <- round(perf_summ$tLM,0)
perf_summ$red <- 0 
perf_summ$red[1:6] <- fred(perf_summ$bsm[1:6],perf_summ$bsm[1])
perf_summ$red[7:12] <- fred(perf_summ$bsm[7:12],perf_summ$bsm[7])
perf_summ$red[13:18] <- fred(perf_summ$bsm[13:18],perf_summ$bsm[13])
perf_summ$bsm <- perf_summ$bsm*100
perf_summ$red <- round(perf_summ$red,1)
table9 <- perf_summ[,c(1,2,3,5,7)]
table9

colnames(table9) <- c("Prediction time","Model","AUC(6 months)","Brier Score- 6 months","% Reduction Brier Score")
table9$Model <- factor(table9$Model,labels=c("Null","Landmark LOCF1", "Landmark LOCF2","Landmark Mixed","Landmark Wavelet","Joint Model"))
View(table9)
print(xtable::xtable(table9,digits=3), include.rownames=FALSE)
