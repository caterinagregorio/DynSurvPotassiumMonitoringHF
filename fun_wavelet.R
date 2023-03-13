
fitwavelet <- function(x,data,lambda) {
  
  load("soglia significatività.RData")
  source("fun_wavelet.R")
  require(tidyverse)
  require(magrittr)
  require(JMbayes2)
  require(WaveletComp)
  require(splines)
  require(plyr)
  require(doParallel)
  
  
    # I select data from subject x
    my.data <- data %>%
      group_by(id) %>%
      filter(id == x) %>%
      ungroup() %>%
      dplyr::select(time, k, k.trend ,k.detrended, dist_m) %>%
      dplyr::rename(date = time)
    


 
  
  
  skip <- FALSE
  
  tryCatch({
    
    # to scale data
    sd <- sd(my.data$k.detrended-mean(my.data$k.detrended))
    my.data$k.detrended.scaled <- scale(my.data$k.detrended,scale=T,center=T)
    
    time0 <- (max(my.data$date)+1):max(my.data$date)+100
    k.detrended0 <- rep(0,length(time0))
    my.data0 <- data.frame(date=time0,
                           k.detrended=k.detrended0,
                           dist_m=rep(1000,length(time0)))
    my.data1 <- bind_rows(my.data,my.data0)
    
    my.w <- analyze.wavelet(
      my.data1,
      "k.detrended",
      loess.span = 0,
      dt = 1,
      dj = 1 / 20,
      lowerPeriod = 2,
      upperPeriod = 365,
      make.pval = F,
      verbose = F
    )
  },
  error = function(e) {
    skip <<- T
  })
  
  if (skip) {
    my.w <- NULL
    my.data$k.est <- my.data$k.trend 
    my.data$id <-rep(x,nrow(my.data)) 
    my.data$b14 <- 0
    my.data$b30 <- 0
    my.data$b90 <- 0
    my.data$b180 <- 0
    my.data$over180 <- 0
    return(my.data)
  } else{
    max.period <- max(my.w$Period)
    
    
    my.w$Power.pval <-
      # Significant p-value only if power > corresponding thresh_quant 
      ifelse(my.w$Power > thresh_quant[1:nrow(my.w$Power), 1:ncol(my.w$Power)], 0.01, 0.99)
    # Significant p-value can exist only when period greater or equal than distance from last measurament
    dist_matrix <- outer(my.w$Period, my.data1$dist_m, ">=")
    my.w$Power.pval <- ifelse(dist_matrix == F, 0.99, my.w$Power.pval)
    
    invisible(capture.output({
      
      
      r14 <- reconstruct(
        my.w,
        sel.period = c(2:pmin(14, max.period)),
        only.sig = T,
        only.ridge = F,
        plot.rec = F,
        lvl=0.09,
        plot.waves = F,
        rescale = F
      )
      
      my.data$b14 <- r14$series[[3]][1:nrow(my.data)]
      my.data$b14 <- replace(my.data$b14, NaN, 0)
      
      
      if (max.period >= 15) {
        
        
        r30 <- reconstruct(
          my.w,
          sel.period = c(15:pmin(30, max.period)),
          only.sig = T,
          only.ridge = F,
          plot.rec = F,
          lvl=0.09,
          plot.waves = F,
          rescale = F
        )
        my.data$b30 <- r30$series[[3]][1:nrow(my.data)]
        
        my.data$b30 <- replace(my.data$b30, NaN, 0)
        
        
        if (max.period > 30) {
          
          
          r90 <- reconstruct(
            my.w,
            sel.period = c(31:pmin(90, max.period)),
            only.sig = T,
            only.ridge = F,
            lvl=0.09,
            plot.rec = F,
            plot.waves = F,
            rescale = F
          )
          
          my.data$b90 <- r90$series[[3]][1:nrow(my.data)]
          my.data$b90 <- replace(my.data$b90, NaN, 0)
          
          
          if (max.period > 90) {
            
            
            
            r180 <- reconstruct(
              my.w,
              sel.period = c(91:180),
              only.sig = T,
              only.ridge = F,
              lvl=0.09,
              plot.rec = F,
              plot.waves = F,
              rescale = F
            )
            
            my.data$b180 <- r180$series[[3]][1:nrow(my.data)]
            
            my.data$b180 <- replace(my.data$b180, NaN, 0)
            
            
            if (max.period > 180) {
              
              rover180 <- reconstruct(
                my.w,
                sel.period = c(181:365),
                only.sig = T,
                only.ridge = F,
                lvl=0.09,
                plot.rec = F,
                plot.waves = F,
                rescale = F
              )
              
              my.data$over180 <- rover180$series[[3]][1:nrow(my.data)]
              
              my.data$over180 <- replace(my.data$over180, NaN, 0)
            } else  {
              my.data$over180 <- 0
            }
            
            
          }else{
            my.data$b180 <- 0
          }}
        else{
          my.data$b90 <- 0
        }
      }else{
        my.data$b30 <- 0
      }
      
      
      
      
      r <- reconstruct(
        my.w,
        sel.period = c(2:365),
        only.sig = T,
        only.ridge = F,
        lvl = 0.09,
        plot.rec = F,
        plot.waves = F)
      
      
    }))
    
  
    my.data$k.est <-replace_na(r$series[[3]][1:nrow(my.data)], 0) + my.data$k.trend 
    my.data$id <-rep(x,nrow(my.data)) 
    
    
    my.data %<>% mutate_all(replace_na, 0)
    
    my.data %<>% mutate_all(replace_na, 0) %>% 
      mutate(
             k.mis=case_when(dist_m==0~k))
    
    
    
    return(my.data)
  }
}



landmark_data_mixed_wavelet <- function(tLM,w,longdata,longdataexp,survdata){
  source("fun_wavelet.R")
  require(plyr)
  require(doParallel)
  require(tidyverse)
  require(magrittr)
  require(rms)
  require(JMbayes)
  require(data.table)
  require(splines)
  
  # One considers all data collected before the landmark time point for 
  # subjects still at risk at the landmark
  LMlong <- longdata %>% 
    arrange(id,time) %>% 
    filter(time<=tLM & fup_days>tLM)
  
  LMlongexp <- longdataexp %>% 
    arrange(id,time) %>% 
    filter(time<=tLM & fup_days>tLM)
  
  # ID at risk at landmark time
  Ri <- LMlong %>% 
    distinct(id) %>% 
    pull()
  
  LMsurv <- survdata %>% 
    arrange(id) %>% 
    filter(id%in%Ri) %>% 
    mutate(time=tLM)
  
  
  # Set up parallel computing
  cl <- makeCluster(2)
  registerDoParallel(cl)
  opts <- list(preschedule=TRUE)
  clusterSetRNGStream(cl, 123)
  LMwaveletfit <- llply(Ri,
                        fitwavelet,
                        data=LMlongexp,
                        .parallel=T,
                        .paropts = list(.options.snow=opts))
  
  LMlongwave <- do.call(bind_rows,LMwaveletfit)
  LMlongwave %<>% dplyr::mutate(over180_d=case_when(abs(over180)>0~1,
                                                    TRUE~0),
                                b180_d=case_when(abs(b180)>0~1,
                                                 TRUE~0),
                                b90_d=case_when(abs(b90)>0~1,
                                                TRUE~0),
                                b30_d=case_when(abs(b30)>0~1,
                                                TRUE~0),
                                b14_d=case_when(abs(b14)>0~1,
                                                TRUE~0),
                                over180=case_when(abs(over180)>0~over180,
                                                  TRUE~0),
                                b180=case_when(abs(b180)>0~b180,
                                               TRUE~0),
                                b90=case_when(abs(b90)>0~b90,
                                              TRUE~0),
                                b30=case_when(abs(b30)>0~b30,
                                              TRUE~0),
                                b14=case_when(abs(b14)>0~b14,
                                              TRUE~0),
                                over180_d2=case_when(over180>0~1,
                                                     over180<0~2,
                                                     TRUE~0),
                                b180_d2=case_when(b180>0~1,
                                                  b180<0~2,
                                                  TRUE~0),
                                b90_d2=case_when(b90>0~1,
                                                 b90<0~2,
                                                 TRUE~0),
                                b30_d2=case_when(b30>0~1,
                                                 b30<0~2,
                                                 TRUE~0),
                                b14_d2=case_when(b14>0~1,
                                                 b14<0~2,
                                                 TRUE~0)) %>%
    dplyr::rename(time=date)
  LMlongwave$tLM <- tLM
  LMlongwave %<>%dplyr::select(id,time,k.trend,k.est,contains("b14"),contains("b30"),contains("b90"),contains("b180"),contains("over180")) 
  
  
  # Last measurament
  LMsurv$k <- LMlongexp %>% 
    filter(time==tLM) %>% 
    dplyr::pull(k)
  
  
  # Long-term 
  LMsurv$k.trend <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(k.trend)
  
  # Long-term + Short term oscillation
  LMsurv$k.est <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(k.est)
  
  # Oscillations at different frequencies
  LMsurv$b14_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b14_d2)
  
  LMsurv$b30_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b30_d2)
  
  LMsurv$b90_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b90_d2)
  
  LMsurv$b180_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b180_d2)
  
  LMsurv$over180_d2 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(over180_d2)
  
  LMsurv$over180 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(over180)
  
  LMsurv$b180 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b180)
  
  LMsurv$b90 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b90)
  
  LMsurv$b30 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b30)
  
  LMsurv$b14 <- LMlongwave %>% 
    filter(time==tLM) %>% 
    dplyr::pull(b14)
  
  LMsurv %<>%mutate(dht=ifelse(fup_days>w+tLM,0,dht),
                    fup_days_w=pmin(w+tLM,fup_days)) 
  
  return(datasurv=LMsurv)
  
  
}


