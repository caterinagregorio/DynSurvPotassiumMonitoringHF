
# tab_long <- function(fit,digits){
#   
#     s <- summary(fit)$`CoefTable-Long`
#     table <- rownames_to_column( as.data.frame(s), "Variable")
#     table$lower <- s[,5]
#     table$upper <- s[,6]
#     table$ci <- paste(round(table$lower,digits),round(table$upper,2),sep=";")
#     table$Value <- round(table$Value,digits)
#     table$Std.Error <- round(table$Std.Dev,digits)
#     table$`p-value` <- round(table$P,4)
#     table <- table %>% select(Variable,Value,Std.Error,ci,`p-value`)
#   
#   return(table)
# }


tab_long <- function(fit,digits,joint=F){
  
  if(joint){
    s <- summary(fit)$Outcome1[-nrow(summary(fit)$Outcome1),]
    table <- rownames_to_column( as.data.frame(s), "Variable")
    table$lower <- s[,3]
    table$upper <- s[,4]
    table$ci <- paste(round(table$lower,digits),round(table$upper,2),sep=";")
    table$Value <- round(table$Mean,digits)
    table$Std.Error <- round(table$StDev,digits)
    table$`p-value` <- round(table$P,4)
    table <- table %>% select(Variable,Value,Std.Error,ci,`p-value`)
    table$`p-value` <- ifelse(table$`p-value` <0.01,"<0.01",as.character(round(table$`p-value`,2))) 
  }else{
    s <- summary(fit)
    table <- rownames_to_column( as.data.frame(s$tTable), "Variable")
    table$lower <- table[,2]-1.96*table[,3]
    table$upper <-table[,2]+1.96*table[,3]
    table$ci <- paste(round(table$lower,digits),round(table$upper,2),sep=";")
    table$Value <- round(table$Value,digits)
    table$Std.Error <- round(table$Std.Error,digits)
    table$`p-value` <- round(table$`p-value`,4)
    table <- table %>% select(Variable,Value,Std.Error,ci,`p-value`)
    table$`p-value` <- ifelse(table$`p-value` <0.01,"<0.01",as.character(round(table$`p-value`,2))) 
  }
  
  
  return(table)
}



tab_surv <- function(fit,digits,joint=T){
  
  if(joint){
  s <- summary(fit)$Survival
  table <- rownames_to_column( as.data.frame(s), "Variable")
  table$lower <- s[,3]
  table$upper <- s[,4]
  table$ci <- paste(round(table$lower,digits),round(table$upper,2),sep=";")
  table$Value <- round(table$Mean,digits)
  table$Std.Error <- round(table$StDev,digits)
  table$`p-value` <- round(table$P,4)
  table <- table %>% select(Variable,Value,Std.Error,ci,`p-value`)
  }
  
  else{
    s <- summary(fit)
    table <- s$coefficients[,c(1,4,6)]
    table <- data.frame(table)
    table$lower <- round(table[,1]-1.96*table[,2],digits)
    table$upper <- round(table[,1]+1.96*table[,2],digits)
    table$Value<- round(table$coef,digits)
    table$Std.Error  <- round(table$robust.se,digits)
    table$`p-value` <- round(table$Pr...z..,4)
    table$ci <- paste(table$lower,table$upper,sep = ";")
    table <- table[,-c(3,4,5)]
    table <- rownames_to_column( as.data.frame(table), "Variable")
    table <- table %>% select(Variable,Value,Std.Error,ci,`p-value`)
  }
  
  table$`p-value` <- ifelse(table$`p-value` <0.01,"<0.01",as.character(round(table$`p-value`,2))) 
  return(table)
}


ggrelrisk <- function(fit,nbasis,b,ref,lab="Biomarker",slope=F,joint=T,var,terms=NULL,mean=0,rect=T){
  
  if(joint){
    s<- summary(fit)
    if(slope)index <- grepl("slope",rownames(s$Survival))
    else{index <- grepl("value",rownames(s$Survival))}
    coef <- s$Survival[index,1]
    lower <- s$Survival[index,3]
    upper <- s$Survival[index,4]
    
    
      relrisk <- coef*tv(b,ord = 4,knots =  c(2,3,4,6,7))
      relrisk_lower <- lower*tv(b,ord = 4,knots =  c(2,3,4,6,7))
      relrisk_upper <- upper*tv(b,ord = 4,knots =  c(2,3,4,6,7))
      
      ref_index <- min(which(round(b,1)==ref))
      
      
      
      dat <- data.frame(b=b,
                        hr=relrisk-relrisk[ref_index],
                        lower=relrisk_lower-relrisk_lower[ref_index],
                        upper=relrisk_upper-relrisk_upper[ref_index])
      
    
    
  }else{
    
    s <- summary(fit)
    
    coef <- s$coefficients[terms,1]
    lower <- s$coefficients[terms,1]+s$coefficients[terms,3]*1.96
    upper <- s$coefficients[terms,1]-s$coefficients[terms,3]*1.96
    
    
    relrisk <- coef*tv(b,ord = 4,knots =  c(2,3,4,6,7))
    relrisk_lower <- lower*tv(b,ord = 4,knots =  c(2,3,4,6,7))
    relrisk_upper <- upper*tv(b,ord = 4,knots =  c(2,3,4,6,7))
    
    ref_index <- min(which(round(b,1)==ref))
    
    
    
    dat <- data.frame(b=b,
                      hr=relrisk-relrisk[ref_index],
                      lower=relrisk_lower-relrisk_lower[ref_index],
                      upper=relrisk_upper-relrisk_upper[ref_index])
    
  }
  
 
  
  
  # dat <- data.frame(b=b,
  #                   hr=relrisk,
  #                   lower=relrisk_lower,
  #                   upper=relrisk_upper)
  
  
 
  
  # lower_bound <- min(dat$b[dat$lower<1])
  # upper_bound <- max(dat$b[dat$lower<1])
  # 
  max <- max(dat$upper)
  min <- min(dat$lower)
  
  range_lower <- min(dat$b[which(dat$lower<0)])
  range_upper <- max(dat$b[which(dat$lower<0)])

  
 if(rect){
   gg <- ggplot(dat)+
     geom_rect(aes(xmin=range_lower,xmax=range_upper,ymin=-Inf,ymax=Inf),fill="lightblue",alpha=0.1)+
     # geom_vline(aes(xintercept=ref),linetype="dashed")+
     geom_ribbon(aes(b,ymin=lower,ymax=upper),alpha=0.3)+
     geom_line(aes(b,hr),size=1)+
     scale_x_continuous(lab,limits = c(min(dat$b),max(dat$b)),breaks = c(2.5,3.5,4.5,5.5,6.5,7.5))+
     scale_y_continuous("Relative Risk",limits = c(-0.5,1))+
     theme_bw()
   
 }
    else{
      gg <- ggplot(dat)+
        # geom_rect(aes(xmin=range_lower,xmax=range_upper,ymin=-Inf,ymax=Inf),fill="lightblue",alpha=0.1)+
        # geom_vline(aes(xintercept=ref),linetype="dashed")+
        geom_ribbon(aes(b,ymin=lower,ymax=upper),alpha=0.3)+
        geom_line(aes(b,hr),size=1)+
        scale_x_continuous(lab,limits = c(min(dat$b),max(dat$b)),breaks = c(2.5,3.5,4.5,5.5,6.5,7.5))+
        scale_y_continuous("Relative Risk",limits = c(-0.75,1))+
        theme_bw()
    }
  
  
  
  return(gg)
  
}



gg_periodogram <- function(id,data){
  require(WaveletComp)
  
  datasel <- data[data$id==id,]
  datasel$date <- datasel$time
  my.w <- analyze.wavelet(datasel, 
                          "k.detrended",
                          loess.span = 0,
                          dt = 1, 
                          dj=1/20,
                          lowerPeriod = 2,
                          upperPeriod = 365,
                          make.pval = F,
                          verbose=F)
  
  im <- wt.image(my.w)
  
  
  datapower <- data.frame(time=rep(my.w$axis.1,each=151),
                          period=rep(my.w$axis.2,times=length(my.w$axis.1)),
                          power=c(my.w$Power))
  
  
  
  newcol <- colorRampPalette(c("white","#225560","#7EA8BE","#FCC05F","#BB4430"))
  ncols <- 1000
  col <- newcol(ncols)#apply the function to get 100 colours
  
  gg <- ggplot(datapower)+
    geom_tile(aes(time/365,period,fill=power))+
    scale_fill_gradientn("Power",colours = col)+
    scale_x_continuous("Time (Years)",breaks = seq(0,max(data$time/365),by=1))+
    scale_y_continuous("Period",breaks = log2(c(2,4,7,14,30,90,180,365)),labels =c(2,4,7,14,30,90,180,365) )+
    # scale_height_continuous(guide="none")+
    ggtitle("Periodogram")+
    theme_classic()+
    theme(text = element_text(face="bold",
                              size=16),
          legend.position = NULL)
  
  return(list(plot=gg,data=datasel))
  
}


plot_osc_all <- function(data){
  
  id <- unique(data$id)
  data <- fitwavelet(id,data,deriv = F)
  
  
  gg<- ggplot(data)+
    geom_line(aes(date/365,b14+k.trend,group=as.factor(id),col="2-14 days"),size=1)+
    geom_line(aes(date/365,b30+k.trend,group=as.factor(id),col="15-30 days"),size=1)+
    geom_line(aes(date/365,b90+k.trend,group=as.factor(id),col="31-90 days"),size=1)+
    geom_line(aes(date/365,b180+k.trend,group=as.factor(id),col="91-180 days"),size=1)+
    geom_line(aes(date/365,over180+k.trend,group=as.factor(id),col="181-365 days"),size=1)+
    geom_line(aes(date/365,k.trend,group=as.factor(id)),size=1,linetype="dashed",color="black")+
    geom_point(aes(date/365,k.mis),shape=8)+
    scale_y_continuous("mmol/L")+
    scale_color_manual("Short-term oscillations",values=c("2-14 days"="#EF476F",
                                                          "15-30 days"="#FFD166",
                                                          "31-90 days"="#06D6A0",
                                                          "91-180 days"="#118AB2",
                                                          "181-365 days"="#073B4C"))+
    scale_linetype("",labels=c("","Current mean value"))+
    scale_x_continuous("Time (Years)",breaks = seq(0,max(data$date/365),by=1))+
    theme_light()+
    theme(text = element_text(face = "bold",
                              size=16),
          legend.position = "bottom")
  return(gg)
  
  
}


plot_osc <- function(data){
  
  id <- unique(data$id)
  data <- fitwavelet(id,data,deriv = F)
  
  
  gg<- ggplot(data)+
    geom_line(aes(date/365,b14,group=as.factor(id),col="2-14 days"),size=1,alpha=1)+
    geom_line(aes(date/365,b30,group=as.factor(id),col="15-30 days"),size=1,alpha=0.6)+
    geom_line(aes(date/365,b90,group=as.factor(id),col="31-90 days"),size=1,alpha=0.6)+
    geom_line(aes(date/365,b180,group=as.factor(id),col="91-180 days"),size=1,alpha=0.6)+
    geom_line(aes(date/365,over180,group=as.factor(id),col="181-365 days"),size=1,alpha=0.6)+
    scale_color_manual("Short-term oscillations",values=c("2-14 days"="#EF476F",
                                                          "15-30 days"="#FFD166",
                                                          "31-90 days"="#06D6A0",
                                                          "91-180 days"="#118AB2",
                                                          "181-365 days"="#073B4C"))+
    scale_x_continuous("Time (Years)",breaks = seq(0,max(data$date/365),by=1))+
    scale_y_continuous("")+
    theme_light()+
    theme(text = element_text(face = "bold",
                              size=16),
          legend.position = "bottom")+
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  return(gg)
  
  
}



plot_pot <- function(data){
  
  id <- unique(data$id)
  data <- fitwavelet(id,data,deriv = F)
  
  
  gg<- ggplot(data)+
    geom_line(aes(date/365,k.trend,group=as.factor(id)),size=1,linetype="dashed",color="black")+
    geom_point(aes(date/365,k.mis),shape=8)+
    scale_x_continuous("Time (Years)",breaks = seq(0,max(data$date/365),by=1))+
    scale_y_continuous("mmol/L")+
    theme_light()+
    theme(text = element_text(face = "bold",
                              size=16),
          legend.position = "bottom")
  
  return(gg)
  
  
}



