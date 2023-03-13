library(nlme)
library(splines)
library(tidyverse)
library(magrittr)

# load data 
load("data.RData")
source("fun_wavelet.R")
source("utils.R")

set.seed(123)
sampleid <-  index_visit %>% pull(id) %>% sample(50)
potassio_sample <- potassio %>% filter(id%in%sampleid) 


##################################
## LONGITUDINAL TRAJECTORIES    ##
##################################

# FIG 1

dht_lab <- c("Censored","Death")
names(dht_lab) <- c(0,1)

potassio_sample %<>% mutate(time_rev=(time-fup_days)/365) 

fig1 <- ggplot(potassio_sample)+
  geom_point(aes(time_rev,k),alpha=0.2,col="grey30")+
  geom_line(aes(time_rev,k,group=as.factor(id)),alpha=0.3,col="grey30",size=1)+
  geom_smooth(aes(time_rev,k,group=as.factor(dht)),col="red",size=1,method = "loess",se=F)+
  facet_grid(~dht,labeller = labeller(dht=dht_lab))+
  scale_y_continuous("mmol/L",limits = c(2,8))+
  scale_x_continuous("Time until end of observation period (Years)")+
  theme_light()+
  theme(text = element_text(face = "bold",
                            size=16))

fig1
ggsave(plot = fig1,filename = "fig/Fig1.jpg",device = "jpeg",dpi=300,
       width = 7,height = 4.5)

potassio_sample %<>%left_join(index_visit[,c(1,58)]) 
save(potassio_sample,file="pot_sample.RData")


###################################################
############# PREPARATION LANDMARK DATASET ########
####################################################


fp <- lme(k ~ gender + age_stand + irc + nyha_34 + tv(time, ord=4, knots = c(0,218, 
                                                                             752, 1447,4040)),
          data = potassio,
          random = ~tv(time, ord=4, knots = c(0,218, 
                                              752, 1447,4040)) | id,
          control = lmeControl(opt = 'optim'))



potassioext$k.trend <- predict(fp,newdata = potassioext)
potassioext$k.detrended <- potassioext$k-potassioext$k.trend


datawavelm <- do.call("rbind",lapply(tLM,function(x)landmark_data_mixed_wavelet(tLM=x,
                                                                                w=width,
                                                                                longdata=potassio,
                                                                                longdataexp=potassioext,
                                                                                survdata=index_visit))) 

datawavelm$iperk=ifelse(datawavelm$k>5.5,1,0)
datawavelm$ipok=ifelse(datawavelm$k<3.5,1,0)
datawavelm$Tstart <- datawavelm$time
datawavelm$tLM <- datawavelm$time



######################
# EXAMPLE PERIOOGRAM #
######################

gg61 <- gg_periodogram( 1164807,data=potassioext)
gg60 <- plot_pot( gg61$data)
gg62 <- plot_osc(gg61$data)
fig6 <- egg::ggarrange(gg60,gg61$plot,gg62,ncol=1,nrow=3)
fig6

ggsave(plot = fig6,filename = "fig/Fig6.jpg",device = "jpeg",dpi=300,width = 9,height = 6)



