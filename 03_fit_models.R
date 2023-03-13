library(nlme)
library(splines)
library(dynpred)
library(survival)
library(JMbayes2)
library(JMbayes)
library(magrittr)
library(tidyverse)
library(flextable)
library(xtable)
library(survival)

# load data 
load("data.RData")
load("landmark_data.RData")
load("pot_sample.RData")
source("fun.R")
source("utils.R")
# potassio %<>% select(-fup_days,-dht)
potassio_sel <- potassio %>% select(id,time,k)
index_visit$age_stand <- scale(index_visit$age)[,1]
potassio$age_stand <- scale(potassio$age)[,1]
labs_var <- c("Age", "Sex, Male", "NYHA III-IV (vs. I-II)", "HFrEF",">3 comorbidities")
potassio %<>% arrange(id,time) 

#######################################################################
############## LOCF  ###########################################
#######################################################################

LMdata$age_stand <- scale(LMdata$age)[,1]
LMdata$k_stand <- scale(LMdata$k)[,1]
LMdata$iperk <- ifelse(LMdata$k>5,1,0)
LMdata$ipok <- ifelse(LMdata$k<3.5,1,0)
                                                
# landmark Cox model with locf marker value
fitLOCF <-
  coxph(
    Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3 +
      tv(k,ord = 4,knots = c(2,3,4,6,7)) + strata(time) + cluster(id),
    data = LMdata
  )


summary(fitLOCF)

# table landmark-LOCF model
table2 <- tab_surv(fitLOCF,joint=F,digits = 2)
table2$Variable <- c(labs_var,paste0("s(k)"))
table2 %>% flextable()
print(xtable(table2,digits=2), include.rownames=FALSE)

fig2 <- ggrelrisk(fitLOCF,lab = "Last measurament of potassium (mmol/L)",joint = F,var="k",mean=0,ref=4,terms=6,b=seq(2.5,7.5,by=0.001))
fig2
ggsave(fig2,filename = "fig/Fig2.jpeg",device = "jpeg",dpi=300)

LMdata %<>% mutate(iperk=ifelse(k>5.5,1,0),
                    ipok=ifelse(k<3.5,1,0)) 


fitLOCF2 <-
  coxph(
    Surv(Tstart, fup_days, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3 +
      iperk +ipok+ strata(time) + cluster(id),
    data = LMdata)

summary(fitLOCF2)

table3 <- tab_surv(fitLOCF2,joint=F,digits = 2)
table3 %>% flextable()
table3$Variable <- c(labs_var,"k>5 mmol/L","k<3.5 mmol/L")
table3 %>% flextable()
print(xtable(table3,digits=2), include.rownames=FALSE)

#######################################################################
###################### MIXED LANDMARK #################################
#######################################################################

# table results linear mixed effect model
table4 <- tab_long(fp,2,joint = F)
table4$Variable <- c("Intercept","Sex, Male","Age","CKD","NYHA III-IV (vs. I-II)","s(t)")
table4 %>% flextable()
print(xtable(table4,digits=2), include.rownames=FALSE)


potassio_sample_long <- potassio_sample %>%
  group_by(id) %>%
  tidyr::expand(time=seq(0,4040)) %>%
  left_join(index_visit,by="id") %>%
  filter(time<=fup_days)

potassio_sample_long$pred<- predict(fp,newdata = potassio_sample_long)
potassio_sample_long %<>% mutate(time_rev=(time-fup_days)/365)

dht_lab <- c("Censored","Death")
names(dht_lab) <- c(0,1)

# FIG 3
fig3 <- ggplot()+
  geom_point(data=potassio_sample,aes(time_rev,k),alpha=0.2,col="grey30")+
  geom_line(data=potassio_sample_long,aes(time_rev,pred,group=as.factor(id)),alpha=0.5,col="darkblue",size=1)+
  facet_grid(~dht,labeller = labeller(dht=dht_lab))+
  scale_y_continuous(expression(paste(hat("m"),"(k)")),limits = c(2,8))+
  scale_x_continuous("Time until end of observation period (Years)")+
  theme_light()+
  theme(text = element_text(face = "bold",
                            size=16))
fig3
ggsave(plot = fig3,filename = "fig/Fig3.jpg",device = "jpeg",dpi=300,
       width = 7,height = 4.5)


fitMIXED <-
  coxph(
    Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34 + lvef_classi +  somma_comor_nc3 +
      tv(k.trend,ord = 4,knots = c(2,3,4,6,7))+ strata(time) + cluster(id),
    data = LMdata,
    method = "breslow"
  )


summary(fitMIXED)

table5 <- tab_surv(fitMIXED,joint=F,digits = 2)
table5 %>% flextable()
table5$Variable <- c(labs_var,"s(m(k))")
table5 %>% flextable()
print(xtable(table5,digits=2), include.rownames=FALSE)

# FIG 4
fig4 <-ggrelrisk(fitMIXED,lab = "Mean potassium value at landmark point (mmol/L)",joint = F,var="k.trend",mean=0,ref=4,terms=6,b=seq(2.5,7.5,by=0.001),rect = F)
fig4
ggsave(fig4,filename = "fig/Fig4.jpeg",device = "jpeg",dpi=300)

#######################################################################
###################### JOINT ###########################################
#######################################################################

CoxFit <-
  coxph(
    Surv(fup_days, dht) ~ age_stand + gender + nyha_34  + lvef_classi +  somma_comor_nc3,
    data = index_visit,
    x = TRUE
  )


# the joint model  with current value structure and slope
fitJoint1 <-jm(CoxFit,
               fp,
               time_var = "time",
               functional_forms = ~ value(tv(k,ord = 4,knots = c(2,3,4,6,7))),
               control = list(n_iter = 50000,
                              n_burnin = 10000)) 

summary(fitJoint1)

# table results longitudinal process
table6 <- tab_long(fitJoint1,2,joint = T)
table6$Variable <- c("Intercept","Sex, Male","Age","CKD","NYHA III-IV (vs. I-II)","s(t)")
table6 %>% flextable()
print(xtable(table6,digits=2), include.rownames=FALSE)


# FIG 5
fig5 <- ggrelrisk(fitJoint1,joint=T,ref=4,lab = "Current mean potassium value (mmol/L)",slope = F,b=seq(2.5,7.5,by=0.1),rect=F)
fig5
ggsave(fig5,filename = "fig/Fig5.jpeg",device = "jpeg",dpi=300)

# table results survival process
table7 <- tab_surv(fitJoint1,2)
table7$Variable <- c(labs_var,"s(m)")
table7 %>% flextable()
print(xtable(table7,digits=2), include.rownames=FALSE)

#######################################################################
############## WAVELET ###########################################
#######################################################################

fitWAVE1 <-
  coxph(
    Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34   + lvef_classi +  somma_comor_nc3 +
      tv(k.trend,ord = 4,knots = c(2,3,4,6,7))+
      tv(b14,ord = 2, knots = c(-2,0,2)) +
      tv(b30,ord = 2, knots = c(-2,0,2)) +
      tv(b90,ord = 2, knots = c(-2,0,2)) +
      tv(b180,ord = 2, knots = c(-2,0,2)) +
      tv(over180,ord = 2, knots = c(-2,0,2)) +
      strata(time) + cluster(id),
    data = LMdata,
    method = "breslow"
  )

fitWAVE2 <-
  coxph(
    Surv(Tstart, fup_days_w, dht) ~ age_stand + gender + nyha_34   + lvef_classi +  somma_comor_nc3 +
      tv(k.trend,ord = 4,knots = c(2,3,4,6,7))+
      b14_d2 +
      b30_d2 +
      b90_d2 +
      b180_d2 +
      over180_d2 +
      strata(time) + cluster(id),
    data = LMdata,
    method = "breslow"
  )

AIC(fitWAVE1)
AIC(fitWAVE2)


# table results survival process
table8 <- tab_surv(fitWAVE2,2,joint = F)
table8$Variable <- c(labs_var,"s(m)","014UP","014DOWN","O30UP","030DOWN","O90UP","090DOWN","0180UP","0180DOWN","0365UP","0365DOWN")
table8 %>% flextable()
print(xtable(table8,digits=2), include.rownames=FALSE)

