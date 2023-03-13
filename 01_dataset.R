library(tidyverse)
library(magrittr)
library(gtsummary)
library(xtable)
library(survival)
library(flextable)


# load data 
load("data//data.RData")

###########
# TABLE 1 #
###########

potassio_id <- potassio %>% distinct(id,.keep_all = T) %>% select(-k) 
index_visit %<>% left_join(potassio_id) 

table1 <- index_visit %>% 
  ungroup() %>% 
  distinct(id,.keep_all = T) %>% 
  select(age,gender,lvef,lvef_classi,nyha,irc,somma_comor_nc3,n_pot) %>% 
  tbl_summary(label = list(n_pot~"Number of potassium measurements")) %>%
  as_hux_table() %>% 
  as_tibble() 

table1 %>% flextable()
print(xtable(table1), include.rownames=FALSE)

potassio%<>% group_by(id) %>%
  mutate(dist=time-lag(time))

median_dist <- potassio %>% group_by(id) %>%
  mutate(dist=time-lag(time)) %>%  
  summarize(median_dist=median(dist,na.rm = T))


summary(median_dist$median_dist)
summary(potassio$dist)

##########
# KM     #
##########

fit <- survfit(Surv(fup_months,dht)~1,data=index_visit)
table2 <- as.data.frame(summary(fit,times = c(12,24,36,48,60))[c(1,2,3,4,5,6,7,14,15)])
table2



#################
# POTASSIUM     #
#################

# Number of measurements
potassio %>% distinct(id,.keep_all = T) %>% pull(n_pot) %>% summary() 

# Distance betweeen measurements
potassio %>% group_by(id) %>% mutate(dist=time-lag(time)) %>% pull(dist) %>% summary()

