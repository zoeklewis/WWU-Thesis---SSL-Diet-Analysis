library(tidyr)
library(tidyverse)
library(ggplot2)
library(car)

final16dat<- read.csv("/Users/zoeklewis/Documents/Thesis/Data /WWU_Thesis_SSL-Diet-Analysis/Counts to RRA /SSL_RRA.csv")

winter <- c("12", "1", "2")
spring <- c("3","4","5")
summer <-c("6","7", "8")

chinookdat <- final16dat  %>% filter(species == "ONCTSH")   %>% select(species,sample,day,month,year,site,PRCNT100)

chinookdatwinter <- chinookdat  %>% filter(month %in% winter ) %>% add_column(method = "winter")

winteraverage <- mean(chinookdatwinter$PRCNT100)

winteraverage


chinookdatspring <- chinookdat %>% filter(month %in% spring) %>% add_column(method = "spring")

springaverage <- mean(chinookdatspring$PRCNT100)


chinookdatsummer <- chinookdat %>% filter(month %in% summer) %>% add_column(method = "summer")


chinookdatseasons <- bind_rows(chinookdatwinter, chinookdatspring, chinookdatsummer) %>% select(species,sample,day,month,year,site,PRCNT100, method)

