library(psych)
library(psychTools)
library(mice)
library(careless)
library(dplyr)
library(magrittr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data(spi)

dat <- spi

#### basic clean-up ####

careless_corruption <- function(dat,
                                pct_careless =.2,
                                pct_types = c(random = 1/3, longstring = 1/3, mixed = 1/3),
                                careless_levels = c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4)
                                )
  {
  #removing useless columns
  dat <- dat %>%
    select(-c(age:ER))
  
  #### careless screening ####
  #careless respondents will be removed with a high cut-off
  
  flagged <- matrix(nrow=nrow(dat),ncol=0)%>%data.frame()
  #mahalanobis distance
  flagged$mahad <- dat %>%
    mahad(flag=T,confidence=0.95,plot=F) %$%
    flagged
  
  #standard deviation
  flagged$sd <- dat %>%
    irv %>%
    {.<1}
  
  #longstring
  #max
  flagged$longMax <- dat %>%
    longstring %>%
    {.>=5}
  
  #avg
  flagged$longAvg <- dat %>%
    longstring(avg=T) %$%
    avgstr %>%
    {.>=1.4}
  
  #summing all flags
  sumFlagged <- flagged %>%
    apply(1,sum,na.rm=T)
  
  #keeping only those with sum=0
  dat <- dat[sumFlagged==0,]
  
  #ordering columns
  orderColumns <- spi.keys %>%
    unlist %>%
    gsub("-","",.)
  
  dat <- dat[,orderColumns]
  
  #### INSERTING CARELESSNES ####
  
  #creating carelessnes
  carelessInjection <- inject_careless(data=dat,
                                       pct_careless = pct_careless,
                                       pct_types = pct_types,
                                       careless_levels = careless_levels)
  
  #careless dataset
  carelessData <- carelessInjection$data_corrupted
  
  #labels
  carelessLabels <- carelessInjection[["labels"]]

  return(list(
    carelessData = carelessData,
    carelessLabels = carelessLabels
  ))
  }








