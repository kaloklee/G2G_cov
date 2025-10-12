#install.packages("discSurv")
#install.packages("eha")

library(discSurv)
library(eha)
library(tidyverse)
data("scania")
head(scania)

Scania_Person <- scania %>%
  mutate(exit = ceiling(exit),
         birthdate = floor(birthdate),
         spell = exit - enter) %>% #spell refers to the observed duration of a person
  mutate(enter = enter - 50,
         exit = exit - 50)

head(Scania_Person)

set.seed(123)
Scania_Person_Train <- sample_frac(Scania_Person, 0.8)
Scania_Person_Test <- Scania_Person[!Scania_Person$id %in% Scania_Person_Train$id,]

#convert the training set
Scania_PersonPeriod_Train <- dataLong(dataShort = Scania_Person_Train, 
                                      timeColumn = "spell", 
                                      eventColumn = "event",
                                      timeAsFactor = F) %>%
  as_tibble() %>%
  mutate(enter = timeInt - 1,
         age = enter + 50) %>%
  select(-obj, -event, -exit) %>%
  rename(event = y,
         exit = timeInt) %>%
  mutate(year = age + birthdate) %>%
  select(id, enter, exit, event, everything()) %>%
  left_join(logrye, by = "year") #joined with the `logrye` data for a variable on yearly food prices

head(Scania_PersonPeriod_Train, 50)
names(Scania_PersonPeriod_Train)
source("~/R/G2G/G2G_background_functions _al.r")
ptm <- proc.time()
abc<-G2G_covariate_MLE(Surv(exit,event) ~ sex + immigrant + foodprices, data=Scania_PersonPeriod_Train, subject="id") 
proc.time() - ptm


s