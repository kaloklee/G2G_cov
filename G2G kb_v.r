################Set directory#########################


setwd("~/R/G2G/")
#source("G2G_background_functions.r")

#######################################################

########################################################
#################Import data############################

kb <- read.csv("kb_covars.csv")

########################################################


# Find the index of the last instance of each person
last_instance_index <- tapply(seq_len(nrow(kb)), kb$id, tail, n = 1)
print(last_instance_index)
#initialize column
kb$Event <- 0

# Create the Event column
kb$Event[last_instance_index] <- ifelse(kb$censor[last_instance_index] == 0,1,0)


View(kb)


kb2<-kb

################Data manipulation#######################

id_df <- aggregate( censor ~ id, data = kb, FUN = max)

library(tidyverse)

kb2 <- kb %>%
  select (-censor) %>%
  inner_join(id_df, by="id")

kb2$status=(kb2$censor+1) %% 2

head(kb2,50)
min(kb2$status)
##########################################################
##############Estimation##################################
kb2$anyp=kb2$anyp

source("G2G_background_functions _al.r")
ptm <- proc.time()
test_sol<-G2G_varying_MLE(Surv(week,status) ~ coupon + anyp, data=kb2, subject="id") 
loop_over_i<-proc.time() - ptm




model_data$X

solution=optim(par=c(2,2,rep(0,dim(X)[3])),fn=G2G_varying,
               y = data_df$time,
               id = which(data_df$status==0),
               X=X,
               method="L-BFGS-B",
               lower=c(1e-5,1e-5,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf),
               hessian = TRUE)
#standard error and 95% CI
solution$par_stderr<-sqrt(diag(solve(solution$hessian)))
solution$par_upper<-solution$par+1.96*solution$par_stderr
solution$par_lower<-solution$par-1.96*solution$par_stderr

#showing the estimated parameters are within the 95% CI
solution$par
solution$par_upper
solution$par_lower
