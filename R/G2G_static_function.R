
library(tidyverse)
library(survival)

#####################Veteran dataset in the Survival######################
head(veteran)


#Have a 2-covariates example:

#create a data frame to store these key variables:
#time: the duration of each subject (uncensored and censored) 
#status: whether the event happens (1) or not (0)
#X's: the remaining columns are the static covariates

#Please note that I divide age and karno by 10
#I don't know why, and maybe it's the estimation code. 
#But it is required to downscale the covariates to avoid numerical error.  
#You should try without the division

data_df<-data.frame(time = veteran$time,
                    status = veteran$status,
                    x1=veteran$age/10,  
                    x2=veteran$karno/10) #


#################Estimation function#################

#Create a logdiffexp function to avoid computation error
logdiffexp <- function (a,b) {
  
  c = pmax(a,b);
  return (c + log(exp(a-c)-exp(b-c))) ;
  
}

#model likelihood function
G2G_static_MLE <- function(par,df) {
  #par: parameters

  r=par[1];
  alpha=par[2];
  beta=par[-c(1:2)];
  
  #y: duration for each person 
  #status: status
  #X: all the covariates in a matrix format
  y = df[,1];
  status = df[,2];
  X = as.matrix(df[, -c(1:2)]);

  #uncensored piece of likelihood 
  uncen=y[which(status==1)];
  C_u = exp(X[which(status==1),] %*% beta); # exp(X*b) -> Perhaps there is a trick to downscale X here (so the user does not)
  LL_uncen=logdiffexp( -r*log(1+C_u*(uncen-1)/alpha), 
                       -r*log(1+C_u*(uncen)/alpha) );

  #censored piece of likelihood  
  cen=y[which(status==0)];
  C_c = exp(X[which(status==0),] %*% beta); # exp(X*b) -> investigate ways to downscale to avoid numercial error
  LL_cen = -r*log(1+C_c*(cen)/alpha);
  
  return (-sum(LL_uncen)-sum(LL_cen));
}

#######Estimation Function Complete##############################################



#we use the standard 'optim' routine for estimation

solution=optim(par=c(0.5,0.5,rep(0,dim(data_df)[2]-2)),
               fn=G2G_static_MLE,
               df = data_df,
               method="L-BFGS-B",
               lower=c(1e-5,1e-5,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf),
               hessian = TRUE)

#standard error and 95% CI
solution$par_stderr<-sqrt(diag(solve(solution$hessian)))
solution$par_upper<-solution$par+1.96*solution$par_stderr
solution$par_lower<-solution$par-1.96*solution$par_stderr

#Show the estimated parameters are within the 95% CI
solution$par
solution$par_lower
solution$par_upper


#Simulate from the estimated gamma distribution visualize the density of the p 

ran <- data.frame(xx = rgamma(10000,solution$par[1],solution$par[2]))
ran$pp <- 1-exp(-ran$xx)
ggplot(ran, aes(x=pp)) + 
  geom_density()


