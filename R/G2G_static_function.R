
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
                    x1=veteran$age,  
                    x2=veteran$karno) #


#################Estimation function#################

# Numerically stable calculation of log(1 - exp(x))
# Following the algorithm of Mächler 2012
# Mächler 2012: https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
# Returns -Inf when x == 0 and NaN when x > 0
log1mexp <- function(x) {
  ifelse(
    x > -0.6931472,  # approx log(2)
    log(-expm1(x)),
    log1p(-exp(x))
  )
}

# Numerically stable calculation of log(exp(a) - exp(b))
# Again following the algorithm of Mächler 2012
# Returns -Inf when a == b (including a == b == -Inf, since log(0) = -Inf)
# Returns NaN when a < b or a == Inf
logdiffexp <- function(a, b) {
  ifelse(
    (a < Inf) & (a > b),
    a + log1mexp(b - a),
    ifelse((a < Inf) & (a == b), -Inf, NaN)
  )
}

#model likelihood function
G2G_static_MLE <- function(par,df) {
  #par: parameters
  
  r=exp(par[1]);
  alpha=exp(par[2]);
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

solution=optim(par=c(log(0.5),log(0.5),rep(0,dim(data_df)[2]-2)),
               fn=G2G_static_MLE,
               df = data_df,
               method="BFGS",
               hessian = TRUE)


#standard error and 95% CI
solution$par_stderr<-sqrt(diag(solve(solution$hessian)))
solution$par_upper<-solution$par+1.96*solution$par_stderr
solution$par_lower<-solution$par-1.96*solution$par_stderr

# Indices of parameters that are on log scale
exp_idx <- 1:2  # always first two

# Exponentiate only the relevant intervals and estimates
solution$par[exp_idx]        <- exp(solution$par[exp_idx])
solution$par_upper[exp_idx]  <- exp(solution$par_upper[exp_idx])
solution$par_lower[exp_idx]  <- exp(solution$par_lower[exp_idx])

#Show the estimated parameters are within the 95% CI
solution$par
solution$par_lower
solution$par_upper


#Simulate from the estimated gamma distribution visualize the density of the p 

ran <- data.frame(xx = rgamma(10000,solution$par[1],solution$par[2]))
ran$pp <- 1-exp(-ran$xx)
ggplot(ran, aes(x=pp)) + 
  geom_density()


