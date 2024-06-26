#################Background functions#################


#Create a logdiffexp function to avoid computation error
logdiffexp <- function (a,b) {
  
  c = pmax(a,b);
  return (c + log(exp(a-c)-exp(b-c))) ;
  
}



#model log-likelihood function
G2G_varying_LL <- function(par,data_df) {
  #par: parameters
  #data_df: a data.frame with these columns:
  #                       id, id of each subject
  #                       time, duration observed for each subject
  #                       status, 0 if the event has not yet occurred, 1 if the event occurs
  #                       name_not_sure, the covariates matrix, 
  #par[1] = mean of BG
  #par[2] = polarization of BG
  #r=par[1]*(1/par[2]-1); 
  #alpha=(1-par[1])*(1/par[2]-1);
  
  r = par[1];
  alpha = par[2];
  coeff=par[-(1:2)];
  
  #X = scale(as.matrix(data_df[,-(1:3)]));
  X = as.matrix(data_df[,-(1:3)]);
  
  #print(length(coeff)==dim(X)[2])
  
  data_df$Ct = exp(X %*% coeff);
  
  data_df$cumsumCt <- ave(data_df$Ct, data_df$id, FUN = cumsum);
  
  #gather unique id information
  
#  print(names(data_df))
  
  id_df <- aggregate(cbind(time, status) ~ id, data = data_df, FUN = max)
  
  #uncensored piece of likelihood
  
  LL_uncen = 0;
  
  uncen_time = unique(id_df[which(id_df$status==1),"time"]);
  
  for (t in uncen_time) {
    
    id_need = id_df[which(id_df$status == 1 & id_df$time == t),"id"];
    
    if (t == 1) {
      
      cumsumCt_b = data_df[which(data_df$id %in% id_need & data_df$time ==1),
                           "cumsumCt"]; 
      
      LL_uncen = LL_uncen+sum(logdiffexp(0,-r*log(1+cumsumCt_b/alpha)));
      
    }
    
    else { 
      
      cumsumCt_a = data_df[which(data_df$id %in% id_need & data_df$time == (t-1)),
                           "cumsumCt"]; 
      cumsumCt_b = data_df[which(data_df$id %in% id_need & data_df$time == (t)),
                           "cumsumCt"]; 
      
      
      LL_uncen= LL_uncen + sum(logdiffexp(-r*log(1+cumsumCt_a/alpha), 
                                      -r*log(1+cumsumCt_b/alpha)));  
      
    }
    
  }
  
  #censored piece of likelihood  
  
  cen_time = unique(id_df[which(id_df$status==0),"time"]);
  
  LL_cen = 0;
  
  for (t in cen_time) {
    
    id_need = id_df[which(id_df$status == 0 & id_df$time == t),"id"];
    
    cumsumCt_b = data_df[which(data_df$id %in% id_need & data_df$time == (t)),
                         "cumsumCt"]; 
    
    LL_cen= LL_cen + sum((-r)*log(1+cumsumCt_b/alpha));
    
  }  
  
  #  print ( -(LL_uncen) )
  #  print( -(LL_cen) )
  #  print(par)
  
  return ( -(LL_uncen)-(LL_cen) );
  
}


#########Command to run MLE########################

#### Data preparation

G2G_varying_MLE <- function(fo, data, subject) {
  
  #fo: the formula in this format: Surv(A,B)
  #data: the data frame
  #id: a text field for the subject
  
  #the dependent variables: time and status
  time_name = all.vars(as.formula(fo))[1];
  status_name = all.vars(as.formula(fo))[2];
  
  df_temp<-data[,c(all.vars(as.formula(fo))[-c(1:2)])];
  
  X<-model.matrix( ~ ., data = df_temp);
  
  model_data<-data.frame(id = unlist(data[,c(subject)]), 
                         time = unlist(data[,c(time_name)]),
                         status = unlist(data[,c(status_name)]),
                         X=X[,-1]);
  #print(head(model_data))
  return( G2G_varying_optim(model_data) ); 
   
}

G2G_varying_optim <- function(model_data) {
  
  nvar=dim(model_data)[2]-1;
  
  solution=optim(par=c(0.5,.05,rep(0,nvar-2)),
                 fn=G2G_varying_LL,
                 data_df = model_data,
                 method="L-BFGS-B",
                 lower=c(.001,.001,rep(-5,nvar-2)), 
                 upper=c(Inf,Inf,rep(5,nvar-2)),
                 control = list(maxit=1000),
                 hessian = TRUE);
  solution$par_stderr<-sqrt(diag(solve(solution$hessian)))
  solution$par_upper<-solution$par+1.96*solution$par_stderr
  solution$par_lower<-solution$par-1.96*solution$par_stderr
  
  return (solution);
  
}

