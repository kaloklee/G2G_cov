#' Log-Likelihood for G2G Model with Time-Varying Covariates
#'
#' Internal function that computes the negative log-likelihood for the G2G model
#' with time-varying covariates. Used by the optimization routine.
#'
#' @param par Numeric vector of parameters. First two elements are r and alpha
#'   (shape and rate parameters), remaining elements are covariate coefficients.
#' @param data_df Data frame with columns:
#'   \describe{
#'     \item{id}{Subject identifier}
#'     \item{time}{Time point for each observation}
#'     \item{status}{Event indicator (0 = no occurence, 1 = event occurred)}
#'     \item{...}{Time-varying covariates (remaining columns)}
#'   }
#' @return Scalar value of the negative log-likelihood
#' @keywords internal
G2G_varying_LL <- function(par,data_df) {
  #par[1] = mean of BG
  #par[2] = polarization of BG
  #r=par[1]*(1/par[2]-1); 
  #alpha=(1-par[1])*(1/par[2]-1);
  
  r = exp(par[1]);
  alpha = exp(par[2]);
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

#' Fit G2G Model with Time-Varying Covariates
#'
#' Fits a Grassia(II)-Gamma (G2G) survival model with time-varying covariates using
#' maximum likelihood estimation.
#'
#' @param fo Formula object specifying the model in the format 
#'   \code{Surv(time, status) ~ x1 + x2 + ...} where time is the observation
#'   time, status is the event indicator, and x1, x2, etc. are time-varying
#'   covariates. Note: The formula should use variable names, not column indices.
#' @param data Data frame containing all variables specified in the formula.
#'   Data should be in long format with one row per time point per subject 
#'   (i.e. period-person format)
#' @param subject Character string specifying the name of the subject ID column
#'   in the data frame. This column identifies which observations belong to
#'   the same subject over time.
#'
#' @return A list containing the optimization results from \code{\link[stats]{optim}}
#'   with additional elements:
#'   \item{par}{Parameter estimates. First two are r and alpha (shape and rate
#'     parameters), remaining are covariate coefficients.}
#'   \item{par_stderr}{Standard errors of parameter estimates}
#'   \item{par_upper}{Upper bounds of 95\% confidence intervals}
#'   \item{par_lower}{Lower bounds of 95\% confidence intervals}
#'   \item{value}{Negative log-likelihood at the optimum}
#'   \item{convergence}{Convergence code (0 indicates successful convergence)}
#'   \item{hessian}{Hessian matrix at the optimum}
#'
#' @export
#' @importFrom stats aggregate optim as.formula ave model.matrix
#'
G2G_varying_MLE <- function(fo, data, subject) {
  
  #fo: the formula in this format: Surv(A,B)
  #data: the data frame
  #id: a text field for the subject
  
  #the dependent variables: time and status
  time_name = all.vars(as.formula(fo))[1];
  status_name = all.vars(as.formula(fo))[2];
  
  #df_temp<-data[,c(all.vars(as.formula(fo))[-c(1:2)])];
  df_temp <- data[, c(all.vars(as.formula(fo))[-c(1:2)]), drop = FALSE];
  X<-model.matrix( ~ ., data = df_temp);
  
  model_data<-data.frame(id = unlist(data[,c(subject)]), 
                         time = unlist(data[,c(time_name)]),
                         status = unlist(data[,c(status_name)]),
                         X=X[,-1]);
  #print(head(model_data))
  return( G2G_varying_optim(model_data) ); 
   
}

#' Internal Optimization for Time-Varying Model
#'
#' Internal wrapper function that performs the actual optimization for the
#' G2G model with time-varying covariates using L-BFGS-B method.
#'
#' @param model_data Data frame prepared by \code{G2G_varying_MLE} containing
#'   id, time, status, and covariate columns
#' @return List containing optimization results from \code{\link[stats]{optim}}
#'   with additional elements for standard errors and confidence intervals:
#'   \item{par}{Parameter estimates}
#'   \item{par_stderr}{Standard errors}
#'   \item{par_upper}{Upper bounds of 95\% confidence intervals}
#'   \item{par_lower}{Lower bounds of 95\% confidence intervals}
#' @keywords internal
G2G_varying_optim <- function(model_data) {
  
  nvar=dim(model_data)[2]-1;
  
  solution=optim(par=c(log(0.5),log(0.5),rep(0,nvar-2)),
                 fn=G2G_varying_LL,
                 data_df = model_data,
                 method="BFGS",
                 control = list(maxit=1000),
                 hessian = TRUE);
  solution$par_stderr<-sqrt(diag(solve(solution$hessian)))
  solution$par_upper<-solution$par+1.96*solution$par_stderr
  solution$par_lower<-solution$par-1.96*solution$par_stderr

  # Exponentiate first two parameters
  exp_idx <- 1:2
  solution$par[exp_idx] <- exp(solution$par[exp_idx])
  solution$par_upper[exp_idx] <- exp(solution$par_upper[exp_idx])
  solution$par_lower[exp_idx] <- exp(solution$par_lower[exp_idx])
  
  return (solution);
  
}

