library(tidyverse)
library(gmrl)
library(splines2)

##################################################
# Spline approximation and new format of the data 
# to fit the GFLMRL model 
##################################################


func_spline = function(data, n, norder, t_int, type){
  
  if(type == 'simple') {nbasis = 3}  # After CV, you may change this value based on the CV results
  if(type == 'complex') {nbasis = 7}
  if(type == 'constant' | type == 'linear') {nbasis = 2} 
  if(type == 'constant' | type == 'linear') {nbasis = 5} # For hypothesis testing
  
  norder = norder # 4 if cubic spline when using the cubic B-spline

  # Natural cubic spline version
  valueBasis <- naturalSpline(x = t_int,
                              df = nbasis, # The # of basis functions
                              intercept = T,  
                              Boundary.knots = c(min(t_int), max(t_int)), # Interval of \mathcal{S}
                              derivs = 0                 # 0th order -> Just the basis as it is
   )
   
  
  ################################
  # Estimate the function using spline and do integration
  #int_\ss Z(s)B(s) ds denote ad covMat
  
  id_col <- c("Y", "delta", "x1", "x2")
  temp3 <- data[, !(colnames(data) %in% id_col)] %>% as.matrix()  #  Z_i(s)
  deltat = (diff(t_int)[1])
  covMat <- (temp3%*%valueBasis)*deltat  # covMat: n \times nbasis dimension
  
  # Now, augment the new covariate to be used in the GFLMRL model
  Wmat = cbind(data$x1, data$x2, covMat) # (x1, x2, \int_\ss B(s) Z(s) ds )
  lst1 = list(Wmat = Wmat, valueBasis = valueBasis)
  return(lst1)
}


#######################################
# Reorganize and correct the functions 
# in GMRL package
#######################################

est_beta_func = function(time,delta,vars,wt=1,model='exp',se.cal=FALSE){
  data = data.frame(obs_time=time,delta=delta,pi=wt,vars)
  data_order = data[order(data$obs_time),]
  z_mat = data_order[,4:ncol(data_order)] 
  
  if(model=='exp'){
    ## estimate initial value from the Cox PH model
    tmp_fit = survival::coxph(Surv(obs_time,delta)~as.matrix(z_mat),data=data_order)  
    beta_ini = -summary(tmp_fit)$coef[,1]
    est = newton(data_order,z_mat,beta_ini,model=model)
    return(list(estimate=est))
  }
  
  if(model=='add'){
    ## estimate inital value from a linear regression
    tmp_fit = stats::lm(obs_time~ as.matrix(z_mat),data=data_order) # 여기 바꿔줌
    beta_ini = summary(tmp_fit)$coef[-1,1]
    est = newton(data_order,z_mat,beta_ini,model=model)
    return(list(estimate=est))
  }
}

 

score_exp = function(data_order,z_mat,beta0){
  n = nrow(data_order)
  col = ncol(z_mat)
  sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi))))) # original
  sn[is.na(sn)] = 0
  bn = apply(apply(apply(data_order$pi*exp(-z_mat%*%beta0),2,rev),2,cumsum),2,rev)/
    rev(cumsum(rev(data_order$pi)))
  bn[is.na(bn)] = 0
  hat_tmp = sn[-n]*bn[-1]
  hat_m0 = rev(cumsum(rev(c(hat_tmp,0)*diff(c(data_order$obs_time,data_order$obs_time[n])))))/sn
  z_bar = apply(apply(apply(data_order$pi*z_mat,2,rev),2,cumsum),2,rev)/rev(cumsum(rev(data_order$pi)))
  exp_part = exp(-z_mat%*%beta0)
  zpart = diff(c(0,data_order$obs_time))*z_bar
  part1 = as.numeric(hat_m0*(data_order$delta))*as.matrix(z_mat-z_bar) 
  part1[is.na(part1)] = 0
  part2 = (z_mat*data_order$obs_time-apply(zpart,2,cumsum))
  part2 = part2*c(exp_part)
  part2[is.na(part2)] = 0
  return(colMeans(data_order$pi*(part1-part2)))
}

score_add = function(data_order,z_mat,beta0){
  n = nrow(data_order)
  col = ncol(z_mat)
  sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))))
  sn[is.na(sn)] = 0
  part1 = c(rev(cumsum(rev(sn[-1]*diff(data_order$obs_time)))),0)
  part2_tmp = sn*data_order$pi*as.matrix(z_mat)%*%beta0*data_order$delta/rev(cumsum(rev(data_order$pi)))
  part2_tmp[is.na(part2_tmp)] = 0
  part2 = rev(cumsum(rev(part2_tmp)))
  hat_m0 = (part1-part2)/sn
  hat_m0[is.na(hat_m0)] = 0
  z_bar = apply(apply(apply(data_order$pi*z_mat,2,rev),2,cumsum),2,rev)/rev(cumsum(rev(data_order$pi)))
  z_bar[is.na(z_bar)] = 0
  u_tmp1 = data_order$pi*(hat_m0+as.matrix(z_mat)%*%beta0)*data_order$delta 
  u_tmp = t(u_tmp1)%*%as.matrix(z_mat-z_bar)  
  return(u_tmp/n)
}


newton = function(data_order,z_mat,beta_ini,model='exp',iter_max=500,tol=0.0001){
  if(model=='exp'){
    for(j in 1:iter_max){
      z_mat = as.matrix(z_mat)
      newbeta = beta_ini - solve(J_mat(data_order,z_mat,beta_ini,model = 'exp'))%*%score_exp(data_order,z_mat,beta_ini)
      if(prod(abs(newbeta-beta_ini)<tol)==T){break}
      beta_ini = newbeta
    }
    return(newbeta)
  }
  
  if(model=='add'){
    for(j in 1:iter_max){
      z_mat = as.matrix(z_mat)
      newbeta = beta_ini - solve(J_mat(data_order,z_mat,beta_ini,model = 'add'))%*%t(score_add(data_order,z_mat,beta_ini))
      if(prod(abs(newbeta-beta_ini)<tol)==T){break}
      beta_ini = newbeta
    }
    return(newbeta)
  }
}
 

##########################################
# Function for the parallel bootstrapping
##########################################

boot_se_ci_paral = function(B, data, t_Interval, model, norder, type){
  
  mc = 4; cores = 2
  n = nrow(data)
  
  registerDoParallel(mc, cores = cores)

  res = foreach(i = 1 : B, .combine = cbind) %dopar% {
    unlist(fit_fda_mrl(data=data[sample(n, n, replace = TRUE), ], 
                norder = norder, t_Interval = t_Interval, 
                model = model, type = type)) %>% as.vector()
  }
  res_coef = res[1:(2+length(t_Interval)),] # alpha1, alpha2, beta1, beta2,.. ,beta(101)
  
  # Standard errors of \hat \alpha and \hat \beta(s)
  se = apply(res_coef, 1, sd)
  se_alp = se[1:2]
  se_beta = se[-c(1:2)]
  
  boot_alp = res_coef[1:2,]
  boot_beta = res_coef[-c(1:2),]
  
  # 95% Pecentile bootstrap confidence intervals for \alpha and each beta(s) 
  ci_alp_lb = apply(boot_alp,1,function(x){quantile(x,0.025,na.rm = T)})
  ci_alp_ub = apply(boot_alp,1,function(x){quantile(x,0.975,na.rm = T)})
  
  ci_beta_lb = apply(boot_beta,1,function(x){quantile(x,0.025,na.rm = T)})
  ci_beta_ub = apply(boot_beta,1,function(x){quantile(x,0.975,na.rm = T)})
  
  ci_alp = c(rbind(ci_alp_lb, ci_alp_ub)) 
  ci_beta = c(rbind(ci_beta_lb, ci_beta_ub))  
  
  return(list(boot_est = res_coef,  
              se_alp = se_alp, se_beta = se_beta, 
              ci_alp = ci_alp, ci_beta = ci_beta))
} 


##########################################
# Main function to fit the GFLMRL model 
##########################################

fn_gflmrl = function(data, norder=4, t_Interval, model = 'exp', type){
  
  Y = data$Y
  delta = data$delta
  
  spl_fit = func_spline(data=data, n=nrow(data), norder=norder, t_int=t_Interval, type=type)
  Wmat_new = spl_fit$Wmat
  valueBasis = spl_fit$valueBasis

  est = est_beta_func(time = Y, delta = delta, vars = Wmat_new,
                      wt=1, model = model, se.cal = F)$estimate    
  
  est_alpha = est[1:2]  # \hat \alpha 
  est_gamma = est[-c(1:2)]  
  est_psi = valueBasis %*% est_gamma   # \hat \beta(s)

  est_all = c(est_alpha, est_psi)
  res_all = list(est_all = est_all)
  return(res_all)
  
}


#############################################################
# Functions to make the table for the simulation performance 
#############################################################

tb_res_alpha = function(alest, alse, alci, alpha0){
  
  bias_alp = (colMeans(alest)-alpha0) 
  ese_alp = apply(alest,2,sd)  # standard deviation of estimates
  ase_alp = colMeans(alse) # Estimates of standard errors of estimators
  
  # CP calculation - Percentile bootstrap version
  ci_L = t(alci[,c(1,3)]) <= alpha0
  ci_U = t(alci[,c(2,4)]) >= alpha0
  cp1 = rowMeans(ci_L*ci_U)
  
  # CP calculation - Wald type version
  ci_L = t(alest) - 1.96*t(alse) <= alpha0
  ci_U = t(alest) + 1.96*t(alse) >= alpha0
  cp2 = rowMeans(ci_L*ci_U)
  
  tb1 = t(round(rbind(bias_alp,ese_alp,ase_alp,cp1,cp2), 3))
  rownames(tb1) = NULL
  tb1
   
}

tb_res_beta = function(best, bse, bci, beta0){
  
  bias_beta = (colMeans(best)-beta0) 
  ese_beta = apply(best,2,sd)  # standard deviation of estimates
  ase_beta = colMeans(bse) # Estimates of standard errors of estimators
  
  # CP calculation - Percentile bootstrap version
  ci_L = t(bci[,seq(1, 2*length(beta0), by = 2)]) <= beta0
  ci_U = t(bci[,seq(2, 2*length(beta0), by = 2)]) >= beta0
  cp1 = rowMeans(ci_L*ci_U)
  
  # CP calculation - Wald type version
  ci_L = t(best) - 1.96*t(bse) <= beta0
  ci_U = t(best) + 1.96*t(bse) >= beta0
  cp2 = rowMeans(ci_L*ci_U)
  
  tb1 = t(round(rbind(bias_beta,ese_beta,ase_beta,cp1,cp2), 3))
  rownames(tb1) = NULL
  tb1
}


# Function to calculate the measure `MSPE`

fn_cnorm = function(est_vec, true_vec, t_int, K){
  
  lt <- length(t_int)         # Number of grid points of \mathcal{S}
  dt <- diff(t_int)[1]        # Step size
  diff <- est_vec - true_vec     # difference between estimated linear predictor and true linear predictor
  
  Phi <- matrix(0, nrow = lt, ncol = K)
  Phi[, 1] <- 1
  for (k in 2:K) {Phi[, k] <- cos((k - 1) * pi * t_int)}
  
  b_diag <- 1 / (3 * 1:K)
  cov <- diag(b_diag)  # covariance matrix
  ms <- Phi %*% cov %*% t(Phi)
  res <- as.numeric(dt^2 * t(diff) %*% ms %*% diff)
  return(sqrt(res))
} 
