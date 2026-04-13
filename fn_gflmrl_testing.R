
############################################
# Testing the statistical significance 
# or homogeneity effect of \beta(s)
############################################
################################################
fit_mrl_constant_test = function(data, norder=4, t_Interval,
                                 model = model, test_range
){
  
  Y = data$Y
  delta = data$delta
  
  # Testing Time range index
  idx_tmin = which(abs(test_range[1] - t_Interval) < 1e-8)
  idx_tmax = which(abs(test_range[2] - t_Interval) < 1e-8)
  
  id_col <- c("Y", "delta", "x1", "x2")
  temp3 <- data[, !(colnames(data) %in% id_col)] %>% as.matrix()  #  Z_i(s)!!
  
  temp3_test <- temp3[,idx_tmin:idx_tmax]
  s_test <- t_Interval[idx_tmin:idx_tmax]
  
  deltat = (diff(t_Interval)[1])
  
  w1 = apply(temp3_test,1,sum)*deltat  
  Wmat_new = cbind(data$x1, data$x2, w1)
  
  est = est_beta_func(time = Y, delta = delta, vars = Wmat_new,
                      wt=1, model = model, se.cal = F)$estimate    # If data generation process v = 2, some seed fails to estimate...
  
  est_alpha = est[1:2]  # \hat \alpha 
  est_gamma = est[-c(1:2)] # hat c0, hat c1  
  
  est_all = c(est_alpha, est_gamma)
  return(est_all)
  
}

func_testing_constant = function(data, t_int,  
                                 boot_coef, est_coef,  
                                 alpha_nom = 0.05, 
                                 c_null, 
                                 test_range){
  
  n = nrow(data)
  s_len = length(est_coef)
  est_beta = est_coef
  B = ncol(boot_coef)
  
  # Testing Time range index
  idx_tmin = which(abs(test_range[1] - t_int) < 1e-8)
  idx_tmax = which(abs(test_range[2] - t_int) < 1e-8)
  
  ################################
  # Type I error rate Calculation 
  ################################
  
  # hat c 
  if(c_null == 'known'){c_hat = 0} # constant c of interest by researcher - This is for statistical significance
  if(c_null == 'unknown'){
    fit_cons = fit_mrl_constant_test(data=data, norder=4, t_Interval=t_int,
                                     model = model, test_range = range(t_int))
    c_hat = fit_cons[3] 
  }
  
  # Compute T_{1n} (Test stat)
  arg2 <- abs(est_beta - c_hat)
  arg2_range = arg2[idx_tmin:idx_tmax]
  t1n = sqrt(n)*max(arg2_range)
  
  # Compute T_{1n}^b (Bootstrapped test stat)
  ee = rep(est_beta, B)
  emat = matrix(ee, ncol = B, byrow=F)
  arg3 = abs(boot_coef[-c(1,2),] - emat)
  arg3_test = arg3[idx_tmin:idx_tmax, ] # Test range of s of interest
  t1b_boot = sqrt(n)*apply(arg3_test, 2, max)
  
  # Compute T_{1n}^b (1-\alpha) : Empirical alpha 
  t1_emp = quantile(t1b_boot, prob = 1- alpha_nom)
  
  # Test - Rejection rule 
  rej <- 1*I(t1n > t1_emp)
  
  # p-value 
  pval_raw <- mean(t1b_boot>=t1n)
  pval <- (1 + sum(t1b_boot >= t1n)) / (B + 1)
  
  return(list(rej =rej, t1n = t1n, t1_emp = t1_emp, chat = c_hat, pval = pval))
}


############################################
# Testing the linearity of \beta(s)
############################################
fit_mrl_linear_test = function(data, norder=4, t_Interval,
                               model = model, test_range
){
  Y = data$Y ;  delta = data$delta ; 
  
  # Testing Time range index
  idx_tmin = which(abs(test_range[1] - t_Interval) < 1e-8)
  idx_tmax = which(abs(test_range[2] - t_Interval) < 1e-8) 
  
  id_col <- c("Y", "delta", "x1", "x2")
  temp3 <- data[, !(colnames(data) %in% id_col)] %>% as.matrix()   
  
  temp3_test <- temp3[,idx_tmin:idx_tmax]  
  s_test <- t_Interval[idx_tmin:idx_tmax]
  
  deltat = (diff(t_Interval)[1])
  
  w1 = apply(temp3_test,1,sum)*deltat   
  w2 = ((temp3_test) %*% (s_test))*deltat
  Wmat_new = cbind(data$x1, data$x2, w1, w2) # (x1, x2,\int_a^b  Z(s) ds, \int_a^b s Z(s) ds )
  
  est = est_beta_func(time = Y, delta = delta, vars = Wmat_new,
                      wt=1, model = model, se.cal = F)$estimate    
  est_alpha = est[1:2]   
  est_gamma = est[-c(1:2)] # hat c0, hat c1
  
  est_all = c(est_alpha, est_gamma)
  return(est_all)
  
}


func_testing_linear = function(data, t_int,  
                               boot_coef, est_coef,  
                               model = model,
                               alpha_nom = 0.05, 
                               c_null, 
                               test_range
){
  
  n = nrow(data)
  s_len = length(est_coef)
  est_beta = est_coef
  B = ncol(boot_coef)
  
  # Testing Time range index
  idx_tmin = which(abs(test_range[1] - t_int) < 1e-8)
  idx_tmax = which(abs(test_range[2] - t_int) < 1e-8)
  
  ################################
  # Type I error rate and Power Calculation 
  ################################
  
  if(c_null == 'known'){  # under the null 
    c_hat = t_int -0.5 ; chat0 = -0.5 ; chat1 = 1  }
  if(c_null == 'unknown'){ 
    fit_linear = fit_mrl_linear_test(data=data, norder=4, t_Interval=t_int,
                                     model = model, test_range = range(t_int))
    
    chat0 = fit_linear[3] ; chat1 = fit_linear[4]
    c_hat = chat0 + chat1*t_int
  }
  
  # Compute T_{1n} (Test stat)
  arg2 <- abs(est_beta - c_hat)
  arg2_range = arg2[idx_tmin:idx_tmax]
  t1n = sqrt(n)*max(arg2_range)

  # Compute T_{1n}^b (Bootstrapped test stat)
  ee = rep(est_beta, B)
  emat = matrix(ee, ncol = B, byrow=F)
  arg3 = abs(boot_coef[-c(1,2),] - emat)
  arg3_test = arg3[idx_tmin:idx_tmax, ] # Test range of s of interest
  t1b_boot = sqrt(n)*apply(arg3_test, 2, max)
  
  # Compute T_{1n}^b (1-\alpha) : Empirical alpha 
  t1_emp = quantile(t1b_boot, prob = 1- alpha_nom)
  
  # Test - Rejection rule 
  rej <- 1*I(t1n > t1_emp)
  
  # p-value 
  pval_raw <- mean(t1b_boot>=t1n)
  pval <- (1 + sum(t1b_boot >= t1n)) / (B + 1)
   
  return(list(rej =rej, t1n = t1n, t1n_emp = t1_emp, chat0 = chat0, chat1 = chat1, pval=pval))
}
