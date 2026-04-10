
############################################
# Testing the statistical significance 
# or homogeneity effect of \beta(s)
############################################

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
  if(c_null == 'unknown'){c_hat <- mean(est_beta[idx_tmin:idx_tmax])}
  
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
  return(list(rej =rej, t1n = t1n, t1_emp = t1_emp, chat = c_hat))
}


############################################
# Testing the linearity of \beta(s)
############################################

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
    dat_lm = data.frame(s_pt = t_int[idx_tmin:idx_tmax] ,
                        beta_pt = est_beta[idx_tmin:idx_tmax])
    lm_fit = lm(beta_pt~s_pt, data=dat_lm)$coef
    
    chat0 = as.numeric(lm_fit[1]) ; chat1 = as.numeric(lm_fit[2])
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
   
  return(list(rej =rej, t1n = t1n, t1n_emp = t1_emp, chat0 = chat0, chat1 = chat1))
}
