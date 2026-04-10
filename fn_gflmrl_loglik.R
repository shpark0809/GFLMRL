library(tidyverse)
library(gmrl)
library(splines2)
library(pracma)
library(BB)

func_spline_tune = function(data, n, norder, t_int, nbasis){
  
  norder = norder # cubic spline  
  
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

# Estimate the m_0(t) given the estimates par 
fn_hatm0 = function(data_order, z_mat, par, model){
  
  n = nrow(data_order)
  col = ncol(z_mat)
  z_mat = as.matrix(z_mat)
  
  if(model == 'exp'){
    sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))))
    sn[is.na(sn)] = 0
    bn = apply(apply(apply(data_order$pi*exp(-z_mat%*%par),2,rev),2,cumsum),2,rev)/
      rev(cumsum(rev(data_order$pi)))
    bn[is.na(bn)] = 0
    hat_tmp = sn[-n]*bn[-1]
    hat_m0 = rev(cumsum(rev(c(hat_tmp,0)*diff(c(data_order$obs_time,data_order$obs_time[n])))))/sn
  }
  
  if(model == 'add'){
    sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))))
    sn[is.na(sn)] = 0
    part1 = c(rev(cumsum(rev(sn[-1]*diff(data_order$obs_time)))),0)
    part2_tmp = sn*data_order$pi*as.matrix(z_mat)%*%par*data_order$delta/rev(cumsum(rev(data_order$pi)))
    part2_tmp[is.na(part2_tmp)] = 0
    part2 = rev(cumsum(rev(part2_tmp)))
    hat_m0 = (part1-part2)/sn
    hat_m0[is.na(hat_m0)] = 0
  }
  return(hat_m0)
}


fn_gflmrl_tune = function(data, norder=4, t_Interval, model = model, df){
  
  Y = data$Y ; delta = data$delta ; n = nrow(data)
  
  spl_fit = func_spline_tune(data=data, n=nrow(data), 
                             norder=norder,
                             t_int=t_Interval, nbasis = df)
  Wmat_new = spl_fit$Wmat
  valueBasis = spl_fit$valueBasis
  
  est = est_beta_func(time = Y, delta = delta, vars = Wmat_new,
                      wt=1, model = model, se.cal = F)$estimate    # 
  
  est_alpha = est[1:2]  # \hat \alpha 
  est_gamma = est[-c(1:2)]  
  est_psi = valueBasis %*% est_gamma   # \hat \beta(s)
  est_all = c(est_alpha, est_psi)
  
  res_all = list(est_all = est_all, est_gamma = est_gamma,
                 est_alpha = est_alpha, Wmat_new = Wmat_new)
  return(res_all)
}



fn_loglik_cv = function(train_data, test_data, df, model=model){
  
  dt_train = train_data[order(train_data$Y),]
  dt_order = dt_train
  
  fit_dfi = fn_gflmrl_tune(data = dt_order, norder=4, 
                             t_Interval, model = model, 
                             df = df)
  
  est_dfi = (fit_dfi$est_all)  # est of alpha and beta(s) from nbasis = df
  est_gamma = fit_dfi$est_gamma
  est_alpha = fit_dfi$est_alpha
  est_dfi_sp = c(est_alpha, est_gamma)
  Wmat_new = fit_dfi$Wmat_new
  
  exp_axmat = as.matrix(Wmat_new[,1:2]) %*% est_alpha
  exp_iza = as.matrix(Wmat_new[,-c(1:2)]) %*% est_gamma
  
  tms = sort((dt_order$Y))
  m0_all = NULL
  data2 = data.frame(obs_time=dt_order$Y,delta=dt_order$delta,pi=1,vars=Wmat_new)
  data_order2 = data2[order(data2$obs_time),]
  z_mat2 = data_order2[,4:ncol(data_order2)] 
  m0_all = fn_hatm0(data_order2, z_mat2, par=est_dfi_sp, model=model)
  
  if(model == 'exp'){m0_all = log(pmax(m0_all, 1e-5))}
  if(model == 'add'){m0_all = pmax(m0_all, 1e-5)}
  
  ################################################################
  # Part 1: Compute the derivative of baseline MRL m_0(t) 
  # Since our m_0(t) estimates could have explosions, 
  # Do smoothing first and then do derivative
  
  fit_m0 <- smooth.spline(data_order2$obs_time, m0_all, spar = 0.9)  # Higher spar, the smoother 
  m0_fun  <- function(x) predict(fit_m0, x)$y
  m0dot_fun <- function(x) predict(fit_m0, x, deriv = 1)$y
   
  
  ## Approximate these values for test data
  test_order = test_data[order(test_data$Y),]
  spl_test = func_spline_tune(data= test_order, n=nrow(test_order), 
                              norder=4, t_int=t_Interval, nbasis =df)
  Wmat_new_test = spl_test$Wmat
  
  eta_train <- as.numeric(as.matrix(Wmat_new) %*% est_dfi_sp)
  eta_test <- as.numeric(as.matrix(Wmat_new_test) %*% est_dfi_sp)
  
  m0_test    <- m0_fun(test_order$Y)
  m0dot_test <- m0dot_fun(test_order$Y)
  
  ############################################# 
  # For each link function, 
  # Part 2: Compute \int_0^t exp(-m_0(u)) du
  # Part 3: Compute log-likelihood 
  
  if(model == 'exp'){
    u_grid <- seq(0, max(dt_order$Y), length.out = 10000)
    f_grid <- exp(-m0_fun(u_grid))
    cum_grid <- cumtrapz(u_grid, f_grid)
    intg_val <- approx(u_grid, cum_grid, xout = dt_order$Y, rule = 2)$y
    I_fun <- approxfun(dt_order$Y, intg_val, rule = 2)  
    I_test <- I_fun(test_order$Y)  
    
    haz = m0dot_test + exp(-m0_test - eta_test )
    haz = pmax(haz, 1e-5) # This is required due to some negative haz coming from numerical studies. 
    arg1 = exp( (m0_fun(0)) - m0_test)
    arg2 = exp(-exp(-eta_test) * I_test)
    suv = arg1*arg2
    suv = pmin(1, pmax(suv, 1e-5))
    loglik = sum(log(suv) + (test_order$delta)*log(haz))
  }
  
  if(model == 'add'){
    
    u_grid <- seq(0, max(dt_order$Y), length.out = 10000)
    m0_grid <- m0_fun(u_grid)
    intg2 <- sapply(seq_along(dt_order$Y), function(i) {
      idx <- u_grid <= (dt_order$Y)[i]
      trapz(u_grid[idx], 1/(m0_grid[idx] + eta_train[i]))
    })
    intg_val <- exp(-intg2)
    I_fun <- approxfun(dt_order$Y, intg_val, rule = 2)  
    I_test <- I_fun(test_order$Y)  

    arg1_haz = pmax(m0dot_test + 1, 1e-5)
    arg2_haz = pmax(m0_test + eta_test, 1e-5)
    haz = arg1_haz / arg2_haz 
    arg1 = (m0_fun(0) + eta_test) / (m0_test + eta_test)
    arg2 = I_test
    suv = arg1*arg2
    suv = pmin(1, pmax(suv, 1e-5))
    loglik = sum(log(suv) + (test_order$delta)*log(haz))
  }
  return(loglik)
}
