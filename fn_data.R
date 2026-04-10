###############################################
# Sample Data generation Code 
# under the GFLMRL model 
###############################################

library(tidyverse)
 
#################################
# Functional covariate generation 
#################################

simBX <- function(t, v, nk, type){
  
  # t: Interval of integral (compact space \mathcal{S})
  # nk: Large K ( = 50) 
  # type: Level of the complexity of \beta(s)
  
  #########################
  # Functional covariate 
  kk <- 1*pi*(seq(1, nk, 1) -1); 
  phik <- sqrt(1)*cos(t%*%t(kk))  
  phik[, 1] <- rep(1, length(t));   
  uk <- runif(nk, -sqrt(1), sqrt(1));
  
  gk <- ((-1)^(seq(2, (nk+1), 1)))*(seq(1, nk, 1)^(-v/2)); # \xi_k
  XK <- gk*uk*t(phik); # xi_k U_{ik} \phi_k(s) 
  X <- (colSums(XK)); #  Z_i(s)
  
  #########################
  # Coefficient of Functional covariate : \beta(s)  
  bk <- 0.8*((-1)^(seq(1, nk, 1)))*(seq(1, nk, 1)^(-2));  
  if(type =='constant'){ B0 <- rep(-0.5, length(t)) } # \psi(s)
  if(type =='linear'){ B0 <- t-0.5 }
  if(type == 'simple'){ B0 <- colSums(bk*t(phik)) }  
  if(type == 'complex'){ B0 <- 1*colSums(bk*t(phik)) + sin(3*pi*t)  + cos(pi*t/4) }   
  if(type == 'complex2'){ 
    B0 <- 1*colSums(bk*t(phik)) + 2*sin(3*pi*t)  + 0.5*cos(pi*t/4); 
    # Stretching the hill and valleys of complex beta(s)
  }   
  alp <- X*B0 # Z_i(s) \times \beta(s) -> 1 \times K vector -> integrand
  
  ###########################
  # Resulting object 
  # B0 : \beta(s) -> Coefficient function  
  # alp : Z_i(s) \psi(s)  
  # X : Z_i(s) -> Functional covariate Z(s)  
  ###########################
  return(list(X = X, alp = alp, B0 = B0));
}



#######################
# Data generation code
#######################

simudata = function(n, # Sample size  
                    alpha0, # scalar covariate coefficient
                    model = c('exp', 'add'),  # Link function of GFLMRL
                    type = c('constant', 'linear','simple', 'complex', 'complex2'), # Complexity of \beta(s)
                    t_int = seq(0,1,0.01), # \mathcal{S} : Interval for functional covariate
                    v = c(1,1.5,2,2.5), 
                    K = 20, # The number of basis
                    lam = 0.3, # Calibration factor to adjust the censoring rate 
                    D1 = D1, D2 = D2 # Hall-Wellner Family to generate baseline MRL m_0(t)
){
  
  #####################
  # Generate Covariates
  ######################
  
  # Scalar covariate 
  if(model == 'exp'){x1 = runif(n,0,1); x2 = rbinom(n,1,0.5); xx = cbind(x1, x2)  }  
  if(model == 'add'){x1 = runif(n,-1,-0.5); x2 = runif(n,-1,-0.5); xx = cbind(x1, x2)}
  
  # Functional covariate 
  m = length(t_int)
  nk = K
  Z <- Zpsi <- matrix(NA, n, m);
  
  for(i in 1:n){
    fc <- simBX(t=t_int, v=v, nk=K, type=type);
    Z[i, ] <- fc$X; 
    Zpsi[i, ] <- fc$alp  
  }
  
  dt <- (diff(t_int)[1]);   
  int_zs <- apply(Zpsi, 1, sum)*dt;  # \int_0^1 Z_i(s) \psi(s) ds (Functional part)
  eta <- (xx %*% alpha0) + int_zs # (Z %*% beta0)
  
  if(model == 'exp'){
    a1 = D1 * exp(eta) ; b1 = D2 * exp(eta) 
    U = runif(n,0,1)
    Time = (b1/a1) * ( (1-U)^(-a1/(a1+1)) - 1) 
  }
  
  if(model == 'add'){
    U = runif(n,0,1)
    t1 = ((D2+eta)/D1) * ((1-U)^(-D1/(D1+1)) -1)
    t2 = -(D2/D1) - eta * log( (1-U)* ((eta/(D2+eta))^((D1+1)/D1)) )
    S0 = (eta/ (D2+eta))^(-(D1+1)/D1) # Suv val at left limit and right limit of break point of T
    SU = 1-U
    Time = ifelse(SU>=S0, t1, t2)
  }
  
  C = rexp(n,lam) # Censoring Time
  
  delta = 1*I(Time<=C)
  Y = pmin(Time,C)
  
  dat = data.frame(Y = Y, delta= delta, xx, Z)
  lst = list(dat = dat, psi0 = fc$B0)
  return(lst)
}



 