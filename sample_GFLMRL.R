
rm(list=ls())
 
source('fn_data.R')
source('fn_gflmrl.R')
source('fn_gflmrl_testing.R')

n = 500
B = 500  # The number of bootstrapping

alpha0 = c(-1,-1)
t_Interval = seq(0,1,0.01)
K = 20 
v = c(1,1.5,2,2.5)[1]

model = c('exp', 'add')[1]  # Link function
type = c('constant', 'linear', 'simple', 'complex')[4]
cen = c(20, 40)[1]
lam = 0.57 # For exp link, complex beta, censoring rate 20%
# Hall-Wellner family control for the baseline MRL func m_0(t) 
D_case = data.frame(D1 = c(-0.5, -1/3), D2 = c(0.5, 1))[2,]
D1 = D_case$D1 ; D2 = D_case$D2  


test_type = c('TypeI', 'Power')[2] # For the power analysis 
test_int = seq(0.1, 0.9, by = 0.01)  # The range of testing of interest
c_null = c('known', 'unknown')[2]  # When we do not know the parameters in the null, choose c_null = unknown
 
###############
# Sample Code #
###############

set.seed(1)  
sim_set = simudata(n=n, 
                   alpha0 = alpha0,
                   t_int = t_Interval, # \mathcal{S}
                   model = model,
                   type = type,
                   v = v, 
                   K = K, # number of basis to generate Z(s)
                   lam = lam, # for censoring time 
                   D1 = D1, D2 = D2  # Hall-Wellner Family to generate baseline MRL m_0(t)
  )
  
dt = sim_set$dat # Data generation
psi0 = sim_set$psi0 # True beta(s) function -> Functional coefficient
mean(dt$delta)

################
# Fit the GFLMRL
################
fit = fn_gflmrl(data=dt, norder=4, t_Interval=t_Interval, model = model, type = type)
res_est = fit$est_all
  
est_alpha = res_est[1:2]    # \hat \alpha
est_beta = res_est[-c(1:2)] # \hat \beta
 
  
####################### 
# Bootstrapping for SE  
####################### 
boot_fit = boot_se_ci_paral(B=B, data=dt, t_Interval=t_Interval,
                            model = model, norder = 4, type = type) # Parallel in Mac

se_alp = boot_fit$se_alp
se_beta = boot_fit$se_beta

########################
# Hypothesis testings #
########################

# 1. Homogeneity H0: \beta(s) = c
ht_homo = func_testing_constant(data=dt, t_int = t_Interval,
                               boot_coef = boot_fit$boot_est,
                               est_coef = est_beta,
                               alpha_nom = 0.05,
                               c_null = c_null, 
                               test_range = range(test_int)
)

rej_homo_ind = as.numeric(ht_homo$rej)
pval_homo = as.numeric(ht_homo$pval)
c(rej_homo_ind, pval_homo)


# 2. Linearity H0: \beta(s) = c0 + c1 s 
ht_linear = func_testing_linear(data=dt, t_int = t_Interval,
                             boot_coef = boot_fit$boot_est,
                             est_coef = est_beta,
                             model = model,
                             alpha_nom = 0.05,
                             c_null = c_null, 
                             test_range = range(test_int)
)
  
rej_lin_ind = as.numeric(ht_linear$rej)
pval_lin = as.numeric(ht_linear$pval)
c(rej_lin_ind, pval_lin)
