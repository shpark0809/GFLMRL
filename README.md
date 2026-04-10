# GFLMRL : Generalized functional linear mean residual life model with right-censored data

R scripts for implementing the methods developed under the generalized functional linear mean residual life (GFLMRL) model.
The model examines the effects of both scalar and functional predictors on the mean residual life (MRL) function using a quasi-partial score estimation procedure.

## Descriptions 
1. `fn_data.R`: Generates data under the GFLMRL model, including the code for generating the functional covariate.
2. `fn_gflmrl.R`: Main script for fitting the GFLMRL model. It can estimate both scalar and functional coefficients. Different spline approximations can be used by replacing `ns()` with other functions, depending on the application.
3. `fn_gflmrl_lik.R`: Selects the appropriate link function and the optimal number of basis functions using the likelihood criterion.
4. `fn_gflmrl_testing.R`: Performs nonparametric resampling-based hypothesis testing to assess whether the functional coefficient has a specific parametric form.
5. `sample_GFLMRL.R`: Provides an example workflow for generating data and demonstrating all the methods developed in this paper.
   

## Related Paper 
- Title: Generalized functional linear mean residual life model with right-censored data
- Authors: Seohyeon Park, Myeonggyun Lee, Mengling Liu
- Journal: 
- Link:
- Doi:
- How to cite:
 
