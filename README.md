# coprimary_CRT

The R code for implementing the sample size and power calculations under a generalized linear mixed model (GLMM) framework as shown in "Power analysis for cluster randomized trials with continuous co-primary endpoints" by Yang et al.

List of Files:
1) powerSampleCal_varCluster_ttest.R = R program for conducting power/sample size calculation based on the t test
2) EM_standard_K2_function.R = R program to operationalize the EM algorithm for MLMM with K=2 coprimary endpoints
3) EM_standard_K3_function.R = R program to operationalize the EM algorithm for MLMM with K=3 coprimary endpoints
4) Illustration_PowerSampleSize_Calculation.pdf = PDF file for illustrating how to use R program powerSampleCal_varCluster_ttest.R to perform power/sample size calculation proposed in the paper

NOTES:  1) To use powerSampleCal_varCluster_ttest.r, please install and load package "mvtnorm". To use EM_standard_K2_function.R or EM_standard_K3_function.R, please install and load package "nlme" and "numDeriv".
        2) The power/sample size calculation is performed based on the intersenction-union test
        3) For the t test, the critical values are set to the common value t_alpha, which is the (1-alpha)th quantile of the t distribution with df = N-2K
     
