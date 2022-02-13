library(mvtnorm)
########################################################################################################################################################
##Power/Sample Size Calculation based on the t test##
# INPUT
# betas: (beta_1,...,beta_K), the vector of treatment effect for (1st,...,Kth) endpoints
# deltas: (delta_1,...,delta_K), the vector of non-inferiority margins, when delta_1 = ... = delta_K = 0,
#          superiority tests are performed on all endpoints
# vars: (var_1,...,var_K), the vector of marginal variance for (1st,...,Kth) endpoints
# rho01: a K by K dimensional matrix for the correlation parameters (rho0^k) and (rho1^kk')
# For rho01: 
#           the diagonal elements correspond to rho0^k's 
#           the off-diagonal elements correspond to (rho1^kk')'s 
#           For example, rho01[1,1] corresponds to rho0^1, which is the ICC for the first endpoint
#                        rho01[1,2] corresponds to rho1^12, which is the correlation of outcomes between subjects on the 1st and 2nd endpoints    
# rho2: a K by K dimensional matrix for the correlation parameters (rho2^kk')
# For rho2:
#           the diagonal elements are 1
#           the off-diagonal elements correspond to (rho2^kk')'s
#           For example, rho2[1,2] corresponds to rho2^12, which is the correlation of outcomes within same subject on the 1st and 2nd endpoints    
# N: number of clusters
# r: proportion of clusters in the intervention arm
# m: mean cluster size
# K: number of endpoints
# alpha: upper bound of type I error rates over the whole null space
# cv: coefficient of variation for variable cluster size
########################################################################################################################################################
####Function to Calculate Power Given Design Configurations based on the t test (intersection-union test)#######
####Critical values c_1,...,c_K are set to t_alpha, (1-alpha)th quantile of the t distribution with df = N-2K###
####if assuming equal cluster size. i.e. cv =0, the calPower_ttestIU function should be used
calPower_ttestIU_var <- function(betas,deltas,vars,rho01,rho2,N,r,m,K,alpha,cv)
{
  #variance of trt assignment
  sigmaz.square <- r*(1-r) 
  #####function to construct covariance matrix Sigma_E for Y_i########
  constrRiE <- function(rho01,rho2,K,vars)
  { rho0k <- diag(rho01)
    SigmaE_Matrix <- diag((1-rho0k)*vars)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]-rho01[row,col])
        }
      }
    }
   # check for matrix positive definite
    if( min( eigen(SigmaE_Matrix)$values) <= 1e-08 ) 
      {print("Warning: the resulting covariance matrix Sigma_E is not positive definite. Check the input of rho01 and rho2.")}
    return(SigmaE_Matrix)
  }
  
  #####function to construct covariance matrix Sigma_phi for Y_i########
  constrRiP <- function(rho01,K,vars)
  { rho0k <- diag(rho01)
  SigmaP_Matrix <- diag(rho0k*vars)
  for(row in 1:K )
  {
    for(col in 1:K)
    {
      if(row != col){
        SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
      }
    }
  }
   # check for matrix positive definite
  if( min( eigen(SigmaP_Matrix)$values) <= 1e-08 ) 
  {print("Warning: the resulting covariance matrix Sigma_phi is not positive definite. Check the input of rho01 and rho2.")}
  return(SigmaP_Matrix)
  }
  ####Define function to calculate covariance between betas#####
  calCovbetas <- function(vars,rho01,rho2,cv, sigmaz.square, m, K){
    sigmaE <- constrRiE(rho01,rho2,K,vars)
    sigmaP <- constrRiP(rho01,K,vars)
    tmp <- solve(diag(1,K)-cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP) ))
    covMatrix <- 1/(m*sigmaz.square)*(sigmaE+m*sigmaP)%*%tmp
    covMatrix <- (covMatrix +t(covMatrix))/2  # symmerize the off-diagonal
    return(covMatrix)  
  }
  
  
  ####Define function to calculate correlation between test statistics #####
  calCorWks <-  function(vars,rho01,rho2, sigmaz.square, cv, m, K)
  {
    top <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square, m, K)
    wCor <- diag(K)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          wCor[row,col] <- top[row,col]/sqrt(top[row,row]*top[col,col])
        }
      }
    }
    return(wCor)
  }
  #print(cv)
  sigmaks.sq <- diag( calCovbetas(vars,rho01,rho2, cv, sigmaz.square, m, K))
  meanVector <- sqrt(N)*(betas-deltas)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars,rho01,rho2, sigmaz.square,cv, m, K)
  criticalValue <- qt(p=(1-alpha), df=(N-2*K))
  pred.power <- pmvt(lower = rep(criticalValue,K),upper=rep(Inf,K),df = (N-2*K) , sigma = wCor,delta=meanVector)[1]
  return(pred.power)
}
########Function to Calculate Sample Size (Total Number of Clusters) based on the t Test (intersection-union test)################
#### if assuming equal cluster size. i.e. cv =0, the calSampleSize_ttestIU function should be used

calSampleSize_ttestIU_var <- function(betas,deltas,vars, rho01,rho2,m,r,K,alpha,power, cv)
{
  lowerBound <- 1
  upperBound <- 1000
  repeat{
    middle <- floor((lowerBound+upperBound)/2)
    power_temp <- calPower_ttestIU_var(betas,deltas,vars,rho01,rho2,N=middle,r,m,K,alpha,cv)
    if(power_temp < power)
    {
      lowerBound <- middle
    }
    if(power_temp > power)
    {
      upperBound <- middle
    }
    if(power_temp == power)
    {
      return(middle)
      break
    }
    if((lowerBound-upperBound) == -1)
    {
      return(upperBound)
      break
    }
  }
}



####Critical values c_1,...,c_K are set to t_alpha, (1-alpha)th quantile of the t distribution with df = N-2K###
#### equal cluster size
calPower_ttestIU <- function(betas,deltas,vars,rho01,rho2,N,r,m,K,alpha)
{
  #variance of trt assignment
  sigmaz.square <- r*(1-r) 
  
  ####Define function to calculate correlation between betas#####
  calCovbetas <- function(vars,rho01,rho2, sigmaz.square, m, K){
    rho0k <- diag(rho01)
    sigmak.square <-(1+(m-1)*rho0k)*vars/(m*sigmaz.square)
    covMatrix <- diag(sigmak.square)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          covMatrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]+(m-1)*rho01[row,col])/(m*sigmaz.square)
        }
      }
    }
    return(covMatrix)  
  }
  
  
  ####Define function to calculate correlation between test statistics #####
  calCorWks <-  function(vars,rho01,rho2, sigmaz.square, m, K)
  {
    top <- calCovbetas(vars,rho01,rho2, sigmaz.square, m, K)
    wCor <- diag(K)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          wCor[row,col] <- top[row,col]/sqrt(top[row,row]*top[col,col])
        }
      }
    }
    return(wCor)
  }
  sigmaks.sq <- diag(calCovbetas(vars,rho01,rho2, sigmaz.square, m, K))
  meanVector <- sqrt(N)*(betas-deltas)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars,rho01,rho2, sigmaz.square, m, K)
  criticalValue <- qt(p=(1-alpha), df=(N-2*K))
  pred.power <- pmvt(lower = rep(criticalValue,K),upper=rep(Inf,K),df = (N-2*K) , sigma = wCor,delta=meanVector)[1]
  return(pred.power)
}
########Function to Calculate Sample Size (Total Number of Clusters) based on the t Test (intersection-union test)################
# equal cluster size
calSampleSize_ttestIU <- function(betas,deltas,vars, rho01,rho2,m,r,K,alpha,power)
{
  lowerBound <- 1
  upperBound <- 1000
  repeat{
    middle <- floor((lowerBound+upperBound)/2)
    power_temp <- calPower_ttestIU(betas,deltas,vars,rho01,rho2,N=middle,r,m,K,alpha)
    if(power_temp < power)
    {
      lowerBound <- middle
    }
    if(power_temp > power)
    {
      upperBound <- middle
    }
    if(power_temp == power)
    {
      return(middle)
      break
    }
    if((lowerBound-upperBound) == -1)
    {
      return(upperBound)
      break
    }
  }
}
