library(tidyverse)
library(tableone)
library(foreach)
theme_set(theme_grey(base_size = 12))
theme_set(theme_bw(base_size = 16))

# source functions
source("powerSampleCal_varCluster_ttest_SY.R")
source("EM_standard_generalf_K2_function.R")

# read in data
KDD2 <- read.csv("Secondary outcomes_K-DPP trial.csv")
KDD1 <- read.csv("Primary outcome_K-DPP trial.csv")

# cluster info
n <- length(unique(KDD1$cluster))
sum_m <- KDD1 %>% group_by(cluster) %>%  dplyr::summarize(count = n())
m <- sum_m$count
CV <- sd(m)/mean(m)
mu = round(mean(m))

# 2 coprimary outcomes

KDD1$sysbp_c <- KDD1$sysbp2-KDD1$sysbp0
KDD1$diabp_c <- KDD1$diabp2-KDD1$diabp0


myVars <- c("pcs0", "mcs0", "sf6d0", "pf0", "sf0", "sysbp_c", "diabp_c", "idrs_c")
tab3 <- CreateTableOne(vars = myVars, strata = "arms0" , data = KDD1)
print(tab3,  formatOptions = list(big.mark = ","))



# sample size calculation
K=2
formula1 <-as.formula(  "sysbp_c  ~ arms0")
formula2 <-as.formula(  "diabp_c ~ arms0")
KDD.comp <- KDD1[-which(is.na(KDD1$diabp_c)),]
param <- EM.estim(data = KDD.comp, formula1, formula2, maxiter=500, epsilon=1e-4, verbose=FALSE)

rho0 <- diag( param$theta$SigmaPhi) /(diag(param$theta$SigmaE) +diag( param$theta$SigmaPhi))
rho01 <- diag(rho0)

for(row in 1:K )
{
  for(col in 1:K)
  {
    if(row != col){
      rho01[row,col] <- param$theta$SigmaPhi[row,col]/prod(sqrt(diag(param$theta$SigmaE) +diag( param$theta$SigmaPhi)))
    }
  }
}

rho2 <- diag(K)
rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <-
  (param$theta$SigmaPhi[lower.tri(param$theta$SigmaPhi)]+param$theta$SigmaE[upper.tri(param$theta$SigmaE)])/prod(sqrt(diag(param$theta$SigmaE) +diag( param$theta$SigmaPhi)))

mar_var <- (diag(param$theta$SigmaE) +diag( param$theta$SigmaPhi))
mar_sd <- sqrt(mar_var) # 9.8, 13.4

# IUT sample size  
nraw <- calSampleSize_ttestIU_var(betas=c(0.3, 0.3)* mar_sd, deltas=c(0,0), vars= mar_var, rho01, rho2, m=mu, r=0.5, K=2, alpha=0.05, power=0.8, cv=CV)
# 49
nraw <- calSampleSize_ttestIU_var(betas=c(0.5, 0.5)* mar_sd, deltas=c(0,0), vars= mar_var, rho01, rho2, m=mu, r=0.5, K=2, alpha=0.05, power=0.8, cv=CV)
# 19
nraw <- calSampleSize_ttestIU_var(betas=c(0.4, 0.4)* mar_sd, deltas=c(0,0), vars= mar_var, rho01, rho2, m=mu, r=0.5, K=2, alpha=0.05, power=0.8, cv=CV)
# 29

# omnibus sample size
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



calSampleSize_omnibus <- function(beta, vars,rho01,rho2,cv, m,r=0.5, K=2, alpha=0.05, dpower=0.8){
n=2*K+2
sigmaz.square <- r*(1-r) 
power=0
while(power<dpower){
  n <-n+2
omega <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square, m, K)
#Variance Inflation of coprimary outcomes
tau <- n*t(beta) %*% solve(omega) %*% beta
Fscore <- qf(1-0.05, df1=K, df2=n-2*K, ncp=0, lower.tail = TRUE, log.p = FALSE)
power <-1- pf(Fscore, df1=K, df2=n-2*K, ncp=tau, lower.tail = TRUE, log.p = FALSE)
}
return(n)
}
nraw <- calSampleSize_omnibus(beta=matrix(c(0.3, 0.3)* mar_sd, 2,1), vars= mar_var, rho01, rho2,cv=CV, m=17, r=0.5, K=2, alpha=0.05, dpower=0.8)
# 48
nraw <- calSampleSize_omnibus(beta=matrix(c(0.4, 0.4)* mar_sd, 2,1), vars= mar_var, rho01, rho2,cv=CV, m=17, r=0.5, K=2, alpha=0.05, dpower=0.8)
# 30
nraw <- calSampleSize_omnibus(beta=matrix(c(0.5, 0.5)* mar_sd, 2,1), vars= mar_var, rho01, rho2,cv=CV, m=17, r=0.5, K=2, alpha=0.05, dpower=0.8)
# 22

# homogeneity sample size


calSampleSize_homo <- function(beta, vars,rho01,rho2,cv, m,r=0.5, K=2, alpha=0.05, dpower=0.8){
  n=2*K+2
  sigmaz.square <- r*(1-r) 
  power=0
  L <- matrix(c(1,-1),1, 2)
  
  while(power<dpower){
    n <-n+2
    omega <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square, m, K)
    #Variance Inflation of coprimary outcomes
    tau <- n*t(L%*%beta) %*% solve(L%*%omega%*%t(L)) %*%(L%*% beta)
    Fscore <- qf(1-0.05, df1=K-1, df2=n-2*K-1, ncp=0, lower.tail = TRUE, log.p = FALSE)
    power <-1- pf(Fscore, df1=K-1, df2=n-2*K-1, ncp=tau, lower.tail = TRUE, log.p = FALSE)
    
  }
  return(n)
}
nraw <- calSampleSize_homo(beta=matrix(c(0.3, 0.7)* mar_sd, 2,1), vars= mar_var, rho01, rho2,cv=CV, m=17, r=0.5, K=2, alpha=0.05, dpower=0.8)
# 24
nraw <- calSampleSize_homo(beta=matrix(c(0.35, 0.7)* mar_sd, 2,1), vars= mar_var, rho01, rho2,cv=CV, m=17, r=0.5, K=2, alpha=0.05, dpower=0.8)
# 38
nraw <- calSampleSize_homo(beta=matrix(c(0.4, 0.7)* mar_sd, 2,1), vars= mar_var, rho01, rho2,cv=CV, m=17, r=0.5, K=2, alpha=0.05, dpower=0.8)
# 76



### heat map vary effect size
BETA1 <- seq(0.2, 0.8, 0.02)
BETA2 <- seq(0.2, 0.8, 0.02)

res=foreach(beta1=BETA1,.combine=cbind)%:%foreach(beta2=BETA2,.combine=cbind)%do%  
  
  { betas <- c(beta1, beta2)* mar_sd
    nraw <- calSampleSize_ttestIU_var(betas, deltas=c(0,0), vars= mar_var, rho01, rho2, m=mu, r=0.5, K=2, alpha=0.05, power=0.8, cv=CV)

    row <- c(beta1, beta2, nraw)
    data.frame(row)
  }
res <- t(res)
colnames(res) <- c( "sysbp_c","diabp_c", "n")


#################################
###### heat map vary ICC values #
##################################

sigmaz.square <- r*(1-r) 
#####function to construct covariance matrix Sigma_E for Y_i########
constrRiE <- function(rho01,rho2,K,vars)
{ rho0k <- diag(rho01)
SigmaE_Matrix <- diag((1-rho0k)*vars)
for(row in 1:K)
{
  for(col in 1:K)
  {
    if(row != col){
      SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]-rho01[row,col])
    }
  }
}
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

calCovbetas_eq <- function(vars,rho01,rho2, sigmaz.square,  m, K){
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


calCorWks <-  function(vars,rho01,rho2, sigmaz.square, cv, m, K)
{
  top <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square,  m, K)
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
RHO0 <- seq(0.01, 0.09, 0.001)
RHO1R <- seq(0.1, 1.5, 0.01) # ratio of RHO1/RHO0
RHO2 <- c(0.4, 0.79)   # rho_2 
K <-k<- 2 # number of outcomes
m <- 17
CV <- cv <- 0.19
r <- 0.5 # proportion of treated
n <- 60
set.seed(321)
vars <- c(178.4, 96)

res<-NULL
beta <- matrix(c(0.3*sqrt(vars[1]), 0.3*sqrt(vars[2])), 2,1)



# power omnibus test


res= foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    rho_0k <- c(rho0, 2.4 *rho0)
    #rho1R <- rho0/2
    
    rho01 <- diag(rho_0k)
    rho01[lower.tri(rho01)] <- rho01[upper.tri(rho01)] <- rho1R*rho0
    rho2 <- diag(K)
    rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <- rho_2
    
    omega <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square, m, K)
    #Variance Inflation of coprimary outcomes
    tau <- n*t(beta) %*% solve(omega) %*% beta
    Fscore <- qf(1-0.05, df1=K, df2=n-2*K, ncp=0, lower.tail = TRUE, log.p = FALSE)
    power <-1- pf(Fscore, df1=K, df2=n-2*K, tau, lower.tail = TRUE, log.p = FALSE)
    row <- c(m,cv, rho0,rho1R, rho_2,tau, power)
    data.frame(row)
  }


res <- t(res)
colnames(res) <- c("m","CV", "rho_0","rho1R","rho_2","tau", "power")

res <- as.data.frame(res)
range(res$power)  # 0.76, 1.00
ll <- round(range(res$power),2)[1]

res$rho_2 <- factor(res$rho_2 , levels = c("0.4", "0.79"), 
                    labels = c(expression(paste(rho[2],"=0.4")),expression(paste(rho[2],"=0.79"))))


p <- ggplot(res, aes(rho_0, rho1R, fill= power)) + 
  geom_tile()+
  scale_fill_gradient(low="yellow", high="brown",breaks=seq(ll,1,0.04))+  
  #scale_fill_gradient(low = "#56B1F7",high ="#132B43" )+  
  geom_point(aes(x= 0.05,y= 1.4),colour="black", size=2)+
  labs(x = expression(rho[0]^1), y= expression(rho[1]/rho[0]^1))+
  facet_wrap(~ rho_2, scales = "free_x",labeller = label_parsed)
ggsave(paste("heatmap_omnibus",Sys.Date(),".png", sep= "_"), plot=p, width = 8, height = 6, units =  "in")

# power heatmap homogeneity test
L <- matrix(c(1,-1),1, 2)
#beta <- matrix(c(0.3*sqrt(vars[1]), 0.7*sqrt(vars[2])), 2,1)
beta <- matrix(c(0.35*sqrt(vars[1]), 0.7*sqrt(vars[2])), 2,1)

res= foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- c(178.4, 96)
    rho_0k <- c(rho0, 2.4 *rho0)
    #rho1R <- rho0/2
    
    rho01 <- diag(rho_0k)
    rho01[lower.tri(rho01)] <- rho01[upper.tri(rho01)] <- rho1R*rho0
    rho2 <- diag(K)
    rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <- rho_2
    
    omega <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square, m, K)
    #Variance Inflation of coprimary outcomes
    tau <- n*t(L%*%beta) %*% solve(L%*%omega%*%t(L)) %*%(L%*% beta)
    Fscore <- qf(1-0.05, df1=K-1, df2=n-2*K-1, ncp=0, lower.tail = TRUE, log.p = FALSE)
    power <-1- pf(Fscore, df1=K-1, df2=n-2*K-1, ncp=tau, lower.tail = TRUE, log.p = FALSE)
    row <- c(m,cv, rho0,rho1R, rho_2,tau, power)
    data.frame(row)
  }


res <- t(res)
colnames(res) <- c("m","CV", "rho_0","rho1R","rho_2","tau", "power")

res <- as.data.frame(res)
range(res$power) # 0.25, 0.98
res$rho_2 <- factor(res$rho_2 , levels = c("0.4", "0.79"), 
                    labels = c(expression(paste(rho[2],"=0.4")),expression(paste(rho[2],"=0.79"))))


p <- ggplot(res, aes(rho_0, rho1R, fill= power)) + 
  geom_tile()+
   scale_fill_gradient(low="yellow", high="brown",breaks=seq(0.25,1,0.1))+  
  geom_point(aes(x= 0.05,y= 1.4),colour="black", size=2)+
  labs(x = expression(rho[0]^1), y= expression(rho[1]/rho[0]^1))+
  facet_wrap(~ rho_2, scales = "free_x",labeller = label_parsed)
ggsave(paste("heatmap_homogeneity",Sys.Date(),".png", sep= "_"), plot=p, width = 8, height = 6, units =  "in")


# power heatmap IU test
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


res= foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- c(178.4, 96)
    rho_0k <- c(rho0, 2.4 *rho0)
    #rho1R <- rho0/2
    
    rho01 <- diag(rho_0k)
    rho01[lower.tri(rho01)] <- rho01[upper.tri(rho01)] <- rho1R*rho0
    rho2 <- diag(K)
    rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <- rho_2
    beta =c(0.3*sqrt(vars[1]), 0.3*sqrt(vars[2]))
    power <- calPower_ttestIU_var(betas=beta,deltas=c(0,0),vars,rho01,rho2,N=n,r=0.5,m,K,alpha=0.05, cv=cv)
    row <- c(m,cv, rho0,rho1R, rho_2, power)
    data.frame(row)
  }


res <- t(res)
colnames(res) <- c("m","CV", "rho_0","rho1R","rho_2", "power")

res <- as.data.frame(res)
range(res$power) #0.67, 0.99

res$rho_2 <- factor(res$rho_2 , levels = c("0.4", "0.79"), 
                labels = c(expression(paste(rho[2],"=0.4")),expression(paste(rho[2],"=0.79"))))


p <- ggplot(res, aes(rho_0, rho1R, fill= power)) + 
  geom_tile()+
   scale_fill_gradient(low="yellow", high="brown",breaks=seq(0.66,1,0.05))+ 
  geom_point(aes(x= 0.05,y= 1.4),colour="black", size=2)+
  labs(x = expression(rho[0]^1), y= expression(rho[1]/rho[0]^1))+ 
  facet_wrap(~ rho_2, scales = "free_x",labeller = label_parsed)
p
ggsave(paste("heatmap_IUT",Sys.Date(),".png", sep= "_"), plot=p, width = 8, height = 6, units =  "in" )


