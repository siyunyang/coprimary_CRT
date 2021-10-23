library(nlme)
library(mvtnorm)
library(numDeriv)

# function to perform EM estimation with K=3 outcomes
EM.estim <- function(data, formula1,formula2,formula3, maxiter=500,epsilon=1e-4
                     , verbose=FALSE){
  # data: source data set
  # formula1: fit mixed model with outcome 1 and the treatment arm. e.g formula1 <-as.formula(  "out1 ~ arm")
  # formula2: fit mixed model with outcome 2 and the treatment arm. e.g formula1 <-as.formula(  "out2 ~ arm")
  # formula2: fit mixed model with outcome 3 and the treatment arm. e.g formula1 <-as.formula(  "out3 ~ arm")
  
  
  
  # fit mixed model to initialize parameters
  fm1 <- lme(formula1, random = ~ 1|cluster, data=data)
  fm2 <- lme(formula2, random = ~ 1|cluster, data=data)
  fm3 <- lme(formula3, random = ~ 1|cluster, data=data)
  
  K <- 3
  zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed, fm3$coefficients$fixed))
  beta1 = zeta[1:2]
  beta2 = zeta[3:4]
  beta3 = zeta[5:6]
  
  m <- as.numeric(table(data$cluster))
  
  s2phi1 <- VarCorr(fm1)[1,1]
  s2phi2 <- VarCorr(fm2)[1,1]
  s2phi3 <- VarCorr(fm3)[1,1]
  
  SigmaPhi <- diag(c(s2phi1, s2phi2, s2phi3))
  InvS2Phi <- solve(SigmaPhi)
  
  s2e1 <- VarCorr(fm1)[2,1]
  s2e2 <- VarCorr(fm2)[2,1]
  s2e3 <- VarCorr(fm3)[2,1]
  
  SigmaE <- diag(c(s2e1, s2e2, s2e3))
  InvS2E <- solve(SigmaE)
  
  out1 <-all.vars(formula1)[1]
  out2 <-all.vars(formula2)[1]
  out3 <-all.vars(formula3)[1]
  
  arm <- all.vars(formula1)[2]
  
  Y <- as.matrix(data[,c(out1,out2, out3)])
  ID <- as.numeric(data$cluster)
  n <- length(unique(ID))
  #X <- as.matrix(cbind(1, data[,"arm"])) # design matrix
  facz <- as.factor(data[,arm])
  z <- as.numeric(facz)-1
  X <- as.matrix(cbind(1, z)) # design matrix
  
  ESSphi1 <- matrix(0,n,K)
  ESSphi2 <- array(0,c(K,K,n))
  
  
  #maxiter=500
  #epsilon=1e-4
  delta = 2*epsilon
  max_modi = 20
  
  converge = 0
  
  # log likelihood
  
  loglik = function(theta){
    beta1 = theta[1:2]
    beta2 = theta[3:4]
    beta3 = theta[5:6]
    
    sphi11 = theta[7]
    sphi21 = theta[8]
    sphi22 = theta[9]
    sphi31 = theta[10]
    sphi32 = theta[11]
    sphi33 = theta[12]
    
    se11 = theta[13]
    se21 = theta[14]
    se22 = theta[15]
    se31 = theta[16]
    se32 = theta[17]
    se33 = theta[18]

    SigmaPhi = matrix(c(sphi11,sphi21,sphi31,sphi21,sphi22, sphi32, sphi31, sphi32,sphi33),3,3)
    SigmaE = matrix(c(se11,se21, se31, se21,se22,se32, se31,se32, se33),3,3)
    InvS2Phi <- solve(SigmaPhi)
    InvS2E <- solve(SigmaE)
    
    temp <- 0
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2, Xj%*%beta3)
      obs = c(t(residj))
      tm1 <- (m[j]-1)*log(det(SigmaE))+log(det(SigmaE+m[j]*SigmaPhi))
      InvSS2 <- solve(SigmaE+m[j]*SigmaPhi)-InvS2E
      Invj <- kronecker(diag(nrow=m[j]),InvS2E) + 
        kronecker(matrix(1,m[j],m[j]),InvSS2)/m[j]
      tm2 <- c(t(obs) %*% Invj %*% obs)
      temp <- temp-(tm1+tm2)/2
    }
    temp
  }
  thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
  LLold <- loglik(thetah)
  
  
  niter=1
  while((niter <= maxiter) & (abs(delta) > epsilon)){
    
    # Expectation step
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2, Xj%*%beta3)
      Vj <- solve(InvS2Phi + m[j]*InvS2E)
      Muj <- as.numeric(Vj %*% InvS2E %*% colSums(residj))
      Nujj <- Vj + tcrossprod(Muj)
      ESSphi1[j,] <- Muj
      ESSphi2[,,j] <- Nujj
    }
    
    # Maximization step - phi
    SigmaPhi <- apply(ESSphi2,1:2, sum)/n
    InvS2Phi <- solve(SigmaPhi)
    
    # Maximization step - zeta
    # Simplify the expression analytically, and obtain simple expression!
    XXt <- crossprod(X)
    Vzeta <-solve(kronecker(InvS2E, XXt))
    rzeta1 <- t(X)%*%(Y[,1]-ESSphi1[ID,1])
    rzeta2 <- t(X)%*%(Y[,2]-ESSphi1[ID,2])
    rzeta3 <- t(X)%*%(Y[,3]-ESSphi1[ID,3])
    
    zeta <- Vzeta %*% rbind(InvS2E[1,1]*rzeta1 + InvS2E[1,2]*rzeta2 + InvS2E[1,3]*rzeta3,
                            InvS2E[2,1]*rzeta1 + InvS2E[2,2]*rzeta2 + InvS2E[2,3]*rzeta3,
                            InvS2E[3,1]*rzeta1 + InvS2E[3,2]*rzeta2 + InvS2E[3,3]*rzeta3
    )
    zeta <- c(zeta)
    beta1 = zeta[1:2]
    beta2 = zeta[3:4]
    beta3 = zeta[5:6]
    
    # Maximization step - epsilon
    re <- Y - cbind(X%*%beta1, X%*%beta2, X%*%beta3)
    rss <- crossprod(re) + rowSums(sweep(ESSphi2,3,m,FUN="*"),dims=2) -
      crossprod(ESSphi1,rowsum(re,ID)) - crossprod(rowsum(re,ID),ESSphi1)
    SigmaE <- rss/sum(m)
    # SigmaE <- diag(diag(SigmaE))
    InvS2E <- solve(SigmaE)
    
    # whether the algorithm converges
    # thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),diag(SigmaE))
    thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
    #thetah = c(zeta,diag(SigmaPhi),SigmaPhi[1,2], c(SigmaE[!lower.tri(SigmaE)]))
    
    LLnew <- loglik(thetah)
    delta <- abs(LLnew - LLold)
    LLold <- LLnew
    converge = (abs(delta)<=epsilon)
    niter <- niter + 1
    if(verbose) cat(paste('iter=',niter),'\t',paste('param.error=',epsilon),'\t',paste('loglik=',LLnew),'\n');  
    
    #print(niter)
    #print(zeta)
    #print(SigmaPhi)
    #print(SigmaE)
    #print(LLnew)
  }
  param <- list(theta=list(zeta=zeta,SigmaE=SigmaE,SigmaPhi=SigmaPhi),loglik=LLnew,eps=epsilon,iter=niter)
  return(param) 
}







