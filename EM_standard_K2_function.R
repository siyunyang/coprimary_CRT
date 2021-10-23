library(nlme)
library(mvtnorm)
library(numDeriv)

# function to perform EM estimation with K=2 outcomes
EM.estim <- function(data, formula1, formula2, maxiter=500,epsilon=1e-4
                     , verbose=FALSE){
  # data: source data set
  # formula1: fit mixed model with outcome 1 and the treatment arm. e.g formula1 <-as.formula(  "out1 ~ arm")
  # formula2: fit mixed model with outcome 2 and the treatment arm. e.g formula1 <-as.formula(  "out2 ~ arm")

  
  # fit mixed model to initialize parameters
  fm1 <- lme(formula1, random = ~ 1|cluster, data=data)
  fm2 <- lme(formula2, random = ~ 1|cluster, data=data)
  K <- 2
  zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed))
  beta1 = zeta[1:2]
  beta2 = zeta[3:4]
  
  m <- as.numeric(table(data$cluster))
  
  s2phi1 <- VarCorr(fm1)[1,1]
  s2phi2 <- VarCorr(fm2)[1,1]
  SigmaPhi <- diag(c(s2phi1, s2phi2))
  InvS2Phi <- solve(SigmaPhi)
  
  s2e1 <- VarCorr(fm1)[2,1]
  s2e2 <- VarCorr(fm2)[2,1]
  SigmaE <- diag(c(s2e1, s2e2))
  InvS2E <- solve(SigmaE)
  
  out1 <-all.vars(formula1)[1]
  out2 <-all.vars(formula2)[1]
  arm <- all.vars(formula1)[2]
  
  Y <- as.matrix(data[,c(out1,out2)])
  ID <- as.numeric(data$cluster)
  n <- length(unique(ID))
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
    sphi11 = theta[5]
    sphi12 = theta[6]
    sphi22 = theta[7]
    se11 = theta[8]
    se12 = theta[9]
    se22 = theta[10]
    SigmaPhi = matrix(c(sphi11,sphi12,sphi12,sphi22),2,2)
    SigmaE = matrix(c(se11,se12,se12,se22),2,2)
    InvS2Phi <- solve(SigmaPhi)
    InvS2E <- solve(SigmaE)
    
    temp <- 0
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
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
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
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
    zeta <- Vzeta %*% rbind(InvS2E[1,1]*rzeta1 + InvS2E[1,2]*rzeta2,
                            InvS2E[2,1]*rzeta1 + InvS2E[2,2]*rzeta2)
    zeta <- c(zeta)
    beta1 = zeta[1:2]
    beta2 = zeta[3:4]
    
    # Maximization step - epsilon
    re <- Y - cbind(X%*%beta1, X%*%beta2)
    rss <- crossprod(re) + rowSums(sweep(ESSphi2,3,m,FUN="*"),dims=2) -
      crossprod(ESSphi1,rowsum(re,ID)) - crossprod(rowsum(re,ID),ESSphi1)
    SigmaE <- rss/sum(m)
    # SigmaE <- diag(diag(SigmaE))
    InvS2E <- solve(SigmaE)
    
    # whether the algorithm converges
    # thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),diag(SigmaE))
    thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
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







