
load.lib<-c("rmutil","tidyverse", "foreach", "reshape2","mvtnorm","nlme","numDeriv", "parallel")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE ,repo="http://cran.rstudio.com/")
sapply(load.lib, require, character=TRUE)


############################
### parameter set up ########
#############################
K <- 2 # number of outcomes
M <- c(60) # mean cluster size
CV <- seq(0.01,0.8, 0.01)
RHO0 <- c(0.01, 0.05, 0.1)  # rho_0k=seq(RHO0,  0.1, length= k)
RHO1R <- c(0.5, 0.75, 0.9) # ratio of RHO1/RHO0
RHO2 <- c(0.2, 0.5)   # rho_2 
r <- 0.5 # proportion of treated
set.seed(321)

res<-NULL

#########################################################
#### Model Test #########################
#########################################################
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

calCovbetas_eq <- function(vars,rho01,rho2, sigmaz.square, m, K){
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
res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  

{ #print(c(k,m,rho0, rho_2, cv))
  vars <- rep(1, K)
  rho_0k <- rep(rho0, 2)
  rho01 <- diag(rho_0k)
  rho01[lower.tri(rho01)] <- rho01[upper.tri(rho01)] <- rho1R*rho0
  rho2 <- diag(K)
  rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <- rho_2
  
   top <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square, m, K)
   bottom <- calCovbetas_eq(vars,rho01,rho2, sigmaz.square, m, K)
   #Variance Inflation of coprimary outcomes
   re_cp <- c(top/bottom)[1]
   # Variance Inflation of a single outcome
   alpha <- (1-rho0)/rho0
   lmd <- m/(m+alpha)
   re_sg <- 1-cv^2*lmd*(1-lmd)
   row <- c(m,cv, rho0,rho1R, rho_2,re_cp, 1/re_sg)
   data.frame(row)
}


res <- t(res)
colnames(res) <- c("m","cv", "rho_0","rho1R","rho_2","re_cp", "re_sg")

#########################
###  Fig 1 ##############
#########################
library(ggplot2)
library(gridExtra)
theme_set(theme_grey(base_size = 12))
theme_set(theme_bw(base_size = 16))

res<-as.data.frame(res)
res$grp <- paste(res$rho_0, res$rho_2,sep='_')
res$rho_0 <- as.factor(res$rho_0)
res$rho1R <- as.factor(res$rho1R)
res$rho_2 <- as.factor(res$rho_2)

res1 <- subset(res, rho_2==0.2 & rho1R==0.5)
res1_long <- gather(res1, type, RE, re_cp:re_sg, factor_key=TRUE)
res1_long$type = factor(res1_long$type, levels=c("re_sg", "re_cp"))

p1 <-res1 %>% ggplot( aes(x=cv, y=re_sg, linetype=rho_0, group=rho_0)) +
  geom_line(size=1, color="#999999") +
  ggtitle("a) single outcome") +
  ylab("Variance Inflation")+  xlab("CV")+ labs(linetype = bquote( rho[0] ))+
  ylim(1, 1.18)

p2 <-res1_long %>% ggplot( aes(x=cv, y=RE, colour=type, 
                          group=interaction(rho_0, type))) +
  geom_line(size=1,aes(linetype=rho_0)) +
  ggtitle(expression("b) co-primary outcomes," ~ rho[1]/rho[0] ~ "=0.5"))+
  ylab("Variance Inflation")+  xlab("CV")+ labs(colour = "# of outcomes",linetype= bquote( rho[0] ))+
  scale_color_manual(values=c( "#999999", "black"), guide='none')+
  ylim(1, 1.18)

res2 <- subset(res, rho_2==0.2 & rho1R==0.75)
res2_long <- gather(res2, type, RE, re_cp:re_sg, factor_key=TRUE)
res2_long$type = factor(res2_long$type, levels=c("re_sg", "re_cp"))

p3 <-res2_long %>% ggplot( aes(x=cv, y=RE, colour=type, 
                               group=interaction(rho_0, type))) +
  geom_line(size=1,aes(linetype=rho_0)) +
  ggtitle(expression("c) co-primary outcomes," ~ rho[1]/rho[0] ~ "=0.75"))+
  ylab("Variance Inflation")+  xlab("CV")+ labs(colour = "# of outcomes",linetype= bquote( rho[0] ))+
  scale_color_manual(values=c( "#999999", "black"), guide='none')+
  ylim(1, 1.18)


res3 <- subset(res, rho_2==0.2 & rho1R==0.9)
res3_long <- gather(res3, type, RE, re_cp:re_sg, factor_key=TRUE)
res3_long$type = factor(res3_long$type, levels=c("re_sg", "re_cp"))

p4 <-res3_long %>% ggplot( aes(x=cv, y=RE, colour=type, 
                               group=interaction(rho_0, type))) +
  geom_line(size=1,aes(linetype=rho_0)) +
  ggtitle(expression("d) co-primary outcomes," ~ rho[1]/rho[0] ~ "=0.9"))+
  ylab("Variance Inflation")+  xlab("CV")+ labs(colour = "# of outcomes",linetype= bquote( rho[0] ))+
  scale_color_manual(values=c( "#999999", "black"), guide='none')+
  ylim(1, 1.18)


p <- grid.arrange(p1,p2, p3,p4, ncol=2, nrow=2)
p

ggsave(paste("numerical_m60_rho2_0.2",Sys.Date(),".png", sep= "_"), plot=p, width = 12, height = 8, units =  "in")

###############################
#### vary rho2 = 0.5, Web Fig 2
###############################
res1 <- subset(res, rho_2==0.5 & rho1R==0.5)
res1_long <- gather(res1, type, RE, re_cp:re_sg, factor_key=TRUE)
res1_long$type = factor(res1_long$type, levels=c("re_sg", "re_cp"))

p1 <-res1 %>% ggplot( aes(x=cv, y=re_sg, linetype=rho_0, group=rho_0)) +
  geom_line(size=1, color="#999999") +
  ggtitle("a) single outcome") +
  ylab("Variance Inflation")+  xlab("CV")+ labs(linetype = bquote( rho[0] ))+
  ylim(1, 1.18)

p2 <-res1_long %>% ggplot( aes(x=cv, y=RE, colour=type, 
                               group=interaction(rho_0, type))) +
  geom_line(size=1,aes(linetype=rho_0)) +
  ggtitle(expression("b) co-primary outcomes," ~ rho[1]/rho[0] ~ "=0.5"))+
  ylab("Variance Inflation")+  xlab("CV")+ labs(colour = "# of outcomes",linetype= bquote( rho[0] ))+
  scale_color_manual(values=c( "#999999", "black"), guide='none')+
  ylim(1, 1.18)

res2 <- subset(res, rho_2==0.5 & rho1R==0.75)
res2_long <- gather(res2, type, RE, re_cp:re_sg, factor_key=TRUE)
res2_long$type = factor(res2_long$type, levels=c("re_sg", "re_cp"))

p3 <-res2_long %>% ggplot( aes(x=cv, y=RE, colour=type, 
                               group=interaction(rho_0, type))) +
  geom_line(size=1,aes(linetype=rho_0)) +
  ggtitle(expression("c) co-primary outcomes," ~ rho[1]/rho[0] ~ "=0.75"))+
  ylab("Variance Inflation")+  xlab("CV")+ labs(colour = "# of outcomes",linetype= bquote( rho[0] ))+
  scale_color_manual(values=c( "#999999", "black"), guide='none')+
  ylim(1, 1.18)


res3 <- subset(res, rho_2==0.5 & rho1R==0.9)
res3_long <- gather(res3, type, RE, re_cp:re_sg, factor_key=TRUE)
res3_long$type = factor(res3_long$type, levels=c("re_sg", "re_cp"))

p4 <-res3_long %>% ggplot( aes(x=cv, y=RE, colour=type, 
                               group=interaction(rho_0, type))) +
  geom_line(size=1,aes(linetype=rho_0)) +
  ggtitle(expression("d) co-primary outcomes," ~ rho[1]/rho[0] ~ "=0.9"))+
  ylab("Variance Inflation")+  xlab("CV")+ labs(colour = "# of outcomes",linetype= bquote( rho[0] ))+
  scale_color_manual(values=c( "#999999", "black"), guide='none')+
  ylim(1, 1.18)

p <- grid.arrange(p1,p2, p3,p4, ncol=2, nrow=2)
p

ggsave(paste("numerical_m60_rho2_0.5",Sys.Date(),".png", sep= "_"), plot=p, width = 12, height = 8, units =  "in")




#########################################
#### omnibus test K=2: Fig 2 ######
#########################################
K <- 2 # number of outcomes
M <- c(60) # mean cluster size
CV <- c(0, 0.5, 0.9)
RHO0 <- c(0.1)  # rho_0k=seq(RHO0,  0.1, length= k)
RHO1R <- seq(0.1, 0.9, 0.01) # ratio of RHO1/RHO0

RHO2 <- 0.2  # rho_2 
r <- 0.5 # proportion of treated
n <- 30
set.seed(321)
beta <- matrix(rep(0.3, K), 2,1)

res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, 2)
    #rho1R <- rho0/2
    
    rho01 <- diag(rho_0k)
    rho01[lower.tri(rho01)] <- rho01[upper.tri(rho01)] <- rho1R*rho0
    rho2 <- diag(K)
    rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <- rho_2
    
    omega <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square, m, K)
    #Variance Inflation of coprimary outcomes
    tau <- n*t(beta) %*% solve(omega) %*% beta
    Fscore <- qf(1-0.05, df1=K, df2=n-2*K, ncp=0, lower.tail = TRUE, log.p = FALSE)
    power <-1- pf(Fscore, df1=K, df2=n-2*K, ncp=tau, lower.tail = TRUE, log.p = FALSE)
    row <- c(m,cv, rho0,rho1R, rho_2,tau, power)
    data.frame(row)
  }


res <- t(res)
colnames(res) <- c("m","CV", "rho_0","rho1R","rho_2","tau", "power")

# plot
library(ggplot2)
library(gridExtra)
theme_set(theme_grey(base_size = 12))
theme_set(theme_bw(base_size = 16))

res<-as.data.frame(res)
res$CV <- as.factor(res$CV)


p1 <-res %>% ggplot( aes(x=rho1R, y=tau, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle("omnibus test") +
  labs(x = expression(rho[1]/rho[0]), y= expression(tau))



p2 <-res %>% ggplot( aes(x=rho1R, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("b) Vary" ~ rho[1]))+
  
  labs(x = expression(rho[1]/rho[0]), y= "Power")+
  theme(legend.position = "none")


######## 
RHO0 <-seq(0.05, 0.3, 0.01) 
#RHO1R <- RHO0/2
RHO1R <- 0.5 # ratio of RHO1/RHO0
RHO2 <- 0.2  # rho_2 

res<-NULL
beta <- matrix(rep(0.3, K), 2,1)
res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, 2)
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

res<-as.data.frame(res)
res$CV <- as.factor(res$CV)


p3 <-res %>% ggplot( aes(x=rho_0, y=tau, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle("omnibus test") +
  labs(x = expression(rho[0]), y= expression(tau))



p4 <-res %>% ggplot( aes(x=rho_0, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("a) Vary" ~ rho[0]))+ylim(c(0.27,0.92))+
  labs(x = expression(rho[0]), y= "Power")+
  theme(legend.position = "none")


#####################
RHO0 <-0.1
RHO1R <- 0.5 # ratio of RHO1/RHO0
RHO2 <- seq(0.1,0.9,0.01)   # rho_2 

res<-NULL
beta <- matrix(rep(0.3, K), 2,1)
res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, 2)
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

res<-as.data.frame(res)
res$CV <- as.factor(res$CV)

p5 <-res %>% ggplot( aes(x=rho_2, y=tau, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle("omnibus test") +
  labs(x = expression(rho[0]), y= expression(tau))

p6 <-res %>% ggplot( aes(x=rho_2, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("c) Vary" ~ rho[2]))+
  labs(x = expression(rho[2]), y= "Power")+  theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=12),legend.text = element_text(size=12))

p <- grid.arrange(p4,p2, p6 ,widths=c(1,1,1.4),ncol=3, nrow=1)
p
ggsave(paste("numerical_omnibus",Sys.Date(),".png", sep= "_"), plot=p, width = 12, height = 4, units =  "in")

####################################################
### Homogeneity test K=2, ########## Web Fig 3 ######
####################################################

### vary rho0
K <- 2 # number of outcomes
M <- c(60) # mean cluster size
CV <- c(0,0.5, 0.9)
RHO0 <- seq(0.01, 0.3, 0.01)  # rho_0k=seq(RHO0,  0.1, length= k)
RHO1R <- c(0.5) # ratio of RHO1/RHO0

RHO2 <- c(0.2)   # rho_2 
r <- 0.5 # proportion of treated
n <- 30
set.seed(321)
L <- matrix(c(1,-1),1, 2)
res<-NULL
beta <- matrix(c(0.3,0.7), 2,1)


res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, K)
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

# plot
library(ggplot2)
library(gridExtra)
theme_set(theme_grey(base_size = 12))
theme_set(theme_bw(base_size = 16))

res<-as.data.frame(res)
res$CV <- as.factor(res$CV)

p1 <-res %>% ggplot( aes(x=rho_0, y=tau, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("a) Vary" ~ rho[0]))+
  labs(x = expression(rho[0]), y= expression(Tau))+
  theme(legend.position = "none")

p2 <-res %>% ggplot( aes(x=rho_0, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("a) Vary" ~ rho[0]))+
  labs(x = expression(rho[0]), y= "Power")+
  theme(legend.position = "none")

### vary rho1
RHO0 <- 0.1  # rho_0k=seq(RHO0,  0.1, length= k)
#RHO1R <- RHO0/2
RHO1R <- seq(0.1,0.9,0.01) # ratio of RHO1/RHO0
res<-NULL
beta <- matrix(c(0.3,0.7), 2,1)
res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, K)
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

res<-as.data.frame(res)
res$CV <- as.factor(res$CV)

p3 <-res %>% ggplot( aes(x=rho1R, y=tau, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("b) Vary" ~ rho[1]))+
  labs(x = expression(rho[1]/rho[0]), y= expression(Tau))+
  theme(legend.position = "none")

p4 <-res %>% ggplot( aes(x=rho1R, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("b) Vary" ~ rho[1]))+
  labs(x = expression(rho[1]/rho[0]), y= "Power")+
  theme(legend.position = "none")

### vary rho2
RHO0 <- 0.1  # rho_0k=seq(RHO0,  0.1, length= k)
#RHO1R <- RHO0/2
RHO1R <- 0.5 # ratio of RHO1/RHO0
RHO2 <- seq(0.1,0.9,0.01)   # rho_2 

res<-NULL
beta <- matrix(c(0.3,0.7), 2,1)

res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, K)
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

res<-as.data.frame(res)
res$CV <- as.factor(res$CV)

p5 <-res %>% ggplot( aes(x=rho_2, y=tau, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("c) Vary" ~ rho[2]))+
  labs(x = expression(rho[2]),  y= expression(Tau))+
  theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=12),legend.text = element_text(size=12))

p6 <-res %>% ggplot( aes(x=rho_2, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("c) Vary" ~ rho[2]))+
  labs(x = expression(rho[2]), y= "Power")+
  theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=12),legend.text = element_text(size=12))

p <- grid.arrange(p2,p4, p6 ,widths=c(1,1,1.3),ncol=3, nrow=1)
p


ggsave(paste("numerical_homogeneity",Sys.Date(),".png", sep= "_"), plot=p, width = 12, height = 4, units =  "in")




######################################
############ IUT power, Web Fig 4 ####
######################################
r=0.5
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
  
K <- 2 # number of outcomes
M <- c(60) # mean cluster size
CV <- c(0, 0.5, 0.9)
RHO0 <- c(0.1)  # rho_0k=seq(RHO0,  0.1, length= k)
RHO1R <- seq(0.1, 0.9, 0.01) # ratio of RHO1/RHO0

RHO2 <- 0.2  # rho_2 
r <- 0.5 # proportion of treated
n <- 30
set.seed(321)
beta <- rep(0.3, K)

res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, 2)
    rho01 <- diag(rho_0k)
    rho01[lower.tri(rho01)] <- rho01[upper.tri(rho01)] <- rho1R*rho0
    rho2 <- diag(K)
    rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <- rho_2
    power <- calPower_ttestIU_var(betas=beta,deltas=c(0,0),vars,rho01,rho2,N=n,r=0.5,m,K,alpha=0.05, cv=cv)
    row <- c(m,cv, rho0,rho1R, rho_2, power)
    data.frame(row)
  }


res <- t(res)
colnames(res) <- c("m","CV", "rho_0","rho1R","rho_2", "power")


# plot
library(ggplot2)
library(gridExtra)
theme_set(theme_grey(base_size = 12))
theme_set(theme_bw(base_size = 16))

res<-as.data.frame(res)
res$CV <- as.factor(res$CV)


p2 <-res %>% ggplot( aes(x=rho1R, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("b) Vary" ~ rho[1]))+
  labs(x = expression(rho[1]/rho[0]), y= "Power")+
  theme(legend.position = "none")
######## 
RHO0 <-seq(0.05, 0.3, 0.01) 
#RHO1R <- RHO0/2
RHO1R <- 0.5 # ratio of RHO1/RHO0
RHO2 <- 0.2  # rho_2 

res<-NULL
beta <- rep(0.3, K)
res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, 2)
    rho01 <- diag(rho_0k)
    rho01[lower.tri(rho01)] <- rho01[upper.tri(rho01)] <- rho1R*rho0
    rho2 <- diag(K)
    rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <- rho_2
    
    power <- calPower_ttestIU_var(betas=beta,deltas=c(0,0),vars,rho01,rho2,N=n,r=0.5,m,K,alpha=0.05, cv=cv)
    
    row <- c(m,cv, rho0,rho1R, rho_2, power)
    data.frame(row)
  }

res <- t(res)
colnames(res) <- c("m","CV", "rho_0","rho1R","rho_2", "power")


res<-as.data.frame(res)
res$CV <- as.factor(res$CV)

p1 <-res %>% ggplot( aes(x=rho_0, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("a) Vary" ~ rho[0]))+
  
  labs(x = expression(rho[0]), y= "Power")+
  theme(legend.position = "none")

#####################
RHO0 <- 0.1  # rho_0k=seq(RHO0,  0.1, length= k)
RHO1R <- 0.5 # ratio of RHO1/RHO0
RHO2 <- seq(0.1,0.9,0.01)   # rho_2 
res<-NULL
beta <- rep(0.3, K)
res= foreach(m=M,.combine=cbind)%:%foreach(cv=CV,.combine=cbind)%:%foreach(rho0=RHO0,.combine=cbind)%:%foreach(rho1R=RHO1R,.combine=cbind)%:%foreach(rho_2=RHO2,.combine=cbind)%do%  
  
  { #print(c(k,m,rho0, rho_2, cv))
    vars <- rep(1, K)
    rho_0k <- rep(rho0, 2)
    rho01 <- diag(rho_0k)
    rho01[lower.tri(rho01)] <- rho01[upper.tri(rho01)] <- rho1R*rho0
    rho2 <- diag(K)
    rho2[lower.tri(rho2)] <- rho2[upper.tri(rho2)] <- rho_2
    
    #Variance Inflation of coprimary outcomes
    power <- calPower_ttestIU_var(betas=beta,deltas=c(0,0),vars,rho01,rho2,N=n,r=0.5,m,K,alpha=0.05, cv=cv)
    
    row <- c(m,cv, rho0,rho1R, rho_2, power)
    data.frame(row)
  }


res <- t(res)
colnames(res) <- c("m","CV", "rho_0","rho1R","rho_2", "power")

res<-as.data.frame(res)

res$CV <- as.factor(res$CV)

p3 <-res %>% ggplot( aes(x=rho_2, y=power, linetype=CV, color=CV)) +
  geom_line(size=1.5) +
  ggtitle(expression("c) Vary" ~ rho[2]))+
  
  labs(x = expression(rho[2]), y= "Power")+
  theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=12),legend.text = element_text(size=12))

p <- grid.arrange(p1,p2,p3 ,widths=c(1,1,1.3),ncol=3, nrow=1)
p


ggsave(paste("numerical_IUT",Sys.Date(),".png", sep= "_"), plot=p, width = 12, height = 4, units =  "in")


#################################################
### Cluster size distribution: Web Figure 1 #####
#################################################
CV <- seq(0.2,0.8, 0.2)
m=60
kappa <- 1/CV^2
theta <- m*CV^2
sd <- sqrt(kappa* theta^2)
png(file="cluster_size.png",  width=1200, height=900,   pointsize = 2.2, res=280)
par(mar = c(18,18,4,4)) # Set the margin on all sides to 2
plot(0, 0, xlim = c(0, 200), ylim = c(0, 0.04), type = "n", cex.axis=4,xaxt='n')
title(ylab="Density", line=8, cex.lab=4)
title(xlab="Cluster Size",line=8, cex.lab=4)
axis( side=1, at=c(0,50,100,150,200), labels=c(0,50,100,150,200), cex.axis=4,  padj = 1)
for (i in (1:length(CV)))
  curve(dgamma(x, shape = kappa[i], scale = theta[i]), from = 0, to = 200, col = i, lty=i, lwd = 2, add = TRUE)

cv.f <- factor(CV, levels= seq(0.2,0.8, 0.2),
               labels = c("CV = 0.2", "CV = 0.4", "CV = 0.6","CV = 0.8"))

colfill<-c(1:(1+length(levels(cv.f))))
legend("topright", levels(cv.f), col =colfill, lty =colfill, cex=4, lwd = 2)
dev.off()




