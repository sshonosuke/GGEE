##------------------------------------------------------------##
##             One-shot simulation experiment for             ##
##       estimation, selection and confidence interval        ##
##------------------------------------------------------------##
rm(list=ls())

library(bindata)
source("GGEE-function.R")
logistic <- function(x){ (1+exp(-x))^(-1) }


## settings
TT <- 10    # 10 or 20
n <- 180    # 180 or 270 
ID <- 1:n
N <- n*TT
G <- 3   # number of groups


## covariates 
p <- 2   # number of covariates
CC <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
zz <- mvrnorm(n*TT, c(0, 0), CC)

x1 <- matrix(zz[,1], n, TT)
x2 <- matrix(zz[,2], n, TT)
X <- array(NA, c(n, TT, p+1))
X[,,1] <- 1
X[,,2] <- x1
X[,,3] <- x2

## regression coefficients
Beta <- cbind(c(0, -2, 0), c(-1, 1, -2), c(1, 1, 2))   # true regression coefficients
gg <- sort( rep(1:G, n/G) )
tBeta <- Beta   # true coefficients without intercept 

tMu <- matrix(NA, n, TT)
for(i in 1:n){
  tMu[i,] <- as.vector(X[i,,]%*%Beta[,gg[i]])
}
Prob <- logistic(tMu) 


## true correlation structure (equi-correlation) 
tR <- matrix(0.5, TT, TT)
diag(tR) <- 1


## data generation 
Y <- matrix(NA, n, TT)
for(i in 1:n){
  Y[i,] <- as.vector(rmvbin(1, margprob=Prob[i,], sigma=tR))
}



## fitting (fixed G) 
meth <- c("ID", "EC", "AR", "UN", "crude")
M <- length(meth)
hBeta <- array(NA, c(M, p+1, G))
hg <- matrix(NA, M, n)

# grouped GEE
fit1 <- GGEE(Y, X, G=3, cor="ID")
hBeta[1,,] <- fit1$Beta
hg[1,] <- fit1$g

fit2 <- GGEE(Y, X, G=3, cor="EC")
hBeta[2,,] <- fit2$Beta
hg[2,] <- fit2$g

fit3 <- GGEE(Y, X, G=3, cor="AR")
hBeta[3,,] <- fit3$Beta
hg[3,] <- fit3$g

fit4 <- GGEE(Y, X, G=3, cor="UN")
hBeta[4,,] <- fit4$Beta
hg[4,] <- fit4$g


# crude clustering 
Ind <- fit1$init$cluster
for(k in 1:G){
  yy <- as.vector(Y[Ind==k,])
  xx <- X[Ind==k,,]
  xx2 <- cbind(as.vector(xx[,,2]), as.vector(xx[,,3]))
  hBeta[5,,k] <- coef(glm(yy~xx2, family="binomial"))
}


# index adjustment
Ind <- matrix(NA, G, M)
for(k in 1:G){
  for(m in 1:M){
    Ind[k, m] <- which.min(apply((hBeta[m,,k]-tBeta)^2, 2, sum))
  }
}

# estimation and classification error
SE <- matrix(NA, M, G)
CE <- rep(NA, M)
dimnames(SE)[[1]] <- names(CE) <- meth
for(m in 1:M){
  SE[m, ] <- apply((hBeta[m,,]-tBeta[,Ind[,m]])^2, 2, sum)
  CE[m] <- mean(ifelse(Ind[hg[m,], m]!=gg, 1, 0))
}


# result 
apply(SE, 1, mean)    # squared error
CE[-M]*100    # classification error 




## selection of G 
CV(Y, X, Gm=7, cor="ID")$opt
CV(Y, X, Gm=7, cor="EC")$opt
CV(Y, X, Gm=7, cor="AR")$opt
CV(Y, X, Gm=7, cor="UN")$opt




## confidence intervals (EX working correlation)
AV1 <- GGEE.AV(Y, X, hBeta=fit2$Beta, hg=fit2$g, hR=fit2$Cor)
AV2 <- GGEE.AV.boot(Y, X, G=3, B=100, cor="EC")
SD1 <- matrix(NA, p+1, G)
SD2 <- matrix(NA, p+1, G)
for(g in 1:G){ 
  SD1[,g] <- sqrt(diag(AV1[g,,])) 
  SD2[,g] <- sqrt(diag(AV2[g,,]))
}
fit2$Beta+1.96*SD1  # upper bound of confidence interval based on analytical formula   
fit2$Beta-1.96*SD1  # lower bound of confidence interval based on analytical formula   
fit2$Beta+1.96*SD2  # upper bound of confidence interval based on cluster bootstrap  
fit2$Beta-1.96*SD2  # lower bound of confidence interval based on cluster bootstrap  



## confidence intervals (AR working correlation)
AV1 <- GGEE.AV(Y, X, hBeta=fit3$Beta, hg=fit3$g, hR=fit3$Cor)
AV2 <- GGEE.AV.boot(Y, X, G=3, B=100, cor="AR")
SD1 <- matrix(NA, p+1, G)
SD2 <- matrix(NA, p+1, G)
for(g in 1:G){ 
  SD1[,g] <- sqrt(diag(AV1[g,,])) 
  SD2[,g] <- sqrt(diag(AV2[g,,]))
}
fit3$Beta+1.96*SD1  # upper bound of confidence interval based on analytical formula   
fit3$Beta-1.96*SD1  # lower bound of confidence interval based on analytical formula   
fit3$Beta+1.96*SD2  # upper bound of confidence interval based on cluster bootstrap  
fit3$Beta-1.96*SD2  # lower bound of confidence interval based on cluster bootstrap  



## confidence intervals (US working correlation)
AV1 <- GGEE.AV(Y, X, hBeta=fit4$Beta, hg=fit4$g, hR=fit4$Cor)
AV2 <- GGEE.AV.boot(Y, X, G=3, B=100, cor="UN")
SD1 <- matrix(NA, p+1, G)
SD2 <- matrix(NA, p+1, G)
for(g in 1:G){ 
  SD1[,g] <- sqrt(diag(AV1[g,,])) 
  SD2[,g] <- sqrt(diag(AV2[g,,]))
}
fit4$Beta+1.96*SD1  # upper bound of confidence interval based on analytical formula   
fit4$Beta-1.96*SD1  # lower bound of confidence interval based on analytical formula   
fit4$Beta+1.96*SD2  # upper bound of confidence interval based on cluster bootstrap  
fit4$Beta-1.96*SD2  # lower bound of confidence interval based on cluster bootstrap  

