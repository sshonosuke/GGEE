##----------------------------------------------------------------##
##    One-shot simulation experiment for prediction performance   ##
##----------------------------------------------------------------##
rm(list=ls())

library(bindata)
library(lme4)
library(glmertree)
library(PRROC)
source("GGEE-function.R")
source("Competitors-function.R")
logistic <- function(x){ (1+exp(-x))^(-1) }


## settings
scenario <- 1     # 1, 2 or 3
TT <- 10    # 10 or 20
n <- 180   
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

new.zz <- mvrnorm(n, c(0, 0), CC)
new.X <- cbind(1, new.zz)
new.x1 <- new.zz[,1]
new.x2 <- new.zz[,2]


## regression coefficients
if(scenario==1){
  Beta <- cbind(c(0, -2, 0), c(-1, 1, -2), c(1, 1, 2))   
  gg <- sample(1:G, n, replace=T)
  iBeta <- Beta[,gg]
}

if(scenario==2){
  Beta <- cbind(c(0, -2, 0), c(-1, 1, -2), c(1, 1, 2))   
  gg <- sort( rep(1:G, n/G) )
  iBeta <- Beta[, gg] + matrix(runif(n*(p+1), -0.5, 0.5), (p+1), n)  
}

if(scenario==3){
  iBeta <- rbind(runif(n, -0.2, 0.2), runif(n, -2, 2), runif(n, 0, 2)) 
}


tMu <- matrix(NA, n, TT)
new.tMu <- c()
for(i in 1:n){
  tMu[i,] <- as.vector(X[i,,]%*%iBeta[,i])
  new.tMu[i] <- as.vector(new.X[i,]%*%iBeta[,i])
}
Prob <- logistic(tMu) 
new.Prob <- logistic(new.tMu) 


## true correlation structure (equi-correlation) 
tR <- matrix(0.5, TT, TT)
diag(tR) <- 1


## data generation 
Y <- matrix(NA, n, TT)
for(i in 1:n){
  Y[i,] <- as.vector(rmvbin(1, margprob=Prob[i,], sigma=tR))
}



## fitting 
# GGEE-EX
sel <- CV(Y, X, Gm=7, cor="EC")
hG <- sel$opt
fit <- GGEE(Y, X, G=hG, cor="EC")
hg1 <- fit$g 
hBeta1 <- t(fit$Beta[,hg1])
pp1 <- logistic( apply(new.X*hBeta1, 1, sum) )

# GGEE-US
fit <- GGEE(Y, X, G=hG, cor="UN")
hg2 <- fit$g 
hBeta2 <- t(fit$Beta[,hg2])
pp2 <- logistic( apply(new.X*hBeta2, 1, sum) )

# random coefficient (RC) 
yy <- as.vector(t(Y))
xx <- matrix(NA, N, p+1)
for(k in 1:(p+1)){
  xx[,k] <- as.vector(t(X[,,k]))
}

ID <- sort(rep(1:n, TT))
fit.re <- glmer(yy~xx[,-1]+(xx[,-1]|ID), family="binomial")
hBeta3 <- as.matrix(coef(fit.re)$ID) 
pp3 <- logistic( apply(new.X*hBeta3, 1, sum) )


# growth mixture model (GMM)
fit.GMM <- GMM(yy, xx, ID, G=hG)
mpp <- logistic( new.X%*%fit.GMM$Beta )
pp4 <- apply(fit.GMM$ZP*mpp, 1, sum)


# tree based model (GLMMT) 
dd <- data.frame(y=as.vector(Y), x1=as.vector(x1), x2=as.vector(x2), ID=ID)
tree.fit <- glmertree(y~x1+x2|ID|(x1+x2), data=dd, family="binomial")
newdata <- data.frame(x1=new.x1, x2=new.x2, ID=1:n)
pp5 <- predict(tree.fit, newdata=newdata, type="response")


# pair-wise lasso (PWL)
fit.PWL <- PWL(yy, xx, ID)
pp6 <- logistic( apply(new.X*fit.PWL, 1, sum) )


# Summary 
Pred <- cbind(pp1, pp2, pp3, pp4, pp5, pp6)
sqrt( apply((Pred-new.Prob)^2, 2, mean) )    # RMSE

