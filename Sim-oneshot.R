##----------------------------------------------------------------##
##        One-shot simulation study with binary response          ##
##----------------------------------------------------------------##
set.seed(1)

library(bindata)
source("GGEE-function.R")

logistic <- function(x){ (1+exp(-x))^(-1) }

# settings
TT <- 10      # the number of repeated measurements 
n <- 180      # the number of subjects 
ID <- 1:n
N <- n*TT
G <- 3        # the number of groups


# covariate
CC <- matrix(c(1, 0.4, 0.4, 1), 2, 2)     # variance-covariance matrix between two covariates
zz <- mvrnorm(n*TT, c(0, 0), CC)
x1 <- matrix(zz[,1], n, TT)
x2 <- matrix(zz[,2], n, TT)
X <- array(NA, c(n, TT, 3))
X[,,1] <- 1
X[,,2] <- x1
X[,,3] <- x2

# true regression coefficients
Beta <- cbind(c(0, -2, 0), c(0, 1, -2), c(0, 1, 2))  

# true grouping
gg <- sort( rep(1:G, n/G) )
 
# true marginal probability
tMu <- matrix(NA, n, TT)
for(i in 1:n){
  tMu[i,] <- as.vector(X[i,,]%*%Beta[,gg[i]])
}
Prob <- logistic(tMu) 

# true correlation structure (exchangeable correlation) 
tR <- matrix(0.5, TT, TT)
diag(tR) <- 1


# data generation 
Y <- matrix(NA, n, TT)
for(i in 1:n){
  Y[i,] <- as.vector(rmvbin(1, margprob=Prob[i,], sigma=tR))
}

# selection of the number of groups 
hG <- CV(Y, X, Gm=7, num=5, cor="EC", print=T)$opt

# grouped GEE with the estimated number of groups
fit <- GR.GEE(Y, X, G=hG, cor="EC")
fit$Beta
fit$g

# asymptotic variance 
AV <- GR.GEE.AV(Y, X, fit$Beta, fit$g, fit$Cor)
AV[1,,]
AV[2,,]
AV[3,,]

