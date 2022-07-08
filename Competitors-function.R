library(glmnet)
library(ape) 



##----------------------------------------------------------------##
##      Growth curve mixture for binary longitudinal data         ##
##----------------------------------------------------------------##
GMM <- function(yy, xx, ID, G=3, maxit=100){
  # preparation 
  p <- dim(xx)[2]
  n <- max(ID)
  # initial values
  Beta.init <- matrix(NA, n, p)
  for(i in 1:n){
    Beta.init[i,] <- coef(glm(yy[ID==i]~xx[ID==i,-1]))
  }
  SSE <- c()
  init.fit <- list()
  for(j in 1:20){
    KM <- kmeans(Beta.init, G)
    SSE[j] <- sum(KM$withinss)
    init.fit[[j]] <- KM
  }
  init.fit <- init.fit[[which.min(SSE)]]
  hBeta <- t(init.fit$centers)
  hPP <- rep(1/G, G)
  
  # replications
  for(j in 1:maxit){
    hBeta.old <- hBeta
    # update classification 
    prob <- logistic(xx%*%hBeta)
    LL <- dbinom(yy, 1, prob=prob, log=T)
    ZP <- matrix(NA, n, G)
    mLL <- matrix(NA, n, G)
    Z <- c()
    for(i in 1:n){
      mLL <- apply(LL[ID==i,], 2, sum)+log(hPP)
      ZP[i,] <- exp(mLL - max(mLL)) / sum(exp(mLL - max(mLL)))
      Z[i] <- sample(1:G, 1, prob=ZP[i,])
    }
    # update parameters 
    for(k in 1:G){
      ww <- ZP[ID, k]
      hBeta[,k] <- coef( glm(yy~xx[,-1], family="quasibinomial", weights=ww) )
    }
    # convergence 
    dd <- sum(abs(hBeta-hBeta.old))/sum(abs(hBeta.old))
    if(dd<0.0001){ break }
  }
  
  # BIC
  prob <- logistic(xx%*%hBeta)
  LL <- dbinom(yy, 1, prob=prob, log=T)
  mLL <- matrix(NA, n, G)
  for(i in 1:n){
    mLL[i,] <- apply(LL[ID==i,], 2, sum)+log(hPP)
  }
  ML <- apply(mLL, 1, max)
  loglike <- ML + log(apply(exp(mLL - ML), 1, sum))
  N <- length(ID)
  BIC <- (-2)*sum(loglike) + log(N)*(p*G+G)
  
  # output
  return(list(Beta=hBeta, ZP=ZP, BIC=BIC))
}





##----------------------------------------------------------------##
##      Pair-wise penalization for binary longitudinal data       ##
##----------------------------------------------------------------##
PWL <- function(yy, xx, ID){
  n <- max(ID)
  p <- dim(xx)[2]
  logistic <- function(x){ (1+exp(-x))^(-1) }
  
  # initial fit 
  Beta.init <- matrix(NA, n, p)
  for(i in 1:n){
    Beta.init[i,] <- coef(glm(yy[ID==i]~xx[ID==i,-1], family="binomial"))
  }
  
  # minimum spaning tree
  MST <- mst(dist(Beta.init))
  H <- c()
  for(i in 1:n){
    for(j in 1:n){
      if(i<j){
        if(MST[i,j]==1){
          h <- rep(0, n)
          h[i] <- 1
          h[j] <- -1
          H <- rbind(H, h)
        }
      }
    }
  }
  
  # design matrix
  XX <- matrix(0, dim(xx)[1], n*p)
  for(i in 1:n){
    sub <- (p*(i-1)+1):(p*i)
    XX[ID==i, sub] <- xx[ID==i,]
  }
  XX <- as(XX, "sparseMatrix")
  
  # transformed design matrix
  HH <- rbind(H, rep(1/n, n))
  HH.mat <- matrix(0, p*n, p*n)
  for(i in 1:n){
    if(i<n){
      sub1 <- (1+p*(i-1)):(i*p)
      sel <- which(HH[i,]==1)
      sub2 <- (1+p*(sel-1)):(sel*p)
      diag(HH.mat[sub1, sub2]) <- 1
      sel <- which(HH[i,]==(-1))
      sub2 <- (1+p*(sel-1)):(sel*p)
      diag(HH.mat[sub1, sub2]) <- (-1)
    }
    if(i==n){
      for(k in 1:p){
        sub1 <- (i-1)*p + k  
        sub2 <- p*(1:n) - (p-k)
        HH.mat[sub1, sub2] <- 1/n
      }
    }
  }
  HH.mat <- as(HH.mat, "sparseMatrix")
  invHH <- solve(HH.mat)
  XH <- XX%*%invHH
  
  # fitting
  pen <- rep(1, n*p)
  pen[((n-1)*p+1):(n*p)] <- 0
  fit <- glmnet(x=XH, y=yy, family="binomial", nlambda=30, intercept=F, penalty.factor=pen)
  
  # BIC
  BIC <- c()
  L <- dim(fit$beta)[2]
  for(j in 1:L){
    mu <- logistic( as.vector(XH%*%as.vector(fit$beta[,j])) )
    BIC[j] <- (-2)*sum(dbinom(yy, 1, mu, log=T)) + log(N)*fit$df[j]
  }
  opt <- which.min(BIC)
  
  # estimates
  hBeta <- as.vector(invHH%*%as.vector(fit$beta[,opt]))
  hBeta <- t(matrix(hBeta, p, n))
  return(hBeta)
}


