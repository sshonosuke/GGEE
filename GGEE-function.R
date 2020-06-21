##----------------------------------------------------------------##
##   Grouped GEE analysis for balanced binary longitudinal data   ##
##----------------------------------------------------------------##

##    This code includes the following 4 main functions      ##
# GEE: computing GEE estimates given working correlation matrix
# GR.GEE: computing grouped GEE estimates given the number of groups
# GR.GEE.AV: asymptotic variance-covariance matrix of grouped GEE estimates
# CV: selection of the number of groups via cross-validation 

library(MASS)
library(flexmix)
logit <- function(x){ 1/(1+exp(-x)) }



##   function for GEE (given correlation matrix)   ##
# YY: (n,TT)-matrix of binary response 
# (n: the number of subjects; TT: the number of repeated measurements)
# XX: (n,TT,p)-array of covariates including an intercept
# (p: the number of covariates)
# RR: correlation matrix
# init.beta: initial values of regression coefficients
# maxit: the maximum number of iterations 
GEE <- function(YY, XX, RR, init.beta=NULL, maxit=100){
  th <- 0.0001
  ep <- 10^(-5)
  n <- dim(YY)[1]
  TT <- dim(YY)[2]
  p <- dim(XX)[3]
  if(is.null(init.beta)){ hB <- rep(0, p) }
  else{ hB <- init.beta }
  
  # Iteration
  for(k in 1:maxit){
    hB0 <- hB
    mat1 <- 0
    mat2 <- 0
    for(i in 1:n){
      Mu <- logit(as.vector(XX[i,,]%*%hB))
      S <- YY[i,]-Mu
      vv <- Mu*(1-Mu)
      vv[vv<th] <- th
      A <- diag(vv)
      V <- sqrt(A)%*%RR%*%sqrt(A) + th*diag(TT) 
      D <- A%*%XX[i,,]
      mat1 <- mat1+t(D)%*%ginv(V)%*%D
      mat2 <- mat2+t(D)%*%ginv(V)%*%S
    }
    hB <- hB+as.vector( ginv(mat1)%*%mat2 )
    dd <- sum( abs(hB-hB0) ) / sum( abs(hB0)+0.0001 )
    if( dd<ep ){ break }
  }
  return(hB)
}






##  Grouped GEE estimation   ##
# Y: (n,TT)-matrix of binary response 
# (n: the number of subjects; TT: the number of repeated measurements)
# X: (n,TT,p)-array of covariates including an intercept
# (p: the number of covariates)
# G: the number of groups
# maxit: the maximum number of iterations 
# cor: working correlation structure (following 4 options)
# "ID": independent; "EC": exchangeable; "AR": AR(1); "UN": unstructured
# init: method for computing initial values (following 2 options) 
# "sep-reg": separate estimation in each subject 
# "mix-reg": mixture of regressions

GR.GEE <- function(Y, X, G=2, maxit=1000, cor="EC", init="sep-reg"){
  th <- 0.0001
  ep <- 10^(-5)
  n <- dim(Y)[1]
  TT <- dim(Y)[2]
  p <- dim(X)[3]
  N <- n*TT
  
  ## correlation matrix
  # independent
  if(cor=="ID"){
    hR <- diag(rep(1, TT))
    halpha <- NULL
  }
  # exchangeable
  if(cor=="EC"){
    halpha <- 0
    IM <- matrix(NA, TT, TT)
    for(t1 in 1:TT){
      for(t2 in 1:TT){  IM[t1, t2] <- ifelse(abs(t1-t2)>0, 1, 0) }
    }
    Rc <- function(alpha){  alpha^IM  }
    hR <- Rc(halpha)
  }
  # AR(1)
  if(cor=="AR"){
    halpha <- 0
    IM <- matrix(NA, TT, TT)
    for(t1 in 1:TT){
      for(t2 in 1:TT){  IM[t1, t2] <- abs(t1-t2) }
    }
    Rc <- function(alpha){  alpha^IM  }
    hR <- Rc(halpha)
  }
  # unstructured
  if(cor=="UN"){
    hR <- diag(TT)
    halpha <- NULL
  }
  
  ## initial fit
  # separate estimation 
  if(init=="sep-reg"){
    Beta.init <- matrix(NA, n, p)
    for(i in 1:n){
      Beta.init[i,] <- coef(glm(Y[i,]~X[i,,-1]))
    }
    
    SSE <- c()
    init.fit <- list()
    for(j in 1:20){
      KM <- kmeans(Beta.init, G)
      SSE[j] <- sum(KM$withinss)
      init.fit[[j]] <- KM
    }
    init.fit <- init.fit[[which.min(SSE)]]
    hg <- init.fit$cluster
    hBeta <- t(init.fit$centers)
  }
  
  # mixture of regressions
  if(init=="mix-reg"){
    yy <- as.vector(t(Y))
    xx <- matrix(NA, N, p)
    for(k in 1:p){
      xx[,k] <- as.vector(t(X[,,k]))
    }
    fit.set <- list()
    LL <- c()
    for(j in 1:30){
      fit.set[[j]] <- flexmix(cbind(yy,1-yy)~xx[,-1], k=G, model=FLXMRglm(family="binomial"))
      LL[j] <- fit.set[[j]]@logLik
    }
    fit <- fit.set[[which.max(LL)]]
    hBeta <- parameters(fit)
    hp <- matrix(NA, n, G)
    for(k in 1:G){
      hp[,k] <- apply(matrix(fit@posterior$scaled[,k], TT, n), 2, mean)
    }
    hg <- apply(hp, 1, which.max)
  }
  
  ## Iterations
  for(k in 1:maxit){
    hBeta0 <- hBeta
    hg0 <- hg
    halpha0 <- halpha
    
    # update Beta
    for(j in 1:G){
      Ind <- (1:n)[hg==j]
      if( length(Ind)>1 ){
        hBeta[,j] <- GEE(YY=Y[Ind,], XX=X[Ind,,], RR=hR)
      }else{
        hBeta[,j] <- 0
      }
    }
    
    # update grouping assignment
    invR <- ginv(hR)
    QQ <- matrix(NA, n, G)
    for(j in 1:G){
      for(i in 1:n){
        Mu <- logit(as.vector(X[i,,]%*%hBeta[,j]))
        vv <- Mu*(1-Mu)
        vv[vv<th] <- th
        res <- Y[i,]-Mu
        QQ[i,j] <- t(res)%*%invR%*%res
      }
    }
    hg <- apply(QQ, 1, which.min)
    
    # update correlation parameter
    if(cor=="EC"){
      val <- 0
      for(i in 1:n){
        Mu <- logit(as.vector(X[i,,]%*%hBeta[,hg[i]]))
        vv <- Mu*(1-Mu)
        vv[vv<th] <- th
        invA.sq <- diag(1/sqrt(vv))
        u <- invA.sq%*%(Y[i,]-Mu)
        mat <- u%*%t(u)
        diag(mat) <- 0
        val <- val + sum(mat)
      }
      halpha <- val/(n*TT*(TT-1))
      halpha <- min(0.95, halpha)
      halpha <- max(0, halpha)
      hR <- Rc(halpha)
    }
    if(cor=="AR"){
      mat <- 0
      for(i in 1:n){
        Mu <- logit(as.vector(X[i,,]%*%hBeta[,hg[i]]))
        vv <- Mu*(1-Mu)
        vv[ vv<th ] <- th
        invA.sq <- diag(1/sqrt(vv))
        u <- invA.sq%*%(Y[i,]-Mu)
        mat <- mat + u%*%t(u)
      }
      opt <- function(a){  sum((Rc(a) - mat/n)^2)  }
      halpha <- optim(par=halpha, fn=opt, method="L-BFGS-B", lower=0, upper=0.95)$par
      hR <- Rc(halpha)
    }
    if(cor=="UN"){
      mat <- 0
      for(i in 1:n){
        Mu <- logit(as.vector(X[i,,]%*%hBeta[,hg[i]]))
        vv <- Mu*(1-Mu)
        vv[vv<th] <- th
        invA.sq <- diag(1/sqrt(vv))
        u <- invA.sq%*%(Y[i,]-Mu)
        mat <- mat + u%*%t(u)
      }
      hV <- mat/n
      DD <- diag( 1/sqrt(diag(hV)) )
      hR <- DD%*%hV%*%DD
    }
    
    # convergence check
    dd <- sum( abs(hBeta-hBeta0) ) / sum( abs(hBeta0) )
    if( dd<ep ){ break }
  }
  
  # Summary
  Res <- list(Beta=hBeta, g=hg, Cor=hR, itr=k)
  return(Res)
}








##  Asymptotic variance covariance matrix   ##
# Y: (n,TT)-matrix of binary response 
# (n: the number of subjects; TT: the number of repeated measurements)
# X: (n,TT,p)-array of covariates including an intercept
# (p: the number of covariates)
# hBeta: (p,G)-matrix of estimated regression coefficients
# hg: estimated grouping assignment 
# hR: estimated working correlation matrix

GR.GEE.AV <- function(Y, X, hBeta, hg, hR){
  th <- 0.0001
  G <- max(hg)
  n <- dim(Y)[1]
  p <- dim(X)[3]
  
  HH <- array(NA, c(n, p, p))
  MM <- array(NA, c(n, p, p))
  for(i in 1:n){
    Mu <- logit(as.vector(X[i,,]%*%hBeta[,hg[i]]))
    resid <- Y[i,]-Mu
    vv <- Mu*(1-Mu)
    vv[vv<th] <- th
    A <- diag(vv)
    V <- sqrt(A)%*%hR%*%sqrt(A)
    IV <- ginv(V)
    D <- A%*%X[i,,]
    HH[i,,] <- t(D)%*%IV%*%D
    MM[i,,] <- t(D)%*%IV%*%resid%*%t(resid)%*%IV%*%D
  }
  
  VV <- array(0, c(G, p, p))
  for(g in 1:G){
    if( sum(hg==g)>1 ){
      sH <- apply(HH[hg==g,,], c(2,3), sum)
      IH <- ginv(sH)
      sM <- apply(MM[hg==g,,], c(2,3), sum)
    }else{
      sH <- HH[hg==g,,]
      IH <- ginv(sH)
      sM <- MM[hg==g,,]
    }
    VV[g,,] <- IH%*%sM%*%IH
  }
  return(VV)
}







##  Cross-Validation for selecting the number of groups   ##
# Y: (n,TT)-matrix of binary response 
# (n: the number of subjects; TT: the number of repeated measurements)
# X: (n,TT,p)-array of covariates including an intercept
# (p: the number of covariates)
# Gm: the maximum number of group 
# num: the number of iterations to compute the criterion 
# cor: working correlations (4 options; see 'GR.GEE' function )
# print: criterion for each G is printed if 'T' 
# maxit: the maximum number of iterations

CV <- function(Y, X, Gm=5, num=10, cor="EC", print=F, maxit=30, init="sep-reg"){
  n <- dim(Y)[1]
  TT <- dim(Y)[2]
  p <- dim(X)[3]
  sM <- round(n/3)
  
  CV.val <- c()
  for(k in 2:Gm){
    val <- c()
    for(s in 1:num){
      Ind <- apply(matrix(sample(1:n, 2*sM), sM, 2), 2, sort)
      om <- (1:n)[-as.vector(Ind)]
      # fitting
      fit1 <- GR.GEE(Y[Ind[,1],], X[Ind[,1],,], G=k, cor=cor, maxit=maxit, init=init)
      hR1 <- fit1$Cor
      hBeta1 <- fit1$Beta
      fit2 <- GR.GEE(Y[Ind[,2],], X[Ind[,2],,], G=k, cor=cor, maxit=maxit, init=init)
      hR2 <- fit2$Cor
      hBeta2 <- fit2$Beta
      # grouping 
      nn <- length(om)
      group <- matrix(NA, nn, 2)
      for(h in 1:nn){
        Mu1 <- logit( X[om[h],,]%*%hBeta1 )
        resid1 <- Y[om[h],]-Mu1
        Mu2 <- logit( X[om[h],,]%*%hBeta2 )
        resid2 <- Y[om[h],]-Mu2
        group[h, 1] <- which.min(diag(t(resid1)%*%ginv(hR1)%*%resid1))
        group[h, 2] <- which.min(diag(t(resid2)%*%ginv(hR2)%*%resid2))
      }
      Ind.G1 <- Ind.G2 <- matrix(NA, nn, nn)
      for(h in 1:nn){
        Ind.G1[h,] <- ifelse(group[h,1]==group[,1], 1, 0)
        Ind.G2[h,] <- ifelse(group[h,2]==group[,2], 1, 0)
      }
      val[s] <- sum(ifelse(Ind.G1+Ind.G2==1, 1, 0))
    }
    CV.val[k] <- mean(val)
    if(print){ print(CV.val[k]) }
  }
  
  # summary
  opt.G <- which.min(CV.val)
  Res <- list(opt=opt.G, CV=CV.val)
  return(Res)
}




