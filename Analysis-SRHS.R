##------------------------------------------------------------##
##            Analysis of health status dataset               ##
##------------------------------------------------------------##
rm(list_ls())


# load dataset
library(LMest)
data(data_SRHS_long)
?data_SRHS_long

data <- data_SRHS_long
ID <- unique(data$id)
n <- max(ID)
TT <- length(data$id)/n


# response 
YY <- matrix(NA, n, TT)
for(i in 1:n){
  YY[i,] <- data$srhs[ID==i]
}
Y <- ifelse(YY>=3, 1, 0)   # binary response 


# covariate
Time <- rep(1:TT, n)
XX <- cbind(1, data$gender-1, ifelse(data$race==2, 1, 0), ifelse(data$race==3, 1, 0),
           ifelse(data$education==4, 1, 0), ifelse(data$education==5, 1, 0),
           data$age, data$age^2, ifelse(Time==2, 1, 0), ifelse(Time==3, 1, 0),
           ifelse(Time==4, 1, 0), ifelse(Time==5, 1, 0), ifelse(Time==6, 1, 0),
           ifelse(Time==7, 1, 0), ifelse(Time==8, 1, 0))

p <- dim(XX)[2]
X <- array(NA, c(n, TT, p))
for(k in 1:p){
  X[,,k] <- t(matrix(XX[,k], TT, n))
}
dimnames(X)[[1]] <- 1:n
dimnames(X)[[2]] <- paste0("time", 1:TT)
dimnames(X)[[3]] <- c("Intercept", "Gender", "Black", "Other", "someCollege", 
                   "College", "Age", "Age2", paste0("time-effect", 2:TT))



## analysis
source("GGEE-function.R")
set.seed(1)

cv <- CV(Y, X, Gm=10, num=10, cor="UN", print=T)
hG <- cv$opt
fit <- GGEE(Y, X, G=hG, maxit=1000, cor="UN", print=T)


#  coefficients and standard errors 
hg <- fit$g
hBeta <- fit$Beta
hC <- fit$Cor
AV <- GGEE.AV(Y, X, hBeta, hg, hC)
SD <- sqrt(apply(AV, 1, diag))


# time effect 
hg <- fit$g
int <- fit$Beta[1,]
TE <- rbind(0, fit$Beta[-(1:8),])
TE <- t(t(TE) + int)

matplot(TE, xaxt="n", type="l", col=1:hG, lty=1, ylab="Time Effect", xlab="Time")
for(k in 1:hG){
  points(1:TT, TE[,k], pch=8, col=k, cex=1)
}
axis(1, c(1,3,5,7,9), c(1,3,5,7,9))
legend("bottomright", paste("Group", 1:hG), col=1:hG, lty=1, ncol=3)



# group-wise age effect
xx <- seq(40, 70, by=5)
trend <- outer(xx, hBeta[7,]) + outer(xx^2, hBeta[8,])
trend.sd <- sqrt(outer(xx^2, (SD^2)[7,]) + outer(xx^4, (SD^2)[8,]))

matplot(trend, xaxt="n", type="l", col=1:hG, lty=1, main="Age effect", ylab="Effect", xlab="Age", ylim=c(-6, 6), lwd=2)
axis(1, at=1:length(xx), labels=xx)
for(g in 1:hG){
  points((trend+1.96*trend.sd)[,g], col=g, type="l", lty=2)
  points((trend-1.96*trend.sd)[,g], col=g, type="l", lty=2)
}
legend("bottomright", paste("Group", 1:hG), col=1:hG, lty=1, ncol=3)




# group-wise probability
Prob <- matrix(NA, TT, hG)
Prob.sd <- matrix(NA, TT, hG)
for(g in 1:hG){
  subY <- Y[hg==g,]
  Prob[,g] <- apply(subY, 2, mean)
  Prob.sd[,g] <- sqrt(Prob[,g]*(1-Prob[,g])/dim(subY)[1])
}

matplot(Prob, xaxt="n", type="l", col=1:hG, lty=1, main="Average probablity", 
        ylab="Probability", xlab="Time", ylim=c(-0.1, 1), lwd=2)
for(k in 1:hG){
  points(1:TT, Prob[,k], pch=8, col=k, cex=1)
}
axis(1, c(1,3,5,7,9), c(1,3,5,7,9))
legend("bottomright", paste("Group", 1:hG), col=1:hG, lty=1, ncol=3)
for(k in 1:hG){
  points(1:TT, Prob[,k]+1.96*Prob.sd[,g], type="l", col=k, lty=2)
  points(1:TT, Prob[,k]-1.96*Prob.sd[,g], type="l", col=k, lty=2)
}


