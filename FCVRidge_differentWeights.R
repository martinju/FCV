# Test script for doing focused type of weighted cross validation, for first ridge, later elastic net and others. 

# Idea: 
# Assign weights to each observation corresponding to the inverse of the distance to covariate value of interest (or other distance measure).
# Then perform leave-one-out cross-validation where the prediction error for each observation being predicted is weighted with these weights
# The model with the smallest test prediction error is chosen

# Experiment: 
# Sample data from linear regression model with several covariates of the same importance
# Define a set of models containing this true model
# Perform leave-one-out weighted cross validation (later 5-fold) and find the best model according to that criterion (AND with equal weighting)
# Predict the mu based on the two models 
# Repeat the procedure several times to see what procedure performs the best.

rm(list = ls())
# Defining a true model ####

beta0 <-rep(0.1,11)  #c(1,-0.1,0.4,0.1,-0.2)
p <- length(beta0)-1
sigma0 = 1

# Defining the simulation experiment #### 

noSim <- 1000
n <- 100
#Defining the ridge Leave-one-out crossvalidation 

# Defining the cross validation function ####
LOOCVRidge <- function(X,Y,focusX){
  #LOOCV does not need a fold partition
  n <- dim(X)[1]
  p <- dim(X)[2]-1

  lambdas <- seq(10^-3,10^4, length.out = 2000)
  noLambda <- length(lambdas)

  CVRegularMat <- matrix(NA,nrow=noLambda,ncol=n)
  CVFocusWeightedMat1 <- matrix(NA,nrow=noLambda,ncol=n)
  CVFocusWeightedMat2 <- matrix(NA,nrow=noLambda,ncol=n)
  CVFocusWeightedMat3 <- matrix(NA,nrow=noLambda,ncol=n)
  CVFocusWeightedMat4 <- matrix(NA,nrow=noLambda,ncol=n)
  CVFocusWeightedMat5 <- matrix(NA,nrow=noLambda,ncol=n)
  
  focusCentred <- scale(X[,-1], center=focusX[-1],scale=FALSE) 
  weigths1 <- 1/sqrt(rowSums((t(t(X)-focusX))^2))
  
  for(i in 1:length(lambdas)){
      H <- X %*% solve(crossprod(X,X)+lambdas[i]*diag(p+1))%*% t(X) 
      e <- (diag(n) - H) %*% Y
      CVRegularMat[i,] <- (e/(1-diag(H)))^2
      CVFocusWeightedMat1[i,] <- (e/(1-diag(H)))^2*weigths1 
      CVFocusWeightedMat2[i,] <- (e/(1-diag(H)))^2/(1- diag(focusCentred%*%solve(crossprod(focusCentred,focusCentred))%*%t(focusCentred)))^2
      CVFocusWeightedMat3[i,] <- (e/(1-diag(H)))^2*exp(2*diag(focusCentred%*%solve(crossprod(focusCentred,focusCentred))%*%t(focusCentred)))
      CVFocusWeightedMat4[i,] <- (e/(1-diag(H)))^2/(1- diag(focusCentred%*%solve(crossprod(focusCentred,focusCentred))%*%t(focusCentred)))^(log(lambdas[i]/40+0.01))
      CVFocusWeightedMat5[i,] <- (e/(1-diag(H)))^2/(1- diag(focusCentred%*%solve(crossprod(focusCentred,focusCentred))%*%t(focusCentred)))^(log(lambdas[i]/10+0.01))
  }
  retList <- list()
  retList$Reg <- rowMeans(CVRegularMat)
  retList$Focus1 <- rowMeans(CVFocusWeightedMat1)
  retList$Focus2 <- rowMeans(CVFocusWeightedMat2)
  retList$Focus3 <- rowMeans(CVFocusWeightedMat3)
  retList$Focus4 <- rowMeans(CVFocusWeightedMat4)
  retList$Focus5 <- rowMeans(CVFocusWeightedMat5)
  retList$lambda <- lambdas
  return(retList)
}

# The actual simulation ####
predMuRegVec <-  rep(NA,noSim)
noWeights <- 5
predMuFocusVec <- matrix(NA,nrow=noWeights,ncol=noSim)
focusXMat <- matrix(NA,nrow=p+1,ncol=noSim)

for (i in 1:noSim){
  X <- cbind(1,matrix(runif(p*n,min=-1,max=1),ncol=p))
  Y <- X%*%beta0 + rnorm(n,mean=0,sd=sigma0)

  # Defines a random focus paramter
  focusX <- c(1,runif(p,min=-2,max=-1))
  
  ## CV stuff
  thisCVRun <- LOOCVRidge(X,Y,focusX=focusX)
  lambdas <- thisCVRun$lambda  
  predMuRegVec[i] <-  focusX %*% solve(crossprod(X,X)+lambdas[which.min(thisCVRun$Reg)]*diag(p+1))%*% t(X) %*% Y
  
  for(j in 1:noWeights){
    predMuFocusVec[j,i] <-  focusX %*% solve(crossprod(X,X)+lambdas[which.min(eval(parse(text=paste('thisCVRun$Focus',j,sep=''))))]*diag(p+1))%*% t(X) %*% Y
  }
  focusXMat[,i]=focusX
  print(i) 
}

trueMuVec <- t(focusXMat)%*%beta0
orderTrueMuVec <- order(trueMuVec)

mean((predMuRegVec-trueMuVec)^2)
rowMeans((predMuFocusVec-t(trueMuVec[,rep(1,noWeights)]))^2)

### To get sound when simulation is completed
library(beepr)
beep()
#### end ####