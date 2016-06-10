# Test script for doing focused type of weighted cross validation, for first ridge, later elastic net and others. 

# Script showing the effect of different foci

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

noSim <- 3000
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
  CVFocusWeightedMat <- matrix(NA,nrow=noLambda,ncol=n)
  
  focusCentred <- scale(X[,-1], center=focusX[-1],scale=FALSE) 
  for(i in 1:length(lambdas)){
      H <- X %*% solve(crossprod(X,X)+lambdas[i]*diag(p+1))%*% t(X) 
      e <- (diag(n) - H) %*% Y
      CVRegularMat[i,] <- (e/(1-diag(H)))^2
      CVFocusWeightedMat[i,] <- (e/(1-diag(H)))^2/(1- diag(focusCentred%*%solve(crossprod(focusCentred,focusCentred))%*%t(focusCentred)))^(log(lambdas[i]/40+0.01))
  }
  retList <- list()
  retList$Reg <- rowMeans(CVRegularMat)
  retList$Focus <- rowMeans(CVFocusWeightedMat)
  retList$lambda <- lambdas
  return(retList)
}

# The acutal simulation ####
predMuRegVec <- rep(NA,noSim)
predMuFocusVec <- rep(NA,noSim)

lambdaRegVec <- rep(NA,noSim)
lambdaFocusVec <- rep(NA,noSim)

focusXMat = matrix(NA,ncol=p+1,nrow=noSim)

tuningPar <- matrix(NA,ncol=2,nrow=noSim)

for (i in 1:noSim){
  X <- cbind(1,matrix(runif(p*n,min=-1,max=1),ncol=p))
  Y <- X%*%beta0 + rnorm(n,mean=0,sd=sigma0)
  
  # Defines a random focus paramter
  focusX <- c(1,runif(p,min=-1,max=0))
  
  ## CV stuff
  thisCVRun <- LOOCVRidge(X,Y,focusX=focusX)
  lambdas <- thisCVRun$lambda
  lambdaRegVec[i] <- lambdas[which.min(thisCVRun$Reg)]
  lambdaFocusVec[i] <- lambdas[which.min(thisCVRun$Focus)]
  
  predMuRegVec[i] <-  focusX %*% solve(crossprod(X,X)+lambdaRegVec[i]*diag(p+1))%*% t(X) %*% Y
  predMuFocusVec[i] <-  focusX %*% solve(crossprod(X,X)+lambdaFocusVec[i]*diag(p+1))%*% t(X) %*% Y
  
  focusXMat[i,]=focusX
  print(i) 
}

#Dividing the focus into groups 
trueMuVec <- focusXMat%*%beta0
orderTrueMuVec <- order(trueMuVec)

mean((predMuRegVec-trueMuVec)^2)
mean((predMuFocusVec-trueMuVec)^2)

index <- apply(apply(focusXMat[,-1],1,function(a){((a<=0))}),2,all)
mean((predMuRegVec[index]-trueMuVec[index])^2)
mean((predMuFocusVec[index]-trueMuVec[index])^2)

index <- apply(apply(focusXMat[,-1],1,function(a){((a<=-0.1))}),2,all)
mean((predMuRegVec[index]-trueMuVec[index])^2)
mean((predMuFocusVec[index]-trueMuVec[index])^2)

index <- apply(apply(focusXMat[,-1],1,function(a){((a<=-0.2))}),2,all)
mean((predMuRegVec[index]-trueMuVec[index])^2)
mean((predMuFocusVec[index]-trueMuVec[index])^2)

index <- apply(apply(focusXMat[,-1],1,function(a){((a<=-0.3))}),2,all)
sum(index)
mean((predMuRegVec[index]-trueMuVec[index])^2)
mean((predMuFocusVec[index]-trueMuVec[index])^2)

### To get sound when simulation is completed
library(beepr)
beep()
#### end ####