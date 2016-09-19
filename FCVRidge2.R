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

beta0 <-rep(2,2)   #c(1,-0.1,0.4,0.1,-0.2)
p <- length(beta0)
sigma0 = 1

# Defining the simulation experiment #### 

noSim <- 1000
n <- 50
#Defining the ridge Leave-one-out crossvalidation 
CV <- function(lambda){
  H <- X %*% solve(crossprod(X,X)+lambda*diag(p))%*% t(X) 
  e <- (diag(n) - H) %*% Y
  mean((e/(1-diag(H)))^2)
}

CVfocus <- function(lambda){
  H <- X %*% solve(crossprod(X,X)+lambda*diag(p))%*% t(X) 
  e <- (diag(n) - H) %*% Y
  mean((e/(1-diag(H)))^2*weights/mean(weights))
}

Weight.function <- function(X){
  X.c <- X-focusX
  exp(-0.02*t(X.c)%*%solve(matrix(c(beta0[1],0.95*sqrt(beta0[1]*beta0[2]),0.95*sqrt(beta0[1]*beta0[2]),beta0[2]),2,2))%*%X.c)
}

# The acutal simulation ####
predMuRegVec <- rep(NA,noSim)
predMuFocusVec <- rep(NA,noSim)

lambdaRegVec <- rep(NA,noSim)
lambdaFocusVec <- rep(NA,noSim)

focusXMat = matrix(NA,ncol=p,nrow=noSim)

tuningPar <- matrix(NA,ncol=2,nrow=noSim)

for (i in 1:noSim){
  X <- scale(matrix(runif(p*n,min=-1,max=1),ncol=p),center = TRUE, scale = FALSE)
  Y <- X%*%beta0 + scale(rnorm(n,mean=0,sd=sigma0),center = TRUE, scale = FALSE)
  
  # Defines a random focus paramter
  focusX <- runif(p,min=-1,max=1)
  
  ## CV stuff
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #focusCentred <- scale(X, center=focusX,scale=FALSE)  
  #weights <- exp(-0.05/(abs(focusX%*%beta0)^0.5)*rowMeans(focusCentred^2)) #

  weights <- apply(X,1, Weight.function)
  lambdaRegVec[i] <- optim(par = 1,CV,lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e6))$par
  print(lambdaRegVec[i])
  if((  lambdaRegVec[i] < 1.01)&&(  lambdaRegVec[i]  > 0.999)){print('Error')}
  
  lambdaFocusVec[i] <-optim(par = 1,CVfocus,lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e6))$par
  print(lambdaFocusVec[i])
  if((  lambdaFocusVec[i] < 1.01)&&(  lambdaFocusVec[i] > 0.999)){print('Error focus')}
  
  predMuRegVec[i] <-  focusX %*% solve(crossprod(X,X)+lambdaRegVec[i]*diag(p))%*% t(X) %*% Y
  predMuFocusVec[i] <-  focusX %*% solve(crossprod(X,X)+lambdaFocusVec[i]*diag(p))%*% t(X) %*% Y
  
  focusXMat[i,]=focusX
}

#Dividing the focus into groups 
trueMuVec <- focusXMat%*%beta0
orderTrueMuVec <- order(trueMuVec)

##Check difference in mean prediction error
mean((predMuRegVec-trueMuVec)^2)
mean((predMuFocusVec-trueMuVec)^2)

##Check mean and median of the ratio of the prediction erros, sould be less than 1 (but highly skewed)
mean((predMuFocusVec-trueMuVec)^2<(predMuRegVec-trueMuVec)^2)

#### end ####