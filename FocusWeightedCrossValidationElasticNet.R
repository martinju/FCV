
# Test script for doing focused type of weighted cross validation with Elastic net:

# Idea: 
# Assign weights to each observation corresponding to the inverse of the distance to covariate value of interest (or other distance measure).
# Then perform 5 fold cross validation where the prediction error for each observation being predicted is weighted with these weights
# The model with the smallest test prediction error is chosen

# Experiment: 
# Sample data from linear regression model with several covariates of variable importance
# Define a set of models containing this true model
# Perform 5 fold weighted cross validation and find the best model according to that criterion (AND with equal weighting)
# Predict the mu based on the two models 
# Repeat the procedure several times to see what procedure performs the best.

library(glmnet)

# Defining a true model ####

beta0 <- c(1,-0.1,0.4,0.1,-0.2)
p <- length(beta0)-1
sigma0 = 1

# Defining the simulation experiment #### 

noSim <- 10^3
n <- 100
kFold <- 5
alpha=0 # alpha is the alpha parameter in the elastic net (alpha=1 gives Lasso, alpha=0 gives ridge)

# Defining the cross validation function ####

CVElasticNet <- function(data,kFold,focusX,alpha){
  n <- dim(data)[1]
  obsPerFold <- n/kFold
  ind <- sample(n)-1
  whichFold <- numeric(0)
  for (l in 1:n){
    whichFold[l] <- ind[l]%/%obsPerFold+1
  }
  
  X <- data[,-1]
  covarDist <- sqrt(rowSums((t(t(X)-focusX))^2))
  weight <- 1/covarDist

  fit0 <- glmnet(as.matrix(X0),Y0,family = "gaussian",alpha=alpha,nlambda=1000)
  useLambda <- fit0$lambda
  noLambda <- length(useLambda)
  
  predErrorRegularMat <- matrix(NA,nrow=noLambda,ncol=kFold)
  predErrorFocusWeightedMat <- matrix(NA,nrow=noLambda,ncol=kFold)
  
  for (k in 1:kFold){
    trainData <- data[which(!(whichFold==k)),]
    testData <- data[which(whichFold==k),]
    testWeighttNorm <- weight[which(whichFold==k)]/sum(weight[which(whichFold==k)])
    
    trainY <- trainData[,1]
    trainX <- trainData[,-1]
    testY <- testData[,1]
    testX <- testData[,-1]
    
    fit <- glmnet(as.matrix(trainX),trainY,family = "gaussian",alpha=alpha,lambda=useLambda)
    predY <- predict(fit,as.matrix(testX),type="response")
    
    predErrorRegularMat[,k] <- colMeans((predY-testY)^2)
    predErrorFocusWeightedMat[,k] <- colSums(diag(testWeighttNorm)%*%((predY-testY)^2))
  }
  
  retList <- list()
  retList$Reg <- rowMeans(predErrorRegularMat)
  retList$Focus <- rowMeans(predErrorFocusWeightedMat)
  retList$lambda <- useLambda
  
  return(retList)
}


# The acutal simulation ####

predMuRegVec <- rep(NA,noSim)
predMuFocusVec <- rep(NA,noSim)

bestRegVec <- rep(NA,noSim)
bestFocusVec <- rep(NA,noSim)

focusXMat = matrix(NA,ncol=5,nrow=noSim)

for (i in 1:noSim){
  X <- cbind(1,matrix(runif(p*n,min=-1,max=1),ncol=p))
  Y <- X%*%beta0 + rnorm(n,mean=0,sd=sigma0)

  data = data.frame(Y=Y,X=X)

  # Defines a random focus paramter
  focusX <- c(1,runif(4,min=-1,max=1))
  
  ## CV stuff
  thisCVRun <- CVElasticNet(data,kFold=kFold,focusX=focusX,alpha=alpha)
  lambdas <- thisCVRun$lambda
  bestRegVec[i] <- lambdas[which.min(thisCVRun$Reg)]
  bestFocusVec[i] <- lambdas[which.min(thisCVRun$Focus)]

  modElasticNetBestReg <- glmnet(as.matrix(X),Y,family="gaussian",alpha=alpha,lambda=bestRegVec[i])
  modElasticNetBestFocus <- glmnet(as.matrix(X),Y,family="gaussian",alpha=alpha,lambda=bestFocusVec[i])

  predMuRegVec[i] <- predict(modElasticNetBestReg,t(focusX),type="response")
  predMuFocusVec[i] <- predict(modElasticNetBestFocus,t(focusX),type="response")
  
  focusXMat[i,]=focusX
  print(i) 
}

trueMuVec <- focusXMat%*%beta0

orderTrueMuVec <- order(trueMuVec)

mean((predMuRegVec-trueMu)^2)
mean((predMuFocusVec-trueMu)^2)

plot(1:noSim,trueMuVec[orderTrueMuVec],type='l')
lines(1:noSim,predMuRegVec[orderTrueMuVec],col=2)
lines(1:noSim,predMuFocusVec[orderTrueMuVec],col=3)

plot(trueMuVec[orderTrueMuVec],predMuRegVec[orderTrueMuVec],col=2,type='l')
lines(trueMuVec[orderTrueMuVec],predMuFocusVec[orderTrueMuVec],col=3)



# mean(bestRegVec)
# mean(bestFocusVec)

# summary(predMuRegVec)
# summary(predMuFocusVec)

# end ####



