
# Test script for doing focused type of weighted cross validation with regression:

# Idea: 
# Assign weights to each observation corresponding to the inverse of the distance to covariate value of interest (or other distance measure).
# Then perform 5 fold cross validation where the prediction error for each observation being predicted is weighted with these weights
# The model with the smallest test prediction error is chosen

# Experiment: 
# Sample data from linear regression model with several covariates of variable importance
# Define a set of models containing this true model
# Perform 5 fold weighted cross validation and find the best model according to that criterion (AND with equal weighting)
# Predict the mu based on the two models and also on the winning AIC and BIC model for comparison
# Repeat the procedure several times to see what procedure performs the best.


# Defining a true model ####

beta0 <- c(1,-0.5,0.5,0.1,0)
p <- length(beta0)-1
sigma0 = 0.2

# Defining the simulation experiment #### 

noSim <- 10^3
n <- 100
kFold <- 5

# Defining the candidate models ####

helpList <- rep(list(0:1),p)
modMat <- cbind(1,expand.grid(helpList))
noMod <- dim(modMat)[1]

# Defines the focus ####

focusX <- c(1,0.7,-0.5,-0.3,-0.2)

# Defining the cross validation function ####

CV <- function(data,kFold,focusX){
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
  
  predErrorRegularMat <- matrix(NA,nrow=noMod,ncol=kFold)
  predErrorFocusWeightedMat <- matrix(NA,nrow=noMod,ncol=kFold)
  
  for (k in 1:kFold){
    trainData <- data[which(!(whichFold==k)),]
    testData <- data[which(whichFold==k),]
    testWeighttNorm <- weight[which(whichFold==k)]/sum(weight[which(whichFold==k)])
    
    trainY <- trainData[,1]
    trainX <- trainData[,-1]
    testY <- testData[,1]
    testX <- testData[,-1]
    
    for (l in 1:noMod){
      theseCovar <- c(which(modMat[l,]==1))
      
      modTrainX <- trainX[,theseCovar]
      modTestX <- testX[,theseCovar]
      
      modLm <- lm(trainY~as.matrix(modTrainX)-1)
      predY <- as.matrix(modTestX)%*%coef(modLm)
      
      predErrorRegularMat[l,k] <- mean((predY-testY)^2)
      predErrorFocusWeightedMat[l,k] <- sum(testWeighttNorm*(predY-testY)^2)
      
    }
  }
  retList <- list()
  retList$Reg <- rowMeans(predErrorRegularMat)
  retList$Focus <- rowMeans(predErrorFocusWeightedMat)
  
  return(retList)
}


# The acutal simulation ####

predMuRegVec <- rep(NA,noSim)
predMuFocusVec <- rep(NA,noSim)
predMuAICVec <- rep(NA,noSim)
predMuBICVec <- rep(NA,noSim)

for (i in 1:noSim){
  X <- cbind(1,matrix(runif(p*n,min=-1,max=1),ncol=p))
  Y <- X%*%beta0 + rnorm(n,mean=0,sd=sigma0)

  data = data.frame(Y=Y,X=X)

  ## CV stuff
  thisCVRun <- CV(data,kFold=kFold,focusX=focusX)
  bestReg <- which.min(thisCVRun$Reg)
  bestFocus <- which.min(thisCVRun$Focus)

  ## AIC and BIC stuff
  AICMod <- rep(NA,noMod)
  BICMod <- rep(NA,noMod)
  for (j in 1:noMod){
    theseCovar <- c(which(modMat[j,]==1))
    
    thisX <- X[,theseCovar]
    
    modLm <- lm(Y~as.matrix(thisX)-1)
    
    AICMod[j] <- AIC(modLm)
    BICMod[j] <- BIC(modLm)
  }
  bestAIC <- which.min(AICMod)
  bestBIC <- which.min(BICMod)
  
  theseCovarBestReg <- c(which(modMat[bestReg,]==1))
  theseCovarBestFocus <- c(which(modMat[bestFocus,]==1))
  theseCovarBestAIC <- c(which(modMat[bestAIC,]==1))
  theseCovarBestBIC <- c(which(modMat[bestBIC,]==1))
  
  XBestReg<- X[,theseCovarBestReg]
  XBestFocus<- X[,theseCovarBestFocus]
  XBestAIC<- X[,theseCovarBestAIC]
  XBestBIC<- X[,theseCovarBestBIC]
  
  modLmBestReg <- lm(Y~as.matrix(XBestReg)-1)
  modLmBestFocus <- lm(Y~as.matrix(XBestFocus)-1)
  modLmBestAIC <- lm(Y~as.matrix(XBestAIC)-1)
  modLmBestBIC <- lm(Y~as.matrix(XBestBIC)-1)
  
  predMuRegVec[i] <- sum(focusX[theseCovarBestReg]*coef(modLmBestReg))
  predMuFocusVec[i] <- sum(focusX[theseCovarBestFocus]*coef(modLmBestFocus))
  predMuAICVec[i] <- sum(focusX[theseCovarBestAIC]*coef(modLmBestAIC))
  predMuBICVec[i] <- sum(focusX[theseCovarBestBIC]*coef(modLmBestBIC))
 
  print(i) 
}

trueMu <- sum(beta0*focusX)

sum((predMuRegVec-trueMu)^2)
sum((predMuFocusVec-trueMu)^2)
sum((predMuAICVec-trueMu)^2)
sum((predMuBICVec-trueMu)^2)

# end ####



