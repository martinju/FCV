## Plot showing the distribution of lambda_x0 in 2 dimensions
## Which shape does sets of identical lambda_x0 have?

rm(list = ls())
# Defining a true model: no intercept, require centering of outcome 
p <- 2
sigma2 <- 2
beta <-rep(1,2)   
n <- 50

#Drawing data matrix and focus from uniform distribution. Oracle, no distributed outcome?
set.seed(111)
X <- scale(matrix(runif(p*n,min=-1,max=1),ncol=p),center = TRUE, scale = FALSE)
set.seed(111)
Y <- X%*%beta + scale(rnorm(n,mean=0,sd=sqrt(sigma2)),center = TRUE, scale = FALSE)

PredictionErrorFocus <- function(x1,x2){
  x0 <- c(x1,x2)
  optim(par = 10,function(lambda){(x0%*%solve(crossprod(X,X)+lambda*diag(p))%*% t(X)%*%Y- x0%*%beta)^2},lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e6))$par
}

betaOLS <- solve(crossprod(X,X))%*% t(X)%*%Y
x1 <- seq(-1,1,length.out = 200)
x2 <- seq(-1,1,length.out = 200)
VecFun <- Vectorize(PredictionErrorFocus)

z <- outer(x1,x2,VecFun)+10^-10

range(z)
hist(z)
hist(log(z))

image(x1,x2,log(z),col=heat.colors(99),breaks=quantile(log(z),probs = exp(seq(log(0.01), log(0.9), length.out = 100))),main = expression(paste("Contours for log ", lambda[x[0]])))
contour(x1,x2,log(z),levels=quantile(log(z),probs = seq(0.2, 0.9, 0.1)),add=T,method='edge')
contour(x1,x2,log(z),levels=c(min(log(z)+2)),add=T,method='edge')
abline(a=0,b=-(betaOLS[1])/betaOLS[2],col='blue',lwd=3)
abline(a=0,b=-beta[1]/beta[2],col='red',lwd=1,lty=3)
text(X,label=round(log(apply(X,1,function(x){PredictionErrorFocus(x[1],x[2])})+10^-14),1))
legend('bottomleft',legend = c('OLS beta','True beta'),lty = c(1,3),lwd = c(3,1),col=c('blue','red') )

#Test curves 
betaOLS
x0 <- c(0.5,-betaOLS[1]/betaOLS[2]*0.5+0.1)
x0
curve(sapply(x,function(lambda){(x0%*%solve(crossprod(X,X)+lambda*diag(p))%*% t(X)%*%Y- x0%*%beta)^2}),0,20)
optim(par = 10,function(lambda){(x0%*%solve(crossprod(X,X)+lambda*diag(p))%*% t(X)%*%Y- x0%*%beta)^2},lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e6))$par

x0 <- c(0.5,-betaOLS[1]/betaOLS[2]*0.5-0.005)
x0
curve(sapply(x,function(lambda){(x0%*%solve(crossprod(X,X)+lambda*diag(p))%*% t(X)%*%Y- x0%*%beta)^2}),0,100)
optim(par = 100,function(lambda){(x0%*%solve(crossprod(X,X)+lambda*diag(p))%*% t(X)%*%Y- x0%*%beta)^2},lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e6))$par

X <- scale(matrix(runif(p*n,min=-1,max=1),ncol=p),center = TRUE, scale = FALSE)
Y <- X%*%beta + scale(rnorm(n,mean=0,sd=sqrt(sigma2)),center = TRUE, scale = FALSE)

PredictionErrorLambdaNull<- function(x1,x2){
  x0 <- c(x1,x2)
  (x0%*%solve(crossprod(X,X)) %*% t(X)%*%Y- x0%*%beta)^2
}

betaOLS <- solve(crossprod(X,X))%*% t(X)%*%Y
x1 <- seq(-2,2,length.out = 600)
x2 <- seq(-2,2,length.out = 600)

NullError <- outer(x1,x2,Vectorize(PredictionErrorLambdaNull))
image(x1,x2,NullError,col=heat.colors(99),breaks=quantile(NullError,probs = exp(seq(log(0.001), log(0.9), length.out = 100))),main = expression(paste("Contours for log ", lambda[x[0]])))
abline(a=0,b=-(betaOLS[1])/betaOLS[2],col='blue',lwd=3)
abline(a=0,b=-beta[1]/beta[2],col='red',lwd=1,lty=3)
abline(a=0,b=-(betaOLS[1]-beta[1])/(betaOLS[2]-beta[2]),col='blue',lwd=3)

