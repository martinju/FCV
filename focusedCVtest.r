

# True model y =  beta0 + beta1x + eps, eps ~ N(0,sigma2)
betaTrue = c(1,0.5)
sigma2  = 0.9

# Defining deterministic weights
x = seq(-2,2,1)
n = length(x)

# Sampling y from the true model, and then keeping these fixed
set.seed(1234)
y = betaTrue[1] + betaTrue[2]*x + rnorm(n,0,sigma2)

# Fitting linear model

lm.mod = lm(y~x)
mod1 = lm.mod$coefficients
mod0 = c(mean(y),0)

y.hat.mod0 = mod0[1]+mod0[2]*x
y.hat.mod1 = mod1[1]+mod1[2]*x


# Plotting the data, and fitted models

plot(x,y)
abline(0,10^5,lwd=1,lty=1)
abline(0,0,lwd=1,lty=1)
abline(betaTrue[1],betaTrue[2],lwd=2,lty=2)
abline(mod0[1],mod0[2],col=2,lwd=2)
abline(mod1[1],mod1[2],col=4,lwd=2)

points(x,y.hat.mod0,col=2)
points(x,y.hat.mod1,col=4)


# Checking which model is best for different x0
x0val = seq(-2,2,0.01)
pred.error.mod0 = ((mod0[1]+ mod0[2]*x0val) - (betaTrue[1]+ betaTrue[2]*x0val))^2
pred.error.mod1 = ((mod1[1]+ mod1[2]*x0val) - (betaTrue[1]+ betaTrue[2]*x0val))^2

best.mod = apply(cbind(pred.error.mod0,pred.error.mod1),MARGIN=1,which.min)-1
first0best = min(which(best.mod==0))
last0best = max(which(best.mod==0))

lines(x0val,rep(0,length(x0val)),col=4,lwd=6)
lines(x0val[first0best:last0best],rep(0,length(first0best:last0best)),col=2,lwd=6)

# Computing CV contributions for each model
X = cbind(1,x)

HatMatrix.mod0 = X[,1]%*%solve(t(X[,1])%*%X[,1])%*%t(X[,1])
HatMatrix.mod1 = X%*%solve(t(X)%*%X)%*%t(X)

h.vec.mod0 = diag(HatMatrix.mod0)
h.vec.mod1 = diag(HatMatrix.mod1)

res.mod0 = y-y.hat.mod0
res.mod1 = y-y.hat.mod1

CV.contrib.mod0 = (res.mod0/(1-h.vec.mod0))^2
CV.contrib.mod1 = (res.mod1/(1-h.vec.mod1))^2

CV.contrib.mod0/CV.contrib.mod1


sum(CV.contrib.mod0)
sum(CV.contrib.mod1)



