

# True model y =  beta0 for x<0, and beta0 + beta1x + eps, eps ~ N(0,sigma2)
betaTrue = c(0,2)
sigma2  = 0.2

y.func <- function(x,betaTrue){
  x0 <- x[x<0]
  x1 <- x[x>=0]
  y0 <- rep(betaTrue[1],length(x0))
  y1 <- betaTrue[1] + betaTrue[2]*x1
  y <- rep(NA,length(x))
  y[x<0] <- y0
  y[x>=0] <- y1
  return(y)
}

### Finding least false versions of this model, by sampling

xx = seq(-2,2,0.001)
yy = y.func(xx,betaTrue)
lm.mod = lm(yy~xx)
mod1.lf = c(1,1)
mod0.lf = c(1,0)

# Defining deterministic weights

# Sampling y from the true model, and then keeping these fixed
# Fitting linear model

# Plotting the data, and fitted models

plot(1,1,type="n",xlim=c(-2,2),ylim=c(-2,3))
abline(0,10^5,lwd=1,lty=1)
abline(0,0,lwd=1,lty=1)

abline(mod0.lf[1],mod0.lf[2],col=2,lwd=2,lty=3)
abline(mod1.lf[1],mod1.lf[2],col=4,lwd=2,lty=3)
lines(seq(-2,2,0.01),y.func(seq(-2,2,0.01),betaTrue),lwd=3,lty=2)


# Checking which LEAST FALSE model is best for different x0
x0val = seq(-2,2,0.01)
pred.error.mod0.lf = ((mod0.lf[1]+ mod0.lf[2]*x0val) - (y.func(x0val,betaTrue)))^2
pred.error.mod1.lf = ((mod1.lf[1]+ mod1.lf[2]*x0val) - (y.func(x0val,betaTrue)))^2

best.mod.lf = apply(cbind(pred.error.mod0.lf,pred.error.mod1.lf),MARGIN=1,which.min)-1


for (i in 1:length(best.mod.lf)){
  lines(c(x0val[i],x0val[i+1]),c(2,2),col=ifelse(best.mod.lf[i]==0,"red","blue"),lwd=6)
}

set.seed(12)

n=10000
x = runif(n,-2,2)
y = y.func(x,betaTrue)+rnorm(n,0,sigma2)
lm.mod = lm(y~x)
mod1 = lm.mod$coefficients
mod0 = c(mean(y),0)
y.hat.mod0 = mod0[1]+mod0[2]*x
y.hat.mod1 = mod1[1]+mod1[2]*x

abline(mod0[1],mod0[2],col=2,lwd=1)
abline(mod1[1],mod1[2],col=4,lwd=1)
pred.error.mod0.x = ((mod0[1]+ mod0[2]*x) - (y))^2
pred.error.mod1.x = ((mod1[1]+ mod1[2]*x) - (y))^2
best.mod = apply(cbind(pred.error.mod0.x,pred.error.mod1.x),MARGIN=1,which.min)-1

points(x[best.mod==0],y[best.mod==0],col="red")
points(x[best.mod==1],y[best.mod==1],col="blue")

binned.x <- cut(x,breaks=seq(-2,2,length.out=n/50+1),include.lowest = T)
levels.x = levels(binned.x)

mean.predmod0 = rep(NA,length(length(levels.x)))
mean.predmod1 = rep(NA,length(length(levels.x)))

for (i in 1:length(levels.x)){
  these <- binned.x==levels.x[i]
  mean.predmod0[i] = mean(((mod0[1]+ mod0[2]*x[these]) - (y[these]))^2)
  mean.predmod1[i] = mean(((mod1[1]+ mod1[2]*x[these]) - (y[these]))^2)
}

lines(seq(-2,2,length.out=n/50),mean.predmod0,col=1,lwd=3)
lines(seq(-2,2,length.out=n/50),mean.predmod1,col=3,lwd=3)


# Checking which model is best for different x0
x0val = seq(-2,2,0.01)
pred.error.mod0 = ((mod0[1]+ mod0[2]*x0val) - (y.func(x0val,betaTrue)))^2
pred.error.mod1 = ((mod1[1]+ mod1[2]*x0val) - (y.func(x0val,betaTrue)))^2

best.mod = apply(cbind(pred.error.mod0,pred.error.mod1),MARGIN=1,which.min)-1

for (i in 1:length(best.mod)){
  lines(c(x0val[i],x0val[i+1]),c(-2,-2),col=ifelse(best.mod[i]==0,"red","blue"),lwd=6)
}




(pred.error.mod0.avgwhere0isbest <- mean(pred.error.mod0[best.mod.lf==0]))
(pred.error.mod1.avgwhere0isbest <- mean(pred.error.mod1[best.mod.lf==0]))

(pred.error.mod0.avgwhere1isbest <- mean(pred.error.mod0[best.mod.lf==1]))
(pred.error.mod1.avgwhere1isbest <- mean(pred.error.mod1[best.mod.lf==1]))



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



mean(CV.contrib.mod0)
mean(CV.contrib.mod1)


plot(log(CV.contrib.mod0[order(x)]))
points(log(CV.contrib.mod1[order(x)]),col=2)

focusVec <- seq(-2,2,0.1)

cv.vec.mod0 = rep(NA,length(focusVec))
cv.vec.mod1 = rep(NA,length(focusVec))


for (i in 1:length(focusVec)){
  focus = focusVec[i]
  weight=1/(focus-x)^2
  weight=weight/sum(weight)

  cv.vec.mod0[i] = mean(CV.contrib.mod0*weight)
  cv.vec.mod1[i] = mean(CV.contrib.mod1*weight)
  print(i)
}

plot(focusVec,cv.vec.mod0,type='l',col="red")
lines(focusVec,cv.vec.mod1,col="blue")




focus=0.4



