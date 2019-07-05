##############################
#### hardy-weinberg model#####
n.o <- 176; n.a<-182; n.b<-60; n.ab<-17; n<-sum(n.o+n.a+n.b+n.ab)
install.packages("rjags")
library(rjags)
install.packages("runjags")
library(runjags)
install.packages("optimx")
library(optimx)

#### generate theta ####
theta1<-runif(1,0,1)
theta2<-runif(1,0,theta1)
theta0<-1-theta1-theta2

hardy<-function(theta) {                   
  theta1<-theta[1]
  theta2<-theta[2]
  theta0<-1-theta1-theta2
  p0<-theta0^2
  p1<-theta1^2+2*theta1*theta0
  p2<-theta2^2+2*theta2*theta0
  p3<-2*theta1*theta2
  
  loglik.hardy<-log((p0^n.o)*(p1^n.a)*(p2^n.b)*(p3^n.ab))
  -loglik.hardy                                        # -(log-likelihood function)
}


theta1<-runif(1,0,1)
theta2<-runif(1,0,theta1)
mle.hardy <- optimx( par = c(theta1,theta2), hardy, hessian=TRUE)
mle.hardy            #theta1=0.2644747, theta2=0.09316877, AIC=492.5353, theta0=0.642344

theta1.mle<-mle.hardy$p1[1]
theta2.mle<-mle.hardy$p2[1]
theta0.mle<-1-theta1.mle-theta2.mle

#### Metropolis Algorithm ####
theta.initial<-NULL
theta.initial<-c(theta1.mle,theta2.mle)
z11<-NULL; z12<-NULL

theta.hastings <- function(M){
  # implementation of the metropolis algorithm
  # M: number of time steps
  theta <- matrix(0,M,2)
  thetastar <- matrix(0,M-1,2)
  theta[1,] <- theta.initial
  for (t in 2:M) {
    z11<-rnorm(1,0,1); z12<-rnorm(1,0,1)
    thetastar[t-1,] <- theta[t-1,]+c(z11,z12)/sqrt(n)
    lik.star<-(((1-thetastar[t-1,1]-thetastar[t-1,2])^2)^n.o)*
      (((thetastar[t-1,1])^2+2*thetastar[t-1,1]*(1-thetastar[t-1,1]-thetastar[t-1,2]))^n.a)*
      (((thetastar[t-1,2])^2+2*thetastar[t-1,2]*(1-thetastar[t-1,2]-thetastar[t-1,2]))^n.b)*
      (2*thetastar[t-1,1]*thetastar[t-1,2])
    lik.old<-(((1-theta[t-1,1]-theta[t-1,2])^2)^n.o)*
      (((theta[t-1,1])^2+2*theta[t-1,1]*(1-theta[t-1,1]-theta[t-1,2]))^n.a)*
      (((theta[t-1,2])^2+2*theta[t-1,2]*(1-theta[t-1,2]-theta[t-1,2]))^n.b)*
      (2*theta[t-1,1]*theta[t-1,2])
    a<-(lik.star*thetastar[t-1,])/(lik.old*theta[t-1,])
    alpha <- min(1,a)
    if (runif(1,0,1) < alpha) {
      theta[t,] <- thetastar[t-1,]
    }
    else {
      theta[t,] <- theta[t - 1,]
    }
  }
  return(theta)
}
final.theta<-rbind(theta[1,],thetastar)

#### M-H algorithm 95% C.I
CI.theta1<-quantile(final.theta[,1], probs=c(0.25,0.75))   
CI.theta2<-quantile(final.theta[,2], probs=c(0.25,0.75))   

x=theta.initial; y=loglik[201]                
a=(y-max.lik)/((x-beta.mle)^2)
J=-2*a
CI1=c(beta.mle-1.96/sqrt(-2*a), beta.mle+1.96/sqrt(-2*a))

out<-nlm(hardy,theta.initial, hessian=TRUE)
jj<-out$hessian
j<-solve(jj)

conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)
theta.initial[1] + c(-1, 1) * crit * sqrt(j[1, 1])

theta.initial[2] + c(-1, 1) * crit * sqrt(j[2, 2])

##### calling M-H algorithm and plotting result
test <- theta.hastings(1000)
plot(test, main="h(a,b | Data)")




#################################
###### JAGS ##############
accidents<-read.csv("Mining-Accidents-CPA.csv", header=T, na.strings="*")
names(accidents)<-c("t","y","h")
head(accidents)
attach(accidents)

# Modeling
modelString="
model
{
  for(i in 1:K){
  y[i]~dpois(mu[i])
  }
  
  for(i in 1:K){
  mu[i]=exp(a-b*h[i])
  h[i]<-ifelse((t[i]-c)>=0, 1, 0)
  }
  
  c~dunif(1851,1962)
  a~dnorm(1,1)
  b~dnorm(1,1)
  
}"
  
writeLines(modelString, "model_accidents.txt")

#Prior
t<-accidents$t; y<-accidents$y; h<-accidents$h
K<-length(y)

#데이터 리스트
dataList=list(y=y, t=t, K=K)

#초기치 설정
initsList=list(a=1,b=1,c=1888)

#JAGS 실행준비
jagsModel=jags.model(file="model_accidents.txt", data=dataList, inits=initsList,
                     n.chains=3, n.adapt=500)
#예비단계 수행
update(jagsModel, n.iter=3000)

#MCMC 표본생성
codaSamples=coda.samples(jagsModel, variable.names=c("a","b","c"),n.iter=5000)

#수렴진단
par(mfrow=c(1,2))
coda::traceplot(codaSamples[,"a"], main="", ylab="a")
acf(codaSamples[,"a"][[1]], plot=T, main="")
coda::traceplot(codaSamples[,"b"], main="", ylab="b")
acf(codaSamples[,"b"][[1]], plot=T, main="")
coda::traceplot(codaSamples[,"c"], main="", ylab="c")
acf(codaSamples[,"c"][[1]], plot=T, main="")

#사후추론
a.Samples=as.matrix(codaSamples[,"a"])
b.Samples=as.matrix(codaSamples[,"b"])
c.Samples=as.matrix(codaSamples[,"c"])
par(mfrow=c(1,1))
plot(density(a.Samples), xlab=expression(a), ylab="posterior density",main="")
plot(density(b.Samples), xlab=expression(b), ylab="posterior density",main="")
plot(density(c.Samples), xlab=expression(c), ylab="posterior density",main="")

hist(a.Samples, xlab="a")
hist(b.Samples, xlab="b")
hist(c.Samples, xlab="c")

### 추정값 및 Bayesian 95% C.I ####
##HPD samles###
HPDsample=function(theta, level=0.95){
  N=length(theta)
  theta.sort=sort(theta)
  M=ceiling(N*level)
  nCI=N-M
  CI.width=rep(0,nCI)
  for( i in 1:nCI) CI.width[i]=theta.sort[i+M]-theta.sort[i]
  index=which.min(CI.width)
  HPD=c(theta.sort[index], theta.sort[index+M])
  return(HPD)
}

MCMCSamples=as.matrix(codaSamples)

para=c("a","b","c")
par(mfrow=c(1,3))
for(i in 1:3){
  plot(density(MCMCSamples[,i]), main="", xlab=para[i])
  HPD=HPDsample(MCMCSamples[,i])
  abline(v=HPD, lty=2)
}

HPD.a=HPDsample(MCMCSamples[,1])
HPD.b=HPDsample(MCMCSamples[,2])
HPD.c=HPDsample(MCMCSamples[,3])

### MLE 추정값 ####
loglik.abc<-function(para){
  a<-para[1]; b<-para[2]; c<-para[3]
  h<-ifelse((t-c)>=0,1,0)
  mu<-exp(a-b*h)
  loglik <- sum(log(dpois(y, mu)))
  return(-loglik)
}

abc.mle<-optim(par=c(mean(a.Samples),mean(b.Samples),mean(c.Samples)),loglik.abc)$par

### Boostrap 신뢰구간
B=1000; p=3; n=length(y)
boots.abc<-matrix(0,B,p)

for(iter in 1:B){                               
  mynum=sample(1:n, n, replace=TRUE)         
  newy=y[mynum]
  #newmu=mu[mynum]                  
  newt=t[mynum]   
  
  loglik.abc<-function(para){
    a<-para[1]; b<-para[2]; c<-para[3]
    h<-ifelse((newt-c)>=0,1,0)
    mu<-exp(a-b*h)
    loglik <- sum(log(dpois(newy, mu)))
    return(-loglik)
  }
  boots.abc[iter,]<-optim(par=c(mean(a.Samples),mean(b.Samples),mean(c.Samples)), 
                          fn=loglik.abc)$par
}

boots.abc.hat<-colMeans(boots.abc)
boots.abc.hat

CI.lower<-NULL; CI.upper<-NULL
for(i in 1:3){
  CI.lower[i]<-sort(boots.abc[,i])[B*0.025+1]
  CI.upper[i]<-sort(boots.abc[,i])[B*0.975]
}
boots.CI<-cbind(CI.lower, boots.abc.hat, CI.upper)