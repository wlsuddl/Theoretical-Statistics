#prior
c<-0.5 
ybar<-c(0.78663,0.76219,0.77715,1.04744,0.8412,0.46245,0.86253,0.69475) 
sigsq0<-c(0.0629,0.0058,0.0071,0.0356,0.0186,0.036,0.1409,0.0962)
K=8
k0=1
mu0<-mean(ybar)
tausq0<-var(ybar)
d<-c*tausq0

######### Gibbs Sampling ##########
######start gibbs#####
Nwarm=1000; Nsim=5000
THETA=matrix(nrow=Nsim, ncol=K)
MTS=matrix(nrow=Nsim, ncol=2)
theta=c(1:K)

#initial values
mu=mu0
tausq=tausq0
sigsq=sigsq0
theta=ybar
for(ns in 1:(Nwarm+Nsim)){
  #generate theta
  for(i in 1:K){
    mi=(ybar[i]/(sigsq[i]^2)+mu/tausq)/(1/tausq+1/(sigsq[i]^2))
    vari=1/(1/tausq+1/(sigsq[i]^2))
    theta[i]=rnorm(1,mi,sqrt(vari))
  }
  #generate mu
  mu=rnorm(1,(K*mean(theta)+mu0*k0)/(K+k0),sqrt(tausq/(K+k0)))
  #generate tausq
  alpha=c
  beta=d
  tausq<-1/rgamma(1,alpha,beta)
  #store
  if(ns>Nwarm){
    THETA[ns-Nwarm,]<-theta
    MTS[ns-Nwarm,]<-c(mu,tausq)
  }
}

#time sequence plot
par(mfrow=c(1,2))
plot(MTS[,1], ylab="mu")
plot(MTS[,2], ylab="tausq")

#auto correlation graph
par(mfrow=c(1,2))
acf(MTS[,1], ylab="ACF mu")
acf(MTS[,2], ylab="ACF tausq")

#posterior density function of mu, sigsq, tausq
par(mfrow=c(1,2))
plot(density(MTS[,1]), xlab=expression(mu), main="")
abline(v=quantile(MTS[,1], c(0.025,0.5,0.975)), lty=c(3,2,3))
plot(density(MTS[,2]), xlab=expression(tau^2), main="")
abline(v=quantile(MTS[,2], c(0.025,0.5,0.975)), lty=c(3,2,3))

#inference
mu.hat=mean(MTS[,1])
tausq.hat=mean(MTS[,2])
theta.hat=apply(THETA,2,mean)

#표본경로그림
#plot of theta.hat_i
par(mfrow=c(1,1))
theta.grid=seq(mu.hat-7*sqrt(tausq.hat), mu.hat+7*sqrt(tausq.hat), length=100)
hist(theta.hat, prob=T, main="")
lines(theta.grid, dnorm(theta.grid, mu.hat, sqrt(tausq.hat)))

ytotal.mean=mean(ybar)
par(mfrow=c(1,1))

plot(c(1:K), ybar,  xlab="i")
points(c(1:K), theta.hat, pch=19)
lines(c(1:K), rep(ytotal.mean, K), lw=4)
for(i in 1:K) lines(c(i,i), c(ybar[i], ytotal.mean))
plot(1:K, abs(ybar-theta.hat), xlab="n_i", ylab="|ybar_i-theta.hat_i|")

#95% 신뢰구간
#HPD 
par(mfrow=c(1,1))
plot(density(THETA[,1]), xlim=c(-3,3),xlab="theta_1", main="")
abline(v=quantile(THETA[,1],c(0.025,0.975)), lty=c(3,3))
plot(density(THETA[,2]), xlab="theta_2", main="")
abline(v=quantile(THETA[,2],c(0.025,0.975)), lty=c(3,3))
plot(density(THETA[,3]), xlab="theta_3", main="")
abline(v=quantile(THETA[,3],c(0.025,0.975)), lty=c(3,3))
plot(density(THETA[,4]), xlim=c(-1,1.5), xlab="theta_4", main="")
abline(v=quantile(THETA[,4],c(0.025,0.975)), lty=c(3,3))
plot(density(THETA[,5]), xlim=c(0,2),xlab="theta_5", main="")
abline(v=quantile(THETA[,5],c(0.025,0.975)), lty=c(3,3))
plot(density(THETA[,6]), xlim=c(-2,2), xlab="theta_6", main="")
abline(v=quantile(THETA[,6],c(0.025,0.975)), lty=c(3,3))
plot(density(THETA[,7]), xlim=c(-8,2), xlab="theta_7", main="")
abline(v=quantile(THETA[,7],c(0.025,0.975)), lty=c(3,3))
plot(density(THETA[,8]), xlim=c(-5,5), xlab="theta_8", main="")
abline(v=quantile(THETA[,8],c(0.025,0.975)), lty=c(3,3))

##### JAGS #####
install.packages("rjags") 
library(rjags) 
install.packages("coda") 
library(coda)

modelString="
model
{
  for(i in 1:K){
  y[i]~dnorm(theta[i], invsigsq[i])
  }
  for(i in 1:K){
  theta[i]~dnorm(mu, invtausq)
  }
  mu~dnorm(mu0, k0*invtausq)
  sigsq0=c(0.0629,0.0058,0.0071,0.0356,0.0186,0.036,0.1409,0.0962)
  invsigsq=1/sigsq0
  invtausq~dgamma(c,d)
  sigsq=1/invsigsq
  tausq=1/invtausq
}
"

writeLines(modelString, "model_ex12_2.txt")
K=8
ybar<-c(0.78663,0.76219,0.77715,1.04744,0.8412,0.46245,0.86253,0.69475)
k0=1
mu0<-mean(ybar)
c<-0.5
tausq0<-var(ybar)
d<-c*tausq0

dataList=list(K=K, y=ybar, c=c,d=d,k0=k0, mu0=mu0)
initsList=list(theta=ybar, mu=mu0, invtausq=1/tausq0)

nChains=3
jagsModel=jags.model(file="model_ex12_2.txt", data=dataList, inits=initsList, n.chains=nChains, n.adapt=100)
update(jagsModel, n.iter=1000)
codaSamples=coda.samples(jagsModel, variable.names=c("theta", "mu", "tausq"), n.iter=30000)
variable.names(codaSamples[[1]])

para=c("mu","tausq")
par(mfrow=c(2,2))
for(i in 1:2){
  traceplot(codaSamples[,para[i]], main="", ylab=para[i])
  acf(codaSamples[,para[i]][[1]], plot=TRUE, main="")
}


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
par(mfrow=c(1,2))
for(i in 1:2){
  plot(density(MCMCSamples[,i]), main="", xlab=para[i])
  HPD=HPDsample(MCMCSamples[,i])
  abline(v=HPD, lty=2)
}

#모수의 주변사후 밀도함수
par(mfrow=c(2,4))
for(i in 3:10){
  plot(density(MCMCSamples[,i]), main="", xlab=para[i])
  HPD=HPDsample(MCMCSamples[,i])
  abline(v=HPD,lty=2)
}

