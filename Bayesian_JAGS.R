install.packages("rjags") 
library(rjags)
install.packages("coda") 
library(coda)

##### Poisson 분포 추정 #####
#모델링 
modelString2=" model
{
  s ~ dpois(n*theta)
  theta ~ dgamma(a,b)
  xnew ~ dpois(theta)
}
"
writeLines(modelString2, "model_ex11_2.txt")

#데이터 리스트
dataList2=list(n=50, s=102, a=2, b=1)

#초기치 설정 
initsList2=list(theta=2.04)

#JAGS 실행준비
jagsModel2=jags.model(file="model_ex11_2.txt", data=dataList2, inits=initsList2,n.chains=3, n.adapt=500)

#예비단계 수행
update(jagsModel2, n.iter=500)

#MCMC 표본생성
codaSamples2=coda.samples(jagsModel2, variable.names=c("theta","xnew"), n.iter=5000)

#수렴진단
par(mfrow=c(1,2))
coda::traceplot(codaSamples2[,"theta"], main="", ylab="theta") 
acf(codaSamples2[,"theta"][[1]], plot=T, main="") 
coda::traceplot(codaSamples2[,"xnew"], main="", ylab="xnew") 
acf(codaSamples2[,"xnew"][[1]], plot=T, main="")

#사후추론 
thetaSamples=as.matrix(codaSamples2[,"theta"]) 
xnewSamples=as.matrix(codaSamples2[,"xnew"])

par(mfrow=c(1,2))
plot(density(thetaSamples), xlab=expression(theta), ylab="posterior density",main="")
hist(xnewSamples, xlab="x_new", ylab="posterior density", main="")

#### 디리슐레 #####
#모델링 
modelString3=" model
{
  Y[1:K] ~ dmulti( pi[], N)
  pi[1] <- a1 * theta + b1
  pi[2] <- a2 * theta + b2
  pi[3] <- a3 * eta + b3
  pi[4] <- a4 * eta + b4
  pi[5] <- c * (1-theta-eta)
  alpha[1:3] ~ ddirch(a[1:3])
  theta <- alpha[1]
  eta <- alpha[2]
  a1<- 1/4
  b1<- 1/8
  a2<- 1/4
  b2<- 1/8
  a3<- 1/4
  b3<- 0
  a4<- 1/4
  b4<- 3/8
  c<- 1/2
}
"
writeLines(modelString3, "model_num11_1.txt")

#데이터 리스트
dataList3=list(N=22, Y=c(14,1,1,1,5), K=5, a=c(1,1,1))

#JAGS 실행준비
jagsModel3=jags.model(file="model_num11_1.txt", data=dataList3, n.chains=3, n.adapt=500)

#예비단계 수행
update(jagsModel3, n.iter=500)

#MCMC 표본 생성
codaSamples3=coda.samples(jagsModel3, variable.names=c("theta", "eta"), n.iter=5000)

#수렴진단
coda::traceplot(codaSamples3[,"theta"], main="", ylab="theta") 
acf(codaSamples3[,"theta"][[1]], main="", plot=T) 
coda::traceplot(codaSamples3[,"eta"], main="", ylab="eta") 
acf(codaSamples3[,"eta"][[1]], main="", plot=T)

#사후추론 
thetaSamples2=as.matrix(codaSamples3[,"theta"]) 
etaSamples=as.matrix(codaSamples3[,"eta"])

plot(density(thetaSamples2), xlab="theta", ylab="posterior density", main="")
plot(density(etaSamples), xlab="eta", ylab="posterior density", main="")


##### Multi Norm #####
#모델링 
modelString4=" model
{
  for(i in 1:nm){
    x[1:2,i] ~ dmnorm(pi, Omega)
  }
  pi<-c(mu1,mu2)
  mu1 ~ dnorm(mu0, invsigsq0)
  mu2 ~ dnorm(mu0, invsigsq0)
  mu0<-0
  invsigsq0<-1/300
}
"
writeLines(modelString4, "model_num11_2.txt")

#데이터 리스트
dataList4=list(nm=30, Omega=matrix(c(1,0.5,0.5,1),ncol=2))

#초기치 설정
initsList4=list(mu1=0.3, mu2=0.58)

#JAGS 실행 준비
jagsModel4=jags.model(file="model_num11_2.txt", data=dataList4, inits=initsList4, n.chains=3, n.adapt=500)

#예비단계 수행 
update(jagsModel4, n.iter=500)

#표본생성
codaSamples4=coda.samples(jagsModel4, variable.names=c("mu1","mu2"), n.iter=5000)
par(mfrow=c(1,2))

#수렴진단
coda::traceplot(codaSamples4[,"mu1"], main="", ylab="theta1") 
acf(codaSamples4[,"mu1"][[1]], plot="T", main="")
coda::traceplot(codaSamples4[,"mu2"], main="", ylab="theta2")
acf(codaSamples4[,"mu2"][[1]], plot="T", main="")

#사후추론 
mu1samples=as.matrix(codaSamples4[,"mu1"]) 
mu2samples=as.matrix(codaSamples4[,"mu2"])
plot(density(mu1samples), xlab="theta1", ylab="posterior density", main="")
plot(density(mu2samples), xlab="theta2", ylab="posterior density", main="")