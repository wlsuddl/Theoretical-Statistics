alzheimer<-read.csv("Alzheimer.csv", header=T, na.strings="*")
names(alzheimer)<-c("y","subject","t","lecithin")
head(alzheimer)
attach(alzheimer)

library(rjags)
library(runjags)
library(optimx)
library(coda)
library(dplyr)
library(MASS)
library(ggplot2)

#변수 정리
y<-alzheimer$y; subject<-alzheimer$subject; t<-alzheimer$t; leci<-alzheimer$lecithin

### plotting
placebo<-NULL; lecithin<-NULL
placebo<-alzheimer[1:125,]
lecithin<-alzheimer[126:235,]
plot(placebo$t+1,placebo$y, pch=8, cex=0.5, xlab="t", ylab="y", ylim=c(0,22))
points(lecithin$t+1,lecithin$y)

### ggplot 
ggplot(alzheimer, aes(x=t,y=y,group=leci))+geom_point(aes(color=as.factor(leci)))+
  ggtitle("Scatter plot : placebo vs lecithin")+theme(plot.title=element_text(size=20))

###### GLM modeling with jags #####
alzh<-read.csv("Alzh.csv", header=T, na.strings="*")
names(alzh)<-c("subject","y0","y1","y2","y3","y4","lecithin")
head(alzh)
attach(alzh)
y<-cbind(alzh$y0,alzh$y1,alzh$y2,alzh$y3,alzh$y4)
subject<-alzh$subject 
leci<-alzh$lecithin

# Modeling
modelString="
model
{
  for(i in 1:n){
  for(j in 1:5){
  y[i,j]~dbinom(p[i,j],30)
  }
  }
  
  
  for(i in 1:n){
  for(j in 1:5){
  p[i,j]=(exp(a+c*leci[i]+(d*leci[i]+b)*
  (j-1)/12))/(1+exp(a+c*leci[i]+(d*leci[i]+b)*(j-1)/12))
  
  }
  }
  
  c~dunif(0,1)
  d~dunif(0,1)
  a~dunif(0,1)
  b~dunif(0,1)
  
}
"
  
writeLines(modelString, "model_lecithin.txt")
  
#Prior
n<-length(index)
  
#데이터 리스트
dataList=list(y=y, n=n,leci=leci)
  
#초기치 설정
initsList=list(a=1,b=1,c=1,d=1)
  
#JAGS 실행준비
jagsModel=jags.model(file="model_lecithin.txt", data=dataList, inits=initsList,
                       n.chains=3, n.adapt=500)
#예비단계 수행
update(jagsModel, n.iter=3000)
  
#MCMC 표본생성
codaSamples=coda.samples(jagsModel, variable.names=c("a","b","c","d"), n.iter=5000)
  
#수렴진단
par(mfrow=c(1,2))
coda::traceplot(codaSamples[,"a"], main="", ylab="a")
acf(codaSamples[,"a"][[1]], plot=T, main="")
coda::traceplot(codaSamples[,"b"], main="", ylab="b")
acf(codaSamples[,"b"][[1]], plot=T, main="")
coda::traceplot(codaSamples[,"c"], main="", ylab="c")
acf(codaSamples[,"c"][[1]], plot=T, main="")
coda::traceplot(codaSamples[,"d"], main="", ylab="d")
acf(codaSamples[,"d"][[1]], plot=T, main="")
  
#사후추론
a.Samples=as.matrix(codaSamples[,"a"])
b.Samples=as.matrix(codaSamples[,"b"])
c.Samples=as.matrix(codaSamples[,"c"])
d.Samples=as.matrix(codaSamples[,"d"])

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
para=c("a","b","c","d")
par(mfrow=c(1,2))                  # about c,d
for(i in 1:2){
  plot(density(MCMCSamples[,i+2]), main="", xlab=para[i+2])
  HPD=HPDsample(MCMCSamples[,i+2])
  abline(v=HPD, lty=2)
}
  
HPD.c=HPDsample(MCMCSamples[,3])
HPD.d=HPDsample(MCMCSamples[,4])
HPD.c; HPD.d
mean(c.Samples); mean(d.Samples)

###### Mixed GLM Modeling with jags #####
alzh<-read.csv("Alzh.csv", header=T, na.strings="*")
names(alzh)<-c("subject","y0","y1","y2","y3","y4","lecithin")
head(alzh)
attach(alzh)
y<-cbind(alzh$y0,alzh$y1,alzh$y2,alzh$y3,alzh$y4)
subject<-alzh$subject 
leci<-alzh$lecithin

# Modeling
modelString="
model
{
  
  for(i in 1:n){
  for(j in 1:5){
  y[i,j]~dbinom(p[i,j],30)
  }
  }
  
  
  for(i in 1:n){
  for(j in 1:5){
  p[i,j]=(exp(alpha[i]+beta[i]*((j-1)/12)))/(1+exp(alpha[i]+beta[i]*((j-1)/12)))
  }
  }
  
  for(i in 1:n){
  alpha[i]=a+delta*u[i]+c*leci[i]
  u[i]~dnorm(0,1)
  beta[i]=b+d*leci[i]
  }
  
  delta~dunif(0,1)
  c~dunif(0,1)
  d~dunif(0,1)
  a~dunif(0,1)
  b~dunif(0,1)
  
}
"
writeLines(modelString, "model_lecithin_mixed.txt")

#Prior
index<-which(t==0) 
n<-length(index)
ui<-rnorm(n,0,1)

#데이터 리스트
dataList=list(y=y, n=n,leci=leci)

#초기치 설정
initsList=list(a=1,b=1,c=1,d=1,delta=1,u=ui)

#JAGS 실행준비
jagsModel=jags.model(file="model_lecithin_mixed.txt", data=dataList, inits=initsList,
                     n.chains=3, n.adapt=500)
#예비단계 수행
update(jagsModel, n.iter=1000)

#MCMC 표본생성
codaSamples=coda.samples(jagsModel, variable.names=c("a","b","c","d","delta"),
                         n.iter=3000)
#수렴진단
par(mfrow=c(1,2))
coda::traceplot(codaSamples[,"c"], main="", ylab="c")
acf(codaSamples[,"c"][[1]], plot=T, main="")
coda::traceplot(codaSamples[,"d"], main="", ylab="d")
acf(codaSamples[,"d"][[1]], plot=T, main="")

#사후추론
MCMCSamples=as.matrix(codaSamples)
a.Samples=as.matrix(codaSamples[,"a"])
b.Samples=as.matrix(codaSamples[,"b"])
c.Samples=as.matrix(codaSamples[,"c"])
d.Samples=as.matrix(codaSamples[,"d"])
delta.Samples=as.matrix(codaSamples[,"delta"])
nn<-length(a.Samples)
ui<-matrix(0,nn,47)
for(i in 1:47){
  ui[,i]<-MCMCSamples[,(i+5)]
}
par(mfrow=c(1,1))
plot(density(a.Samples), xlab=expression(a), ylab="posterior density",main="")
plot(density(b.Samples), xlab=expression(b), ylab="posterior density",main="")
plot(density(c.Samples), xlab=expression(c), ylab="posterior density",main="")
plot(density(d.Samples), xlab=expression(d), ylab="posterior density",main="")

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

para=c("a","b","c","d","delta")
par(mfrow=c(1,2))                  # about c,d
for(i in 3:4){
  plot(density(MCMCSamples[,i]), main="", xlab=para[i])
  HPD=HPDsample(MCMCSamples[,i])
  abline(v=HPD, lty=2)
}
HPD.c=HPDsample(MCMCSamples[,3])
HPD.d=HPDsample(MCMCSamples[,4])
HPD.c; HPD.d
mean(c.Samples); mean(d.Samples)
mean(a.Samples); mean(delta.Samples)
