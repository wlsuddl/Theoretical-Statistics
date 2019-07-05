install.packages("dplyr")
library(dplyr)

getwd()
setwd("/Users/Jinyoung-Kim/Google 드라이브/대학원1학기/이통1/과제6")
bankruptcy<-read.csv("bankruptcy.csv", header=T, sep=",")
bankrup<-bankruptcy[,-c(1,3)]        ## 1,3 열 제거
head(bankrup)
y<-bankrup$Y
X<-bankrup[,-1]
attach(bankrup)

bank.glm<-glm(Y~., data=bankrup, family=binomial(link=logit))

##backward elimination, graph at p 
bank.glm<-glm(Y~.,family=binomial,data=bankrup)
dat1<-bankrup
min.index<-NULL; result.AIC<-NULL; c<-NULL; variable<-NULL

for(i in 1:23){
  min.index[i]<-which.min(abs(summary(bank.glm)$coefficients[-1,3]))
  c<-(min.index+1)                     #절편이 첫번째 변수이므로 
  dat1<-dat1[,-c[i]]
  variable[[i]]<-names(dat1)
  bank.glm<-glm(Y~.,data=dat1,family=binomial)
  result.AIC[i]<-extractAIC(bank.glm)[2]
}

which.min(result.AIC)        ##19번째 결과에 해당하는 변수들이 포함된 모형의 AIC가 최소.                             
result.AIC[which.min(result.AIC)]
variable[[19]]                 ## R5, R6, R16, R18, R22 

result.AIC[24]<-glm(Y~1, data=X, family=binomial(link=logit))$aic          ##only intercept model
plot(seq(1,24),result.AIC[24:1],type="l",xlab="p",ylab="AIC",col="orange",lwd=2)
points(5, result.AIC[19],col="blue",cex=2,pch=20) 

#R18, R5, R16, R6, R22 넣었을 때 AIC
lm=glm(Y~R18+R5+R16+R6+R22, data=dat1, family=binomial(link=logit))
summary(lm)$aic        #동일한 AIC값     
result.AIC

### 최적 예측 모형 구하기
lm=glm(Y~R5+R6+R16+R18+R22, data=dat1, family=binomial(link=logit))


### LRT 검정
p.lm=predict(lm, type='response')
p.full=predict(bank.glm, type='response')
loglik.lm=sum(Y*log(p.lm)+(1-Y)*log(1-p.lm))
loglik.full=sum(Y*log(p.full)+(1-Y)*log(1-p.full))
LRT=-2*(loglik.lm-loglik.full)

1-pchisq(LRT,19)                  # df=19

### 최종 예측 모형의 산점도
finaldata<-select(bankrup,variable[[19]])
x<-select(finaldata,variable[[19]])[,-1]
y<-finaldata$Y
group1<-filter(finaldata,Y==1)
group0<-filter(finaldata,Y==0)
z.group0<-predict(lm)[1:50]                 #no=10, n1=50
z.group1<-predict(lm)[51:100]

#scatterplot of predictors
pairs(finaldata[,2:6], pch = c(1,3)[y+1],col=c("sky blue","orange")[y+1])   #0=bankruptcy(sky blue)

#boxplot
boxplot(data.frame(bankruptcy=z.group0,solvent=z.group1),main="Boxplot of Z")

#kernel density plot
plot(density(z.group1), col=2, ylim=c(0,0.21), main='kernel density plot')
lines(density(z.group0, adjust=TRUE), col=4, lty=2)

legend("topright", legend=c('solvent','bankruptcy'), lty=c(1,2), col=c(2,4))

