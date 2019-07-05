lung<-read.csv("lung.csv", header=T, na.strings="*")
names(lung)<-c("age","smoke","pop","cancer")
head(lung)
attach(lung)

####### 변수 정리 #########
y<-lung$cancer; n<-lung$pop; x1<-lung$age; x2<-lung$smoke
ybar<-y/n; ai<-(37+5*x1)

##### Poisson GLM 적합 ######
as.numeric(smoke)
as.numeric(age)

dat1<-data.frame(lnn=log(n),lna=log(ai),S=as.factor(smoke),yi=y)
h1<-glm(yi~.,data=dat1,family="poisson")

dat2<-data.frame(lnn=log(n),Z=as.factor(age),S=as.factor(smoke),yi=y)
h2<-glm(yi~., data=dat2, family="poisson")

#모형선택
h3 <- glm(yi~offset(lnn)+lna+S,data=dat1, family=poisson)
anova(h3,h1, test="Chisq")

#모형 적합도
deviance.r<-summary(h3)$deviance.resid
di<-h3$residuals             # residuals
D<-h3$deviance
shapiro.test(di)             # 잔차 정규성 검정.
# p-value < 0.05 이므로 정규성 따름.
mu<-fitted(h3)              # fitted values

pearson.r<-(y-mu)/sqrt(mu)     # pearson residuals
Xsquare<- (pearson.r)^2

X.square <- sum(Xsquare)       # 62.41 == (m-p)
par(mfrow=c(2,2))
plot(h3)

# deviance residual Q-Q plot
par(mfrow=c(1,2))
qqnorm(deviance.r, pch=20, col="sky blue", main="deviance residual")
qqline(deviance.r,col="navy",lty=3)

# pearson residual Q-Q plot
qqnorm(pearson.r, pch=20, col="light green", main="pearson residual")
qqline(pearson.r, col="dark green",lty=3)

##### 변수변환,교호작용 고려한 모형에 AIC 판정기준을 이용한 변수선택
sel.model <-glm(yi~offset(lnn)+lna*S,data=dat1, family=poisson)        #interaction
summary(sel.model)
final.model<-step(sel.model)
summary(final.model)

## logistic GLM
dat3<-data.frame(lna=log(ai),S=as.factor(smoke),yi=y,ni=n)
logit<-glm(cbind(yi,ni-yi)~S+lna,data=dat3,family="binomial"(link="logit"))
summary(logit)
