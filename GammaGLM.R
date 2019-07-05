library(MASS)

seed<-read.csv("yield.csv")
yi <- seed[,5]; xi <- seed[,2]

###### fitting gamma GLM
gamma.model <- glm(yi~I(1/xi), family=Gamma)
summary(gamma.model, dispersion=gamma.dispersion(gamma.model))

###### MLE and 95% CI of alpha
#> gamma.shape(gamma.model)
#Alpha: 11.267622
#SE:     2.867251

alpha<- as.numeric(gamma.shape(gamma.model)[1])              
se <- as.numeric(gamma.shape(gamma.model)[2])
alpha; se
(alpha.ci <- c(alpha-1.96*se, alpha+1.96*se))

####### MLE and 95% CI of beta
confint(gamma.model)


###### 블록 효과와 비선형 항을 포함하는 Gamma GLM 모형에 변수선택 기법 적용 #####
###### 최적 모형 선택
B <- as.factor(seed[,3])

## variable selection
gamma.model2<-glm(yi~B+B*I(1/xi)+xi, family=Gamma)
summary(gamma.model2, dispersion=gamma.dispersion(gamma.model2))    # AIC=330.91
gamma.step<-step(gamma.model2,trace=0)
summary(gamma.step, dispersion=gamma.dispersion(gamma.step))        # AIC=323.14

alpha2<- as.numeric(gamma.shape(gamma.step)[1])              
se2 <- as.numeric(gamma.shape(gamma.step)[2])
(alpha2.ci <- c(alpha2-1.96*se2, alpha2+1.96*se2))

####### MLE and 95% CI of beta
confint(gamma.step)

## zi : number of plants harvested variable
zi <-seed[,4]

# 최적모형 찾기
glm(cbind(zi,xi-zi)~B+xi+B*xi,family="binomial"(link="probit"))
glm(cbind(zi,xi-zi)~B+xi+B*xi,family="binomial"(link="logit"))
glm(cbind(zi,xi-zi)~B+xi+B*xi,family="binomial"(link="cloglog"))


