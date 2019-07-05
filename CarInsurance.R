library(ggplot2)

car<-read.table("car.txt", header=T)
head(car)

##변수간 plot 그림으로 주효과와 교호작용 효과 알 수 있다. 
qplot(car$Merit, car$Claims/car$Insured) + 
  geom_line(aes(group = car$Class, color = factor(car$Class)), size=1)+
  labs(title="Merit vs. Claims",x ="Meirt", y = "Claims per Insured")  
qplot(car$Class, car$Claims/car$Insured) + 
  geom_line(aes(group = car$Merit, color = factor(car$Merit)), size=1)+
  labs(title="Class vs. Claims",x ="Class", y = "Claims per Insured")

car$Merit<-factor(car$Merit); car$Class<-factor(car$Class)

####### 포아송 분포 가정 시, 적절한 회귀모형 구축 #######
glm1<-glm(Claims~.*., data=car, family=poisson(link="log"))
step1<-step(glm1, direction="both")

#### 감마 분포 가정 시 #####
glm2<-glm(Cost~., data=car, family=Gamma(link="identity"))
#적정 보험료 = 1인당 평균 보험 청구* 1건당 평균 보험 지급액
lambda<-step1$fitted.values/car$Insured
mu<-glm2$fitted.values/car$Claims
ins <- lambda*mu                # 적정 보험료 계산식.
cbind(Merit=car$Merit, Class=car$Class, Optimal.Premium=ins)

dat4<-data.frame(lambda=lambda, mu=mu)
summary(lm(mu~lambda))

#그룹별 평균사고빈도 와 평균사고심도 산점도
ggplot(data=dat4, aes(x=lambda, y=mu))+geom_point()+geom_smooth(method="lm")+
  annotate("text", x =0.11, y=0.35, label = "mu= 0.24788 + 0.40984*lambda")+ 
  annotate("text", x =0.11, y=0.34, label = "R2.adj=0.1366 ")+
  labs(x="lambda",y="mu")+ggtitle("lambda vs mu")
