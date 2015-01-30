rslt<-NULL
n.data<- 500
sample.size<- 500
sigma.U<-1
sigma.X<-1
sigma.Z<-1
sigma.T<-1

rslt <- numeric(n.data)
for (i in 1:n.data)
{
  U<- rnorm(sample.size, sd=sigma.U)
  e.X<- rnorm(sample.size, sd=sigma.X)
  e.T<- rnorm(sample.size, sd=sigma.T)
  Z<- rnorm(sample.size, sd=sigma.Z)
  X<- exp(U+Z+e.X)
  logt<- X+U+e.T
  
  X.lm<- lm(log(X)~Z)
  X.pred<- exp(predict(X.lm))  
  T.lm<- lm(logt~X.pred)
  rslt[i] <- summary(T.lm)$coef[2,1]
}

summary(rslt)


# removing transformation
rslt <- numeric(n.data)
for (i in 1:n.data)
{
  U<- rnorm(sample.size, sd=sigma.U)
  e.X<- rnorm(sample.size, sd=sigma.X)
  e.T<- rnorm(sample.size, sd=sigma.T)
  Z<- rnorm(sample.size, sd=sigma.Z)
  X<- exp(U+Z+e.X)
  logt<- X+U+e.T
  
  X.lm<- lm(X~Z)
  X.pred<- (predict(X.lm))  
  T.lm<- lm(logt~X.pred)
  rslt[i] <- summary(T.lm)$coef[2,1]
}

summary(rslt)
