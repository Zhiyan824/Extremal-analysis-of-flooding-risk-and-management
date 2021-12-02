
##bond simulation 
##Rain level modeling 
##package 
library(dplyr)
library(tseries)
library(extRemes)
library(fExtremes)
library(ggplot2)

library(extRemes)
library(timeDate);library(timeSeries);library(fBasics);
library(fGarch);library(fExtremes)

##time series
rainlevel_ts <- ts(rainlevel,frequency = 12,start = c(1/1/2006))
plot(rainlevel_ts,xlab='Starting from 2006 loss',ylab='amount of loss')

##fitting
par(mfrow=c(2,2))
rainlevel_ts
fit <- fevd(rainlevel_ts, threshold = 600 , type = 'GP')
profliker(fit, type="parameter")
summary(fit)
plot(fit)


## Triggering procedure

## preparation
install.packages('SMFI5')
install.packages("xlsx")
library(extRemes)
library(SMFI5)
library(xlsx)

## initial setting 
lambda <- 2.38
threshold <- 600
maxTime <- 3

##coupon period, couponPayDate in years and face value K in millions
delta <- 1/4
s <- seq(from = 1, to = round(maxTime/delta), by = 1)
K <- 1000

## from fitting result
scale <- 258.5489066
shape <- -0.1805117
ratio<-  scale/shape
maximum<- threshold-scale/shape

## iteration for price for m=10^5 path
m <- 10^3

## interest rate parameter estimation and simulation
setwd("C:/Users/win/Desktop")
TreasuryBillData <- "C:/Users/win/Desktop/treasurybilldata.xlsx"
treasuryBillFrame <- read.xlsx(TreasuryBillData ,1)

ShiborData <- "C:/Users/win/Desktop/Shibordata.xlsx"
Shiborframe <- read.xlsx(ShiborData,1)

set.seed(10)
outr = est.vasicek(treasuryBillFrame, days = 365) 
alphar = outr$param[1] 
betar = outr$param[2]/100
thetar = outr$param[3]/100
r0 <- 2.28 
rt = sim.vasicek(outr$param[1], outr$param[2], outr$param[3], r0, 365*maxTime, 1/365)

outl = est.vasicek(Shiborframe, days = 365)
alphal = outl$param[1] 
betal = outl$param[2]/100
thetal = outl$param[3]/100
l0 <- 2.425
lt = sim.vasicek(outl$param[1], outl$param[2], outl$param[3], l0, 365*maxTime, 1/365)


## assumption for R coupon rate 
couponRate <- 3.4

##corrlation coefficient
setwd("C:/Users/win/Desktop")
comparePath <- "C:/Users/win/Desktop/comparationdata.xlsx"
compareframe <- read.xlsx(comparePath,1)
a=compareframe[,1]
b=compareframe[,2]
cor<- cor(a,b,method="pearson")


## preparation of prcing function 
Bts <- function(t, s, alphar){
  (1-exp(-alphar*(s-t)))/alphar
}

Ats <- function(t, s, alphar, betar, thetar){
  B = Bts(t, s, alphar)
  exp( (B - (s-t) )*(alphar^2*betar - thetar^2/2)/alphar^2 - thetar^2*B^2/(4*alphar) )
}

EQ2Dts<- function(t, s, alphar, betar, thetar, r0, rt){
  if(t==0)
  {   
    Ats(t, s, alphar, betar, thetar)*exp(-Bts(t, s, alphar)*r0/100)
  } 
  else 
  {
    Ats(t, s, alphar, betar, thetar)*exp(-Bts(t, s, alphar)*rt[round(t*365)]/100)
  }
}

##Appendex realization A.4.
C1ts <- function(t, s, alphar, alphal, betar, betal, thetar, thetal, cor){
  (betar-thetar^2/(2*alphar))*(s-t)+3* thetar^2/(4*alphar^2) + cor*thetal*thetar/(alphal*(alphal-alphar))
  + betal-betar/alphar
}


C2ts<- function(t, s, alphar, alphal, betar, betal, thetar, thetal,cor){
  thetar^2*exp(-2*alphar*(s-t))/(4*alphar^2) + (betar/alphar - thetar^2/alphar^2) *exp(-alphar*(s-t))
  + (cor*thetar*thetal/(alphal*alphar)-betal)*exp(alphal*(s-t))
  - cor*thetar*thetal/(alphar*(alphal-alphar))*exp((alphal-alphar)*(s-t))
  - thetal^2/(4*alphal)*exp(2*alphal*(s-t))
}


pBts<- function(alphal, t, s){
  exp(-alphal*(s-t))
}

pAts = function(t, s, alphar, alphal, betar, betal, thetar, thetal, cor){
  exp(-(C1ts(t, s, alphar, alphal, betar, betal, thetar, thetal, cor) + C2ts(t, s, alphar, alphal, betar, betal, thetar, thetal, cor)))
}

EQ2Dtsi<- function(t, s, alphar, alphal, betar, betal, thetar, thetal, cor, r0, l0, rt){
  if(t==0){
    pAts(t, s, alphar, alphal, betar, betal, thetar, thetal, cor)*
      exp(-Bts(t, s, alphar)*r0+pBts(alphal, t, s)*l0)
    - Ats(t,s,alphar, betar,thetar)*exp(-Bts(t,s,alphar)*r0)
  }
  else {
    pAts(t, s, alphar, alphal, betar, betal, thetar, thetal, cor)*exp(-Bts(t,s,alphar)*rt[round(365*t)]/100
    +pBts(alphal,s,t)*lt[round(365*t)]/100)
    - Ats(t,s,alphar, betar,thetar)*exp(-Bts(t,s,alphar)*rt[round(365*t)]/100)
  }
}


## inter-arrival time generation
waitsGen <- function(lambda, maxTime)
{
  #Input: 
  #     lambda, the intensity of Poisson process $N(t)$
  #     maxTime, the upper bound, i.e., $t \in [0,maxTime]$
  #Output: 
  #     head(waits,-1), sup{k: W_1 + W_2 +\cdots+ W_k<= maxTime} 
  #     the number of samples of inter-arrival times required 
  #     simulation such that the total time <= maxTime
  waits <- NULL
  i <- 1
  while(sum(waits) < maxTime)
  {
    samp <- rexp(n = 1, rate = lambda)
    waits[i] <- samp
    i <- i+1
  }
  return(head(waits, -1))
}

arrivalsGen <- function(lambda, maxTime)
{ 
  # arrivals[i]= W_1+\cdots+W_i, i\le waits
  # the arrival time of all arriving events 
  waits <- waitsGen(lambda, maxTime)
  arrivals <- NULL
  cumSum <- cumsum(waits)
  return(cumSum)
}


ppGen <- function(lambda, maxTime)
{ 
  # pp, the number of events occuring in [0,t]
  # for any timeslot t \in [0, maxTime]. 
  arrivals <- arrivalsGen(lambda, maxTime)
  maxJumps <- length(arrivals)
  pp <- NULL
  for(i in 1 : maxJumps){
    pp[i] <- sum(arrivals <= arrivals[i])
  }
  return(cbind(arrivals,pp))
}

## given trigger value according to every arrival spot time; 
triggerGen <- function(fraction, u, scale, shape, lambda, maxTime, kappa, threshold) {
  ## Input: 
  # fraction, proportion of principal wiped out if trigger
  # u, used to define trigger if severity in those intervals with endpoints specified
  # scale, shape, threshold used to generated distorted severity using GP model based on POT
  # kappa, parameter to define Wang's distortion
  ## Output:
  # data.frame with arrival spot times, severity, trigger process

  arrivals <- arrivalsGen(lambda, maxTime)
  arrivals
  
  maxJumps <- length(arrivals)
  if (maxJumps == 0){ 
    trigger <- 0 
    severity <- 0
    arrival <- 0 
    return(t(as.matrix(c(0,0,0))))
  }else{
  severity = rep(NA, maxJumps)
 
  for (i in 1 : maxJumps){
    
    severity[i] = ((1-pnorm(qnorm(runif(1,min=0,max=1))+kappa))^(-shape)-1)*scale/shape + threshold
    
  }
  
  trigger <- rep(1,maxJumps)
  L <- length(fraction)
  N <- matrix(nrow = L, ncol = maxJumps)
  for(i in 1 : maxJumps){
    for(j in 1:(L-1))
    {
      N[j,i] <- sum(severity[1:i]> u[j] & severity[1:i] <= u[j+1])
      N[j,i]
    }
    N[L,i] <- sum(severity[1:i]> u[L])
    trigger[i] <- fraction %*% N[,i]
  }
  trigger
  
  df <- data.frame(arrivalST = arrivals, severity = severity, triggerresult = trigger)
  return(df)
  }
}

## generate price
priceGen <- function(fraction, u, scale, shape, lambda, maxTime, kappa, 
                     K, alphar, alphal, betar, betal, thetar, 
                     thetal, cor, r0, l0, rt, threshold, delta, couponRate)
{

  ## calling function 

  df <- triggerGen(fraction, u, scale, shape, lambda, maxTime, kappa, threshold)

  arrivals <- df[,1]
  arrivals
  maxJumps <- length(arrivals)
  trigger <- df[,3]
  trigger
  maxJumps
  ## length of arrival=0

  ## the payoff PI(Yt), t = arrivalST in Example 5.2
  arrivalpayoff <- rep(NA, maxJumps)
  if (length(arrivals)==0){
    couponPDatePayoff <- rep(1, couponMaxNumber)
  }else{
  for (i in 1:maxJumps){
    arrivalpayoff[i] <- max(1-trigger[i], 0)
  }
  
  arrivalpayoff
  if (any(trigger >= 1)){
    tau <- arrivals[min(which(trigger >= 1))] 
  }else{
    tau <- maxTime+0.0001
  }
  tau
  couponMaxNumber <- round(min(maxTime,tau)/delta)
  if(couponMaxNumber==0){
    PriceCoupon<- 0
  }else{
    couponPDatePayoff <- rep(1, couponMaxNumber)
    couponPDatePayoff
    ## pricing for the first part of P0
    PriceCouponPart <- rep(0, couponMaxNumber)
    for (i in 2:couponMaxNumber){
      if(any(arrivals <= delta*(i-1))){
        couponPDatePayoff[i] <- 
          arrivalpayoff[max(which(arrivals <= delta*(i-1)))]
      }
    }
    couponPDatePayoff}}
  
  
    for (s in 1:couponMaxNumber){
      PriceCouponPart[s] <-  K*delta*couponRate *couponPDatePayoff[s]*EQ2Dts(0,s*delta,alphar,betar,thetar,r0,rt)
      + K*delta*couponPDatePayoff[s]*EQ2Dtsi(0, delta*s, alphar, alphal, betar, betal, thetar, thetal, cor, r0, l0, rt)##
    } 
    PriceCouponPart
    PriceCoupon <- sum(PriceCouponPart)
  
 
  if (0 <= tau && tau < maxTime){
    ## the final price when tau less than coupon date means principal has been run off 
    ProporPrincPrice <-  K*delta*( tau - floor(tau/delta)*delta )*
      ( ifelse(tau< delta, 1, couponPDatePayoff[floor(tau/delta)]) ) * 
      ( couponRate * EQ2Dts(0, tau, alphar, betar, thetar, r0, rt)
        + EQ2Dtsi(0, tau, alphar, alphal, betar, betal, thetar, thetal, cor, r0, l0, rt) )
  }else{
    ProporPrincPrice <- K*arrivalpayoff[maxJumps] * EQ2Dts(0, maxTime, alphar, betar, thetar, r0, rt)
  }
  Price <- PriceCoupon + ProporPrincPrice
  Price
  return(Price)
}

couponRate<- 0.34
alphar
alphal
betar
betal
thetal
thetar
r0
l0
m <- 100000 ## the time of simulation 
threshold<- 600
cor
u <- c(626.4, 744.2, 849.0, 985.1) #quantile(data, seq(0.6, 0.8, by=0.1))
fraction <- c(0.005, 0.015, 0.15, 0.2)

## kappa V.S. Price
kappa <- seq(from = 0.02, to = 1.5, by = 0.05)
nk <- length(kappa)
pricekappa <- rep(NA, nk)
pricev <- rep(0,m)
for (k in 1:nk){
  print(k)
  for (t in 1:m){
    ## iteration for price for m=10^5 path
    pricev[t] <- priceGen(fraction, u, scale, shape, lambda, maxTime, kappa[k], 
                          K, alphar, alphal, betar, betal, thetar, 
                          thetal, cor, r0, l0, rt, threshold, delta, couponRate)
  }
  pricekappa[k] <- mean(pricev)
}
print(pricekappa)

 # draw chart for kappa V.S. price
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(kappa, pricekappa, ann=F)



fitted.model<-lm(pricekappa~kappa)
summary(fitted.model)
library(ggplot2)
b<-fitted.model$coefficients[[1]]
a<-fitted.model$coefficients[[2]]

fitted.curve<-function(y){
  return((y-b)/a)
}
fitted.curve(1000)

ggplot()+
  geom_point(aes(x=kappa,y=pricekappa),
             size=3,color="black",alpha=0.3)+
  geom_abline(intercept = fitted.model$coefficients[[1]],
              slope = fitted.model$coefficients[[2]],
              size=2,color="black",alpha=0.8)+
  theme_bw()+
  geom_hline(yintercept = 1000,lty="dashed")+
  annotate(geom = "segment",x=fitted.curve(1000),
           xend = fitted.curve(1000),y=1000,yend = -Inf,
           lty="dashed")+
  geom_point(aes(x=fitted.curve(1000),y=1000),size=1,shape=6,
             color="black",alpha=0.9)+
  annotate(geom = "text",x=fitted.curve(1000),y=1000,
           label=round(fitted.curve(1000),2),
           vjust=4,color="black")+
  labs(x="kappa",y="Price")


##sensitivity for shape parameter
m=100000
kappa<- 0.42
shape<- seq(from=-0.5 , to=-0.1 , by= 0.005)
nk <- length(shape)
priceshape <- rep(NA, nk)
pricev <- rep(0,m)
for (k in 1:nk){
  scale<- shape[k]*-1432.311
  for (t in 1:m){
    
    ## iteration for price for m=10^5 path
    pricev[t] <- priceGen(fraction, u, scale, shape[k], lambda, maxTime, kappa, 
                          K, alphar, alphal, betar, betal, thetar, 
                          thetal, cor, r0, l0, rt, threshold, delta, couponRate)
  }
  priceshape[k] <- mean(pricev)
}


fitted.model<-lm(priceshape~shape)
summary(fitted.model)


par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(shape,priceshape,ann=F)
library(ggplot2)
b<-fitted.model$coefficients[[1]]
a<-fitted.model$coefficients[[2]]

fitted.curve<-function(y){
  return((y-b)/a)
}
fitted.curve(1000)

ggplot()+
  geom_point(aes(x=shape,y=priceshape),
             size=3,color="black",alpha=0.3)+
  geom_smooth(aes(x=shape,y=priceshape),color="blue")+
  theme_bw()+
  labs(x="shape parameter",y="Price")

par(pin = c(5,3))

lines(lowess(shape,priceshape)) 



##plot distortion plot

f = function (x,kappa) {
  1- pnorm(qnorm(1-(1+(-0.1805117)*x/258.5489066)^(1/0.1805117) ) - kappa )
  
}

kappa<- 0

x = seq(0.5, 900, length=3000)
y = rep(0, length(x))
for (i in 1:length(x)) {
  y[i] = f(x[i],kappa)
}

plot(x, y, type='l',lwd = 4,pch=15,lty=1,col='red')

kappa<- 1

x = seq(0.5, 900, length=3000)
y = rep(0, length(x))
for (i in 1:length(x)) {
  y[i] = f(x[i],kappa)
}


lines(x, y, col='blue',lwd = 4,pch=16,lty=2)

kappa<- 1.5

x = seq(0.5, 900, length=3000)
y = rep(0, length(x))
for (i in 1:length(x)) {
  y[i] = f(x[i],kappa)
}

lines(x, y, col='darkorchid4',lwd = 4, pch=17,lty=3)

kappa<- 2

x = seq(0.5, 900, length=3000)
y = rep(0, length(x))
for (i in 1:length(x)) {
  y[i] = f(x[i],kappa)
}

lines(x, y,lwd = 4, pch=18,lty=4)
legend(720,1.036,cex=0.8,c("0","1","1.5","2"),col=c("red","blue","darkorchid4","black"),text.col=c("red","blue","darkorchid4","black"),pch=c(15,16,17,18),lty=c(1,2,3,4),title="kappa")
