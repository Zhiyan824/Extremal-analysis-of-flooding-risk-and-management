library(dplyr)
library(extRemes)
library(fExtremes)
library(MASS)
library(scales)
library(nortest)

loss <- read.csv("loss.csv", sep = ",")
loss$billion <- loss$eco_loss/10
plot(loss$Maxpoint, log(loss$billion), ylab = "national economic loss (on log scale)", xlab = "Maximum Point Precipitation")+abline(lm(log(loss$billion)~loss$Maxpoint),col="red")
hist(loss$billion,main="", xlab="Economic loss (billion)")
hist(loss$Maxpoint, main="", xlab="Max point Precipitation (mm)")
cor.test(loss$Maxpoint, log(loss$billion)) ## Duration effect

### economic loss descriptive analysis
mean(loss$billion)
median(loss$billion)
range(loss$billion)
mean(((loss$billion-mean(loss$billion))/sd(loss$billion))^3) # value greater than 0, right-skew
mean(((loss$billion-mean(loss$billion))/sd(loss$billion))^4) # value greater than 3, heavy tail

### max point precipitation
mean(loss$Maxpoint)
median(loss$Maxpoint)
range(loss$Maxpoint)
mean(((loss$Maxpoint-mean(loss$Maxpoint))/sd(loss$Maxpoint))^3) # right skew
mean(((loss$Maxpoint-mean(loss$Maxpoint))/sd(loss$Maxpoint))^4) # heavy tail

mrlplot(loss$billion,xlab = "threshold")
mePlot(loss$eco_loss,doplot=TRUE) 

par(mfrow=c(3,1), mar=c(2,2,2,2))
threshrange.plot(loss$billion, r=c(15,30), type = 'PP', nint = 50)
u <- quantile(loss$billion,0.7)

## PP model
PP <- fevd(loss$billion, threshold = u, type = 'PP')
summary(PP)
location1 <- 84.34
scale1 <- 6.64
shape1 <- -0.29
ci(PP, type = "parameter") 
# confidence interval of three parameters

## GP model
GP <- fevd(loss$billion, threshold = u, type = "GP")
summary(GP)
ci(GP, type = "parameter")
summary(gpdFit(loss$billion, u, type = "mle"))

# Exponential(shape = 0)
Exp <- fevd(loss$billion, threshold = u, type = "Exponential")
summary(Exp)
ci(Exp, type = "parameter")

lr.test(Exp,PP) ## choose PP model, PP is better (p<2.2e-16)
lr.test(GP,PP) ## choose PP model, PP is better (p<2.2e-16)

# diagnostic plots (density plot cannot be plotted)
plot(PP,"qq") # qq plot of the data quantile against the fitted model quantile
plot(fit,"qq2")# of quantile from model-simulated data against the fitted model quantile

# extract exceedance
exceedances <- subset(loss, loss$billion > u)
frequency <- as.data.frame(table(exceedances$Year))

# estimation of lambda
pois <-  fitdistr(frequency$Freq, "Poisson")
lambda <- pois$estimate
lambda

# test whether it fits Poisson distribution (p=0.3626)
ks.test(frequency$Freq, ppois,lambda = mean(frequency$Freq)) 

ES <- lambda * (u + (scale / (1 - shape)))
ES
(EX <- ES/lambda)

# extract exceedance
exceedances <- subset(loss, loss$billion > u)
frequency <- as.data.frame(table(exceedances$Year))

# estimation of lambda
pois <-  fitdistr(frequency$Freq, "Poisson")
lambda <- pois$estimate
lambda

# test whether it fits Poisson distribution (p=0.3626)
ks.test(frequency$Freq, ppois,lambda = mean(frequency$Freq)) 

# Calculation of E(S): expected annual economic loss 

ES <- lambda * (u + (scale / (1 - shape)))
ES
(EX <- ES/lambda)

# VaR and CVaR 

n <- length(loss$billion) # number of all cases
nu <- length(exceedances$billion) # number of exceedances
(n-nu)/n
q <- c(0.99,0.975,0.95,0.9,0.85) # confidence levels

VaR <- (u + (scale / shape)*((((n/nu)*(1-q))^(-shape)) - 1 ))
VaR

# VaR actual spillover rate
E_rate0 <- c()
for (i in 1:length(q)){
  excess <- subset(loss$billion,loss$billion>VaR[i])
  ni <- length(excess)
  E_rate0 <- c(E_rate0,ni/n)
}

CVaR <- (VaR / (1-shape)) + (scale - shape * location) / (1-shape)
CVaR

# CVaR actual spillover rate
E_rate1 <- c()
for (i in 1:length(q)){
  excess1 <- subset(loss$billion, loss$billion > CVaR[i])
  ni1 <- length(excess1)
  E_rate1 <- c(E_rate1,ni1/n)
}

var <- data.frame("confidence level" = percent(q, 0.01), "VaR" = VaR, "significance level" = percent(1-q, 0.01), "actual spillover" = percent(E_rate0, 0.01))
var

cvar <- data.frame("confidence" = percent(q, 0.01), "CVaR" = CVaR, "significance level" = percent(1-q,0.01), "actual spillover" = percent(E_rate1, 0.01))
cvar
