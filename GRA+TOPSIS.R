## MCDM: GRA
## Data Input 
setwd("C:/Users/dell/Desktop")
data <- read.csv("data.csv",header=T) #normalization has been done in excel
rownames(data) = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13")
data <- data
referColNum <- 1
distingCoeff <- 0.5 
y = data[referColNum]

##Grey Relational Coefficients Calculation
nr = nrow(y)
data1 = data[-referColNum]
nc = ncol(data1)

##Grey Relational Coefficients Calculation
diff = data1

for (i in 1:nc) {
  diff[i] = abs(y-data1[i])
}
diff <- diff
mi = min(diff)
mx = max(diff)

relCoeff = (mi + distingCoeff*mx) / (diff + distingCoeff*mx)

##GRA with Equal Weight
relDegree = rep(NA,nc) 
for (i in 1:nc) 
{
  relDegree[i] = mean(relCoeff[,i])  
}

##Ranking
data1_order = data1[order(relDegree, decreasing = TRUE)]
relationalDegree = relDegree[order(relDegree, decreasing = TRUE)]
relDes = rep(NA, nc)
data1_names = names(data1_order)
data1_names
for (i in 1:nc) {
  relDes[i] = paste(names(df)[referColNum], data1_names[i], sep = "~")
}
names(relationalDegree) = relDes
relDes = as.data.frame(relDes, row.names = NULL, optional = FALSE)

## Initial Capital Allocation
## ignore Regional Economic Development Indiactors and Normalization
library(readxl)
data0 <- read_excel("C:/Users/dell/Desktop/data0.xlsx")
dataReduced = data0[-(7:13),]
rownames(dataReduced) = c("X1","X2","X3","X4","X5","X6")
maxrow = apply(dataReduced,1,max)
minrow = apply(dataReduced,1,min)
dataNormalized = (dataReduced - minrow)/(maxrow - minrow) # the larger the better

## Grey Relational Coefficients Calculation
referColNum1 = as.data.frame(apply(dataNormalized,1,max))
distingCoeff = 0.5
nr1 = nrow(dataNormalized)
nc1 = ncol(dataNormalized)

diff1 = dataNormalized

for (i in 1:nc1) {
  diff1[i] = abs(referColNum1-dataNormalized[i])
}

diff1 <- diff1
mi = min(diff1)
mx = max(diff1)

relCoeff1 = (mi + distingCoeff*mx) / (diff1 + distingCoeff*mx)    

## Grey Relational Degree Calculation with Equal Weights
relDegree1 = rep(NA,nc1) 
for (i in 1:nc1) 
{
  relDegree1[i] = mean(relCoeff1[,i])  
}

## GRA Ranking without Economic Development Indicators
dataNormalized_order = dataNormalized[order(relDegree1, decreasing = TRUE)]
relationalDegree1 = relDegree1[order(relDegree1, decreasing = TRUE)]
relDes1 = rep(NA, nc1)
dataNormalized_names = names(dataNormalized_order)
dataNormalized_names
for (i in 1:nc1) 
{
  relDes1[i] = paste(names(df)[referColNum1], dataNormalized_names[i], sep = "~")
}
names(relationalDegree1) = relDes1
relDes1 = as.data.frame(relDes1, row.names = NULL, optional = FALSE)


## Entropy Weight Calculation
## Data input
setwd("C:/Users/dell/Desktop")
data <- read.csv("data.csv",header=T) 
#normalization has been done in excel 
#the larger the better and the smaller the better
rownames(data) = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13")
data2 = data[,-1]
##Coordinate translation normalization (prepare for later entropy weight)
data2 = data +0.01
t = rowSums(data1)
data2 = data2/t
rowSums(data2) 

##Compute the entropy
K = 1/log(length(data1))
a = log(data1)
dataprime = data1*a
datanew = as.data.frame(rowSums(dataprime))
entropy = (-1)*K*datanew

##Calculation of the weight of each criterion
difference = 1-entropy
weight = difference/sum(difference)
weight

# GRA with Entropy weight
relDegree[i] = sum(relCoeff[,i]%*% weight) 
}
relDegree

#GRA with Entropy weight ranking
data2_order = data2[order(relDegree, decreasing = TRUE)]
relationalDegree = relDegree[order(relDegree, decreasing = TRUE)]
relDes = rep(NA, nc)
data2_names = names(data2_order)
data2_names
for (i in 1:nc) 
{
  relDes[i] = paste(names(df)[referColNum], data2_names[i], sep = "~")
}
names(relationalDegree) = relDes
relDes = as.data.frame(relDes, row.names = NULL, optional = FALSE)


## MCDM: TOPSIS 
# Normalization
library(readxl)
matrixR <- read_excel("C:/Users/dell/Desktop/data0.xlsx")
matrixR = t(matrixR)
nrowR = nrow(matrixR)
ncolR = ncol(matrixR)
# Normalization by larger the better and 
#smaller the better for the first 6 and the last 7 column
MAX6 <- apply(matrixR[,1:6], 2, max)
MIN7 <- apply(matrixR[,7:ncolR], 2, min)

matrixV <- matrixR

for (i in 1:nrowR)
{
  matrixV[i,1:6] <- matrixR[i,1:6]/MAX6
  matrixV[i,7:ncolR] <- MIN7/matrixR[i,7:ncolR]
}

# Weighed Decision Matrix VW Calculation
matrixW <- matrix(rep(t(weight), nrowR), nrow=nrowR, by=2)

matrixVW <- matrixV * matrixW

## Determine the PIS and NIS
Astar = apply(matrixVW,2,max) #PIS
Aminus = apply(matrixVW,2,min) #NIS

## Distance Calculation
nr = nrow(matrixVW)
nc = ncol(matrixVW)
Sstar = rep(NA,nr)
for (i in 1:nr) 
{
  Sstar[i] = sqrt(sum((matrixVW[i,]-Astar)^2))
}

Sminus = rep(NA,nr)
for (i in 1:nr) 
{
  Sminus[i] = sqrt(sum((matrixVW[i,]-Aminus)^2))
}

## Closeness coefficient Calculation
Cstar = Sminus/(Sstar + Sminus)

## TOPSIS Ranking
ID = seq(1,19)
combine = as.data.frame(cbind(ID,Cstar))
combine <- combine[order(combine$Cstar,decreasing = TRUE),]



