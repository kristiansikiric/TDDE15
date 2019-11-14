
###### Functions #######
posteriorGP = function(X,y,XStar,hyperParam,sigmaNoise){
  n = length((X))
  K = squeredExpKernel(X,X,hyperParam[1],hyperParam[2])
  L = t(chol(K+(sigmaNoise^2)*diag(n)))
  alpha = solve(t(L),solve(L,y))
  kstar = squeredExpKernel(X,XStar,hyperParam[1], hyperParam[2])
  fbarstar = t(kstar)%*%alpha
  
  v = solve(L, kstar)
  V.f = squeredExpKernel(XStar,XStar,hyperParam[1], hyperParam[2])-t(v)%*%v
  
  
  return(list("Mean" = fbarstar, "Variance" = diag(V.f)))
}

squeredExpKernel= function(x1,x2,sigmaF,l){
  n1 = length(x1)
  n2 = length(x2)
  k = matrix(NA,n1,n2)
  for (i in 1:n2){
    k[,i] = sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(k)
}

GPkernel = function(sigmaf,ell){
  squeredExpKernel= function(x,y){
    n1 = length(x)
    n2 = length(y)
    k = matrix(NA,n1,n2)
    for (i in 1:n2){
      k[,i] = sigmaf^2*exp(-0.5*( (x-y[i])/ell)^2 )
    }
    return(k)
  }
  class(squeredExpKernel) <- "kernel"
  return(squeredExpKernel)
}

plotGP = function(mean,variance,grid,x,y){
  plot(grid,mean,ylim = c(min(mean-1.96*sqrt(variance))
                          ,max(mean+1.96*sqrt(variance))),
       type = "l")
  lines(grid,
        mean+1.96*sqrt(variance), 
        col = rgb(0, 0, 0, 0.3))
  lines(grid,
        mean-1.96*sqrt(variance), 
        col = rgb(0, 0, 0, 0.3))
  lines(grid,
        mean+1.96*sqrt(variance)+sigmaNoise^2,
        col = rgb(0, 0, 0, 1))
  lines(grid,
        mean-1.96*sqrt(variance)+sigmaNoise^2, 
        col = rgb(0, 0, 0, 1))
  points(x,y)
}

GP.periodic.kernel = function(sigmaf, l1,l2,d){
  periodicKernel = function(x1,x2){
    return(  (sigmaF^2)*exp(-(2*sin(pi*abs(x1-x2)/d)^2)/l1^2)*
               exp(-0.5*(abs(x1-x2)^2)/l2^2))
  }
  class(periodicKernel) <- "kernel"
  return(periodicKernel)
}




###########################
###### Assignment 1 #######
###########################
###### 1.2 ######
sigmaF = 1
l = 0.3
hyperParam = c(sigmaF,l)
sigmaNoise = 0.1
X = 0.4
y = 0.719

XStar = seq(-1,1,by=0.01)
GP = posteriorGP(X,y,XStar,hyperParam,sigmaNoise)
plotGP(GP$Mean,GP$Variance,XStar,X,y)

##### 1.3 ######
sigmaF = 1
l = 0.3
hyperParam = c(sigmaF,l)
sigmaNoise = 0.1
X = c(0.4,-0.6)
y = c(0.719,-0.044)

XStar = seq(-1,1,by=0.01)
GP = posteriorGP(X,y,XStar,hyperParam,sigmaNoise)
plotGP(GP$Mean,GP$Variance,XStar,X,y)

##### 1.4 ######
sigmaF = 1
l = 0.3
hyperParam = c(sigmaF,l)
sigmaNoise = 0.1
X = c(-1.0, -0.6, -0.2, 0.4, 0.8)
y = c(0.768, -0.044, -0.940, 0.719, -0.664)

XStar = seq(-1,1,by=0.01)
GP = posteriorGP(X,y,XStar,hyperParam,sigmaNoise)
plotGP(GP$Mean,GP$Variance,XStar,X,y)

##### 1.5 ######
sigmaF = 1
l = 1
hyperParam = c(sigmaF,l)
sigmaNoise = 0.1
X = c(-1.0, -0.6, -0.2, 0.4, 0.8)
y = c(0.768, -0.044, -0.940, 0.719, -0.664)

XStar = seq(-1,1,by=0.01)
GP = posteriorGP(X,y,XStar,hyperParam,sigmaNoise)
plotGP(GP$Mean,GP$Variance,XStar,X,y)

##########################
###### Assignment 2 ######
##########################

###### Innit ######
library(kernlab)
data = read.csv(
  "https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv",
  header=TRUE,
  sep=";")

time = seq(1,2190,by=5)
day = rep(seq(1,365,by=5),6)
###### 2.1 ######
X = c(1,3,4)
XStar = c(2,3,4)
x = 1
xprim = 2

gpK = GPkernel(sigmaf = 20,ell = 0.2)
gpK(x = 1, y = 2)
K = kernelMatrix(gpK,x=X,y=XStar) 

zGrid <- seq(0.01, 1, by = 0.01)
count = 0
covs = rep(0,length(zGrid))
for (z in zGrid){
  count = count + 1
  covs[count] <- GPkernel(sigmaf = 20, ell = 0.2)(0,z)
}
plot(zGrid, covs, type = "l", xlab = "ell")
###### 2.2 ######
sigmaF = 20
l = 0.2
temp = data[time,]$temp
lm.fit = lm(temp ~ time + I(time^2))
sigmaNoise = sd(lm.fit$residuals)

GPfit = gausspr(x = time,
                y = temp,
                kernel = GPkernel,
                kpar = list(sigmaf = sigmaF, ell = l),
                var = sigmaNoise^2)
meanPred = predict(GPfit,time) #Predict on train data
plot(time,temp)
lines(time,meanPred, col = "green")

###### 2.3 ######
var = posteriorGP(X = scale(time),
                  y = scale(temp),
                  XStar = scale(time),
                  hyperParam = c(sigmaF,l),
                  sigmaNoise = sigmaNoise)$Variance
lines(time,meanPred+1.96*sqrt(var), col = "black")
lines(time,meanPred-1.96*sqrt(var), col = "black")

###### 2.4 ######
sigmaF = 20
l = 0.2
temp = data[time,]$temp
lm.fit = lm(temp ~ day + I(day^2))
sigmaNoise = sd(lm.fit$residuals)

GPfit = gausspr(x = day,
                y = temp,
                kernel = GPkernel,
                kpar = list(sigmaf = sigmaF, ell = l),
                var = sigmaNoise^2)
meanPred = predict(GPfit,day) #Predict on train data
lines(time,meanPred, col = "red")

###### 2.5 ######
sigmaf = 10
l1 = 1
l2 = 10
d = 365/sd(time)
temp = data[time,]$temp
lm.fit = lm(temp ~ time + I(time^2))
sigmaNoise = sd(lm.fit$residuals)

GPfit = gausspr(x = time,
                y = temp,
                kernel = GP.periodic.kernel,
                kpar = list(sigmaf = sigmaf, l1 = l1, l2 = l2, d=d),
                var = sigmaNoise^2)
meanPred = predict(GPfit,time) #Predict on train data
#plot(time,temp)
lines(time,meanPred, col = "purple")


##########################
###### Assignment 3 ######
##########################
###### Innit ######
library(AtmRay)
data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv",
                 header=FALSE,
                 sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
set.seed(111); SelectTraining <- sample(1:dim(data)[1], size = 1000,
                         replace = FALSE)
train = data[SelectTraining,]
test = data[-SelectTraining,]
###### 3.1 ######
GP.fit = gausspr(fraud ~ varWave+skewWave, data = train)

x1 <- seq(min(train$varWave),max(train$varWave),length=100)
x2 <- seq(min(train$skewWave),max(train$skewWave),length=100)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))
gridPoints <- data.frame(gridPoints)
names(gridPoints) <- names(data)[1:2]
probPreds <- predict(GP.fit, gridPoints, type="probabilities")

contour(x1,x2,matrix(probPreds[,2],100,byrow = TRUE), 20)

points(train[train$fraud == 1,"varWave"],
       train[train$fraud == 1,"skewWave"],
       col = "blue")

points(train[train$fraud == 0,"varWave"],
       train[train$fraud == 0,"skewWave"],
       col = "red")

con.mat.train = table(predict(GP.fit,train),train$fraud)
acc.train = sum(diag(con.mat.train)/sum(con.mat.train))

###### 3.2 ######
con.mat.test = table(predict(GP.fit,test),test$fraud)
acc.test = sum(diag(con.mat.test)/sum(con.mat.test))

###### 3.3 ######
GP.fit.all = gausspr(fraud ~ ., data = train)
con.mat.test.all = table(predict(GP.fit.all,test),test$fraud)
acc.test.all = sum(diag(con.mat.test.all)/sum(con.mat.test.all))
