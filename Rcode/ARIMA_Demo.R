rm(list = ls())
cat("\014")
Sys.setlocale(category = "LC_ALL", locale = "english")


# Stationarity ------------------------------------------------------------
# White Noise
e <- rnorm(1000, mean = 0, sd = 1)
plot(e,type='l')
hist(e,20)
acf(e)

# Deterministic trend process
TT    <- 100
alpha <- 0
beta  <- 0.5
e     <- rnorm(TT, mean=0, sd=1)
trend <- 1:TT
x     <- alpha+beta*trend+e
plot(x,type='l',xlab="",ylab="")
grid()
acf(x)

# Change in Variance
NN <- 100
TT <- 100
X  <- matrix(NA,TT,NN)

for (i in 1:NN){
  e.part1 <- rnorm(TT/2, mean=0, sd=1)
  e.part2 <- rnorm(TT/2, mean=0, sd=2)
  X[,i] <-  c(e.part1, e.part2)
}

matplot(X,type='l',xlab="",ylab="")
grid()


# AR ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# AR(1)
ar1 <- arima.sim(list(order=c(1,0,0),ar=0.5),n=1000)
plot(ar1,type='l')
acf(ar1)
pacf(ar1)

# AR(2)
ar2 <- arima.sim(list(order=c(2,0,0),ar=c(-0.2,0.35)),n=1000)
plot(ar2,type='l')
acf(ar2)
pacf(ar2)

# MA ----------------------------------------------------------------------------------------------------------------------------------------------------------------
ma1 <- arima.sim(list(order=c(0,0,1),ma=0.5),n=10000)
plot(ma1,type='l')
acf(ma1)
pacf(ma1)

ma2 <- arima.sim(list(order=c(0,0,2),ma=c(0.3,0.5)),n=10000)
plot(ma2,type='l')
acf(ma2)
pacf(ma2)


# ARMA --------------------------------------------------------------------------------------------------------------------------------------------------------------
phi <- 0.9
theta <- 0.5
my.model <- list(order=c(1,0,1),ar=phi,ma=theta)
arma11 <- arima.sim(my.model,n=10000)
plot(arma11,type='l')
acf(arma11)
pacf(arma11)

# Fitting an ARMA ---------------------------------------------------------------------------------------------------------------------------------------------------
library(readxl)
library(zoo)
library(tseries)
library(urca)
library(lmtest)

# Load Data
Raw.Data <- read_excel("UK GDP Quarterly.xls")
Year  <- as.yearqtr(Raw.Data$TIME)
Year1 <- Year[-1]
GDP <- Raw.Data$Value
LGDP <- log(GDP)
DLGDP <- diff(LGDP)

# Data Visualisation
plot(Year, GDP, type='l')
grid()

plot(Year1, DLGDP, type='l',ylab='Growth Rate of GDP', xlab='Year')
grid()

hist(DLGDP,100, main="", xlab="")

boxplot(DLGDP,horizontal=TRUE)

acf(DLGDP, 24, main="")
grid()

pacf(DLGDP, 24, main="")
grid()

# Model Selection ARMA(p,q) by AIC
pLag <- 0:5
qLag <- 0:5
np   <- length(pLag)
nq   <- length(qLag)
IC   <- matrix(NA, np, nq)

for (i in 1:np){
  for (j in 1:nq){
    p <- pLag[i]
    q <- qLag[j]
    
    ifit    <- arima(DLGDP, order=c(p,0,q))
    IC[i,j] <- AIC(ifit)
  }
}

write.table(IC, "clipboard", sep="\t", row.names=FALSE)

idx <- which(IC == min(IC), arr.ind = TRUE)
best.p=idx[1]-1
best.q=idx[2]-1

# Model Diagnostics
fit.best <- arima(DLGDP, order=c(best.p,0,best.q))
coeftest(fit.best)

# Overfitting
overfit <- arima(DLGDP, order=c(3,0,4))
coeftest(overfit)

# Time Series Plot of Residuals 
my.resid <- resid(fit.best)
plot(Year1, my.resid, type='l', xlab="Year", ylab="", main="Residuals")
grid()

# ACF Plot of Residuals
acf(my.resid,50, main="")
grid()

# Test of autocorrelation on residuals by Ljung-Box Test
Box.test(my.resid,type="Ljung-Box",lag=10, fitdf=best.p+best.q)
Box.test(my.resid,type="Ljung-Box",lag=20, fitdf=best.p+best.q)
Box.test(my.resid,type="Ljung-Box",lag=30, fitdf=best.p+best.q)

# Normality Test
qqnorm(my.resid)
qqline(my.resid)

jarque.bera.test(my.resid)


# Seasonal ARIMA ----------------------------------------------------------------------------------------------------------------------------------------------------
library(forecast)
model <- Arima(ts(rnorm(100),freq=4), order=c(1,1,1), seasonal=c(1,1,1),
               fixed=c(phi=0.5, theta=-0.4, Phi=0.3, Theta=-0.2))
sarima <- simulate(model, nsim=1000)
plot(sarima,type='l')
acf(sarima)
pacf(sarima)

