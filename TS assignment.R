library(forecast)
library(ggplot2)
library(gridExtra)

# Part 1
# Q1
n = 1000    # length of time points

phi1 = -0.022
phi2 = -0.405
phi3 = 0.493
psi1 = 0.550
psi2 = 0.276
psi3 = 0.490
psi4 = 0.682
psi5 = 0.222
psi6 = -0.232

ar3 = c(phi1, phi2, phi3)   # parameter vector for autoregressive process
ma6 = c(psi1, psi2, psi3, psi4, psi5, psi6)   # parameter vector for moving average process

ar3ma6_model = list(ar = ar3, ma = ma6)

ar3ma6_series = arima.sim(model = ar3ma6_model, n = n)  # simulation of time series

# plot time series
autoplot(ar3ma6_series, main="Simulated Time Series of ARMA(3,6)", xlab="Time", ylab="y")

# calculate ACF and PACF
ar3ma6_acf = acf(ar3ma6_series, lag.max = 30, plot = FALSE)
ar3ma6_pacf = pacf(ar3ma6_series, lag.max = 30, plot=FALSE)

# plot ACF and PACF
par(mfrow = c(2, 1))
par(mar = c(5,5,3,1))
plot(ar3ma6_acf, ylim=c(-0.35,1), main="")
title("Autocorrelation Function of Original Time Series", line = 0.5)
plot(ar3ma6_pacf, ylim=c(-0.3,1), main="")
title("Partial Autocorrelation Function of Original Time Series", line = 0.5)

# Q2
z1 <- polyroot(c(1,-phi1,-phi2,-phi3))
Mod(z1) > 1    # causal stationary if TRUE

z2 <- polyroot(c(1,psi1,psi2,psi3,psi4,psi5,psi6))
Mod(z2) > 1    # invertibility if TRUE

# Q3
data1 <- read.delim("data1.txt", sep=" ")[,1]    # import data
data1 = ts(data1)   # transform data to time series

# Fit model with fixed p,d,q values
arima_model <-arima(data1, order = c(0,0,3))

# estimated coefficients and innovation variance
coef_1 <- arima_model$coef; coef_1
sigma_sq_eps_1 <- arima_model$sigma2; sigma_sq_eps_1

# simulate for the next 1000 values
model <- list(order = c(0,0,3), ma = c(coef_1[1],coef_1[2],coef_1[3])) #assign coefficients
generated_data <- arima.sim(model, n, sd=sqrt(sigma_sq_eps_1))
plot(seq(0,length(data1)-1),data1,ylim=c(-4500,3000),ylab="y",xlab="Time",
     type="l",main="Original Time Series")
plot(seq(0,length(generated_data)-1),generated_data,type="l",ylim=c(-1000,1000),
     main="Next 1000 Values in Realisation",ylab="y",xlab="Time")

par(mfrow=c(1,2))
par(mar = c(5,4,3,1))
modelacf <- acf(data1, main="")
title("ACF of Original Time Series", line = 1)
genacf <- acf(generated_data, main="")
title("ACF of Simulated Time Series", line = 1)


# Part 2
# Q1
data2 <- read.delim("data2.txt", sep=" ")[,1]   # import data
data2 <- ts(data2)    # transform data to time series

data2diff1 <- data2[-1] - data2[-length(data2)]     # first differences
data2diff2 <- data2diff1[-1] - data2diff1[-length(data2diff1)]    # second differences

# plot original time series, first and second differences
par(mfrow=c(3,3))
par(mar = c(4,4,3,1))
plot(data2, type="l", ylab="y", main="Original Time Series")
plot(data2diff1, type="l", xlab="Time", ylab="y", main="First Difference")
plot(data2diff2, type="l", xlab="Time", ylab="y", main="Second Difference")

# Q2
# plot ACF for original time series, first and second differences
data2acf = acf(data2, lag.max = 20, plot = TRUE, main="")
title("ACF of Original Time Series", line = 1)
data2diff1acf = acf(data2diff1, lag.max = 20, plot = TRUE, main="")
title("ACF of First Difference", line = 1)
data2diff2acf = acf(data2diff2, lag.max = 20, plot = TRUE, main="")
title("ACF of Second Difference", line = 1)

# plot PACF for original time series, first and second differences
data2pacf = pacf(data2, lag.max = 20, plot = TRUE, main="")
title("PACF of Original Time Series", line = 1)
data2diff1pacf = pacf(data2diff1, lag.max = 20, plot = TRUE, main="")
title("PACF of First Difference", line = 1)
data2diff2pacf = pacf(data2diff2, lag.max = 20, plot = TRUE, main="")
title("PACF of Second Difference", line = 1)

AICmat <- matrix(rep(0,16),nrow=4,ncol=4)
colnames(AICmat) <- c("q=0","q=1","q=2","q=3")
rownames(AICmat) <- c("p=0","p=1","p=2","p=3")
BICmat <- matrix(rep(0,16),nrow=4,ncol=4)
colnames(BICmat) <- c("q=0","q=1","q=2","q=3")
rownames(BICmat) <- c("p=0","p=1","p=2","p=3")

for (p in 1:4){
  for (q in 1:4){
    arima_test <-arima(data2, order = c(p-1,2,q-1))
    AICmat[p,q] <- log(arima_test$sigma2) + 2*((p-1)+(q-1))/n
    BICmat[p,q] <- log(arima_test$sigma2) + log(n)*((p-1)+(q-1))/n
  }
}

AICmat
BICmat

# Q3
# Fit model with selected p,d,q values
arima_model2 <-arima(data2, order = c(2,2,1))
arima_model2$coef     # estimated coefficients
arima_model2$sigma2   # estimated innovation variance

# Q4
# 95% forecast interval for next 500 values
forecasted = forecast(arima_model2, h = 500, level=95)
autoplot(forecasted,legend=TRUE)
