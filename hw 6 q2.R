# HW6 Q2
## LOADING PACKAGES AND SETTING SEED 

library(stats)
library(ltsa)
library(dplyr)
library(ggplot2)
set.seed(4857)

## MODEL CREATION 

#models with WN ~ N(0,1)
model_1_1 <- rnorm(500)
model_2_1 <- arima.sim(model = list(
  ar = c(0.3), ma=c(0), sd=1), n=500)
model_3_1 <- arima.sim(model = list(
  ar = c(0.6,-0.2), ma=c(-0.2), sd=1), n=500)

#models with WN ~ t(4)
model_1_2 <- rt(500, 4)
model_2_2 <- arima.sim(n=500, model = list(
  ar = c(0.3), ma=c(0)), rand.gen=function(n,...) rt(n=500,df=4))
model_3_2 <- arima.sim(n=500, model = list(
  ar = c(0.6,-0.2), ma=c(-0.2)), rand.gen = function(n,...) rt(n=500,df=4))

# models with WN ~ (X^2 - 2)
model_1_3 <- (rchisq(500,2)) - 2
model_2_3 <- arima.sim(n=500, model = list(
  ar = c(0.3), ma=c(0)), rand.gen=function(n,...) rchisq(500,2) - 2)
model_3_3 <- arima.sim(n=500, model = list(
  ar = c(0.6,-0.2), ma=c(-0.2)), rand.gen=function(n,...) rchisq(500,2) - 2)

## THEO ACFV FOR EACH MODEL 

# models with WN ~ N(0,1)
theo_1_1 <- tacvfARMA(phi = c(0), theta = c(0), maxLag=20)
theo_2_1 <-tacvfARMA(phi = c(0.3), theta=c(0), maxLag=20)
theo_3_1 <- tacvfARMA(phi = c(0.6,-0.2), theta = c(-0.2),
                       maxLag = 20)

# models with WN ~ t(4)
theo_1_2 <-  tacvfARMA(phi = c(0), theta = c(0), maxLag=20, sigma2 = 2)
theo_2_2 <- tacvfARMA(phi = c(0.3), theta=c(0), maxLag=20, sigma2 = 2)
theo_3_2 <- tacvfARMA(phi = c(0.6,-0.2), theta = c(-0.2),
                      maxLag = 20, sigma2 = 2)
# models with WN ~ (X^2 - 2)
theo_1_3 <-  tacvfARMA(phi = c(0), theta = c(0), maxLag=20, sigma2 = 4)
theo_2_3 <- tacvfARMA(phi = c(0.3), theta=c(0), maxLag=20, sigma2 = 4)
theo_3_3 <- tacvfARMA(phi = c(0.6,-0.2), theta = c(-0.2),
                      maxLag = 20, sigma2 = 4)

acvf_i <- acf(model_1_1, type = "covariance", plot=FALSE, lag.max=20)
acvf_i$acf[21] # lag 1 starts at 1
## SIMULATOR 

sim <- function(model, theo_acvf){
  performance <- c()
  # for loop to gets lags 1-20
  for(h in 2:21) {
    theo_acvf_i <- theo_acvf[h]
    count = 0
    left <- abs(theo_acvf_i - (2 / sqrt(500)))
    right <- abs(theo_acvf_i + (2 / sqrt(500)))
    # for loop for 1000 iterations of testing if sample ACFV is in interval
    for (i in 1:1000) {
      acvf_i <- acf(model, type = "covariance", plot=FALSE)
      if (between(acvf_i$acf[h],left, right) == TRUE) {
        count = count+1
      }
    }
    performance <- append(performance, count / 1000)
  }
  return(performance)
}

## RESULTS OF SIMULATOR 

res_1_1 <- sim(model_1_1, theo_1_1)
res_2_1 <- sim(model_2_1, theo_2_1)
res_3_1 <- sim(model_3_1, theo_3_1)
res_1_2 <- sim(model_1_2, theo_1_2)
res_2_2 <- sim(model_2_2, theo_2_2)
res_3_2 <- sim(model_3_2, theo_3_2)
res_1_3 <- sim(model_1_3, theo_1_3)
res_2_3 <- sim(model_2_3, theo_2_3)
res_3_3 <- sim(model_3_3, theo_3_3)
res_1_1
res_2_1

## BOXPLOTS
# sequence of boxplots, one for each lag
# should be decay in the medium of that boxplot
boxplot(res_1_1, res_2_1, res_3_1, res_1_2, res_2_2, res_3_2,
        res_3_1, res_3_2, res_3_3)
boxplot(res_1_1)

df <- data.frame(matrix(c(res_1_1, res_2_1, res_3_1,
                         res_1_2, res_2_2, res_3_2,
                         res_1_3,res_2_3,res_3_3),
                nrow=9, ncol=20),
  row.names = c("model_1_1", "model_2_1", "model_3_1",
  "model_1_2", "model_2_2", "model_3_2",
  "model_1_3", "model_2_3", "model_3_3"))
view(df)

for(h in 2:21) {
  theo_acvf_i <- theo_3_1[h]
  print(theo_acvf_i)
}
















## can't get specific lag to be outputted 
model_3_1
theo_3_1[2]
acvf_3_1 <- acf(model_3_1, type = "covariance", plot=FALSE, lag.max=20)
left <- abs(theo_3_1[1] - (2 / sqrt(500)))
right <- abs(theo_3_1[1] + (2 / sqrt(500)))
between(acvf_3_1$acf[1], left, right)

