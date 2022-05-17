# HW 6 Q3
# MA(1) Model

## PART I:

library(stats)
ma_model <- arima.sim(model = list(ma=c(-0.6),
                                    sd=1), n=48)

# recursive function to get e_t
e_t <- function(t, ma, theta){
  if (t == 0){
    return(0)
  }
  return(ma[t] - e_t(t-1, ma, theta)*theta)
}


# get sigma^2 for log likelihood function 
sigma_2 <- acf(ma_model, type = "covariance", plot = FALSE)$acf[1]

# creating log likelihood function 
log_like <- function(theta, sigma_2){
  T = 48
  # get sum of squares of e_ts
  ss_e = 0
  for (t in 1:48){
    ss_e = ss_e + e_t(t, ma_model, theta)^2
  }
  # return negative log likelihood as the optime() maximizes 
  return(-1*(-T/2)*(log(2*pi*sigma_2)) - (ss_e / 2*sigma_2))
}

# creating single parameter function to plu into optim 
fr <- function(x){
  return (log_like(x,sig_2))
}

# use optim to estimate theta
ma <- arima.sim(model = list(ma=c(0.8) , sd=1), n=48)
sig_2 = acf(ma,type='covariance', plot=FALSE)$acf[1]
res <- optim(-1,fr, lower= -0.9, upper = 0.9, method='Brent')

## PART II & III: 
# run part 1 many times with new simulated series and collect results
results <- c()
for (i in 1:100){
  ma <- arima.sim(model = list(ma=c(0.8) , sd=1), n=48)
  sig_2 = acf(ma,type='covariance', plot=FALSE)$acf[1]
  res <- optim(par=0.6,fr,upper=0.9,lower=-0.9,method='Brent')$par[1]
  results <- append(results,res)
}
# view sampling distribution
results
hist(results)

# PART IV:
# calculate variance 
var(results)

# calcualte expected variance using the true parametre and n=48
(1-(-0.6))/48



