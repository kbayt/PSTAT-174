## PART 1
#set up: AR Model
ar_model <- arima.sim(model = list(ar=c(0.8),
                                    sd=1), n=48)


# Method of Moments Estimate (358)
# part i: method-of-moments
# from method of moments we get the the sample parameter = sample lag 1
sample_acf <- acf(ar_model, type = "correlation", plot=FALSE)
mom_estimate <- sample_acf$acf[2,1,1]

# part ii. Least Squares Estimation 
# from least squares we get that the sample parameter = sample lag 1
# note that this is the same as method of moment
ls_estimate <- sample_acf$acf[2,1,1]

# part iii. Maximum Likelihood Estimate 
# create function for sum within log likelihood 
sum=0
sum_yt <- function(phi){
  for (i in 2:48){
    sum = sum + (ar_model[i]-(phi*ar_model[i-1]))^2
  }
  return(sum)
}

# get sample acvf(0) to estimate variance
sigma_sq <- acf(ar_model, type="covariance", plot=FALSE)$acf[1]
# function for log likelohood 
ll <- function(phi){
  return(((-48/2)*log(2*pi)) - ((48/2)*log(sigma_sq)) -
           ((1/2)*log(1-(phi^2))) -
           (1/(2*sigma_sq))*((sum_yt(phi))+ ((1 - phi^2)*(ar_model[1]^2))))

}
# use optim to ply in innovation for phi and lof likelihood function 
optim(0.5, ll, lower = -1, upper = 1, method = "Brent")
arima(ar_model)
