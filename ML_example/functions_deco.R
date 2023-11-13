# Functions:

library(scoringRules)

crps_ecdf_list <- function(y, p,x) {
  w <- lapply(p, function(x) c(x[1], diff(x)))
  
  crps0 <- function(y, p, w, x) 2 * sum(w * ((y < x) - p + 0.5 * w) * (x - y))
  mapply(crps0, y = y, p, w = w, x = x)
}

crps_ecdf <- function(y, grid_vals, ecdf) {
  
  x <- lapply(seq_len(length(y)), function(i) grid_vals)
  p <- lapply(seq_len(nrow(ecdf)), function(i) ecdf[i,])#as.list( t(ecdf) )
  w <- lapply(p, function(x) c(x[1], diff(x)))
  
  crps0 <- function(y, p, w, x) 2 * sum(w * ((y < x) - p + 0.5 * w) * (x - y))
  mapply(crps0, y = y, p, w = w, x = x)
}

crps_unc <- function(obs){
  obs_ensemble <- matrix(rep(obs, length(obs)), nrow = length(obs), ncol = length(obs), byrow = TRUE)
  return(mean(crps_sample(obs, dat = obs_ensemble)))
}


gaussian_crps <- function(y, mu, sigma){
  return(mean(scoringRules::crps_norm(y, mean = mu, sd = sigma)))
}

bounds_norm_minmax <- function(y, mu, sigma, epsilon, delta){
  average  <- mean(y)
  a <- min(y)
  b <- max(y)
  n <- length(y)
  tmp  = rep(NA,n)
  for ( i in 1:n){
    integrand1 = function(x) {return((pnorm(x,mean=mu[i],sd=sigma[i]))^2)}
    integrand2 = function(x) {return((1-pnorm(x,mean=mu[i],sd=sigma[i]))^2)}
    tmp[i]= integrate(integrand1,lower=-Inf,upper=a)$value + integrate(integrand2,lower=b,upper=Inf)$value
  }
  if (mean(tmp) <= epsilon){
    return(c(a, b))
  } else {
    while (mean(tmp)>epsilon){
      a = a - delta
      b = b + delta
      for ( i in 1:n){
        integrand1 = function(x) {return((pnorm(x,mean=mu[i],sd=sigma[i]))^2)}
        integrand2 = function(x) {return((1-pnorm(x,mean=mu[i],sd=sigma[i]))^2)}
        tmp[i]=integrate(integrand1,lower=-Inf,upper=a)$value +integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(c(a, b))
  }
}


func_crps_normAB <- function(y, mu, sigma, a, b){
  n = length(y)
  tmp  = rep(NA,n)
  for ( i in 1:n){
    integrand1 = function(x) {return((pnorm(x,mean=mu[i],sd=sigma[i]))^2)}
    integrand2 = function(x) {return((pnorm(x,mean=mu[i],sd=sigma[i]) - 1)^2)}
    tmp[i]=integrate(integrand1,lower=a,upper=y[i])$value + integrate(integrand2,lower=y[i],upper=b)$value
  }
  return(mean(tmp))
}

func_crps_normAB2 <- function(y, mu, sigma, a, b){
  return(mean(scoringRules::crps_cnorm(y, location = mu, scale = sigma, lower = a, upper = b)))
}
