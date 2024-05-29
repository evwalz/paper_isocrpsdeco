# Functions:

library(scoringRules)
library(extraDistr)
library(bmixture)
library(EnvStats)
library(evd)

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




bounds_norm_mix <- function(y,idr_preds, h, epsilon, delta){

  a <- min(y)
  b <- max(y)
  n <- length(y)
  tmp  = rep(NA,n)
  for ( i in 1:n){
    mu <- idr_preds[[i]]$points
    weights <- diff(c(0, idr_preds[[i]]$cdf))
    # extraDistr::pmixnorm(grid_vals[k], mean = mean, alpha = weights, sd = rep(h, length(mean)))
    integrand1 = function(x) {return((extraDistr::pmixnorm(x,mean=mu,sd=rep(h, length(mu)), alpha = weights))^2)}
    integrand2 = function(x) {return((1-extraDistr::pmixnorm(x,mean=mu,sd=rep(h, length(mu)), alpha = weights))^2)}
    tmp[i]= integrate(integrand1,lower=-Inf,upper=a)$value + integrate(integrand2,lower=b,upper=Inf)$value
  }
  if (mean(tmp) <= epsilon){
    return(c(a, b))
  } else {
    while (mean(tmp)>epsilon){
      a = a - delta
      b = b + delta
      for ( i in 1:n){
        mu <- idr_preds[[i]]$points
        weights <- diff(c(0, idr_preds[[i]]$cdf))
        integrand1 = function(x) {return((extraDistr::pmixnorm(x,mean=mu,sd=rep(h, length(mu)), alpha = weights))^2)}
        integrand2 = function(x) {return((1-extraDistr::pmixnorm(x,mean=mu,sd=rep(h, length(mu)), alpha = weights))^2)}
        tmp[i]=integrate(integrand1,lower=-Inf,upper=a)$value +integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(c(a, b))
  }
}

mix_pt <- function(x_mean, h, weights, df){
  return(sum(weights * pt(x_mean / h, df = df)))
}

mix_pt2 <- function(x, mean, h, weights, df){
  n <- length(mean)
  tmp <- 0
  for (j in 1:n){
    tmp <- tmp + weights[j]*pt((x-mean[j]) / h, df = df)
  }
  return(tmp)
}

bounds_t_mix <- function(y,idr_preds, h, df, epsilon, delta){
  #average  <- mean(y)
  #print(average)
  #t <- max(max(y)-average,average-min(y))
  #print(t)
  a <- min(y)
  b <- max(y)
  n <- length(y)
  tmp  = rep(NA,n)
  for ( i in 1:n){
    mu <- idr_preds[[i]]$points
    weights <- diff(c(0, idr_preds[[i]]$cdf))
    
    # extraDistr::pmixnorm(grid_vals[k], mean = mean, alpha = weights, sd = rep(h, length(mean)))
    integrand1 = function(x) {return((mix_pt2(x,mu,h, weights, df))^2)}
    integrand2 = function(x) {return((1-mix_pt2(x,mu,h, weights, df))^2)}
    tmp[i]= integrate(integrand1,lower=-Inf,upper=a)$value + integrate(integrand2,lower=b,upper=Inf)$value
  }
  if (mean(tmp) <= epsilon){
    return(c(a, b))
  } else {
    while (mean(tmp)>epsilon){
      a = a - delta
      b = b + delta
      for ( i in 1:n){
        mu <- idr_preds[[i]]$points
        weights <- diff(c(0, idr_preds[[i]]$cdf))
        integrand1 = function(x) {return((mix_pt2(x,mu,h, weights, df))^2)}
        integrand2 = function(x) {return((1-mix_pt2(x,mu,h, weights, df))^2)}
        tmp[i]=integrate(integrand1,lower=-Inf,upper=a)$value +integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(c(a, b))
  }
}



bounds_norm_mix_cp <- function(y,ens, h, epsilon, delta){
  a <- min(y)
  b <- max(y)
  n <- length(y)
  tmp  = rep(NA,n)
  for ( i in 1:n){
    mean <- sort(unique(ens[i,]))
    colnames(mean) <- NULL
    mean <- unlist(mean)
    ecdf_fun <- ecdf(ens[i,])
    ecdf_vals <- ecdf_fun(mean)
    weights <-  diff(c(0, ecdf_vals))
    sd <- rep(h, length(mean))
    #mu <- mean[i,]
    #sd <- sig[i,]
    #print(i)
    #weights <- diff(c(0, idr_preds[[i]]$cdf))
    # extraDistr::pmixnorm(grid_vals[k], mean = mean, alpha = weights, sd = rep(h, length(mean)))
    integrand1 = function(x) {return((extraDistr::pmixnorm(x,mean=mean,sd=sd, alpha = weights))^2)}
    integrand2 = function(x) {return((1-extraDistr::pmixnorm(x,mean=mean,sd=sd, alpha = weights))^2)}
    tmp[i]= integrate(integrand1,lower=-Inf,upper=a)$value + integrate(integrand2,lower=b,upper=Inf)$value
  }
  if (mean(tmp) <= epsilon){
    return(c(a, b))
  } else {
    while (mean(tmp)>epsilon){
      a = a - delta
      b = b + delta
      for ( i in 1:n){
        mean <- sort(unique(ens[i,]))
        colnames(mean) <- NULL
        mean <- unlist(mean)
        ecdf_fun <- ecdf(ens[i,])
        ecdf_vals <- ecdf_fun(mean)
        weights <-  diff(c(0, ecdf_vals))
        sd <- rep(h, length(mean))
        integrand1 = function(x) {return((extraDistr::pmixnorm(x,mean=mean,sd=sd, alpha = weights))^2)}
        integrand2 = function(x) {return((1-extraDistr::pmixnorm(x,mean=mean,sd=sd, alpha = weights))^2)}
        tmp[i]=integrate(integrand1,lower=-Inf,upper=a)$value +integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(c(a, b))
  }
}

# y_test,ens, h, df, epsilon, delta

bounds_t_mix_cp <- function(y,ens, h, df, epsilon, delta){
  #average  <- mean(y)
  #print(average)
  #t <- max(max(y)-average,average-min(y))
  #print(t)
  a <- min(y)
  b <- max(y)
  n <- length(y)
  
  tmp  = rep(NA,n)
  for ( i in 1:n){
    mean <- sort(unique(ens[i,]))
    colnames(mean) <- NULL
    mean <- unlist(mean)
    ecdf_fun <- ecdf(ens[i,])
    ecdf_vals <- ecdf_fun(mean)
    weights <-  diff(c(0, ecdf_vals))
    #sd <- rep(h, length(mean))
    
    # extraDistr::pmixnorm(grid_vals[k], mean = mean, alpha = weights, sd = rep(h, length(mean)))
    integrand1 = function(x) {return((mix_pt2(x,mean,h, weights, df))^2)}
    integrand2 = function(x) {return((1-mix_pt2(x,mean,h, weights, df))^2)}
    tmp[i]= integrate(integrand1,lower=-Inf,upper=a)$value + integrate(integrand2,lower=b,upper= Inf)$value
  }
  
  if (mean(tmp) <= epsilon){
    return(c(a, b))
  } else {
    while (mean(tmp)>epsilon){
      a = a - delta
      b = b + delta
      for ( i in 1:n){
        mean <- sort(unique(ens[i,]))
        colnames(mean) <- NULL
        mean <- unlist(mean)
        ecdf_fun <- ecdf(ens[i,])
        ecdf_vals <- ecdf_fun(mean)
        weights <-  diff(c(0, ecdf_vals))
        integrand1 = function(x) {return((mix_pt2(x,mean,h, weights, df))^2)}
        integrand2 = function(x) {return((1-mix_pt2(x,mean,h, weights, df))^2)}
        tmp[i]=integrate(integrand1,lower=-Inf,upper=a)$value +integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(c(a, b))
  }
}


bounds_norm_mix_mc <- function(y,mean, sig, epsilon, delta){
  #average  <- mean(y)
  #print(average)
  #t <- max(max(y)-average,average-min(y))
  #print(t)
  a <- min(y)
  b <- max(y)
  n <- length(y)
  tmp  = rep(NA,n)
  for ( i in 1:n){
    mu <- mean[i,]
    sd <- sig[i,]
    #print(i)
    #weights <- diff(c(0, idr_preds[[i]]$cdf))
    # extraDistr::pmixnorm(grid_vals[k], mean = mean, alpha = weights, sd = rep(h, length(mean)))
    integrand1 = function(x) {return((extraDistr::pmixnorm(x,mean=mu,sd=sd, alpha = rep(1, length(mu))/length(mu)))^2)}
    integrand2 = function(x) {return((1-extraDistr::pmixnorm(x,mean=mu,sd=sd, alpha = rep(1, length(mu))/ length(mu)))^2)}
    tmp[i]= integrate(integrand1,lower=-Inf,upper=a)$value + integrate(integrand2,lower=b,upper=Inf)$value
  }
  if (mean(tmp) <= epsilon){
    return(c(a, b))
  } else {
    while (mean(tmp)>epsilon){
      print(mean(tmp))
      a = a - delta
      b = b + delta
      for ( i in 1:n){
        mu <- mean[i,]
        sd <- sig[i,]
        integrand1 = function(x) {return((extraDistr::pmixnorm(x,mean=mu,sd=sd, alpha = rep(1, length(mu))/ length(mu)))^2)}
        integrand2 = function(x) {return((1-extraDistr::pmixnorm(x,mean=mu,sd=sd, alpha = rep(1, length(mu))/ length(mu)))^2)}
        tmp[i]=integrate(integrand1,lower=-Inf,upper=a)$value +integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(c(a, b))
  }
}



