# Interval functions for Precip Case Study

#############
### EMOS: ###
##############

# a = 0 fixed
# paramaters of GEV distribution: loc, scale, shape (=-1*c) (all vectors of length from y)
# epsilon and delta for [a, b]
# return b

bounds_pgev_minmax <- function(y, loc,scale, shape,  epsilon, delta){
  a <- 0
  b <- max(y)
  n <- length(y)
  tmp  = rep(NA,n)
  for ( i in 1:n){
    integrand2 = function(x) {return((1-evd::pgev(x,loc=loc[i],scale = scale[i], shape = shape[i]))^2)}
    tmp[i]= integrate(integrand2,lower=b,upper=Inf)$value
  }
  
  if (mean(tmp) <= epsilon){
    return(b)
  } else {
    while (mean(tmp)>epsilon){
      b = b + delta
      for ( i in 1:n){
        integrand2 = function(x) {return((1-evd::pgev(x,loc=loc[i],scale = scale[i], shape =  shape[i]))^2)}
        tmp[i]=integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(b)
  }
}



#############
### BMA: ####
##############

# a = 0 fixed
# gamma distribution with p0 mass at zero
# epsilon and delta for [a, b]
# return b

# weights: vector of length 3
# p0, mean, variance: matrix of dim: (length(y)) x 3


bounds_bma_minmax <- function(y, weights, p0, mean, variance, epsilon, delta){
  a <- 0
  b <- max(y)
  n <- length(y)
  tmp  = rep(NA,n)
  rate <- mean / variance
  for ( i in 1:n){
    #integrand2 = function(x) {return((1-sum(weights*(p0[i,]+(1-p0[i,])*pgamma(x,shape=mean[i,]*rate[i,],rate=rate[i,]))))^2)}
    integrand2 = function(x) {return((1-((weights[1]*(p0[i,1]+(1-p0[i,1])*pgamma(x,shape=mean[i,1]*rate[i,1],rate=rate[i,1]))) + (weights[2]*(p0[i,2]+(1-p0[i,2])*pgamma(x,shape=mean[i,2]*rate[i,2],rate=rate[i,2]))) + (weights[3]*(p0[i,3]+(1-p0[i,3])*pgamma(x,shape=mean[i,3]*rate[i,3],rate=rate[i,3])))))^2)}
    tmp[i]= integrate(integrand2,lower=b,upper=Inf)$value
  }
  
  if (mean(tmp) <= epsilon){
    return(b)
  } else {
    while (mean(tmp)>epsilon){
      b = b + delta
      for ( i in 1:n){
        integrand2 = function(x) {return((1-((weights[1]*(p0[i,1]+(1-p0[i,1])*pgamma(x,shape=mean[i,1]*rate[i,1],rate=rate[i,1]))) + (weights[2]*(p0[i,2]+(1-p0[i,2])*pgamma(x,shape=mean[i,2]*rate[i,2],rate=rate[i,2]))) + (weights[3]*(p0[i,3]+(1-p0[i,3])*pgamma(x,shape=mean[i,3]*rate[i,3],rate=rate[i,3])))))^2)}
        tmp[i]=integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(b)
  }
}



#############
### HCLR: ###
##############

# a = 0 fixed
# logistic distribution with loc and scale
# epsilon and delta for [a, b]
# return b

boundshclr_minmax <- function(y, loc, scale, epsilon, delta){
  a <- 0
  b <- max(y)
  n <- length(y)
  tmp  = rep(NA,n)
  for ( i in 1:n){
    integrand2 = function(x) {return((1-plogis(x,location=loc[i],scale = scale[i]))^2)}
    tmp[i]= integrate(integrand2,lower=b,upper=Inf)$value
  }
  
  if (mean(tmp) <= epsilon){
    return(b)
  } else {
    while (mean(tmp)>epsilon){
      b = b + delta
      for ( i in 1:n){
        integrand2 = function(x) {return((1-plogis(x,location=loc[i],scale = scale[i]))^2)}
        tmp[i]=integrate(integrand2,lower=b,upper=Inf)$value
      }
    }
    return(b)
  }
}

#-------------------------------------------------------------------------------
# Paras function

paras_EMOS <- function(fit, Data_test){
  gini.md <- function(x,na.rm=FALSE)  {     ## Michael Scheuerer's code
    if(na.rm & any(is.na(x)))  x <- x[!is.na(x)] 
    n <-length(x)
    return(4*sum((1:n)*sort(x,na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)
  }
  S <- fit$s
  A <- fit$a
  B <- fit$B
  C <- fit$c
  D <- fit$d
  SHAPE <- fit$q
  MEAN <- SCALE <- LOC  <-  rep(NaN, dim(Data_test)[1])
  for (i in 1:dim(Data_test)[1]) {
    f = Data_test[i,]
    MEAN[i] <- as.numeric(c(A,B)%*%c(1,f)+S*mean(f==0, na.rm = TRUE))  #location of GEV
    SCALE[i] <- C + D*gini.md(f, na.rm = TRUE)  #scale of GEV
    LOC[i] <- as.numeric(MEAN[i] - SCALE[i]*(gamma(1-SHAPE)-1)/SHAPE)
  }
  return(list(LOC = LOC, SCALE = SCALE, SHAPE = rep(SHAPE, dim(Data_test)[1])))
}

paras_BMA <- function(fit, Data_test){
  inverselogit<-function(x) {
    #
    # copyright 2006-present, University of Washington. All rights reserved.
    # for terms of use, see the LICENSE file
    
    if (is.na(x)) return(NA)
    if (x >= 0) {
      if (-x >= log(.Machine$double.eps)) {
        x <- exp(-x)
        1/(1+x)
      }
      else 1
    }
    else {
      if (x >= log(.Machine$double.xmin)) {
        x <- exp(x)
        x/(1+x)
      }
      else 0
    }
  }
  powfun <- function(x,power) x^power
  powinv <- function(x,power) x^(1/power)
  WEIGHTS <- fit$weights
  coef <- fit$varCoefs
  VAR <- matrix(NA, nrow =dim(Data_test)[1], ncol = 3)
  MEAN <- matrix(NA, nrow =dim(Data_test)[1], ncol = 3)
  PROB0 <- matrix(NA, nrow =dim(Data_test)[1], ncol = 3)
  for (i in 1:dim(Data_test)[1]) {
    f = Data_test[i,]
    
    VAR[i,] <- coef[1] + coef[2]*f
    
    fTrans <- sapply(f, powfun, power = fit$power)
    
    MEAN[i,] <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)
    
    PROB0[i,] <- sapply(apply(rbind( 1, fTrans, f == 0)*fit$prob0coefs,
                              2,sum), inverselogit)
    
  }
  return(list(P0 = PROB0, MEAN = MEAN, VAR = VAR , Weights = WEIGHTS))
}


