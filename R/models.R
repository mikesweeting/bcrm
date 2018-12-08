#-----------------------------------------------------------------------
## HT Gamma
## WinBUGS code for a hyperbolic tangent model with Gamma prior on alpha 
#-----------------------------------------------------------------------
HTGamma <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  pow((exp(d[i]) / (exp(d[i]) + exp(-d[i]))),  alpha)
  }
  alpha ~ dgamma(p1,  p2)
}

#-----------------------------------------------------------------------
## HT Unif
## WinBUGS code for a hyperbolic tangent model with Uniform prior on alpha 
#-----------------------------------------------------------------------
HTUnif <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  pow((exp(d[i]) / (exp(d[i]) + exp(-d[i]))),  alpha)
  }
  alpha ~ dunif(p1,  p2)
}
#-----------------------------------------------------------------------
## HT Lognormal
## WinBUGS code for a hyperbolic tangent model with Lognormal prior on alpha 
#-----------------------------------------------------------------------
HTLogNormal <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  pow((exp(d[i]) / (exp(d[i]) + exp(-d[i]))),  alpha)
  }
  prec <- 1/p2
  alpha ~ dlnorm(p1, prec)
}

#-----------------------------------------------------------------------
## Logistic Gamma
## WinBUGS code for a one-parameter logistic model with Gamma prior on alpha (slope of logistic function) 
## Intercept parameter is fixed at 3.0
#-----------------------------------------------------------------------
LogisticGamma <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  exp(3.0 + alpha * d[i]) / (1 + exp(3.0 + alpha * d[i]))
  }
  alpha ~ dgamma(p1,  p2)
}

#-----------------------------------------------------------------------
## Logistic Uniform
## WinBUGS code for a one-parameter logistic model with Uniform prior on alpha (slope of logistic function) 
## Intercept parameter is fixed at 3.0
#-----------------------------------------------------------------------
LogisticUnif <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  exp(3.0 + alpha * d[i]) / (1 + exp(3.0 + alpha * d[i]))
  }
  alpha ~ dunif(p1,  p2)
}

#-----------------------------------------------------------------------
## Logistic Lognormal
## WinBUGS code for a one-parameter logistic model with Lognormal prior on alpha (slope of logistic function) 
## Intercept parameter is fixed at 3.0
#-----------------------------------------------------------------------
LogisticLogNormal <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  exp(3.0 + alpha * d[i]) / (1 + exp(3.0 + alpha * d[i]))
  }
  prec <- 1/p2
  alpha ~ dlnorm(p1,  prec)
}

#-----------------------------------------------------------------------
## Power Gamma
## WinBUGS code for a power model with Gamma prior on alpha
#-----------------------------------------------------------------------
PowerGamma <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  pow(d[i],  alpha)
  }
  alpha ~ dgamma(p1,  p2)
}

#-----------------------------------------------------------------------
## Power Uniform
## WinBUGS code for a power model with Uniform prior on alpha
#-----------------------------------------------------------------------
PowerUnif <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  pow(d[i],  alpha)
  }
  alpha ~ dunif(p1,  p2)
}
## Power LogNormal
## WinBUGS code for a power model with LogNormal prior on alpha
#-----------------------------------------------------------------------
PowerLogNormal <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    p[i]  <-  pow(d[i],  alpha)
  }
  prec <- 1/p2
  alpha ~ dlnorm(p1,  prec)
}
#-----------------------------------------------------------------------
## 2-parameter Logistic Bivariate Lognormal
## WinBUGS code for a two-parameter logistic model with Bivariate Lognormal prior on alpha[1] (intercept) and alpha[2] (slope of logistic function) 
#-----------------------------------------------------------------------
TwoPLogisticLogNormal <- function(){
  for (i in 1:N1){
    s[i] ~ dbin(p[i],  n[i])    
    logit(p[i])  <-  log.alpha[1] + alpha[2] * d[i]
  }
  alpha[1] <- exp(log.alpha[1])
  alpha[2] <- exp(log.alpha[2])
  Omega[1:2, 1:2] <- inverse(p2[, ])
  log.alpha[1:2] ~ dmnorm(p1[],  Omega[, ])
}

