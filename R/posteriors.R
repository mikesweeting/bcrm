# ----------------------------------------------------------------------
# generates a vector of values and prior distribution
#
#     prior.alpha --> list containing the information for prior
#               [[1]] - the prior distribution type:
#                     1 - gamma: mean=a*b; var=a*b*b
#                     2 - uniform: a+b*unif(0, 1)
#           3 - lognormal: lnorm(a, b),  where mean=a,  var=b
#           4 - log Multivariate normal(a, b),  where a=mean vector,  b=Variance-covariance matrix
#               [[2]] - a: first parameter of the prior distribution
#               [[3]] - b: second parameter of the prior distribution
#    
# ----------------------------------------------------------------------

#' Samples from the specified prior distribution.
#' 
#' A sample of specified size is obtained from the prior distribution.
#' 
#' A vector of size \code{n} is returned from the specified prior distribution.
#' 
#' @param prior.alpha A list of length 3 containing the distributional
#' information for the prior. The first element is a number from 1-4 specifying
#' the type of distribution. Options are \enumerate{ \item Gamma(a, b),  where
#' a=shape,  b=scale: mean=a*b,  variance=a*b*b \item Uniform(a, b),  where a=min, 
#' b=max \item Lognormal(a, b),  where a=mean on the log scale,  b=variance on the
#' log scale \item Bivariate Lognormal(a, b),  where a=mean vector on the log
#' scale,  b=Variance-covariance matrix on the log scale. This prior should be
#' used only in conjunction with a two-parameter logistic model.  } The second
#' and third elements of the list are the parameters a and b,  respectively.
#' @param n The number of samples.
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK),  drawing on code originally developed by J. Jack Lee and Nan
#' Chen,  Department of Biostatistics,  the University of Texas M. D. Anderson
#' Cancer Center
#' @seealso \code{\link{bcrm}},  \code{\link{find.x}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#' @examples
#' 
#' prior.alpha <- list(1, 1, 1)
#' samples.alpha <- getprior(prior.alpha, 2000)
#' hist(samples.alpha)
#' 
#' @export getprior
getprior  <-  function(prior.alpha,  n) {
  type  <-  prior.alpha[[1]]
  a     <-  prior.alpha[[2]]
  b     <-  prior.alpha[[3]]    
  if ( type==1 ) {
    prior <- rgamma(n, a)*b        
  }
  else if ( type==2 ) {
    prior <- sapply(1:length(a), function(i){runif(n, a[i], b[i])})
  }
  else if (type==3) { 
    prior <- rlnorm(n, a, sqrt(b))
  }
  else if (type==4) {
    log.prior <- rmvnorm(n, a, b)
    prior <- exp(log.prior)
  }
  return (prior)
}

# ----------------------------------------------------------------------
#     returns the vector containing sampled data from winbug
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#          "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
#     burnin.itr  --> number of burn-in (adaptive) iterations
#     production.itr --> number of production iterations
# ----------------------------------------------------------------------
Posterior.rjags  <-  function(tox,  notox, sdose, ff,  prior.alpha,  burnin.itr,  production.itr)
{    
  all.patient  <-  tox + notox
  datan  <- all.patient[all.patient!=0]
  datas  <- tox[all.patient!=0]
  datad  <- sdose[all.patient!=0]
  k  <-  length(datan)
  if (k == 1)
  {
    datan  <-  c(datan,  0)
    datas  <-  c(datas,  0)
    datad  <-  c(datad,  0)
  }
  mydata  <-  list(N1 = k,  s = datas, n = datan, d = datad,  p1 = prior.alpha[[2]],  p2 = prior.alpha[[3]])
  model.file <- if (prior.alpha[[1]] == 1)
  {
    if (ff == "ht")
      HTGamma
    else if (ff == "logit1")
      LogisticGamma
    else if (ff == "power")
      PowerGamma
    else stop("Functional form not currently available with specified prior distribution")
  }
  else if (prior.alpha[[1]] == 2)
  {
    if (ff == "ht")
      HTUnif
    else if (ff == "logit1")
      LogisticUnif
    else if (ff == "power")
      PowerUnif
    else stop("Functional form not currently available with specified prior distribution")
  }        
  else if (prior.alpha[[1]] == 3)
  {
    if (ff == "ht")
      HTLogNormal
    else if (ff == "logit1")
      LogisticLogNormal
    else if (ff == "power")
      PowerLogNormal
    else stop("Functional form not currently available with specified prior distribution")
  }        
  else if (prior.alpha[[1]] == 4)
  {
    if (ff == "logit2")
      TwoPLogisticLogNormal            
    else stop("Functional form not currently available with specified prior distribution")
  }    
  path.model <- file.path(tempdir(),  "model.file.txt")
  R2WinBUGS::write.model(model.file, path.model)
  inits.list <- if (ff=="logit2")
    list(list(log.alpha=c(-3, 0), .RNG.seed=sample(1:1e6,  size=1), .RNG.name="base::Wichmann-Hill"), list(log.alpha=c(-3, 0), .RNG.seed=sample(1:1e6,  size=1), .RNG.name="base::Wichmann-Hill"))
  else list(list(alpha=1, .RNG.seed=sample(1:1e6,  size=1), .RNG.name="base::Wichmann-Hill"), list(alpha=1, .RNG.seed=sample(1:1e6,  size=1), .RNG.name="base::Wichmann-Hill"))
  
  jagsobj <- rjags::jags.model(path.model, data=mydata, n.chains=2, quiet=TRUE, inits=inits.list)
  update(jagsobj, n.iter=burnin.itr, progress.bar="none")
  tt <- rjags::jags.samples(jagsobj,  "alpha",  n.iter=production.itr/2, progress.bar="none")
  if(ff=="logit2"){
    t <- cbind(c(tt$alpha[1, , ]), c(tt$alpha[2, , ]))
  } else {   t <-  c(tt$alpha) }
  return(t)
}

# ----------------------------------------------------------------------
#     returns the vector containing sampled data from winbug
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#          "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
#     burnin.itr  --> number of burn-in iterations
#     production.itr --> number of production iterations
# ----------------------------------------------------------------------
Posterior.BRugs  <-  function(tox,  notox, sdose, ff,  prior.alpha,  burnin.itr,  production.itr)
{    
  all.patient  <-  tox + notox
  datan  <- all.patient[all.patient!=0]
  datas  <- tox[all.patient!=0]
  datad  <- sdose[all.patient!=0]
  k  <-  length(datan)
  if (k == 1)
  {
    datan  <-  c(datan,  0)
    datas  <-  c(datas,  0)
    datad  <-  c(datad,  0)
  }
  mydata  <-  list(N1 = k,  s = datas, n = datan, d = datad,  p1 = prior.alpha[[2]],  p2 = prior.alpha[[3]])
  model.file <- if (prior.alpha[[1]] == 1)
  {
    if (ff == "ht")
      HTGamma
    else if (ff == "logit1")
      LogisticGamma
    else if (ff == "power")
      PowerGamma
    else stop("Functional form not currently available with specified prior distribution")
  }
  else if (prior.alpha[[1]] == 2)
  {
    if (ff == "ht")
      HTUnif
    else if (ff == "logit1")
      LogisticUnif
    else if (ff == "power")
      PowerUnif
    else stop("Functional form not currently available with specified prior distribution")
  }        
  else if (prior.alpha[[1]] == 3)
  {
    if (ff == "ht")
      HTLogNormal
    else if (ff == "logit1")
      LogisticLogNormal
    else if (ff == "power")
      PowerLogNormal
    else stop("Functional form not currently available with specified prior distribution")
  }        
  else if (prior.alpha[[1]] == 4)
  {
    if (ff == "logit2")
      TwoPLogisticLogNormal            
    else stop("Functional form not currently available with specified prior distribution")
  }    
  path.model <- file.path(tempdir(),  "model.file.txt")
  path.data <- file.path(tempdir(),  "data.file.txt")
  R2WinBUGS::write.model(model.file, path.model)
  BRugs::bugsData(mydata, path.data)
  BRugsSeed <- sample(1:14, 1)
  BRugs::BRugsFit(path.model, path.data,  inits=NULL,  numChains = 2,  parametersToSave="alpha", 
                  nBurnin = burnin.itr,  nIter = production.itr/2,  BRugsVerbose = FALSE, DIC=F, seed=BRugsSeed)
  if(ff=="logit2"){
    t <- cbind(BRugs::samplesSample("alpha[1]"), BRugs::samplesSample("alpha[2]"))
  } else {   t <-  BRugs::samplesSample("alpha") }
  return(t)
}

# ----------------------------------------------------------------------
#     returns the vector containing sampled data from winbug
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#         "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
#     burnin.itr  --> number of burn-in iterations
#     production.itr --> number of production iterations
#   bugs.directory --> directory that contains the WinBUGS executable,  defaults to C:/Program Files/WinBUGS14/
# ----------------------------------------------------------------------
Posterior.R2WinBUGS  <-  function(tox,  notox, sdose, ff,  prior.alpha,  burnin.itr,  production.itr, bugs.directory)
{    
  all.patient  <-  tox + notox
  datan  <- all.patient[all.patient!=0]
  datas  <- tox[all.patient!=0]
  datad  <- sdose[all.patient!=0]
  k  <-  length(datan)
  if (k == 1)
  {
    datan  <-  c(datan,  0)
    datas  <-  c(datas,  0)
    datad  <-  c(datad,  0)
  }
  mydata  <-  list(N1 = k,  s = datas, n = datan, d = datad,  p1 = prior.alpha[[2]],  p2 = prior.alpha[[3]])
  ## initdat <- list(list(alpha = 0.5), list(alpha=1))
  parameters <- "alpha"
  model.file <- if (prior.alpha[[1]] == 1)
  {
    if (ff == "ht")
      HTGamma
    else if (ff == "logit1")
      LogisticGamma
    else if (ff == "power")
      PowerGamma
    else stop("Functional form not currently available with specified prior distribution")
  }
  else if (prior.alpha[[1]] == 2)
  {
    if (ff == "ht")
      HTUnif
    else if (ff == "logit1")
      LogisticUnif
    else if (ff == "power")
      PowerUnif
    else stop("Functional form not currently available with specified prior distribution")
  }        
  else if (prior.alpha[[1]] == 3)
  {
    if (ff == "ht")
      HTLogNormal
    else if (ff == "logit1")
      LogisticLogNormal
    else if (ff == "power")
      PowerLogNormal
    else stop("Functional form not currently available with specified prior distribution")
  }        
  else if (prior.alpha[[1]] == 4)
  {
    if (ff == "logit2")
      TwoPLogisticLogNormal            
    else stop("Functional form not currently available with specified prior distribution")
  }    
  R2WinBUGS.seed <- sample(1:1e6, 1)
  res <- R2WinBUGS::bugs(mydata,  inits=NULL,  parameters,  model.file, n.chains=2, n.thin=1, n.iter=burnin.itr+production.itr/2, n.burnin=burnin.itr, DIC=F, bugs.directory=bugs.directory, bugs.seed=R2WinBUGS.seed)  
  if(is.null(dim(res$summary))){
    res$summary <- t(as.matrix(res$summary))
  }
  if(any(res$summary[, "Rhat"]>1.01)){ 
    warning("Convergence may not have been achieved: Trace plots shown")
    par(mfrow=c(dim(res$sims.array)[[3]], 1))
    for(i in 1:dim(res$sims.array)[[3]]){
      plot(res$sims.array[, 1, i], type="l")
      lines(res$sims.array[, 2, i], col=2)
    }
  }
  t <- apply(res$sims.array, 3, rbind)
  return(t)
}

# ----------------------------------------------------------------------
#     returns the posterior mean of alpha and actualised values of toxic probabilities at each dose level
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#         "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
# ----------------------------------------------------------------------
Posterior.exact <- function(tox, notox, sdose, ff, prior.alpha){    
  all.patient  <-  tox + notox
  data.tox  <- tox[all.patient!=0]
  data.notox  <- notox[all.patient!=0]
  data.dose  <- sdose[all.patient!=0]
  
  ## Following code fixes bug if no data available (prior only)
  if(length(data.dose)==0){
    data.dose <- sdose[1]
    data.tox <- data.notox <- 0
  }  
  wmodel <- which.f(ff)
  
  prior <- switch(prior.alpha[[1]]
                , "1"=function(alpha, prior.alpha){dgamma(alpha, shape=prior.alpha[[2]], scale=prior.alpha[[3]])}
                , "2"=function(alpha, prior.alpha){dunif(alpha, min=prior.alpha[[2]], max=prior.alpha[[3]])}
                , "3"=function(alpha, prior.alpha){dlnorm(alpha, meanlog=prior.alpha[[2]], sdlog=sqrt(prior.alpha[[3]]))}
                , "4"=function(alpha, prior.alpha){1/(alpha[, 1]*alpha[, 2])*dmvnorm(log(alpha), mean=prior.alpha[[2]], sigma=prior.alpha[[3]])})
  
  if(prior.alpha[[1]]!=4){
    ## Scaling factor to prevent likelihood getting too small
    C <- 1/prod(sapply(1:length(data.dose), function(d){wmodel(data.dose[d], 1)^data.tox[d]*(1-wmodel(data.dose[d], 1))^data.notox[d]}))
    lik <- function(dose, tox, notox, alpha, C){
      l <- rep(1, length(alpha))
      for(d in 1:length(dose)){
        l <- l*wmodel(dose[d], alpha)^tox[d]*(1-wmodel(dose[d], alpha))^notox[d]
      }
      C*l
    }
    int.norm.constant <- function(alpha, dose, tox, notox, prior.alpha){lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
    int.alpha.mean <- function(alpha, dose, tox, notox, prior.alpha){alpha*lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
    int.dose.mean <- function(alpha, new.dose, dose, tox, notox, prior.alpha){wmodel(new.dose, alpha)*lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
    int.dose.sd <- function(alpha, new.dose, dose.mean, dose, tox, notox, prior.alpha){(wmodel(new.dose, alpha)-dose.mean)^2*lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
    
    norm.constant <- integrate(int.norm.constant, ifelse(prior.alpha[[1]]==2, prior.alpha[[2]], 0), ifelse(prior.alpha[[1]]==2, prior.alpha[[3]], Inf), dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]
    alpha.mean <- integrate(int.alpha.mean, ifelse(prior.alpha[[1]]==2, prior.alpha[[2]], 0), ifelse(prior.alpha[[1]]==2, prior.alpha[[3]], Inf), dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
    dose.mean <- sapply(sdose, function(d){integrate(int.dose.mean, ifelse(prior.alpha[[1]]==2, prior.alpha[[2]], 0), ifelse(prior.alpha[[1]]==2, prior.alpha[[3]], Inf), new.dose=d, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant})
    dose.sd <- sapply(1:length(sdose), function(d){sqrt(integrate(int.dose.sd, ifelse(prior.alpha[[1]]==2, prior.alpha[[2]], 0), ifelse(prior.alpha[[1]]==2, prior.alpha[[3]], Inf), new.dose=sdose[d], dose.mean=dose.mean[d], dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant)})
    
    cdf <- function(par){integrate(int.norm.constant, ifelse(prior.alpha[[1]]==2, prior.alpha[[2]], 0), par, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant}
    fn <- function(par, q){abs(cdf(par)-q)}
    alpha.quantiles <- sapply(c(0.025, 0.25, 0.5, 0.75, 0.975), function(q){  max.x <- seq(0, 10*alpha.mean, length=100)[which.min(sapply(seq(0, 10*alpha.mean, length=100), function(i){fn(i, q)}))+1]
    optimize(fn, c(0, max.x), q=q,  tol = .Machine$double.eps^0.50)$minimum
    })
    dose.quantiles <- sapply(sdose, function(d){wmodel(d, sort(alpha.quantiles, TRUE))})
    rownames(dose.quantiles) <- c("2.5%", "25%", "50%", "75%", "97.5%")
  } else {
    ## VECTORISED FORM OF LIK
    lik.2param <- function(dose, tox, notox, alpha){
      l <- rep(1, nrow(alpha))
      for(d in 1:length(dose)){
        l <- l*wmodel(dose[d], alpha)^tox[d]*(1-wmodel(dose[d], alpha))^notox[d]
      }
      l
    }
    joint.posterior <- function(alpha1, alpha2, dose, tox, notox, prior.alpha){lik.2param(dose, tox, notox, cbind(alpha1, alpha2))*prior(cbind(alpha1, alpha2), prior.alpha)}
    joint.posterior.2 <- function(alpha2, alpha1, dose, tox, notox, prior.alpha){lik.2param(dose, tox, notox, cbind(alpha1, alpha2))*prior(cbind(alpha1, alpha2), prior.alpha)}
    marginal.alpha1 <- function(alpha1, dose, tox, notox, prior.alpha){sapply(alpha1, function(a1){integrate(joint.posterior.2, 0, Inf, alpha1=a1, dose=dose, tox=tox, notox=notox, prior.alpha=prior.alpha)$value})}  
    marginal.alpha2 <- function(alpha2, dose, tox, notox, prior.alpha){sapply(alpha2, function(a2){integrate(joint.posterior, 0, Inf, alpha2=a2, dose=dose, tox=tox, notox=notox, prior.alpha=prior.alpha)$value})}
    int.alpha1.mean <- function(alpha1, dose, tox, notox, prior.alpha){alpha1*marginal.alpha1(alpha1, dose, tox, notox, prior.alpha)}
    int.alpha2.mean <- function(alpha2, dose, tox, notox, prior.alpha){alpha2*marginal.alpha2(alpha2, dose, tox, notox, prior.alpha)}
    int.dose.mean.2param <- function(alpha1, alpha2, new.dose, dose, tox, notox, prior.alpha){wmodel(new.dose, cbind(alpha1, alpha2))*joint.posterior(alpha1, alpha2, dose, tox, notox, prior.alpha)}
    int.dose.mean.2 <- function(alpha2, new.dose, dose, tox, notox, prior.alpha){sapply(alpha2, function(a2){integrate(int.dose.mean.2param, 0, Inf, alpha2=a2, new.dose=new.dose, dose=dose, tox=tox, notox=notox, prior.alpha=prior.alpha)$value})}
    int.dose.sd.2param <- function(alpha1, alpha2, new.dose, dose.mean, dose, tox, notox, prior.alpha){(wmodel(new.dose, cbind(alpha1, alpha2))-dose.mean)^2*joint.posterior(alpha1, alpha2, dose, tox, notox, prior.alpha)}
    int.dose.sd.2 <- function(alpha2, new.dose, dose.mean, dose, tox, notox, prior.alpha){sapply(alpha2, function(a2){integrate(int.dose.sd.2param, 0, Inf, alpha2=a2, new.dose=new.dose, dose.mean=dose.mean, dose=dose, tox=tox, notox=notox, prior.alpha=prior.alpha)$value})}
    
    norm.constant <- integrate(marginal.alpha2, 0, Inf, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]
    alpha1.mean <- integrate(int.alpha1.mean, 0, Inf, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
    alpha2.mean <- integrate(int.alpha2.mean, 0, Inf, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
    alpha.mean <- c(alpha1.mean, alpha2.mean)
    dose.mean <- sapply(sdose, function(d){integrate(int.dose.mean.2, 0, Inf, new.dose=d, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant})
    dose.sd <- sapply(1:length(sdose), function(d){sqrt(integrate(int.dose.sd.2, 0, Inf, new.dose=sdose[d], dose.mean=dose.mean[d], dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant)})
    
    dose.quantiles <- NULL
    
    # Using cuhre
    #prior <- function(alpha, prior.alpha){1/(alpha[1]*alpha[2])*dmvnorm(log(alpha), mean=prior.alpha[[2]], sigma=prior.alpha[[3]])}
    #wmodel <- function(dose, alpha) {1/(1+exp(-log(alpha[1])-alpha[2]*dose))}
    ## Scaling factor to prevent likelihood getting too small
    #C <- 1/prod(wmodel(data.dose, c(1, 1))^data.tox*(1-wmodel(data.dose, c(1, 1)))^data.notox)
    #lik <- function(dose, tox, notox, alpha, C){
    #  l <- prod(wmodel(dose, alpha)^tox*(1-wmodel(dose, alpha))^notox)
    #  l
    #}
    #int.norm.constant <- function(alpha, dose, tox, notox, prior.alpha){lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
    #norm.constant <- cuhre(2, 1, int.norm.constant, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, lower=c(0, 0), upper=c(Inf, Inf), flags=list(verbose=0))$value
    #int.dose.mean <- function(alpha, new.dose, dose, tox, notox, prior.alpha){wmodel(new.dose, alpha)*lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
    #dose.mean <- cuhre(2, length(sdose), int.dose.mean, new.dose=sdose, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, lower=c(0, 0), upper=c(500, 500))$value/norm.constant
    
  }
  return(list(alpha.mean=alpha.mean, dose.mean=dose.mean, dose.sd=dose.sd, dose.quantiles=dose.quantiles))
}


# ----------------------------------------------------------------------
#     Cut-down version of Posterior.exact for use with simulations and two-parameter models
#
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#     sdose       --> vector containing the dose level
#     ff          --> the model applied in the study
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#         "logit2" - Two-parameter logistic
#     prior.alpha  --> list of prior distribution information for parameter alpha
#   pointest   --> Which summary estimate of the posterior distribution should be used to choose next dose,  "plugin" (default),  "mean" or a numerical quantile between 0 and 1 (e.g. 0.5). If NULL then escalation based on posterior intervals (Loss function approach) is assumed
# ----------------------------------------------------------------------
Posterior.exact.sim <- function(tox, notox, sdose, ff, prior.alpha, pointest){    
  
  all.patient  <-  tox + notox
  data.tox  <- tox[all.patient!=0]
  data.notox  <- notox[all.patient!=0]
  data.dose  <- sdose[all.patient!=0]
  
  ## Following code fixes bug if no data available (prior only)
  if(length(data.dose)==0){
    data.dose <- sdose[1]
    data.tox <- data.notox <- 0
  }  
  wmodel <- which.f(ff)
  
  prior <- switch(prior.alpha[[1]]
                , "1"=function(alpha, prior.alpha){dgamma(alpha, shape=prior.alpha[[2]], scale=prior.alpha[[3]])}
                , "2"=function(alpha, prior.alpha){dunif(alpha, min=prior.alpha[[2]], max=prior.alpha[[3]])}
                , "3"=function(alpha, prior.alpha){dlnorm(alpha, meanlog=prior.alpha[[2]], sdlog=sqrt(prior.alpha[[3]]))}
                , "4"=function(alpha, prior.alpha){1/(alpha[, 1]*alpha[, 2])*dmvnorm(log(alpha), mean=prior.alpha[[2]], sigma=prior.alpha[[3]])})
  
  if(prior.alpha[[1]]!=4){
    ## Scaling factor to prevent likelihood getting too small
    C <- 1/prod(sapply(1:length(data.dose), function(d){wmodel(data.dose[d], 1)^data.tox[d]*(1-wmodel(data.dose[d], 1))^data.notox[d]}))
    lik <- function(dose, tox, notox, alpha, C){
      l <- rep(1, length(alpha))
      for(d in 1:length(dose)){
        l <- l*wmodel(dose[d], alpha)^tox[d]*(1-wmodel(dose[d], alpha))^notox[d]
      }
      C*l
    }
    int.norm.constant <- function(alpha, dose, tox, notox, prior.alpha){lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
    norm.constant <- integrate(int.norm.constant, ifelse(prior.alpha[[1]]==2, prior.alpha[[2]], 0), ifelse(prior.alpha[[1]]==2, prior.alpha[[3]], Inf), dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]
    if(pointest=="plugin"){
      int.alpha.mean <- function(alpha, dose, tox, notox, prior.alpha){alpha*lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
      alpha.mean <- integrate(int.alpha.mean, ifelse(prior.alpha[[1]]==2, prior.alpha[[2]], 0), ifelse(prior.alpha[[1]]==2, prior.alpha[[3]], Inf), dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
      dose.mean <- NULL
    } else if(pointest=="mean"){
      int.dose.mean <- function(alpha, new.dose, dose, tox, notox, prior.alpha){wmodel(new.dose, alpha)*lik(dose, tox, notox, alpha, C)*prior(alpha, prior.alpha)}
      alpha.mean <- NULL
      dose.mean <- sapply(sdose, function(d){integrate(int.dose.mean, ifelse(prior.alpha[[1]]==2, prior.alpha[[2]], 0), ifelse(prior.alpha[[1]]==2, prior.alpha[[3]], Inf), new.dose=d, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant})
    } 
  } else {
    ## VECTORISED FORM OF LIK
    lik.2param <- function(dose, tox, notox, alpha){
      l <- rep(1, nrow(alpha))
      for(d in 1:length(dose)){
        l <- l*wmodel(dose[d], alpha)^tox[d]*(1-wmodel(dose[d], alpha))^notox[d]
      }
      l
    }
    joint.posterior <- function(alpha1, alpha2, dose, tox, notox, prior.alpha){lik.2param(dose, tox, notox, cbind(alpha1, alpha2))*prior(cbind(alpha1, alpha2), prior.alpha)}
    joint.posterior.2 <- function(alpha2, alpha1, dose, tox, notox, prior.alpha){lik.2param(dose, tox, notox, cbind(alpha1, alpha2))*prior(cbind(alpha1, alpha2), prior.alpha)}
    marginal.alpha1 <- function(alpha1, dose, tox, notox, prior.alpha){sapply(alpha1, function(a1){integrate(joint.posterior.2, 0, Inf, alpha1=a1, dose=dose, tox=tox, notox=notox, prior.alpha=prior.alpha)$value})}  
    marginal.alpha2 <- function(alpha2, dose, tox, notox, prior.alpha){sapply(alpha2, function(a2){integrate(joint.posterior, 0, Inf, alpha2=a2, dose=dose, tox=tox, notox=notox, prior.alpha=prior.alpha)$value})}
    int.alpha1.mean <- function(alpha1, dose, tox, notox, prior.alpha){alpha1*marginal.alpha1(alpha1, dose, tox, notox, prior.alpha)}
    int.alpha2.mean <- function(alpha2, dose, tox, notox, prior.alpha){alpha2*marginal.alpha2(alpha2, dose, tox, notox, prior.alpha)}
    int.dose.mean.2param <- function(alpha1, alpha2, new.dose, dose, tox, notox, prior.alpha){wmodel(new.dose, cbind(alpha1, alpha2))*joint.posterior(alpha1, alpha2, dose, tox, notox, prior.alpha)}
    int.dose.mean.2 <- function(alpha2, new.dose, dose, tox, notox, prior.alpha){sapply(alpha2, function(a2){integrate(int.dose.mean.2param, 0, Inf, alpha2=a2, new.dose=new.dose, dose=dose, tox=tox, notox=notox, prior.alpha=prior.alpha)$value})}
    int.dose.sd <- function(alpha1, alpha2, new.dose, dose.mean, dose, tox, notox, prior.alpha){(wmodel(new.dose, cbind(alpha1, alpha2))-dose.mean)^2*joint.posterior(alpha1, alpha2, dose, tox, notox, prior.alpha)}
    int.dose.sd.2 <- function(alpha2, new.dose, dose.mean, dose, tox, notox, prior.alpha){sapply(alpha2, function(a2){integrate(int.dose.sd, 0, Inf, alpha2=a2, new.dose=new.dose, dose.mean=dose.mean, dose=dose, tox=tox, notox=notox, prior.alpha=prior.alpha)$value})}
    
    norm.constant <- integrate(marginal.alpha2, 0, Inf, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]
    if(pointest=="plugin"){
      alpha1.mean <- integrate(int.alpha1.mean, 0, Inf, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
      alpha2.mean <- integrate(int.alpha2.mean, 0, Inf, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant
      alpha.mean <- c(alpha1.mean, alpha2.mean)
      dose.mean <- NULL
    } else if(pointest=="mean"){
      alpha.mean <- NULL
      dose.mean <- sapply(sdose, function(d){integrate(int.dose.mean.2, 0, Inf, new.dose=d, dose=data.dose, tox=data.tox, notox=data.notox, prior.alpha=prior.alpha, rel.tol=.Machine$double.eps^0.5)[[1]]/norm.constant})
    }
  }
  return(list(alpha.mean=alpha.mean, dose.mean=dose.mean))
}