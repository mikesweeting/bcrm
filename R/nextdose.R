# ----------------------------------------------------------------------
#     find the new dose and posteria mean of alpha
#
#     samples.alpha  --> sampled data of variable(s) alpha
#     sdose     --> vector of allowable dose values
#     ff       --> integer value:
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#             "logit2" - Two-parameter Logistic
#     target.tox    --> desired probability of event (toxicity)
#  constrain    --> TRUE or FALSE
#  first       --> TRUE or FALSE
#  pointest    --> "plugin",  "mean" or a quantile
#  current    --> current dose level (of last treated cohort) - 0 if this is first cohort
# quantiles --> quantiles of the posterior distribution to calculate for each dose
# only.below --> TRUE or FALSE; if TRUE, only choose largest dose that has estimate below
#                target (if no dose below target, choose lowest dose level unless stopping
#                criteria satisfied); otherwise, choose dose level closest to target
#
# ----------------------------------------------------------------------
nextdose <- function(samples.alpha, sdose, ff, target.tox, constrain, first, pointest, current, tox.cutpoints, loss, quantiles, only.below){
  f <- which.f(ff)
  k <- length(sdose)
  samples.sdose <- sapply(sdose, function(x){f(x, samples.alpha)})
  mean <- apply(samples.sdose, 2, mean)
  median <- apply(samples.sdose, 2, median)
  sd <- apply(samples.sdose, 2, sd)
  quantiles <- apply(samples.sdose, 2, quantile, quantiles)
  
  if(!is.null(loss)){
    probs <- sapply(1:k, function(i){hist(samples.sdose[, i],  c(0, tox.cutpoints, 1),  plot = FALSE)$counts/dim(samples.sdose)[1]})
    est <- apply(loss*probs, 2, sum)
    target <- 0
    # Next dose
    ndose <- if(!constrain){
      if(!first){which.min(abs(est-target))}else{min(current+1,k)}
    } else {which.min(abs(est[1:min(current+1, k)]-target))}
    
  } else {
    probs <- NULL
    if(pointest=="plugin"){
      mean.alpha <- apply(as.matrix(samples.alpha), 2, mean)
      if(length(mean.alpha)>1) mean.alpha <- t(as.matrix(mean.alpha))
    } 
    if(is.numeric(pointest)){
      mtd <- find.x(ff, ptox=target.tox, alpha=samples.alpha)
      target <- quantile(mtd, pointest)
    } else {
      target <- target.tox
    }
    est <- if(pointest=="plugin"){ f(sdose, mean.alpha)
    } else {if(pointest=="mean"){apply(samples.sdose, 2, mean)
    } else { sdose }}
  
  # Next dose
  est.below<-est[which(est<=target)]
  ndose <- if(only.below==FALSE){if(!constrain){
    if(!first){which.min(abs(est-target))}else{
      if(length(current)!=0){min(current+1,k)}else{which.min(abs(est-target))}
      }
  } else {which.min(abs(est[1:min(current+1, k)]-target))}
  }else{
    if(length(est.below)==0){1}else{
    if(!constrain){
      if(!first){which.min(abs(est.below-target))}else{
        if(length(current)!=0){min(current+1,k)}else{which.min(abs(est.below-target))}
        }
    } else {which.min(abs(est.below[1:min(current+1, k)]-target))}}
    
  }
  }
  return(list(ndose=ndose, est=est, mean=mean, sd=sd, quantiles=quantiles, target=target, probs=probs))
}


# ----------------------------------------------------------------------
#     find the new dose and posteria mean of alpha
#
#     alpha     --> posterior mean of alpha and posterior mean p(DLT) at each dose level
#     sdose     --> vector of allowable dose values
#     ff       --> integer value:
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#         "logit2" - Two-parameter Logistic
#     target.tox    --> desired probability of event (toxicity)
#  constrain    --> TRUE or FALSE
#  first       --> TRUE or FALSE
#  pointest    --> "mean" or "plugin"
#  current    --> current dose level (of last treated cohort) - 0 if this is first cohort
# only.below --> TRUE or FALSE; if TRUE, only choose largest dose that has estimate below
#                target (if no dose below target, choose lowest dose level unless stopping
#                criteria satisfied); otherwise, choose dose level closest to target
#
# ----------------------------------------------------------------------
nextdose.exact <- function(alpha, sdose, ff, target.tox, constrain, first, pointest, current, only.below){
  f <- which.f(ff)
  k <- length(sdose)
  est <- if(pointest=="mean") alpha$dose.mean
  else if(pointest=="plugin") f(sdose, alpha$alpha.mean)
  else stop("Quantile estimation not available for exact computation,  please use pointest='mean' or 'plugin'")
  mean <- alpha$dose.mean
  sd <- alpha$dose.sd
  quantiles <- alpha$dose.quantiles
  est.below<-est[which(est<=target.tox)]
  ndose <- if(only.below==FALSE){if(!constrain){
    if(!first){which.min(abs(est-target.tox))}else{
      if(length(current)!=0){min(current+1,k)}else{which.min(abs(est-target.tox))}
      }
  } else {if(!first){which.min(abs(est[1:min(current+1, k)]-target.tox))}else{min(current+1,k)}}
  }else{
    if(length(est.below)==0){1}else{
      if(!constrain){
        if(!first){which.min(abs(est.below-target.tox))}else{
          if(length(current)!=0){min(current+1,k)}else{which.min(abs(est.below-target.tox))}
        }
      } else {which.min(abs(est.below[1:min(current+1, k)]-target.tox))}}
    
  }
  return(list(ndose=ndose, est=est, mean=mean, sd=sd, quantiles=quantiles))
}

# ----------------------------------------------------------------------
#     Cut-down version of nextdose.exact for use with simulations
#
#     alpha     --> posterior mean of alpha and posterior mean p(DLT) at each dose level
#     sdose     --> vector of allowable dose values
#     ff       --> integer value:
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#         "logit2" - Two-parameter Logistic
#     target.tox    --> desired probability of event (toxicity)
#  constrain    --> TRUE or FALSE
#  first       --> TRUE or FALSE
#  pointest    --> "mean" or "plugin"
#  current    --> current dose level (of last treated cohort) - 0 if this is first cohort
# only.below --> TRUE or FALSE; if TRUE, only choose largest dose that has estimate below
#                target (if no dose below target, choose lowest dose level unless stopping
#                criteria satisfied); otherwise, choose dose level closest to target
#
# ----------------------------------------------------------------------
nextdose.exact.sim <- function(alpha, sdose, ff, target.tox, constrain, first, pointest, current, only.below){
  f <- which.f(ff)
  k <- length(sdose)
  est <- if(pointest=="mean") alpha$dose.mean
  else if(pointest=="plugin") f(sdose, alpha$alpha.mean)
  else stop("Quantile estimation not available for exact computation,  please use pointest='mean' or 'plugin'")
  est.below<-est[which(est<=target.tox)]
  ndose <- if(only.below==FALSE){if(!constrain){
    if(!first){which.min(abs(est-target.tox))}else{
      if(length(current)!=0){min(current+1,k)}else{which.min(abs(est-target.tox))}
    }
  } else {if(!first){which.min(abs(est[1:min(current+1, k)]-target.tox))}else{min(current+1,k)}}
  }else{
    if(length(est.below)==0){1}else{
      if(!constrain){
        if(!first){which.min(abs(est.below-target.tox))}else{
          if(length(current)!=0){min(current+1,k)}else{which.min(abs(est.below-target.tox))}
        }
      } else {which.min(abs(est.below[1:min(current+1, k)]-target.tox))}}
    
  }
  return(list(ndose=ndose, est=est))
}


# -----------------------------------------------------------------------
#     Return the calculation models
#     ff       --> integer value:
#                    "ht" - hyperbolic tangent
#                    "logit1" - logistic
#                    "power" - Power
#         "logit2" - Two-parameter logistic
# -----------------------------------------------------------------------
which.f  <-  function( ff ) {
  return( switch(ff, 
                 ht=function(dose, alpha) {((tanh(dose)+1)/2)^alpha}, 
                 logit1=function(dose, alpha) {1/(1+exp(-3-alpha*dose))}, 
                 power=function(dose, alpha) {dose^alpha}, 
                 logit2=function(dose, alpha) {1/(1+exp(-log(alpha[, 1])-alpha[, 2]*dose))} ))
}

# ----------------------------------------------------------------------
#     Find the dose corresponding to a certain p(DLT)
#
#     ff         --> functional form for the probability density
#     ptox       --> desired probability of DLT
#     alpha      --> the parameter(s) of the model
# ----------------------------------------------------------------------


#' Obtain samples from the maximum tolerated dose (MTD) distribution.
#' 
#' Given a posterior (or prior) sample of the parameters,  this function inverts
#' the given functional form to obtain samples from the MTD distribution.
#' 
#' Given a posterior (or prior) sample of the parameters,  this function inverts
#' the given functional form to obtain samples from the MTD distribution or any
#' other targeted quantile.
#' 
#' @param ff A string indicating the functional form of the dose-response
#' curve. Options are \describe{ \item{ht}{ 1-parameter hyperbolic tangent}
#' \item{logit1}{ 1-parameter logistic} \item{power}{ 1-parameter power}
#' \item{logit2}{ 2-parameter logistic} }
#' @param ptox The required probability of DLT. For example,  if the MTD
#' distribution is sought then set \code{ptox} to the target toxicity level.
#' @param alpha A sample from the posterior (or prior) distribution of the
#' parameter(s).
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK),  drawing on code originally developed by J. Jack Lee and Nan
#' Chen,  Department of Biostatistics,  the University of Texas M. D. Anderson
#' Cancer Center
#' @seealso \code{\link{bcrm}},  \code{\link{getprior}},  \code{\link{Posterior.exact}},  \code{\link{Posterior.BRugs}},  \code{\link{Posterior.R2WinBUGS}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#' @examples
#' 
#' ## Dose-escalation cancer trial example as described in Neuenschwander et al 2008.
#' ## Pre-defined doses
#' dose <- c(1, 2.5, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250)
#' ## Pre-specified probabilities of toxicity
#' ## [dose levels 11-15 not specified in the paper,  and are for illustration only]
#' p.tox0 <- c(0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050,
#'   0.100, 0.170, 0.300, 0.400, 0.500, 0.650, 0.800, 0.900)
#' ## Data from the first 5 cohorts of 18 patients
#' tox <- c(0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0)
#' notox <- c(3, 4, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' ## Target toxicity level
#' target.tox <- 0.30
#' 
#' ## Prior distribution for the MTD given a lognormal(0, 1.34^2) distribution for alpha
#' ## and a power model functional form
#' prior.alpha <- list(3, 0, 1.34^2)
#' ff <- "power"
#' samples.alpha <- getprior(prior.alpha, 2000)
#' mtd <- find.x(ff, target.tox, alpha=samples.alpha)
#' hist(mtd)
#' 
#' ## Standardised doses
#' sdose <- find.x(ff, p.tox0, alpha=1)
#' 
#' ## Posterior distribution of the MTD (on standardised dose scale) using data 
#' ## from the cancer trial described in Neuenschwander et al 2008.
#' ## Using rjags
#' \dontrun{
#' posterior.samples <- Posterior.rjags(tox, notox, sdose, ff, prior.alpha
#'   , burnin.itr=2000, production.itr=2000)
#' posterior.mtd <- find.x(ff, target.tox, alpha=posterior.samples)
#' hist(posterior.mtd)
#' }
#' 
#' 
#' @export find.x
find.x  <-  function( ff,  ptox,  alpha ) {
  if ( ff == "ht" ) {
    x <- atanh(2*ptox^(1/alpha)-1)
  }
  else if ( ff == "logit1" )
    x  <-  (qlogis(ptox)-3)/alpha
  else if ( ff == "power")
    x  <-  exp(log(ptox)/alpha)
  else if (ff=="logit2"){
    if(is.vector(alpha))  alpha <- matrix(alpha, ncol=2)
    x  <-  (qlogis(ptox)-log(alpha[, 1]))/alpha[, 2]
  }
  return( x )
}

