#' Bayesian Continual Reassessment Method for Phase I Dose-Escalation Trials
#' 
#' Implements a wide variety of Bayesian CRM designs. The program can run
#' interactively, allowing the user to enter outcomes after each cohort has
#' been recruited, or via simulation to assess operating characteristics.
#' 
#' \tabular{ll}{ Package: \tab bcrm\cr Type: \tab Package\cr Version: \tab
#' 0.5.1\cr Date: \tab 2019-04-03\cr License: \tab GPL (>= 2) \cr LazyLoad:
#' \tab yes\cr }
#' 
#' @name bcrm-package
#' @docType package
#' @author Michael Sweeting \email{michael.sweeting@@leicester.ac.uk}
#' 
#' Graham Wheeler \email{graham.wheeler@@ucl.ac.uk}
#' 
#' Maintainer: Graham Wheeler \email{graham.wheeler@@ucl.ac.uk}
#' @references Sweeting M., Mander A., Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
NULL





#' Calculate posterior distribution of CRM model parameter(s)
#' 
#' Given the prior, functional form of the dose-toxicity model, and data, this
#' function returns the posterior distribution (via either MCMC samples, or
#' posterior summaries) of the CRM model parameter(s)
#' 
#' \code{Posterior.exact} produces posterior summary statistics of the CRM
#' model parameter(s), and probabilities of toxicity at the dose levels using
#' exact Bayesian computation (integration) of the posterior distribution. If
#' \code{Posterior.BRugs} or \code{Posterior.R2WinBUGS} is specified, then
#' posterior samples of the CRM model parameter(s) is returned by the function.
#' 
#' @aliases Posterior.exact Posterior.rjags Posterior.BRugs Posterior.R2WinBUGS
#' @param tox A vector of length \code{k} listing the number of patients who
#' have experienced the outcome (toxicity) at each dose level \code{1,...,k}.
#' If missing, then it is assumed that no data have thus far been collected.
#' @param notox A vector of length \code{k} listing the number of patients who
#' have not experienced the outcome (toxicity) at each dose level
#' \code{1,...,k}. If missing, then it is assumed that no data have thus far
#' been collected.
#' @param sdose A vector of length \code{k} listing the standardised doses used
#' in the CRM model.
#' @param ff A string indicating the functional form of the dose-response
#' curve. Options are \describe{ \item{ht}{ 1-parameter hyperbolic tangent}
#' \item{logit1}{ 1-parameter logistic} \item{power}{ 1-parameter power}
#' \item{logit2}{ 2-parameter logistic} }
#' @param prior.alpha A list of length 3 containing the distributional
#' information for the prior. The first element is a number from 1-4 specifying
#' the type of distribution. Options are \enumerate{ \item Gamma(a,b), where
#' a=shape, b=scale: mean=a*b, variance=a*b*b \item Uniform(a,b), where a=min,
#' b=max \item Lognormal(a,b), where a=mean on the log scale, b=standard
#' deviation on the log scale \item Bivariate Lognormal(a,b), where a=mean
#' vector on the log scale, b=Variance-covariance matrix on the log scale. This
#' prior should be used only in conjunction with a two-parameter logistic
#' model.  } The second and third elements of the list are the parameters a and
#' b, respectively.
#' @param burnin.itr Number of burn-in iterations (default 2000).
#' @param production.itr Number of production iterations (default 2000).
#' @param bugs.directory Directory that contains the WinBUGS executable if
#' \code{method="R2WinBUGS"}. Defaults to "C:/Program Files/WinBUGS14/".
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge, UK)
#' @seealso \code{\link{bcrm}}
#' @references Sweeting M., Mander A., Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#' @examples
#' 
#' ## Dose-escalation cancer trial example as described in Neuenschwander et al 2008.
#' ## Pre-defined doses
#' dose<-c(1,2.5,5,10,15,20,25,30,40,50,75,100,150,200,250)
#' ## Pre-specified probabilities of toxicity
#' ## [dose levels 11-15 not specified in the paper, and are for illustration only]
#' p.tox0<-c(0.010,0.015,0.020,0.025,0.030,0.040,0.050,0.100,0.170,0.300,0.400,0.500,0.650
#'   ,0.800,0.900)
#' ## Data from the first 5 cohorts of 18 patients
#' tox<-c(0,0,0,0,0,0,2,0,0,0,0,0,0,0,0)
#' notox<-c(3,4,5,4,0,0,0,0,0,0,0,0,0,0,0)
#' ## Target toxicity level
#' target.tox<-0.30
#' ## Lognormal prior
#' prior.alpha<-list(3,0,1.34^2)
#' ## Power functional form
#' ff<-"power"
#' ## Standardised doses
#' sdose<-find.x(ff,p.tox0,alpha=1)
#' 
#' ## Posterior distribution of the model parameter using exact computation
#' post.exact<-Posterior.exact(tox,notox,sdose,ff,prior.alpha)
#' print(post.exact)
#' 
#' ## Posterior distribution of the model parameter using rjags
#' post.rjags<-Posterior.rjags(tox,notox,sdose,ff,prior.alpha
#'   ,burnin.itr=2000,production.itr=2000)
#' print(mean(post.rjags))
#' hist(post.rjags)
#' 
#' ## Posterior distribution of the model parameter using BRugs (Windows and i386 Linux only)
#' if(Sys.info()["sysname"] %in% c("Windows","Linux")){
#' 	post.BRugs<-Posterior.BRugs(tox,notox,sdose,ff,prior.alpha
#' 	  ,burnin.itr=2000,production.itr=2000)
#' 	print(mean(post.BRugs))
#' 	hist(post.BRugs)
#' 	}
#' 