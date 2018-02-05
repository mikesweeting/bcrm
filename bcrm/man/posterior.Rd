\name{Posterior}
\alias{Posterior.exact}
\alias{Posterior.rjags}
\alias{Posterior.BRugs}
\alias{Posterior.R2WinBUGS}
\title{
Calculate posterior distribution of CRM model parameter(s) 
}
\description{
Given the prior, functional form of the dose-toxicity model, and data, this function returns the posterior distribution (via either MCMC samples, or posterior summaries) of the CRM model parameter(s)
}
\usage{
Posterior.exact(tox,notox,sdose,ff,prior.alpha)
Posterior.rjags(tox,notox,sdose,ff,prior.alpha,burnin.itr,production.itr)
Posterior.BRugs(tox,notox,sdose,ff,prior.alpha,burnin.itr,production.itr)
Posterior.R2WinBUGS(tox,notox,sdose,ff,prior.alpha,burnin.itr,production.itr
  ,bugs.directory)
}

\arguments{
   \item{tox}{A vector of length \code{k} listing the number of patients who have experienced the outcome (toxicity) at each dose level \code{1,...,k}. If missing, then it is assumed that no data have thus far been collected.}
   \item{notox}{A vector of length \code{k} listing the number of patients who have not experienced the outcome (toxicity) at each dose level \code{1,...,k}. If missing, then it is assumed that no data have thus far been collected.}
   \item{sdose}{A vector of length \code{k} listing the standardised doses used in the CRM model.}
   \item{ff}{
A string indicating the functional form of the dose-response curve. Options are
    \describe{
     \item{ht}{ 1-parameter hyperbolic tangent}
      \item{logit1}{ 1-parameter logistic}
      \item{power}{ 1-parameter power}   
  \item{logit2}{ 2-parameter logistic}
    } 
}
 \item{prior.alpha}{
A list of length 3 containing the distributional information for the prior. The first element is a number from 1-4 specifying the type of distribution. Options are
    \enumerate{
      \item Gamma(a,b), where a=shape, b=scale: mean=a*b, variance=a*b*b
      \item Uniform(a,b), where a=min, b=max
      \item Lognormal(a,b), where a=mean on the log scale, b=standard deviation on the log scale
	\item Bivariate Lognormal(a,b), where a=mean vector on the log scale, b=Variance-covariance matrix on the log scale. This prior should be used only in conjunction with a two-parameter logistic model.
    } 
The second and third elements of the list are the parameters a and b, respectively.
}
  \item{burnin.itr}{
Number of burn-in iterations (default 2000).
}
  \item{production.itr}{
Number of production iterations (default 2000).
}
  \item{bugs.directory}{
Directory that contains the WinBUGS executable if \code{method="R2WinBUGS"}. Defaults to "C:/Program Files/WinBUGS14/".
}

}
\details{
\code{Posterior.exact} produces posterior summary statistics of the CRM model parameter(s), and probabilities of toxicity at the dose levels using exact Bayesian computation (integration) of the posterior distribution. If \code{Posterior.BRugs} or \code{Posterior.R2WinBUGS} is specified, then posterior samples of the CRM model parameter(s) is returned by the function.
}

\author{
Michael Sweeting \email{mjs212@medschl.cam.ac.uk} (University of Cambridge, UK)
}
\references{
 Sweeting M., Mander A., Sabin T. \pkg{bcrm}: Bayesian Continual Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal of Statistical Software} (2013) 54: 1--26. \url{http://www.jstatsoft.org/article/view/v054i13}
}
\seealso{
\code{\link{bcrm}}
}
\examples{
## Dose-escalation cancer trial example as described in Neuenschwander et al 2008.
## Pre-defined doses
dose<-c(1,2.5,5,10,15,20,25,30,40,50,75,100,150,200,250)
## Pre-specified probabilities of toxicity
## [dose levels 11-15 not specified in the paper, and are for illustration only]
p.tox0<-c(0.010,0.015,0.020,0.025,0.030,0.040,0.050,0.100,0.170,0.300,0.400,0.500,0.650
  ,0.800,0.900)
## Data from the first 5 cohorts of 18 patients
tox<-c(0,0,0,0,0,0,2,0,0,0,0,0,0,0,0)
notox<-c(3,4,5,4,0,0,0,0,0,0,0,0,0,0,0)
## Target toxicity level
target.tox<-0.30
## Lognormal prior
prior.alpha<-list(3,0,1.34^2)
## Power functional form
ff<-"power"
## Standardised doses
sdose<-find.x(ff,p.tox0,alpha=1)

## Posterior distribution of the model parameter using exact computation
post.exact<-Posterior.exact(tox,notox,sdose,ff,prior.alpha)
print(post.exact)

## Posterior distribution of the model parameter using rjags
post.rjags<-Posterior.rjags(tox,notox,sdose,ff,prior.alpha
  ,burnin.itr=2000,production.itr=2000)
print(mean(post.rjags))
hist(post.rjags)

## Posterior distribution of the model parameter using BRugs (Windows and i386 Linux only)
if(Sys.info()["sysname"] \%in\% c("Windows","Linux"))\{
	post.BRugs<-Posterior.BRugs(tox,notox,sdose,ff,prior.alpha
	  ,burnin.itr=2000,production.itr=2000)
	print(mean(post.BRugs))
	hist(post.BRugs)
	\}
}
