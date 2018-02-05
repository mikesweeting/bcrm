\name{plot.threep3}
\alias{plot.threep3}
\title{
Plot the operating characteristics from a standard 3+3 trial
}
\description{
Plots of the operating characteristics obtained from a standard 3+3 trial, using \code{\link{threep3}}
}
\usage{
\method{plot}{threep3}(x, file = NULL, \ldots)
}
\arguments{
  \item{x}{
	An object of class "threep3", as returned by \code{\link{threep3}}.
}
  \item{file}{
File name where the operating characteristic plot is stored, in a pdf format.
}
 \item{...}{Further arguments passed to or from other methods}
}
\details{
This function plots the sample size distribution, the experimentation distribution, the recommended dose distribution and the percentage of subjects who experience the toxicity outcome (dose-limiting toxicity) for the standard 3+3 trial.
}
\author{
Michael Sweeting \email{mjs212@medschl.cam.ac.uk} (University of Cambridge, UK) 
}
\references{
 Sweeting M., Mander A., Sabin T. \pkg{bcrm}: Bayesian Continual Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal of Statistical Software} (2013) 54: 1--26. \url{http://www.jstatsoft.org/article/view/v054i13}
}
\seealso{
\code{\link{threep3}}
}
