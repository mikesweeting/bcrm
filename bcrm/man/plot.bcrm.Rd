\name{plot.bcrm}
\alias{plot.bcrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the estimated dose-toxicity curve
}
\description{
The estimated dose-toxicity curve using the Bayesian continuous reassessment method is plotted for the patients thus far recruited into the trial
}
\usage{
\method{plot}{bcrm}(x, file = NULL, each = FALSE, trajectory = FALSE, \ldots)
}
\arguments{
  \item{x}{
	An object of class "bcrm", as returned by \code{\link{bcrm}}
}
  \item{each}{
	Should posterior summaries be plotted after each recruited cohort? Defaults to FALSE.
}
  \item{trajectory}{
	Should the sequential dose trajectory of the recruited patients be plotted, along with the observed toxicities? Defaults to FALSE.
}
  \item{file}{
File name where the dose-response plots are stored, in a pdf format. The program will amend the current sample size to the end of the file name.
}
 \item{...}{Further arguments passed to or from other methods}
}
\details{
The estimated 2.5\%, 25\%, 50\%, 75\%, 97.5\% quantiles of the probability of toxicity are plotted for each dose. Additionally, a histogram of the number of toxicities and non-toxicities is plotted at each experimented dose. 

If \code{trajectory = TRUE} then the sequential dose trajectory and observed toxicities are plotted.
}
\author{
Michael Sweeting \email{mjs212@medschl.cam.ac.uk} (University of Cambridge, UK) 
}
\references{
 Sweeting M., Mander A., Sabin T. \pkg{bcrm}: Bayesian Continual Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal of Statistical Software} (2013) 54: 1--26. 
 \url{http://www.jstatsoft.org/article/view/v054i13}
}
\seealso{
\code{\link{bcrm}}
}
