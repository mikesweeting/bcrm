\name{print.bcrm}
\alias{print.bcrm}
\alias{print.bcrm.sim}
\title{
Print information regarding a trial conducted using the Bayesian continuous reassessment method
}
\description{
Print method for a trial or series of trials conducted using a \code{\link{bcrm}} model. 
}
\usage{
\method{print}{bcrm}(x, \ldots)
\method{print}{bcrm.sim}(x, tox.cutpoints = NULL, trajectories = FALSE, threep3 = FALSE, \ldots)
}
\arguments{
  \item{x}{
An object of class "bcrm" or "bcrm.sim" as returned by \code{\link{bcrm}}
}
 \item{tox.cutpoints}{
An optional argument passed to \code{print.bcrm.sim} specifying the cutpoints of toxicity for which the operating characteristics are to be categorised. Defaults to \code{seq(from=0,to=1,by=0.2)}
}
 \item{trajectories}{
Should the individual simulation dose and outcome trajectories be returned? Defaults to FALSE. 
}
 \item{threep3}{
Should operating characteristics of a standard 3+3 rule-based design be displayed alongside those from the \code{bcrm} design? Defaults to \code{FALSE}.
}
 \item{...}{Further arguments passed to or from other methods}
}
\details{
If a single trial is conducted, then the \code{\link{print}} function currently produces summary information about the design used, the data observed, current posterior estimates of toxicity, and the next recommended dose level. 
If a simulation study is conducted, then the following operating characteristics are printed:
\describe{
	\item{Experimentation proportion}{Proportion of patients recruited to each dose, and to each true region of toxicity, across the simulated trials}
	\item{Recommendation proportion}{Proportion of trials that recommend each of the dose levels as the final maximum tolerated dose (i.e. with toxicity "closest" to the target toxicity level), and the associated regions of true toxicity for the recommended MTDs}
}
If \code{trajectories = TRUE} then the dose level administered and outcome observed are returned as matrices for every patient (column) in every simulation (row). 
If \code{threep3 = TRUE} then the operating characteristics of the standard 3+3 design are displayed alongside those of the \code{bcrm} design (see \code{\link{threep3}} for more details).
}
\value{
The following two components are returned from \code{print.bcrm.sim}:
\item{exp}{A matrix with number of rows equal to the number of doses, and number of columns equal to the number of simulations. Gives the experimentation proportions for each dose within each simulation.}
\item{rec}{A vector with length equal to the number of simulations, giving the recommended MTD for each simulation.}
}

\author{
Michael Sweeting \email{mjs212@medschl.cam.ac.uk} (University of Cambridge, UK)
}
\references{
 Sweeting M., Mander A., Sabin T. \pkg{bcrm}: Bayesian Continual Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal of Statistical Software} (2013) 54: 1--26. \url{http://www.jstatsoft.org/article/view/v054i13}
}
\seealso{
\code{\link{bcrm}}, \code{\link{threep3}}
}
