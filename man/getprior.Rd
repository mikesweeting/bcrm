\name{getprior}
\alias{getprior}
\title{
Samples from the specified prior distribution.
}
\description{
A sample of specified size is obtained from the prior distribution.
}
\usage{
getprior(prior.alpha, n)
}
\arguments{
  \item{prior.alpha}{
	A list of length 3 containing the distributional information for the prior. The first element is a number from 1-4 specifying the type of distribution. Options are
    \enumerate{
      \item Gamma(a,b), where a=shape, b=scale: mean=a*b, variance=a*b*b
      \item Uniform(a,b), where a=min, b=max
      \item Lognormal(a,b), where a=mean on the log scale, b=variance on the log scale
	\item Bivariate Lognormal(a,b), where a=mean vector on the log scale, b=Variance-covariance matrix on the log scale. This prior should be used only in conjunction with a two-parameter logistic model.
    } 
The second and third elements of the list are the parameters a and b, respectively.
}
 \item{n}{The number of samples.}
}
\details{
A vector of size \code{n} is returned from the specified prior distribution.
}
\author{
Michael Sweeting \email{mjs212@medschl.cam.ac.uk} (University of Cambridge, UK), drawing on code originally developed by J. Jack Lee and Nan Chen, Department of Biostatistics, the University of Texas M. D. Anderson Cancer Center
}
\references{
 Sweeting M., Mander A., Sabin T. \pkg{bcrm}: Bayesian Continual Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal of Statistical Software} (2013) 54: 1--26. \url{http://www.jstatsoft.org/article/view/v054i13}
}
\seealso{
\code{\link{bcrm}}, \code{\link{find.x}}
}
\examples{
prior.alpha<-list(1,1,1)
samples.alpha<-getprior(prior.alpha,2000)
hist(samples.alpha)
}
