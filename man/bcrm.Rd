% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcrm.R
\name{bcrm}
\alias{bcrm}
\title{Bayesian Continual Reassessment Method for Phase I Dose-Escalation Trials}
\usage{
bcrm(
  stop = list(nmax = NULL, nmtd = NULL, precision = NULL, nmin = NULL, safety = NULL),
  data = NULL,
  p.tox0 = NULL,
  sdose = NULL,
  dose = NULL,
  ff,
  prior.alpha,
  cohort = 3,
  target.tox,
  constrain = TRUE,
  only.below = FALSE,
  sdose.calculate = "mean",
  pointest = "plugin",
  tox.cutpoints = NULL,
  loss = NULL,
  start = NULL,
  simulate = FALSE,
  nsims = 1,
  truep = NULL,
  threep3 = FALSE,
  threep3.start = 1,
  threep3.esc.only = FALSE,
  method = "exact",
  burnin.itr = 2000,
  production.itr = 2000,
  bugs.directory = "c:/Program Files/WinBUGS14/",
  plot = FALSE,
  seed = NULL,
  quietly = 10,
  file = NULL,
  N,
  tox,
  notox
)
}
\arguments{
\item{stop}{A list of stopping rules for the trial. One or more of the
following options should be specified \describe{ \item{list("nmax")}{ The
maximum sample size of the trial} \item{list("safety")}{Stops the trial if
the posterior probability of the lowest dose being above the target toxicity
level is greater than this number} \item{list("nmtd")}{ The maximum number
to be treated at final maximum tolerated dose (MTD) estimate,  \emph{i.e.} if
the next recommended dose has already been administered to \code{nmtd}
patients then the trial will stop} \item{list("precision")}{ A vector of the
lower and upper percentage points that the MTD 95\% credible intervals for
the risk of toxicity should lie within} \item{list("nmin")}{ The minimum
sample size of the trial. To be used in conjunction with \code{nmtd} or
\code{precision}} }}

\item{data}{A named data frame giving information about dose and toxicity
from previously recruited patients. If missing,  then it is assumed that no
data have thus far been collected. Contains the following variables:
\describe{ \item{list("patient")}{ Recruited patient numbers, 
\code{1, ..., n}} \item{list("dose")}{ Dose levels of recruited patients, 
ranging from \code{1, ..., k}} \item{list("tox")}{ An indicator variable for
each patient (1=toxicity,  0=no toxicity)} }}

\item{p.tox0}{A vector of length \code{k} listing the prior probabilities of
experiencing the outcome at each dose level \code{1, ...k} (CRM "skeleton").
The standardised dose levels are formed from these probabilities using the
inverse of the functional form,  with a plug-in estimate for the prior mean
or median of alpha,  as specified in \code{ff},  \code{prior.alpha} and
\code{sdose.calculate}. Alternatively standardised dose levels can be given
directly using \code{sdose}.}

\item{sdose}{A vector of length \code{k} listing the standardised doses to
be used in the CRM model. Only required if \code{p.tox0} is missing.}

\item{dose}{Optional vector of length \code{k} of actual doses for plotting
purposes}

\item{ff}{A string indicating the functional form of the dose-response
curve. Options are \describe{ \item{ht}{ 1-parameter hyperbolic tangent}
\item{logit1}{ 1-parameter logistic} \item{power}{ 1-parameter power}
\item{logit2}{ 2-parameter logistic} }}

\item{prior.alpha}{A list of length 3 containing the distributional
information for the prior. The first element is a number from 1-4 specifying
the type of distribution. Options are \enumerate{ \item Gamma(a, b),  where
a=shape,  b=scale: mean=a*b,  variance=a*b*b \item Uniform(a, b),  where a=min, 
b=max \item Lognormal(a, b),  where a=mean on the log scale,  b=variance on the
log scale \item Bivariate Lognormal(a, b),  where a=mean vector on the log
scale,  b=Variance-covariance matrix on the log scale. This prior should be
used only in conjunction with a two-parameter logistic model.  } The second
and third elements of the list are the parameters a and b,  respectively.}

\item{cohort}{The size of each cohort of patients that are sequentially
recruited to the trial. Defaults to 3}

\item{target.tox}{The target toxicity probability. Defaults to 1/3.}

\item{constrain}{Should a dose-skipping constraint be placed on the
escalation procedure, as imposed by a modified CRM? Defaults to TRUE.}

\item{only.below}{Should the dose chosen for the next patient be that with
the largest risk of DLT that is not greater than the target (TRUE), or the
dose that is closest to the target (FALSE); defaults to FALSE. If TRUE, and no
dose has risk of DLT below or equal to target, the dose chosen will be the lowest
dose level UNLESS one of the pre-defined stopping criteria is satisfied.
NOTE: \code{only.below} applies when \code{loss = NULL}.}

\item{sdose.calculate}{What plug-in estimate of the prior alpha should be
used to calculate the standardised doses? Options are \code{"mean"}
(default) or \code{"median"}. Only required if \code{sdose} is missing.}

\item{pointest}{Which summary estimate of the posterior distribution should
be used to choose the next dose. Options are \code{"plugin"} (default) where
the posterior mean of the model parameter(s) is plugged into the function
form to obtain estimates of toxicity,  or \code{"mean"} where the posterior
mean probabilities of toxicity are directly used. Alternatively,  a number
between 0 and 1 can be specified representing the quantile of the maximum
tolerated dose (MTD) posterior distribution (e.g. 0.5 specifies the
posterior median). This produces an Escalation With Overdose Control (EWOC)
design if the quantile is less than 0.5 (see details). Currently,  EWOC
designs must be fit using MCMC methods.}

\item{tox.cutpoints}{A vector of cutpoints for toxicity intervals if these
are to be used to choose next dose. Defaults to NULL. For example
Underdosing [0, 0.2],  Target dosing (0.2,  0.35],  Excessive toxicity (0.35, 
0.60],  Unacceptable toxicity (0.60,  1.00] set
\code{tox.cutpoints=c(0.2, 0.35, 0.60)}.}

\item{loss}{A vector of length \code{length(tox.cutpoints)+1} specifying the
losses associated with each toxicity interval,  e.g. Underdosing = 1,  Target
dosing = 0,  Excessive toxicity = 1,  Unacceptable toxicity = 2.}

\item{start}{Dose level to be used at the beginning of the trial. Required if
\code{constrain=TRUE}.}

\item{simulate}{Should a simulation be conducted to assess operating
characteristics? Defaults to \code{TRUE}. If \code{FALSE},  a single CRM
trial is run interactively,  allowing the user to input outcomes after each
cohort is recruited.}

\item{nsims}{Number of simulations to perform if \code{simulate==TRUE}
(defaults to 1).}

\item{truep}{A vector of length k giving the true probabilities of the
outcome (toxicity) at each dose level \code{1, ..., k} in order to simulate
data. Only required if \code{simulate=TRUE}}

\item{threep3}{Should operating characteristics of a standard 3+3 rule-based
design be calculated,  for comparison with \code{bcrm} design? Defaults to
\code{FALSE}. Only used in a simulation study,  i.e. when
\code{simulate=TRUE}.}

\item{threep3.start}{Starting dose level for when \code{threep3 = TRUE}. Defaults to 1, i.e. the lowest dose level}

\item{threep3.esc.only}{Whether to forbid de-escalation of doses when \code{threep3 = TRUE}. Defaults to \code{FALSE}}

\item{method}{Optimisation method: options are \code{"exact"} (the default), 
\code{"rjags"},  \code{"BRugs"},  or \code{"R2WinBUGS"}.}

\item{burnin.itr}{Number of burn-in iterations (default 2000).}

\item{production.itr}{Number of production iterations (default 2000).}

\item{bugs.directory}{Directory that contains the WinBUGS executable if
\code{method="R2WinBUGS"}. Defaults to "C:/Program Files/WinBUGS14/".}

\item{plot}{Should the dose-response curve be plotted after each cohort has
been entered? Defaults to FALSE.}

\item{seed}{Integer defining the state of the random number generator to
allow reproducible results. The default is to not specify a seed.}

\item{quietly}{How often to send a message back indicating how many
simulated trials have been performed. Defaults to \code{quietly = 10}
and the function reports back after every 10th simulation.}

\item{file}{File name where the dose-response plots are stored,  in a pdf
format. The program will amend the current sample size to the end of the
file name.}

\item{N}{Final sample size (deprecated). To be replaced with \code{stop} in
future versions.}

\item{tox}{(Deprecated). A vector of length \code{k} listing the number of
patients who have experienced the outcome (toxicity) at each dose level
\code{1, ..., k}.}

\item{notox}{(Deprecated). A vector of length \code{k} listing the number of
patients who have not experienced the outcome (toxicity) at each dose level
\code{1, ..., k}.}
}
\value{
\code{bcrm} returns an object of class "bcrm" or "bcrm.sim"; the
latter occurring when a simulation has been conducted
(\code{simulate=TRUE}). The function \code{\link{print}} (i.e.
\code{\link{print.bcrm}} or \code{\link{print.bcrm.sim}}) can be used to
obtain summary information about the design used,  the data observed,  current
posterior estimates of toxicity,  and the next recommended dose level.

An object of class "bcrm" is a list with the following components:
\item{dose}{Range of doses} \item{sdose}{Standardised doses} \item{tox}{A
vector of length \code{k} listing the number of patients who have
experienced the outcome (toxicity) at each dose level \code{1, ..., k}}
\item{notox}{A vector of length \code{k} listing the number of patients who
have not experienced the outcome (toxicity) at each dose level
\code{1, ..., k}} \item{ndose}{A list of lists containing for each cohort the
components \code{ndose},  the dose level recommended for the next patient, 
\code{est},  the estimated probabilities of toxicity using the chosen metric
(e.g. plugin,  mean,  quantile),  \code{mean},  the posterior mean probability
of toxicity at each dose,  \code{sd},  the posterior standard deviation for
probability of toxicity at each dose,  \code{quantiles},  the posterior
quantiles for probability of toxicity at each dose. This information is only
provided for cohorts recruited subsequent to any data given using \code{tox}
and \code{notox}. The first component relates to the prior information.}
\item{constrain}{Whether a constrained CRM design was used} \item{start}{The
starting dose for the latest run of the model if \code{constrain=TRUE}}
\item{target.tox}{The target toxicity level} \item{ff}{The functional form of
the dose-toxicity model;  "ht" = Hyperbolic tangent,  "logit1" = 1-parameter
logistic,  "power" = Power,  "logit2" = 2-parameter logistic} \item{method}{The calculation
method used} \item{pointest}{The summary estimate used to choose the next
dose,  see \code{pointest}} \item{prior.alpha}{Information about the prior
used for \eqn{alpha},  see \code{prior.alpha}} \item{data}{A data frame with
variables `patient',  `dose' and `tox' listing the dose levels and outcomes
of all patients in the trial}

An object of class "bcrm.sim" is a list of length \code{nsims}. Each
component is itself a list with components similar to those obtained from a
"bcrm" object. The print function,  \code{\link{print.bcrm.sim}} should be
used to obtain operating characteristics from the simulation.
}
\description{
Implements a wide variety of Bayesian CRM designs,  including 1-parameter, 
2-parameter and Escalation With Overdose Control (EWOC) designs. The program
can run interactively,  allowing the user to enter outcomes after each cohort
has been recruited,  or via simulation to assess operating characteristics.
}
\details{
\code{bcrm} implements a Bayesian continual reassessment method (CRM)
(O'Quigley \emph{et al.},  1990); an adaptive design in which cohorts of
patients are sequentially recruited into a Phase I trial. A binary toxicity
outcome is assumed (e.g. Dose Limiting Toxicity / No Dose Limiting
Toxicity). The current cohort are given a dose "closest" to the specified
target toxicity level,  as estimated from the posterior distributions of
toxicity at each dose level from the patients thus far recruited. If
\code{pointest="mean"} then the posterior mean probability of toxicity is
used to choose the next dose. If \code{pointest="plugin"},  however,  the
posterior mean of the model parameter(s) is plugged-into the functional form
of the dose-toxicity model. To implement an EWOC design (Babb \emph{et al.}, 
1998),  \code{pointest} should be a quantile,  \emph{q},  between 0 and 0.5.
The posterior distribution of the MTD (the dose in which the probability of
toxicity is equal to the target toxicity) is then calculated and the next
patient is given dose closest to the \emph{q}th quantile of the MTD
distribution.

Alternatively,  escalation can be based on intervals of toxicity from the
posterior distribution using a loss function,  see Neuenschwander \emph{et
al.},  2008. To implement this approach,  the user should specify the
cutpoints of the toxicity intervals using \code{tox.cutpoints} and the
associated losses using \code{loss}.

The possible choice of dose-toxicity model can be specified using \code{ff}, 
and includes the 1-parameter hyperbolic tangent,  logistic or power "working
models",  and the 2-parameter logistic model as follows: \describe{
\item{Hyperbolic
Tangent}{\deqn{p(Tox|d^*)=\left[(tanh(d^*)+1)/2\right]^\alpha}{p(Tox|d*)=[(tanh(d*)+1)/2]^\alpha}}
\item{Logistic (1-parameter)}{\deqn{p(Tox|d^*)=\frac{\exp(3+\alpha
d^*)}{1+\exp(3+\alpha d^*)}}{p(Tox|d*)=exp(3+\alpha d*)/(1+exp(3+\alpha
d*))}} \item{Power}{\deqn{p(Tox|d^*)={d^*}^\alpha}{p(Tox|d*)=d*^\alpha}}
\item{Logistic
(2-parameter)}{\deqn{p(Tox|d^*)=\frac{\exp(\log(\alpha_1)+\alpha_2
d^*)}{1+\exp(\log(\alpha_1)+\alpha_2 d^*)}}{
p(Tox|d*)=exp(log(\alpha_1)+\alpha_2 d*)/(1+exp(log(\alpha_1)+\alpha_2
d*))}} } where \eqn{\alpha>0} is the single positive-valued parameter for
the 1-parameter models,  and \eqn{\log(\alpha_1)}{log(\alpha_1)} and
\eqn{\alpha_2>0} are the intercept and slope parameters of the 2-parameter
model.

The standardised doses,  \eqn{d^*}{d*},  are specified by the user using
\code{sdose},  or alternatively the prior probability of toxicity at each
dose level is specified using \code{p.tox0}. If the latter is used,  then the
standardised doses are calculated using the inverse of the functional form
and a plug-in estimate of the prior mean or median,  as specified in
\code{sdose.calculate},  as follows \deqn{d^* = f^{-1}(\code{p.tox0}, \alpha=
a)}{d* = f^{-1}(p.tox0, \alpha= a)} where \eqn{f^{-1}} is the the inverse of
the chosen functional form,  and the parameter(s) of the model are set equal
to \eqn{a},  either the prior mean or median of \eqn{\alpha}.

Data that have already been accrued can be entered using the \code{data}
argument. A constrained CRM design can be implemented using
\code{constrain=TRUE},  in which case dose-skipping is prohibited (i.e. the
next cohort can only be dosed up to one dose level above the current
cohort). If a constrained model is used then the starting dose must be
specified using \code{start}. Alternatively, if data have already been
accrued,  then the dose level of the last recruited patient determines the
constraint for the next patient.

The prior is set using \code{prior.alpha}. For example
\code{prior.alpha=list(1, 1, 1)} specifies a Gamma prior with shape and scale
parameters both equal to one (\emph{i.e.} an Exponential(1) distribution), 
whilst \code{prior.alpha=list(2, 0, 10)} specifies a Uniform(0, 10) prior.

To specify a fixed maximum sample size of size \code{m} use
\code{stop=list(nmax=m)}. Alternatively,  the trial can stop after \code{m2}
patients have been treated at the current MTD estimate,  by setting
\code{stop=list(nmtd=m2)}.

To implement a safety constraint as specified in Zohar and Chevret (2001)
specify \code{stop=list(safety=p)},  where the trial is stopped if the
posterior probability that the lowest dose is greater than the target
toxicity probability is greater than \code{p}.

To stop the trial when the MTD estimate is within a certain level of
precision,  use \code{stop=list(precision=c(l, u))},  where \code{l} and
\code{u} are the lower and upper percentage points that the MTD 95\%
credible intervals for the risk of toxicity should lie within. Finally,  to
prevent the trial stopping too early using these rules,  the argument
\code{stop=list(nmin=m3)} can be used to ensure the sample size is greater
than or equal to \code{m3}. Stopping rules can be used on their own or in
combination.

The trial can be run interactively using \code{simulate=FALSE},  where the
user enters the outcomes for each new cohort,  or as a simulation study when
\code{simulate=TRUE}.

The default calculations use exact methods (\code{method="exact"}) to
calculate the mean and quantiles for the posterior distributions. There are
three choices for MCMC calculations: \code{method="rjags"}, 
\code{method="BRugs"} or \code{method="R2WinBUGS"}. The first uses the JAGS
software,  the second uses OpenBUGS,  whilst the latter uses WinBUGS. To
implement these methods,  users require one or more of these packages to be
installed on their system.

A simulated \code{bcrm} design can be compared with the standard 3+3
rule-based method,  see \code{\link{threep3}} for more details.
}
\note{
Currently,  the re-parameterisation of the two-parameter model proposed
by (Babb \emph{et al.},  1998) is not implemented. Therefore,  users wishing
to implement an EWOC design should check whether their choice of prior for
the model parameter(s) translates to a sensible prior for the MTD
distribution before they implement the design. For example \preformatted{
prior.alpha <- list(1, 1, 1); ff <- "ht"; target.tox <- 0.2
samples.alpha <- getprior(prior.alpha, 2000)
mtd <- find.x(ff, target.tox, alpha=samples.alpha) hist(mtd) }

One-parameter models are designed as working models only,  and should not be
used with an escalation strategy based on intervals of the posterior
probabilities of toxicity.
}
\examples{

## Dose-escalation cancer trial example as described in Neuenschwander et al 2008.
## Pre-defined doses
dose <- c(1, 2.5, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250)
## Pre-specified probabilities of toxicity
## [dose levels 11-15 not specified in the paper,  and are for illustration only]
p.tox0 <- c(0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050,
  0.100, 0.170, 0.300, 0.400, 0.500, 0.650, 0.800, 0.900)
## Data from the first 5 cohorts of 18 patients
data <- data.frame(patient=1:18, dose=rep(c(1:4, 7), c(3, 4, 5, 4, 2)), tox=rep(0:1, c(16, 2)))
## Target toxicity level
target.tox <- 0.30

## A 1-parameter power model is used,  with standardised doses calculated using 
## the plug-in prior median
## Prior for alpha is lognormal with mean 0 (on log scale) 
## and variance 1.34^2 (on log scale)
## The recommended dose for the next cohort if posterior mean is used
\dontrun{
Power.LN.bcrm <- bcrm(stop=list(nmax=18), data=data, p.tox0=p.tox0, dose=dose
  , ff="power", prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=FALSE
  , sdose.calculate="median", pointest="mean")
print(Power.LN.bcrm)
plot(Power.LN.bcrm)
}

## Simulate 10 replicate trials of size 36 (cohort size 3) using this design 
## with constraint (i.e. no dose-skipping) and starting at lowest dose
## True probabilities of toxicity are set to pre-specified probabilities (p.tox0) 
\dontrun{
Power.LN.bcrm.sim <- bcrm(stop=list(nmax=36), p.tox0=p.tox0, dose=dose, ff="power"
  , prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=TRUE
  , sdose.calculate="median", pointest="mean", start=1, simulate=TRUE, nsims=10, truep=p.tox0)
print(Power.LN.bcrm.sim)
plot(Power.LN.bcrm.sim)
}

## Comparing this CRM design with the standard 3+3 design 
## (only considering the first 12 dose levels)
\dontrun{
Power.LN.bcrm.compare.sim <- bcrm(stop=list(nmax=36), p.tox0=p.tox0[1:12], dose=dose[1:12]
  , ff="power", prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=TRUE
  , sdose.calculate="median", pointest="mean", start=1, simulate=TRUE, nsims=50
  , truep=p.tox0[1:12], threep3=TRUE)
print(Power.LN.bcrm.compare.sim, threep3=TRUE)
plot(Power.LN.bcrm.compare.sim, threep3=TRUE)
}

## A 2-parameter model,  using priors as specified in Neuenschwander et al 2008.
## Posterior mean used to choose the next dose
## Standardised doses using reference dose,  250mg
sdose <- log(dose/250)
## Bivariate lognormal prior for two parameters
mu <- c(2.15, 0.52)
Sigma <- rbind(c(0.84^2, 0.134), c(0.134, 0.80^2))
## Using rjags (requires JAGS to be installed)
\dontrun{
TwoPLogistic.mean.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
  , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
  , constrain=FALSE, pointest="mean", method="rjags")
print(TwoPLogistic.mean.bcrm)
plot(TwoPLogistic.mean.bcrm)
}

## A 2-parameter model,  using an EWOC design with feasibility bound (MTD quantile) 
## of 0.25 to choose the next dose
## Using rjags (requires JAGS to be installed)
\dontrun{
TwoPLogistic.EWOC0.25.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose, dose=dose
    , ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox, constrain=FALSE
    , pointest=0.25, method="rjags")
print(TwoPLogistic.EWOC0.25.bcrm)
plot(TwoPLogistic.EWOC0.25.bcrm)
}

## A 2-parameter model,  using a loss function based on intervals of toxicity to choose
## the next dose
## Using rjags (requires JAGS to be installed)
\dontrun{
## Toxicity cut-points
tox.cutpoints <- c(0.2, 0.35, 0.6)
## Losses associated with toxicity intervals 
## [0, 0.2]=1,  (0.2, 0.35]=0,  (0.35, 0.6]=1,  (0.6, 1]=2
loss <- c(1, 0, 1, 2)
TwoPLogistic.tox.intervals.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
  , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
  , constrain=FALSE, tox.cutpoints=tox.cutpoints, loss=loss, method="rjags")
print(TwoPLogistic.tox.intervals.bcrm)
plot(TwoPLogistic.tox.intervals.bcrm)
## Greater loss associated with overdosing and unacceptable toxicity
## [0, 0.2]=1,  (0.2, 0.35]=0,  (0.35, 0.6]=2,  (0.6, 1]=4
loss2 <- c(1, 0, 2, 4)
TwoPLogistic.tox.intervals.2.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
  , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
  , constrain=FALSE, tox.cutpoints=tox.cutpoints, loss=loss2, method="rjags")
print(TwoPLogistic.tox.intervals.2.bcrm)
plot(TwoPLogistic.tox.intervals.2.bcrm)
}


}
\references{
Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
of Statistical Software} (2013) 54: 1--26.
\doi{10.18637/jss.v054.i13}

O'Quigley J.,  Pepe M.,  Fisher L. Continual reassessment method: a practical
design for phase I clinical trials in cancer. \emph{Biometrics} (1990) 46:
33--48.

Babb J.,  Rogatko A.,  Zacks S. Cancer phase I clinical trials: efficient dose
escalation with overdose control. \emph{Statistics in Medicine} (1998) 17:
1103--1120.

Neuenschwander B.,  Branson M.,  Gsponer T. Critical aspects of the Bayesian
approach to phase I cancer trials. \emph{Statistics in Medicine} (2008) 27:
2420--2439.

Zohar S.,  Chevret S. The continual reassessment method: comparison of
Bayesian stopping rules for dose-ranging studies. \emph{Statistics in
Medicine} (2001) 20: 2827--2843.
}
\seealso{
\code{\link{print.bcrm}},  \code{\link{print.bcrm.sim}}, 
\code{\link{plot.bcrm}},  \code{\link{plot.bcrm.sim}},  \code{\link{threep3}}
}
\author{
Michael Sweeting (University of Leicester, UK; \email{michael.sweeting@leicester.ac.uk})
and Graham Wheeler (University College London, UK; \email{graham.wheeler@ucl.ac.uk}),
drawing on code originally developed by J. Jack Lee and Nan Chen,
Department of Biostatistics, University of Texas M.D. Anderson Cancer Center.
}
