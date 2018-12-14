## Bayesian CRM - extending original code of J. Jack Lee and Nan Chen,  Department of Biostatistics,  the University of Texas M. D. Anderson Cancer Center
## Now using exact inference,  rjags,  R2WinBUGS or BRugs 

### MJS & GMW 10/10/18

if(getRversion() >= "2.15.1") globalVariables(c("N1", "pow", "d", "alpha", "p2", "logit <- ", "log.alpha", "inverse", "patient", "toxicity", "est", "target.tox", "q2.5", "q97.5", "q50", "q25", "q75", "ndose", "tox.cutpoints", "xmin", "ymin", "xmax", "ymax", "Loss", "Outcome", "traj", "Statistic", "..density..", "Toxicity", "weight", "Method", "..count..", "rec", "obs", "obs.rounded", "count", "n"))

# ----------------------------------------------------------------------
#   bcrm. Conduct a Bayesian CRM trial simulating outcomes from a true model
#   stop       --> A list of stopping rules containing information on when the trial should be terminated
#      nmax = Maximum sample size
#      nmtd = Maximum number to be treated at final MTD estimate      
#      precision = Lower and Upper risks of toxicity that must contain the 95% posterior interval (2.5%ile and 97.5%ile) of the recommended dose in order to stop the trial
#      nmin = Minimum sample size required in the trial
#     safety = Posterior probability that DLT rate of all doses is higher than the TTL - then trial is stopped for safety
#  data    --> A named data frame giving information about dose and toxicity from previously recruited patients. Contains the following variables
#      "patient" - recruited patient number (1, ..., n)
#      "dose" - dose level 
#      "tox" - Indicator variable (1=toxicity,  0=no toxicity)  
#  p.tox0     --> prior probabilities of response at dose levels
#     sdose      --> standardised dose levels (given if p.tox0 is missing)
#    dose    --> optional vector of dose labels (for printing and plotting purposes)
#  ff       --> functional form of dose-response curve
#                    "ht" - hyperbolic tangent
#                    "logit1" - 1-parameter logistic
#                    "power" - Power
#          "logit2" - Two-parameter logistic     
#     prior.alpha--> list containing the information to generate the
#                   prior distribution for parameter (vector) alpha
#                   [[1]] - the prior distribution type:
#                       1 - gamma: mean=a*b; var=a*b*b
#                       2 - uniform: unif(a, b)
#          3 - lognormal: lnorm(a, b),  where a=mean on the log scale,  b=var on the log scale
#             4 - log Multivariate normal(a, b),  where a=mean vector on log scale,  b=Variance-covariance matrix (log scale) (for multiparameter functional forms)
#                   [[2]] - a: first parameter (scalar/vector) of the prior distribution
#                   [[3]] - b: second parameter (scalar/vector) of the prior distribution
#     cohort     --> cohort size - default = 3
#     target.tox --> target toxicity 
#   constrain --> TRUE (default) or FALSE
#    sdose.calculate -> What plug-in estimate of the prior alpha should be used to calculate the standardised doses? Options are "mean" (default) or "median"
#   pointest   --> Which summary estimate of the posterior distribution should be used to choose next dose,  "plugin" (default),  "mean" or a numerical quantile between 0 and 1 (e.g. 0.5). If NULL then escalation based on posterior intervals (Loss function approach) is assumed
#  tox.cutpoints --> Specify the cutpoints for toxicity intervals if these are to be used to choose next dose,  e.g. Underdosing [0, 0.2],  Target dosing (0.2,  0.35],  Excessive toxicity (0.35,  0.60],  Unacceptable toxicity (0.60,  1.00]
#  loss --> The losses associated with each toxicity interval,  e.g. Underdosing = 1,  Target dosing =0,  Excessive toxicity=1,  Unacceptable toxicity=2
#     start      --> Starting dose level. Required if constrain=TRUE and no previous data is provided
#    simulate -->  Perform a simulation to assess operating characteristics (Default=TRUE). If FALSE,  a single CRM trial is run interactively,  allowing the user to input outcomes after each cohort is recruited?
#    nsims  --> No. of simulations to perform if simulate==T (defaults to 1)
#   truep       --> True probabilities of response at dose levels to simulate data. Only should be specified if simulate=TRUE
#    threep3     --> If TRUE (default is FALSE) then operating characteristics of the simulated design are compared against a standard rule-based 3+3 design
#  method      --> Optimisation method: options are "exact" (the default),  "rjags",  "BRugs",  or "R2WinBUGS"
#     burnin.itr --> No. of burn-in iterations
#     production.itr --> No. of production iterations
#    bugs.directory --> directory that contains the WinBUGS executable if R2WinBUGS is being used,  defaults to "C:/Program Files/WinBUGS14/"
#  plot       --> Should dose response curve be plotted after each cohort is entered? Defaults to FALSE
#   file       --> name of the file where the dose-response plots are stored,  in a pdf format
#   seed    --> (Optional) Initial seed value (integer).
# quietly --> Suppress output 
#  N --> Final sample size (Deprecated). To be replaced with stop in future versions
#     tox         --> number of successes (toxicities) (Deprecated)
#     notox       --> number of failed patients (no-toxicities) (Deprecated)
# ----------------------------------------------------------------------


#' Bayesian Continual Reassessment Method for Phase I Dose-Escalation Trials
#' 
#' Implements a wide variety of Bayesian CRM designs,  including 1-parameter, 
#' 2-parameter and Escalation With Overdose Control (EWOC) designs. The program
#' can run interactively,  allowing the user to enter outcomes after each cohort
#' has been recruited,  or via simulation to assess operating characteristics.
#' 
#' \code{bcrm} implements a Bayesian continual reassessment method (CRM)
#' (O'Quigley \emph{et al.},  1990); an adaptive design in which cohorts of
#' patients are sequentially recruited into a Phase I trial. A binary toxicity
#' outcome is assumed (e.g. Dose Limiting Toxicity / No Dose Limiting
#' Toxicity). The current cohort are given a dose "closest" to the specified
#' target toxicity level,  as estimated from the posterior distributions of
#' toxicity at each dose level from the patients thus far recruited. If
#' \code{pointest="mean"} then the posterior mean probability of toxicity is
#' used to choose the next dose. If \code{pointest="plugin"},  however,  the
#' posterior mean of the model parameter(s) is plugged-into the functional form
#' of the dose-toxicity model. To implement an EWOC design (Babb \emph{et al.}, 
#' 1998),  \code{pointest} should be a quantile,  \emph{q},  between 0 and 0.5.
#' The posterior distribution of the MTD (the dose in which the probability of
#' toxicity is equal to the target toxicity) is then calculated and the next
#' patient is given dose closest to the \emph{q}th quantile of the MTD
#' distribution.
#' 
#' Alternatively,  escalation can be based on intervals of toxicity from the
#' posterior distribution using a loss function,  see Neuenschwander \emph{et
#' al.},  2008. To implement this approach,  the user should specify the
#' cutpoints of the toxicity intervals using \code{tox.cutpoints} and the
#' associated losses using \code{loss}.
#' 
#' The possible choice of dose-toxicity model can be specified using \code{ff}, 
#' and includes the 1-parameter hyberbolic tangent,  logistic or power "working
#' models",  and the 2-parameter logistic model as follows: \describe{
#' \item{Hyperbolic
#' Tangent}{\deqn{p(Tox|d^*)=\left[(tanh(d^*)+1)/2\right]^\alpha}{p(Tox|d*)=[(tanh(d*)+1)/2]^\alpha}}
#' \item{Logistic (1-parameter)}{\deqn{p(Tox|d^*)=\frac{\exp(3+\alpha
#' d^*)}{1+\exp(3+\alpha d^*)}}{p(Tox|d*)=exp(3+\alpha d*)/(1+exp(3+\alpha
#' d*))}} \item{Power}{\deqn{p(Tox|d^*)={d^*}^\alpha}{p(Tox|d*)=d*^\alpha}}
#' \item{Logistic
#' (2-parameter)}{\deqn{p(Tox|d^*)=\frac{\exp(\log(\alpha_1)+\alpha_2
#' d^*)}{1+\exp(\log(\alpha_1)+\alpha_2 d^*)}}{
#' p(Tox|d*)=exp(log(\alpha_1)+\alpha_2 d*)/(1+exp(log(\alpha_1)+\alpha_2
#' d*))}} } where \eqn{\alpha>0} is the single positive-valued parameter for
#' the 1-parameter models,  and \eqn{\log(\alpha_1)}{log(\alpha_1)} and
#' \eqn{\alpha_2>0} are the intercept and slope parameters of the 2-parameter
#' model.
#' 
#' The standardised doses,  \eqn{d^*}{d*},  are specified by the user using
#' \code{sdose},  or alternatively the prior probability of toxicity at each
#' dose level is specified using \code{p.tox0}. If the latter is used,  then the
#' standardised doses are calculated using the inverse of the functional form
#' and a plug-in estimate of the prior mean or median,  as specified in
#' \code{sdose.calculate},  as follows \deqn{d^* = f^{-1}(\code{p.tox0}, \alpha=
#' a)}{d* = f^{-1}(p.tox0, \alpha= a)} where \eqn{f^{-1}} is the the inverse of
#' the chosen functional form,  and the parameter(s) of the model are set equal
#' to \eqn{a},  either the prior mean or median of \eqn{\alpha}.
#' 
#' Data that have already been accrued can be entered using the \code{data}
#' argument. A constrained CRM design can be implemented using
#' \code{constrain=TRUE},  in which case dose-skipping is prohibited (i.e. the
#' next cohort can only be dosed up to one dose level above the current
#' cohort). If a constrained model is used then the starting dose must be
#' specified using \code{start}. Alternatively,  if data have already been
#' accrued,  then the dose level of the last recruited patient determines the
#' constraint for the next patient.
#' 
#' The prior is set using \code{prior.alpha}. For example
#' \code{prior.alpha=list(1, 1, 1)} specifies a Gamma prior with shape and scale
#' parameters both equal to one (\emph{i.e.} an Exponential(1) distribution), 
#' whilst \code{prior.alpha=list(2, 0, 10)} specifies a Uniform(0, 10) prior.
#' 
#' To specify a fixed maximum sample size of size \code{m} use
#' \code{stop=list(nmax=m)}. Alternatively,  the trial can stop after \code{m2}
#' patients have been treated at the current MTD estimate,  by setting
#' \code{stop=list(nmtd=m2)}.
#' 
#' To implement a safety constraint as specified in Zohar and Chevret (2001)
#' specify \code{stop=list(safety=p)},  where the trial is stopped if the
#' posterior probability that the lowest dose is greater than the target
#' toxicity probability is greater than \code{p}.
#' 
#' To stop the trial when the MTD estimate is within a certain level of
#' precision,  use \code{stop=list(precision=c(l, u))},  where \code{l} and
#' \code{u} are the lower and upper percentage points that the MTD 95\%
#' credible intervals for the risk of toxicity should lie within. Finally,  to
#' prevent the trial stopping too early using these rules,  the argument
#' \code{stop=list(nmin=m3)} can be used to ensure the sample size is greater
#' than or equal to \code{m3}. Stopping rules can be used on their own or in
#' combination.
#' 
#' The trial can be run interactively using \code{simulate=FALSE},  where the
#' user enters the outcomes for each new cohort,  or as a simulation study when
#' \code{simulate=TRUE}.
#' 
#' The default calculations use exact methods (\code{method="exact"}) to
#' calculate the mean and quantiles for the posterior distributions. There are
#' three choices for MCMC calculations: \code{method="rjags"}, 
#' \code{method="BRugs"} or \code{method="R2WinBUGS"}. The first uses the JAGS
#' software,  the second uses OpenBUGS,  whilst the latter uses WinBUGS. To
#' implement these methods,  users require one or more of these packages to be
#' installed on their system.
#' 
#' A simulated \code{bcrm} design can be compared with the standard 3+3
#' rule-based method,  see \code{\link{threep3}} for more details.
#' 
#' @param stop A list of stopping rules for the trial. One or more of the
#' following options should be specified \describe{ \item{list("nmax")}{ The
#' maximum sample size of the trial} \item{list("safety")}{Stops the trial if
#' the posterior probability of the lowest dose being above the target toxicity
#' level is greater than this number} \item{list("nmtd")}{ The maximum number
#' to be treated at final maximum tolerated dose (MTD) estimate,  \emph{i.e.} if
#' the next recommended dose has already been administered to \code{nmtd}
#' patients then the trial will stop} \item{list("precision")}{ A vector of the
#' lower and upper percentage points that the MTD 95\% credible intervals for
#' the risk of toxicity should lie within} \item{list("nmin")}{ The minimum
#' sample size of the trial. To be used in conjunction with \code{nmtd} or
#' \code{precision}} }
#' @param data A named data frame giving information about dose and toxicity
#' from previously recruited patients. If missing,  then it is assumed that no
#' data have thus far been collected. Contains the following variables:
#' \describe{ \item{list("patient")}{ Recruited patient numbers, 
#' \code{1, ..., n}} \item{list("dose")}{ Dose levels of recruited patients, 
#' ranging from \code{1, ..., k}} \item{list("tox")}{ An indicator variable for
#' each patient (1=toxicity,  0=no toxicity)} }
#' @param p.tox0 A vector of length \code{k} listing the prior probabilities of
#' experiencing the outcome at each dose level \code{1, ...k} (CRM "skeleton").
#' The standardised dose levels are formed from these probabilities using the
#' inverse of the functional form,  with a plug-in estimate for the prior mean
#' or median of alpha,  as specified in \code{ff},  \code{prior.alpha} and
#' \code{sdose.calculate}. Alternatively standardised dose levels can be given
#' directly using \code{sdose}.
#' @param sdose A vector of length \code{k} listing the standardised doses to
#' be used in the CRM model. Only required if \code{p.tox0} is missing.
#' @param dose Optional vector of length \code{k} of actual doses for plotting
#' purposes
#' @param ff A string indicating the functional form of the dose-response
#' curve. Options are \describe{ \item{ht}{ 1-parameter hyperbolic tangent}
#' \item{logit1}{ 1-parameter logistic} \item{power}{ 1-parameter power}
#' \item{logit2}{ 2-parameter logistic} }
#' @param prior.alpha A list of length 3 containing the distributional
#' information for the prior. The first element is a number from 1-4 specifying
#' the type of distribution. Options are \enumerate{ \item Gamma(a, b),  where
#' a=shape,  b=scale: mean=a*b,  variance=a*b*b \item Uniform(a, b),  where a=min, 
#' b=max \item Lognormal(a, b),  where a=mean on the log scale,  b=variance on the
#' log scale \item Bivariate Lognormal(a, b),  where a=mean vector on the log
#' scale,  b=Variance-covariance matrix on the log scale. This prior should be
#' used only in conjunction with a two-parameter logistic model.  } The second
#' and third elements of the list are the parameters a and b,  respectively.
#' @param cohort The size of each cohort of patients that are sequentially
#' recruited to the trial. Defaults to 3
#' @param target.tox The target toxicity probability. Defaults to 1/3.
#' @param constrain Should a dose-skipping constraint be placed on the
#' escalation procedure,  as imposed by a modified CRM? Defaults to TRUE.
#' @param sdose.calculate What plug-in estimate of the prior alpha should be
#' used to calculate the standardised doses? Options are \code{"mean"}
#' (default) or \code{"median"}. Only required if \code{sdose} is missing.
#' @param pointest Which summary estimate of the posterior distribution should
#' be used to choose the next dose. Options are \code{"plugin"} (default) where
#' the posterior mean of the model parameter(s) is plugged into the function
#' form to obtain estimates of toxicity,  or \code{"mean"} where the posterior
#' mean probabilities of toxicity are directly used. Alternatively,  a number
#' between 0 and 1 can be specified representing the quantile of the maximum
#' tolerated dose (MTD) posterior distribution (e.g. 0.5 specifies the
#' posterior median). This produces an Escalation With Overdose Control (EWOC)
#' design if the quantile is less than 0.5 (see details). Currently,  EWOC
#' designs must be fit using MCMC methods.
#' @param tox.cutpoints A vector of cutpoints for toxicity intervals if these
#' are to be used to choose next dose. Defaults to NULL. For example
#' Underdosing [0, 0.2],  Target dosing (0.2,  0.35],  Excessive toxicity (0.35, 
#' 0.60],  Unacceptable toxicity (0.60,  1.00] set
#' \code{tox.cutpoints=c(0.2, 0.35, 0.60)}.
#' @param loss A vector of length \code{length(tox.cutpoints)+1} specifying the
#' losses associated with each toxicity interval,  e.g. Underdosing = 1,  Target
#' dosing =0,  Excessive toxicity=1,  Unacceptable toxicity=2
#' @param start Dose level used at the beginning of the trial. Required if
#' \code{constrain=TRUE}.
#' @param simulate Should a simulation be conducted to assess operating
#' characteristics? Defaults to \code{TRUE}. If \code{FALSE},  a single CRM
#' trial is run interactively,  allowing the user to input outcomes after each
#' cohort is recruited.
#' @param nsims Number of simulations to perform if \code{simulate==TRUE}
#' (defaults to 1).
#' @param truep A vector of length k giving the true probabilities of the
#' outcome (toxicity) at each dose level \code{1, ..., k} in order to simulate
#' data. Only required if \code{simulate=TRUE}
#' @param threep3 Should operating characteristics of a standard 3+3 rule-based
#' design be calculated,  for comparison with \code{bcrm} design? Defaults to
#' \code{FALSE}. Only used in a simulation study,  i.e. when
#' \code{simulate=TRUE}.
#' @param method Optimisation method: options are \code{"exact"} (the default), 
#' \code{"rjags"},  \code{"BRugs"},  or \code{"R2WinBUGS"}.
#' @param burnin.itr Number of burn-in iterations (default 2000).
#' @param production.itr Number of production iterations (default 2000).
#' @param bugs.directory Directory that contains the WinBUGS executable if
#' \code{method="R2WinBUGS"}. Defaults to "C:/Program Files/WinBUGS14/".
#' @param plot Should the dose-response curve be plotted after each cohort has
#' been entered? Defaults to FALSE.
#' @param seed Integer defining the state of the random number generator to
#' allow reproducible results. The default is to not specify a seed.
#' @param quietly Should the simulation number output be suppressed when
#' running bcrm?  Defaults to FALSE.
#' @param file File name where the dose-response plots are stored,  in a pdf
#' format. The program will amend the current sample size to the end of the
#' file name.
#' @param N Final sample size (deprecated). To be replaced with \code{stop} in
#' future versions.
#' @param tox (Deprecated). A vector of length \code{k} listing the number of
#' patients who have experienced the outcome (toxicity) at each dose level
#' \code{1, ..., k}.
#' @param notox (Deprecated). A vector of length \code{k} listing the number of
#' patients who have not experienced the outcome (toxicity) at each dose level
#' \code{1, ..., k}.
#' @return \code{bcrm} returns an object of class "bcrm" or "bcrm.sim"; the
#' latter occurring when a simulation has been conducted
#' (\code{simulate=TRUE}). The function \code{\link{print}} (i.e.
#' \code{\link{print.bcrm}} or \code{\link{print.bcrm.sim}}) can be used to
#' obtain summary information about the design used,  the data observed,  current
#' posterior estimates of toxicity,  and the next recommended dose level.
#' 
#' An object of class "bcrm" is a list with the following components:
#' \item{dose}{Range of doses} \item{sdose}{Standardised doses} \item{tox}{A
#' vector of length \code{k} listing the number of patients who have
#' experienced the outcome (toxicity) at each dose level \code{1, ..., k}}
#' \item{notox}{A vector of length \code{k} listing the number of patients who
#' have not experienced the outcome (toxicity) at each dose level
#' \code{1, ..., k}} \item{ndose}{A list of lists containing for each cohort the
#' components \code{ndose},  the dose level recommended for the next patient, 
#' \code{est},  the estimated probabilities of toxicity using the chosen metric
#' (e.g. plugin,  mean,  quantile),  \code{mean},  the posterior mean probability
#' of toxicity at each dose,  \code{sd},  the posterior standard deviation for
#' probability of toxicity at each dose,  \code{quantiles},  the posterior
#' quantiles for probability of toxicity at each dose. This information is only
#' provided for cohorts recruited subsequent to any data given using \code{tox}
#' and \code{notox}. The first component relates to the prior information.}
#' \item{constrain}{Whether a constrained CRM design was used} \item{start}{The
#' starting dose for the latest run of the model if \code{constrain=TRUE}}
#' \item{target.tox}{The target toxicity level} \item{ff}{A number from 1-4
#' identifying the functional form,  1 = Hyperbolic tangent,  2 = 1-parameter
#' logistic,  3 = Power,  4 = 2-parameter logistic} \item{method}{The calculation
#' method used} \item{pointest}{The summary estimate used to choose the next
#' dose,  see \code{pointest}} \item{prior.alpha}{Information about the prior
#' used for \eqn{alpha},  see \code{prior.alpha}} \item{data}{A data frame with
#' variables `patient',  `dose' and `tox' listing the dose levels and outcomes
#' of all patients in the trial}
#' 
#' An object of class "bcrm.sim" is a list of length \code{nsims}. Each
#' component is itself a list with components similar to those obtained from a
#' "bcrm" object. The print function,  \code{\link{print.bcrm.sim}} should be
#' used to obtain operating characteristics from the simulation.
#' @note Currently,  the re-parameterisation of the two-parameter model proposed
#' by (Babb \emph{et al.},  1998) is not implemented. Therefore,  users wishing
#' to implement an EWOC design should check whether their choice of prior for
#' the model parameter(s) translates to a sensible prior for the MTD
#' distribution before they implement the design. For example \preformatted{
#' prior.alpha <- list(1, 1, 1) ff <- "ht" target.tox <- 0.2
#' samples.alpha <- getprior(prior.alpha, 2000)
#' mtd <- find.x(ff, target.tox, alpha=samples.alpha) hist(mtd) }
#' 
#' One-parameter models are designed as working models only,  and should not be
#' used with an escalation strategy based on intervals of the posterior
#' probabilities of toxicity.
#' @author Michael Sweeting \email{michael.sweeting@@leicester.ac.uk}
#' (University of Leicester,  UK) and Graham Wheeler
#' (graham.wheeler@@ucl.ac.uk),  drawing on code originally developed by J. Jack
#' Lee and Nan Chen,  Department of Biostatistics,  the University of Texas M. D.
#' Anderson Cancer Center
#' @seealso \code{\link{print.bcrm}},  \code{\link{print.bcrm.sim}}, 
#' \code{\link{plot.bcrm}},  \code{\link{plot.bcrm.sim}},  \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#' 
#' O'Quigley J.,  Pepe M.,  Fisher L. Continual reassessment method: a practical
#' design for phase I clinical trials in cancer. \emph{Biometrics} (1990) 46:
#' 33--48.
#' 
#' Babb J.,  Rogatko A.,  Zacks S. Cancer phase I clinical trials: efficient dose
#' escalation with overdose control. \emph{Statistics in Medicine} (1998) 17:
#' 1103--1120.
#' 
#' Neuenschwander B.,  Branson M.,  Gsponer T. Critical aspects of the Bayesian
#' approach to phase I cancer trials. \emph{Statistics in Medicine} (2008) 27:
#' 2420--2439.
#' 
#' Zohar S.,  Chevret S. The continual reassessment method: comparison of
#' Bayesian stopping rules for dose-ranging studies. \emph{Statistics in
#' Medicine} (2001) 20: 2827--2843.
#' @examples
#' 
#' ## Dose-escalation cancer trial example as described in Neuenschwander et al 2008.
#' ## Pre-defined doses
#' dose <- c(1, 2.5, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250)
#' ## Pre-specified probabilities of toxicity
#' ## [dose levels 11-15 not specified in the paper,  and are for illustration only]
#' p.tox0 <- c(0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.100, 0.170, 0.300, 0.400, 0.500, 0.650
#'   , 0.800, 0.900)
#' ## Data from the first 5 cohorts of 18 patients
#' data <- data.frame(patient=1:18, dose=rep(c(1:4, 7), c(3, 4, 5, 4, 2)), tox=rep(0:1, c(16, 2)))
#' ## Target toxicity level
#' target.tox <- 0.30
#' 
#' ## A 1-parameter power model is used,  with standardised doses calculated using 
#' ## the plug-in prior median
#' ## Prior for alpha is lognormal with mean 0 (on log scale) 
#' ## and standard deviation 1.34 (on log scale)
#' ## The recommended dose for the next cohort if posterior mean is used
#' Power.LN.bcrm <- bcrm(stop=list(nmax=18), data=data, p.tox0=p.tox0, dose=dose
#'   , ff="power", prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=FALSE
#'   , sdose.calculate="median", pointest="mean")
#' print(Power.LN.bcrm)
#' plot(Power.LN.bcrm)
#' 
#' ## Simulate 10 replicate trials of size 36 (cohort size 3) using this design 
#' ## with constraint (i.e. no dose-skipping) and starting at lowest dose
#' ## True probabilities of toxicity are set to pre-specified probabilities (p.tox0) 
#' Power.LN.bcrm.sim <- bcrm(stop=list(nmax=36), p.tox0=p.tox0, dose=dose, ff="power"
#'   , prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=TRUE
#'   , sdose.calculate="median", pointest="mean", start=1, simulate=TRUE, nsims=10, truep=p.tox0)
#' print(Power.LN.bcrm.sim)
#' plot(Power.LN.bcrm.sim)
#' 
#' ## Comparing this CRM design with the standard 3+3 design 
#' ## (only considering the first 12 dose levels)
#' \dontrun{
#' Power.LN.bcrm.compare.sim <- bcrm(stop=list(nmax=36), p.tox0=p.tox0[1:12], dose=dose[1:12]
#'   , ff="power", prior.alpha=list(3, 0, 1.34^2), target.tox=target.tox, constrain=TRUE
#'   , sdose.calculate="median", pointest="mean", start=1, simulate=TRUE, nsims=50
#'   , truep=p.tox0[1:12], threep3=TRUE)
#' print(Power.LN.bcrm.compare.sim, threep3=TRUE)
#' plot(Power.LN.bcrm.compare.sim, threep3=TRUE)
#' }
#' 
#' ## A 2-parameter model,  using priors as specified in Neuenschwander et al 2008.
#' ## Posterior mean used to choose the next dose
#' ## Standardised doses using reference dose,  250mg
#' sdose <- log(dose/250)
#' ## Bivariate lognormal prior for two parameters
#' mu <- c(2.15, 0.52)
#' Sigma <- rbind(c(0.84^2, 0.134), c(0.134, 0.80^2))
#' ## Using rjags (requires JAGS to be installed)
#' TwoPLogistic.mean.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
#'   , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
#'   , constrain=FALSE, pointest="mean", method="rjags")
#' print(TwoPLogistic.mean.bcrm)
#' plot(TwoPLogistic.mean.bcrm)
#' 
#' ## A 2-parameter model,  using an EWOC design with feasibility bound (MTD quantile) 
#' ## of 0.25 to choose the next dose
#' ## Using rjags (requires JAGS to be installed)
#' \dontrun{
#' TwoPLogistic.EWOC0.25.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose, dose=dose
#'     , ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox, constrain=FALSE
#'     , pointest=0.25, method="rjags")
#' print(TwoPLogistic.EWOC0.25.bcrm)
#' plot(TwoPLogistic.EWOC0.25.bcrm)
#' }
#' 
#' ## A 2-parameter model,  using a loss function based on intervals of toxicity to choose
#' ## the next dose
#' ## Using rjags (requires JAGS to be installed)
#' \dontrun{
#' ## Toxicity cut-points
#' tox.cutpoints <- c(0.2, 0.35, 0.6)
#' ## Losses associated with toxicity intervals 
#' ## [0, 0.2]=1,  (0.2, 0.35]=0,  (0.35, 0.6]=1,  (0.6, 1]=2
#' loss <- c(1, 0, 1, 2)
#' TwoPLogistic.tox.intervals.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
#'   , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
#'   , constrain=FALSE, tox.cutpoints=tox.cutpoints, loss=loss, method="rjags")
#' print(TwoPLogistic.tox.intervals.bcrm)
#' plot(TwoPLogistic.tox.intervals.bcrm)
#' ## Greater loss associated with overdosing and unacceptable toxicity
#' ## [0, 0.2]=1,  (0.2, 0.35]=0,  (0.35, 0.6]=2,  (0.6, 1]=4
#' loss2 <- c(1, 0, 2, 4)
#' TwoPLogistic.tox.intervals.2.bcrm <- bcrm(stop=list(nmax=18), data=data, sdose=sdose
#'   , dose=dose, ff="logit2", prior.alpha=list(4, mu, Sigma), target.tox=target.tox
#'   , constrain=FALSE, tox.cutpoints=tox.cutpoints, loss=loss2, method="rjags")
#' print(TwoPLogistic.tox.intervals.2.bcrm)
#' plot(TwoPLogistic.tox.intervals.2.bcrm)
#' }
#' 
#' 
#' @export bcrm
bcrm <- function(stop=list(nmax=NULL, nmtd=NULL, precision=NULL,
                           nmin=NULL, safety=NULL),
                 data=NULL, p.tox0=NULL, sdose=NULL, dose=NULL, ff,
                 prior.alpha, cohort=3, target.tox, constrain=TRUE,
                 sdose.calculate="mean", pointest="plugin",
                 tox.cutpoints=NULL, loss=NULL,
                 start=NULL, simulate=FALSE, nsims=1, truep=NULL,
                 threep3=FALSE, method="exact", burnin.itr=2000,
                 production.itr=2000,
                 bugs.directory="c:/Program Files/WinBUGS14/",
                 plot=FALSE, seed=NULL,  quietly=FALSE,  file=NULL,
                 N, tox, notox){
  
  # Checks of argument inputs  
  if(missing(N) & is.null(stop$nmax) & is.null(stop$nmtd) & is.null(stop$precision))
    stop("At least one stopping rule must be provided using the stop argument")
  if(!missing(N)){
    stop$nmax <- N
    warning("N is deprecated and users should now use the stop argument to specify the maximum sample size")
  }
  if(!missing(tox) | !missing(notox)){
    stop("tox and nontox arguments are deprecated and users should now use the data argument to specify previous data,  see ?bcrm")
  }
  if(!(length(stop$precision) %in% c(0, 2))) stop("stop$precision must be a vector of length two")
  if(!is.null(stop$nmax) & !is.null(stop$nmin)) {if(stop$nmin>stop$nmax) stop("stop$nmin must be less than stop$nmax")}
  if(!is.null(stop$safety)){
    if(stop$safety<=0 | stop$safety>=1) stop("stop$safety must be a probability between 0 and 1")
  }
  if(missing(p.tox0) & missing(sdose)) stop("Either p.tox0 or sdose must be specified")
  if(!missing(p.tox0) & !missing(sdose)) stop("Only one of p.tox0 and sdose must be specified")
  if(sdose.calculate!="mean" & sdose.calculate!="median") stop("sdose.calculate must be either `mean' or `median'")
  if((is.character(pointest) & pointest!="mean" & pointest!="plugin") | is.numeric(pointest) & (pointest<0 | pointest>1)) stop("pointest must be either `plugin',  `mean' or an EWOC feasibility quantile between 0 and 1")
  if(is.numeric(pointest) & method=="exact") stop("EWOC design must be fitted using MCMC methods")
  if(!is.null(tox.cutpoints) & method=="exact") stop("Escalation based on toxicity intervals must be fit using MCMC. Please specify either method=`rjags',  method='BRugs' or method='R2WinBUGS'")
  
  if(simulate & is.null(truep)) stop("truep must be specified if simulating data")
  if(simulate){ 
    plot <- FALSE
  }
  if(!(method %in% c("exact", "rjags", "BRugs", "R2WinBUGS"))) stop("method must be either `exact',  `rjags',  `BRugs' or `R2WinBUGS'")
  ## Check to see if ff is one of "ht", "logit1", "power", "logit2"
  if((!ff %in% c("ht", "logit1", "power", "logit2"))) stop("ff must be one of `ht',  `logit1',  `power' or `logit2'")
  
  if(ff=="logit2" & method=="exact") warning("Exact method slow for 2-parameter model,  suggest using rjags (MCMC)")
  if(constrain & is.null(start) & is.null(data)) stop("A starting dose level must be specified using `start' if constrain==TRUE")
  if((!is.null(tox.cutpoints) & is.null(loss)) | (is.null(tox.cutpoints) & !is.null(loss))) stop("Both tox.cutpoints and loss must be specified to conduct escalation based on toxicity intervals")
  if(!is.null(tox.cutpoints) & length(loss)!=length(tox.cutpoints)+1) stop("The number of losses must be one more than the number of cutpoints")
  if(!is.null(tox.cutpoints)) pointest <- NULL
  if(!is.null(tox.cutpoints) & ff!="logit2") warning("One-parameter models are designed as working models only,  and should not be used with an escalation strategy based on intervals of the posterior probabilities of toxicity")
  
  if(ff=="logit2" & (length(prior.alpha[[2]])<2 | length(prior.alpha[[3]])<2)) stop("second and third components of `prior.alpha' must be vectors of size 2")
  
  # Set seed if specified
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(missing(sdose)){
    alpha.prior.plug <- if(prior.alpha[[1]]==1){
      ifelse(sdose.calculate=="mean", prior.alpha[[2]]*prior.alpha[[3]], median(getprior(prior.alpha,  10000)))
    } else if(prior.alpha[[1]]==2){
      0.5*(prior.alpha[[2]]+prior.alpha[[3]])
    } else if(prior.alpha[[1]]==3){
      ifelse(sdose.calculate=="mean", exp(prior.alpha[[2]]+prior.alpha[[3]]/2), exp(prior.alpha[[2]]))
    } else if(prior.alpha[[1]]==4){
      if(sdose.calculate=="mean"){exp(prior.alpha[[2]]+diag(prior.alpha[[3]])/2)} else {exp(prior.alpha[[2]])}
    }     
    sdose <- find.x(ff, p.tox0, alpha=alpha.prior.plug)
  }
  # Checking that length of truep in simulation study is same length as sdose
  if(simulate & length(truep)!=length(sdose)) stop("Length of truep must be the same as length of sdose or p.tox0")  
  # Checking that length of dose (if specified) is same length as sdose
  if(length(dose)>0 & length(dose)!=length(sdose)) stop("Length of dose must be the same as length of sdose or p.tox0")
  # Check data contains the correct variables
  if(!is.null(data)){
    if(any(!(c("patient", "dose", "tox") %in% names(data)))) stop("data must have variables named 'patient',  'dose' and 'tox'")
    data <- data[order(data$patient), ]
    if(any(data$patient != 1:dim(data)[1])) stop("'patient' variable in data must be an ascending vector of positive integers")
    if(any(!(data$tox %in% c(0, 1)))) stop("'tox' variable in data must be a vector of zeros (no toxicity) and ones (toxicity)")
    if(any(!(data$dose %in% 1:length(sdose)))) stop(paste("'dose' variable in data must contain the dose levels (1 to ", length(sdose), ")", sep=""))   
    if(!is.null(start)) warning("start no longer needs to be specified   if data is given; using last recruited patient as current dose")
    start <- as.numeric(data$dose[1])
  }
  
  ## Cannot calculate quantiles if method=="exact" so stop$precision and stop$safety cannot be used
  if(method=="exact" & (!is.null(stop$precision) | !is.null(stop$safety))){ stop("exact method cannot be used with stop$precision or stop$safety,  please use MCMC instead")}
  
  ## Allowing access to fast exact computation if method="exact" & simulate=TRUE & stopping rules do not depend on posterior quantiles
  if(method=="exact" & simulate & is.null(stop$precision)){ method <- "exact.sim" }
  ## If method=="exact" and a two-parameter model is fitted,  only relevant escalation posterior quantities are calculated (apart from trial end)
  if(method=="exact" & ff=="logit2"){ 
    method <- "exact.sim"
    if(plot){
      plot <- FALSE
      warning("Plot function not available for exact computation of 2-parameter model")
    }
  }
  
  k <- length(sdose)
  
  # Set up
  
  if(is.null(stop$safety)){
    quantiles <- c(0.025, 0.25, 0.50, 0.75, 0.975)
  } else {
    quantiles <- sort(unique(c(0.025, 0.25, 0.50, 0.75, 0.975, 1-stop$safety)))
  }
  
  
  if (simulate){
    new.tox <- new.notox <- rep(0, k)
    current <- start - 1
    
    alpha <- switch(method
                    , rjags=getprior(prior.alpha,  10000)
                    , BRugs=getprior(prior.alpha,  10000)
                    , R2WinBUGS=getprior(prior.alpha, 10000)
                    , exact=Posterior.exact(new.tox, new.notox, sdose, ff, prior.alpha)
                    , exact.sim=Posterior.exact.sim(new.tox, new.notox, sdose, ff, prior.alpha, pointest)
    )
    prior.ndose <- switch(method
                          , rjags=nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss, quantiles)
                          , BRugs=nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss, quantiles)
                          , R2WinBUGS=nextdose(alpha, sdose, ff, target.tox, constrain, pointest, current, tox.cutpoints, loss, quantiles)
                          , exact=nextdose.exact(alpha, sdose, ff, target.tox, constrain, pointest, current)
                          , exact.sim=nextdose.exact.sim(alpha, sdose, ff, target.tox, constrain, pointest, current)
    )
    ndose <- prior.ndose
    
    out <- lapply(1:nsims, simFun)
    class(out) <- "bcrm.sim"
    if (threep3){
      if(length(truep) > 9){
        message("Warning: Calculation of all 3+3 designs may take a long time,  continue?  ")
        yn  <-  readline()
        if (yn!="y" & yn=="Y") return(out)
      }
      message("  Calculating operating characteristics of a standard 3+3 trial for comparison...")
      out[[1]]$threep3 <- threep3(truep, start=start, dose=dose) 
    }
  } else { # not simulating
    new.tox <- as.numeric(xtabs(tox ~ factor(dose, levels=1:k), data=data))
    new.notox <- as.numeric(xtabs((1 - tox) ~ factor(dose,levels=1:k), data=data))
    current <- as.numeric(data$dose[dim(data)[1]])
    alpha <-switch(method
                   ,rjags=Posterior.rjags(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr)
                   ,BRugs=Posterior.BRugs(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr)
                   ,R2WinBUGS=Posterior.R2WinBUGS(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr,bugs.directory)
                   ,exact=Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha)
                   ,exact.sim=Posterior.exact.sim(new.tox,new.notox,sdose,ff,prior.alpha,pointest))
    
    ncurrent <- sum(new.tox + new.notox)
    
    newdata <- data
    
    prior.ndose<-switch(method
                        ,rjags=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss,quantiles)
                        ,BRugs=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss,quantiles)
                        ,R2WinBUGS=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss,quantiles)
                        ,exact=nextdose.exact(alpha,sdose,ff,target.tox,constrain,pointest,current)
                        ,exact.sim=nextdose.exact.sim(alpha,sdose,ff,target.tox,constrain,pointest,current)
    )
    ndose <- prior.ndose
    
    stopped <- stop.check(stop, ncurrent, ndose, new.tox, new.notox, simulate)
    ndose <- stopped$ndose

    results <- list(dose=dose, sdose=sdose, tox=new.tox, notox=new.notox, ndose=list(ndose), constrain=constrain, start=start, target.tox=target.tox, ff=ff, method=method, pointest=pointest, tox.cutpoints=tox.cutpoints, loss=loss, prior.alpha=prior.alpha, data=data)
    class(results) <- "bcrm"
    
    if(plot){
      plot(results, file)
    }
    
    while(!stopped$stop){
      current <- ndose[[1]]
      ncurrent <- ncurrent + cohort
      
      interact<-crm.interactive(new.tox, new.notox, ncurrent,
                                cohort, ndose, dose)
      if(interact$bk == TRUE){
        results$data<-newdata
        return(results)
      }
      y <- interact$y
      new.tox <- interact$tox
      new.notox <- interact$notox
      current <- interact$ans
      
      currentdata <- data.frame(patient = (ncurrent - cohort+1):ncurrent,
                                dose = rep(current, cohort),
                                tox = y)
      newdata <- rbind(newdata, currentdata)
      
      alpha<-switch(method
                    ,rjags=Posterior.rjags(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr)
                    ,BRugs=Posterior.BRugs(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr)
                    ,R2WinBUGS=Posterior.R2WinBUGS(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr,bugs.directory)
                    ,exact=Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha)
                    ,exact.sim=Posterior.exact.sim(new.tox,new.notox,sdose,ff,prior.alpha,pointest)
      )
      ndose<-switch(method
                    ,rjags=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss,quantiles)
                    ,BRugs=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss,quantiles)
                    ,R2WinBUGS=nextdose(alpha,sdose,ff,target.tox,constrain,pointest,current,tox.cutpoints,loss,quantiles)
                    ,exact=nextdose.exact(alpha,sdose,ff,target.tox,constrain,pointest,current)
                    ,exact.sim=nextdose.exact.sim(alpha,sdose,ff,target.tox,constrain,pointest,current)
      )
      
      stopped<-stop.check(stop,ncurrent,ndose,new.tox,new.notox,simulate)
      ndose<-stopped$ndose ## update ndose in case no doses are deemed safe

      results$tox<-new.tox
      results$notox<-new.notox
      results$ndose[[length(results$ndose)+1]]<-ndose
      results$data<-newdata

      if(plot){
        plot(results,file)
      }
    }

    out <- results
  }
  
  return(out)
}



# ----------------------------------------------------------------------
#   CRM INTERACTIVE. Conduct a CRM trial interactively allowing user to specify outcomes after each cohort has been recruited
#     tox         --> number of successes (toxicities)
#     notox       --> number of failed patients (no-toxicities)
#  ncurrent --> Current no. of patients in trial
#  cohort   --> Cohort size
#    ndose    --> Proposed dose level for next cohort
#  dose     --> Dose labels for each dose
# ----------------------------------------------------------------------
crm.interactive <- function(tox, notox, ncurrent, cohort, ndose, dose){
  k  <-  length(tox)
  repeat {
    # y.j is toxicity information from the treatment of patient j in the current
    # cohort of patients at dose level ndose
    y <- vector()
    
    if(is.null(dose)){
      message("\n\n RECOMMENDED DOSE LEVEL FOR PATIENTS ", ncurrent-cohort+1, " to ", ncurrent,  "IS:",  ndose[[1]])
      ans  <-  get.dose.level(k, ncurrent, cohort)     
    } else {
      message("\n\n RECOMMENDED DOSE FOR PATIENTS ", ncurrent-cohort+1, " to ", ncurrent,  "IS:",  dose[ndose[[1]]])
      ans  <-  get.dose(dose, ncurrent, cohort)     
    }
    if (ans==-2)
      ans  <-  ifelse(is.null(dose), ndose[[1]], dose[ndose[[1]]])
    if (ans==0) {
      message("\n\n EXIT AND RETURN THE RESULTS SO FAR?")            
      message("\n DO YOU REALLY WANT TO EXIT ? (Y/N)  ")
      yn  <-  readline()
      if (yn=="y" || yn=="Y") break
      else                    next
    }  
    
    for(j in (ncurrent-cohort+1):ncurrent){
      message("ENTER TOXICITY DATA FOR PATIENT", j, "(1=TOX,  0=NO TOX): ")
      y  <-  c(y, get.answer())               
    }
    # give the user a last chance to modify the treatment and outcome
    message("\n\t\t ENTERED VALUES:")
    if(is.null(dose)){
      message("\n DOSE LEVEL ...",  ans)
    } else {
      message("\n DOSE ...",  ans)
    }
    message("\n TOXICITIES ....", y)
    message("\n PRESS `RETURN' IF OK,  OR ANY OTHER KEY TO ENTER NEW VALUES ")
    key  <-  readline()
    if (nchar(key)==0) break
  }
  if (ans==0)  
    return(list(tox=tox,  notox=notox,  y=y,  bk=TRUE,  ans=ans))
  if(!is.null(dose)){
    ans <- which(dose==ans)
  }
  notox[ans] <- notox[ans]+(cohort-sum(y))
  tox[ans] <- tox[ans]+sum(y)
  return(list(tox=tox,  notox=notox,  y=y,  bk=FALSE,  ans=ans))
}


# ----------------------------------------------------------------------
#     User inputs toxicity info from the treatment at a dose level:
#   Must be either 0 and 1
# ----------------------------------------------------------------------
get.answer  <-  function() {   
  repeat {
    ans  <-  as.numeric(readline())
    if (is.na(ans)) next
    if (ans != floor(ans)) next
    if (ans>=0 & ans<=1) return(ans)
  }
}

# ----------------------------------------------------------------------
#     ask the user to input an integer number <= n.
#     the number is checked to belong to [-1, n] and also to be
#     an integer.  `ENTER' returns to the caller with no action, 
#
#     n - biggest number to be accepted
#    ncurrent - patient no. for current patient
#     cohort - cohort size
# ----------------------------------------------------------------------
get.dose.level  <-  function( n , ncurrent, cohort) {
  repeat {
    message("\n\n ENTER DOSE LEVEL BETWEEN 1 AND ",  n, " FOR PATIENTS ", ncurrent-cohort+1, " to ", ncurrent)
    message("\n (`RETURN' TO ACCEPT RECOMMENDATION,  0 TO EXIT AND RETURN CURRENT RESULTS)  ")
    
    ans  <-  readline()
    if ( nchar(ans)==0 ) return( -2 )
    ans  <-  as.integer(ans)
    if (is.na(ans)) next
    if ( -1<=ans && ans<=n ) return( ans )
  }
}

# ----------------------------------------------------------------------
#     ask the user to input a dose from those given
#     `ENTER' returns to the caller with no action, 
#
#     dose  - dose labels
#    ncurrent - patient no. for current patient
#     cohort - cohort size
# ----------------------------------------------------------------------
get.dose  <-  function( dose , ncurrent, cohort) {
  repeat {
    message("\n\n ENTER DOSE FOR PATIENTS ", ncurrent-cohort+1, " to ", ncurrent)
    message("\n POSSIBLE CHOICES ARE ", dose)
    message("\n (`RETURN' TO ACCEPT RECOMMENDATION,  0 TO EXIT AND RETURN CURRENT RESULTS)  ")
    
    ans  <-  readline()
    if ( nchar(ans)==0 ) return( -2 )
    if (is.na(ans)) next
    if ( ans %in% dose | ans==0) return( ans )
  }
}




#-----------------------------------------------------------------------
#    Plot function for an object of class bcrm
# -----------------------------------


#' Plot the estimated dose-toxicity curve
#' 
#' The estimated dose-toxicity curve using the Bayesian continuous reassessment
#' method is plotted for the patients thus far recruited into the trial
#' 
#' The estimated 2.5\%,  25\%,  50\%,  75\%,  97.5\% quantiles of the probability
#' of toxicity are plotted for each dose. Additionally,  a histogram of the
#' number of toxicities and non-toxicities is plotted at each experimented
#' dose.
#' 
#' If \code{trajectory = TRUE} then the sequential dose trajectory and observed
#' toxicities are plotted.
#' 
#' @param x An object of class "bcrm",  as returned by \code{\link{bcrm}}
#' @param each Should posterior summaries be plotted after each recruited
#' cohort? Defaults to FALSE.
#' @param trajectory Should the sequential dose trajectory of the recruited
#' patients be plotted,  along with the observed toxicities? Defaults to FALSE.
#' @param file File name where the dose-response plots are stored,  in a pdf
#' format. The program will amend the current sample size to the end of the
#' file name.
#' @param ... Further arguments passed to or from other methods
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{bcrm}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
plot.bcrm <- function(x, file=NULL, each=FALSE, trajectory=FALSE, ...){
  dose <- if(is.null(x$dose)) x$sdose else x$dose
  dose.label <- if(is.null(x$dose)) "Standardised dose" else "Dose"
  f <- which.f(x$ff)
  if(trajectory){
    df <- x$data
    df$Toxicity <- ifelse(df$tox==1, "Yes", "No")
    cols <-  c("No" = "black", "Yes" = "red")
    a <- ggplot()+geom_point(aes(x=patient, y=dose, shape=Toxicity, colour=Toxicity), size=3, data=df)+
      scale_colour_manual(values=cols)+
      xlab("Patient")+ylab("Dose Level")
    print(a)
    if(!is.null(file))
      ggsave(paste(file, ".pdf", sep=""), ...)
  } else {
    if(!each){
      if(x$method %in% c("exact", "exact.sim") & x$ff=="logit2"){
        df <- data.frame(dose=dose, target.tox=x$target.tox,
                         est=x$ndose[[length(x$ndose)]]$est)
      } else {
        df <- data.frame(dose=dose, target.tox=x$target.tox,
                         est=x$ndose[[length(x$ndose)]]$est,
                         mean=x$ndose[[length(x$ndose)]]$mean,
                         q2.5=x$ndose[[length(x$ndose)]]$quantiles["2.5%", ],
                         q25=x$ndose[[length(x$ndose)]]$quantiles["25%", ],
                         q50=x$ndose[[length(x$ndose)]]$quantiles["50%", ],
                         q75=x$ndose[[length(x$ndose)]]$quantiles["75%", ],
                         q97.5=x$ndose[[length(x$ndose)]]$quantiles["97.5%", ])
      }
      df2 <- data.frame(dose=factor(c(rep(dose, x$tox), rep(dose, x$notox)), levels=dose),
                        Outcome=factor(c(rep("DLT", sum(x$tox)), rep("No DLT", sum(x$notox))), levels=c("DLT", "No DLT")))
    } else { # each=TRUE
      if(x$method %in% c("exact", "exact.sim") & x$ff=="logit2"){
        df <- data.frame(dose=rep(dose, length(x$ndose)), target.tox=x$target.tox, cohort=rep(0:(length(x$ndose)-1), each=length(dose)), est=c(sapply(1:length(x$ndose), function(i){x$ndose[[i]]$est})), 
                         ndose=rep(sapply(1:length(x$ndose), function(i){dose[x$ndose[[i]]$ndose]}), each=length(dose)))
      } else {
        df <- data.frame(dose=rep(dose, length(x$ndose)), target.tox=x$target.tox, cohort=rep(0:(length(x$ndose)-1), each=length(dose)), est=c(sapply(1:length(x$ndose), function(i){x$ndose[[i]]$est})), mean=c(sapply(1:length(x$ndose), function(i){x$ndose[[i]]$mean})), q2.5=c(sapply(1:length(x$ndose), function(i){x$ndose[[i]]$quantiles["2.5%", ]})), q25=c(sapply(1:length(x$ndose), function(i){x$ndose[[i]]$quantiles["25%", ]})), 
                         q50=c(sapply(1:length(x$ndose), function(i){x$ndose[[i]]$quantiles["50%", ]})), q75=c(sapply(1:length(x$ndose), function(i){x$ndose[[i]]$quantiles["75%", ]})), q97.5=c(sapply(1:length(x$ndose), function(i){x$ndose[[i]]$quantiles["97.5%", ]})), 
                         ndose=rep(sapply(1:length(x$ndose), function(i){dose[x$ndose[[i]]$ndose]}), each=length(dose)))
      }
      df2 <- data.frame()
    }
    a <- if(is.null(x$loss)){
      if(!each){
        if(x$method %in% c("exact", "exact.sim") & x$ff=="logit2"){
          ggplot() +
            geom_point(aes(x=dose, y=est), data=df) +
            geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2) +
            xlab(dose.label) +
            ylab("Probability of DLT") +
            ylim(0, 1) +
            ggtitle("Posterior point estimates \n Diamond shows next recommended dose") +
            geom_point(aes(x=dose, y=est), data=df[x$ndose[[length(x$ndose)]][[1]], ],
                       size=4, col=4, shape=9)
        } else {
          ggplot(df) +
            geom_errorbar(aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red") +
            geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75),
                            fill="red", shape=3)+
            geom_hline(aes(yintercept=target.tox), col=4, linetype=2) +
            xlab(dose.label) +
            ylab("Probability of DLT") +
            ylim(0, 1) +
            ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \nPoints = p(DLT) estimates; Diamond = recommended dose") +
            geom_point(aes(x=dose, y=q50), size = 2) +
            geom_point(aes(x=dose, y=q50), data=df[x$ndose[[length(x$ndose)]][[1]], ],
                       size=4, col=4, shape=9)
        }
      } else {
        if(x$method %in% c("exact", "exact.sim") & x$ff=="logit2"){
          ggplot() +
            geom_point(aes(x=dose, y=est), data=df) +
            geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2)  +
            xlab(dose.label) +
            ylab("Probability of DLT") +
            ylim(0, 1) +
            ggtitle("Posterior point estimates \n Diamond shows next recommended dose") +
            geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9) +
            facet_wrap(~ cohort) 
        } else {
          ggplot() +
            geom_errorbar(aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red", data=df) +
            geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
            geom_hline(aes(yintercept=target.tox), data=df, col=4, linetype=2) +
            xlab(dose.label) +
            ylab("Probability of DLT") +
            ylim(0, 1) +
            ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose") +
            geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9)+
            facet_wrap(~ cohort) 
        }
      }
    } else { 
      df.intervals <- data.frame(cohort=rep(0:(length(x$ndose)-1), each=length(x$loss)), xmin=min(dose), xmax=max(dose), ymin=c(0, tox.cutpoints), ymax=c(x$tox.cutpoints, 1), Loss=x$loss)
      if(!each){
        ggplot() +
          geom_errorbar(aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red", data=df) +
          geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
          geom_point(aes(x=dose, y=q50), data=df[x$ndose[[length(x$ndose)]][[1]], ],
                     size=4, col=4, shape=9)+
          geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=Loss),
                    data=df.intervals, alpha=0.3) +
          scale_fill_gradient(breaks=sort(unique(df.intervals$Loss)),
                              high="red", low="#99ccff", guide="legend")+
          xlab(dose.label)+ylab("Probability of DLT") +
          ylim(0, 1)+
          ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose")
      } else {
        ggplot() +
          geom_errorbar(aes(x=dose, ymin=q2.5, ymax=q97.5), colour="red", data=df) +
          geom_pointrange(aes(x=dose, y=q50, ymin=q25, ymax=q75), data=df, fill="red") +
          geom_point(aes(x=ndose, y=q50), data=df[df$dose==df$ndose, ], size=4, col=4, shape=9) +
          geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=Loss), data=df.intervals, alpha=0.3) +
          xlab(dose.label) +
          ylab("Probability of DLT") +
          ylim(0, 1) +
          ggtitle("Posterior p(DLT) quantiles: 2.5%,  25%,  50%,  75%,  97.5% \n Diamond shows next recommended dose")+
          facet_wrap(~ cohort)
      }
    }
    b <- if(nrow(df2)!=0) {
      ggplot() +
        geom_bar(aes(x=dose, fill=Outcome), data=df2) +
        xlab(dose.label) +
        ylab("Number") +
        scale_fill_hue(limits=c("DLT", "No DLT"))
    } else { NULL  }
    
    if(!is.null(file))
      pdf(paste(file, sum(x$tox+x$notox), ".pdf", sep=""), ...)
    if(!each){
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(2, 1)))
      vplayout <- function(x, y)  viewport(layout.pos.row=x, layout.pos.col=y)
      print(a, vp=vplayout(1, 1))  
      if(!is.null(b)) print(b, vp=vplayout(2, 1))  
    } else {
      print(a)
    } 
    if(!is.null(file))
      dev.off()
  }
}

#-----------------------------------------------------------------------
#    Plot function for an object of class bcrm.sim
#    threep3     --> If TRUE (default is FALSE) then operating characteristics of the simulated design are compared against a standard rule-based 3+3 design
# -----------------------------------


#' Plot the operating characteristics from the simulated trials
#' 
#' Plots of the operating characteristics obtained from a CRM simulation.
#' 
#' This function plots the sample size distribution (if variable),  the
#' experimentation distribution,  the recommended dose distribution and the
#' percentage of subjects who experience the toxicity outcome (dose-limiting
#' toxicity). If \code{trajectories = TRUE} then summary statistics of
#' administered dose levels for each patient are plotted instead. If
#' \code{threep3 = TRUE} then the operating characteristics of the standard 3+3
#' design are plotted alongside those of the \code{bcrm} design (see
#' \code{\link{threep3}} for more details).
#' 
#' @param x An object of class "bcrm.sim",  as returned by \code{\link{bcrm}}
#' when conducting a simulation.
#' @param trajectories Should a summary plot of the trajectories of
#' administered dose levels be plotted? Defaults to FALSE.
#' @param file File name where the operating characteristic plot is stored,  in
#' a pdf format.
#' @param threep3 Should operating characteristics of a standard 3+3 rule-based
#' design be plotted alongside the \code{bcrm} design? Defaults to
#' \code{FALSE}.
#' @param ... Further arguments passed to or from other methods
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{print.bcrm.sim}},  \code{\link{bcrm}}, 
#' \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
plot.bcrm.sim <- function(x, trajectories=FALSE, file=NULL, threep3=FALSE, ...){
  dose <- if(is.null(x[[1]]$dose)) x[[1]]$sdose else x[[1]]$dose
  dose.label <- if(is.null(x[[1]]$dose)) "Standardised dose" else "Dose"
  
  if(trajectories){
    ## sample size
    n <- sapply(x, function(i){dim(i$data)[1]})
    traj.mat <- matrix(NA, nrow=length(x), ncol=max(n))
    tox.mat <- matrix(NA, nrow=length(x), ncol=max(n))
    for(i in 1:length(x)){
      traj.mat[i, 1:n[i]] <- x[[i]]$data$dose
      tox.mat[i, 1:n[i]] <- x[[i]]$data$tox
    }
    traj.df <- data.frame(patient=rep(1:max(n), each=5), Statistic=factor(rep(c("Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum"), max(n)), levels=c("Minimum", "Lower Quartile", "Median", "Upper Quartile", "Maximum")), traj=c(apply(traj.mat, 2, quantile, na.rm=T)))
    df <- data.frame(patient=rep(1:max(n), each=length(x)), sim=rep(1:length(n), max(n)), traj=c(traj.mat), Toxicity=ifelse(c(tox.mat)==1, "Yes", "No"))
    lt <- c("Median" = 1, "Lower Quartile" = 2, "Upper Quartile" = 2,  "Minimum" = 4, "Maximum"=4)
    cols <-  c("No" = "black", "Yes" = "red")
    if(length(x)>1){
      b <- ggplot()+geom_step(aes(x=patient, y=traj, group=Statistic, linetype=Statistic), size=1.2, colour="blue", data=traj.df)+
        scale_linetype_manual(values=lt)+
        xlab("Patient")+ylab("Dose Level")
      if(!is.null(file))
        pdf(paste(file, ".pdf", sep=""), ...)
      
      print(b)
      if(!is.null(file))
        dev.off()
    } else {
      a <- ggplot()+scale_colour_manual(values=cols)+
        xlab("Patient")+ylab("Dose Level")+
        geom_point(aes(x=patient, y=traj, col=Toxicity), data=df)
      print(a)
      if(!is.null(file))
        ggsave(paste(file, ".pdf", sep=""), ...)
    }
  } else {
    if(threep3 & is.null(x[[1]]$threep3)){
      message("Calculating 3+3 operating characteristics....")
      x[[1]]$threep3 <- threep3(x[[1]]$truep, x[[1]]$start)
    }
    # sample size
    n <- sapply(x, function(i){dim(i$data)[1]})
    df.n <- data.frame(n)
    if(!threep3){
      a <- if(min(n)!=max(n)){
        ggplot()+geom_histogram(aes(x=n, y=100*..density..), data=df.n, binwidth=1)+xlab("Sample size")+ylab("Percent")+ggtitle("Sample size")
      } else { NULL }
    } else {
      n.threep3 <- x[[1]]$threep3$ssize
      df.n.threep3 <- data.frame(n=c(n, n.threep3), weight=c(rep(1/length(n), length(n)), x[[1]]$threep$prob), Method=rep(c("CRM", "3+3"), c(length(n), length(n.threep3))))
      a <- ggplot()+stat_bin(aes(x=n, y=100*..density.., weight=weight, fill=Method), data=df.n.threep3, binwidth=1, position="dodge")+xlab("Sample size")+ylab("Percent")+ggtitle("Sample size")
    }
    
    # experimentation
    exp <- rep(dose, apply(sapply(x, function(i){(i$tox+i$notox)}), 1, sum))
    if(!threep3){
      df.exp <- data.frame(exp=factor(exp))
      b <- ggplot()+geom_bar(aes(x=exp, y=100*..count../sum(..count..)), data=df.exp)+xlab(dose.label)+ylab("Percent")+ggtitle("Experimentation")
    } else {
      exp.threep3 <- rep(dose, 10000*x[[1]]$threep3$exp)
      df.exp.threep3 <- data.frame(exp=factor(c(exp, exp.threep3)), Method=rep(c("CRM", "3+3"), c(length(exp), length(exp.threep3))), weight=c(rep(1/length(exp), length(exp)), rep(1/length(exp.threep3), length(exp.threep3))))
      b <- ggplot()+geom_bar(aes(x=exp, y=100*..count.., weight=weight, fill=Method), data=df.exp.threep3, position="dodge")+xlab(dose.label)+ylab("Percent")+ggtitle("Experimentation")
    }
    
    # recommendation
    rec <-   sapply(x, function(i){ifelse(i$ndose[[length(i$ndose)]]$ndose==0, 0, dose[i$ndose[[length(i$ndose)]]$ndose])})
    if(!threep3){
      df.rec <- data.frame(rec=factor(rec))
      c <- ggplot()+geom_bar(aes(x=rec, y=100*..count../sum(count)), data=df.rec)+xlab(dose.label)+ylab("Percent")+ggtitle("Recommendation")
    } else {
      rec.threep3 <- dose[x[[1]]$threep3$mtd]
      df.rec.threep3 <- data.frame(rec=factor(c(rec, rec.threep3)), weight=c(rep(1/length(rec), length(rec)), x[[1]]$threep$prob[x[[1]]$threep3$mtd!=0]), Method=rep(c("CRM", "3+3"), c(length(rec), length(rec.threep3))))
      c <- ggplot()+geom_bar(aes(x=rec, y=100*..count.., weight=weight, fill=Method), data=df.rec.threep3, position="dodge")+xlab(dose.label)+ylab("Percent")+ggtitle("Recommendation")
    }
    
    # observed DLTs
    obs <- sapply(x, function(i){100*sum(i$tox)/sum(i$tox+i$notox)})
    if(!threep3){
      bw <- max(diff(range(obs))/30, 1)
      # old version 0.4.6
      #df.obs <- data.frame(obs=obs)
      #d <- ggplot()+geom_histogram(aes(x=obs, y=100*..count../sum(..count..)), data=df.obs, binwidth=bw)+
      #  xlab("Percentage of subjects with DLTs")+ylab("Percent of trials")+ggtitle("DLTs")
      # new version >= 0.4.7 
      df.obs <- data.frame(obs=bw*round(obs/bw))
      d <- ggplot()+geom_bar(aes(x=obs, y=100*..count../sum(..count..)), data=df.obs)+
        xlab("Percentage of subjects with DLTs")+ylab("Percent of trials")+ggtitle("DLTs")
    } else {
      obs.threep3 <- 100*x[[1]]$threep3$dlt.no/x[[1]]$threep3$ssize
      range <- range(c(obs, obs.threep3))
      bw <- diff(range)/30
      
      df.obs.threep3 <- data.frame(obs=c(obs, obs.threep3), weight=c(rep(1/length(obs), length(obs)), x[[1]]$threep$prob), Method=rep(c("CRM", "3+3"), c(length(obs), length(obs.threep3))))
      df.obs.threep3 <- subset(df.obs.threep3, df.obs.threep3$weight>0)
      
      df.obs.threep3$obs.rounded <- bw*round(c(obs, obs.threep3)/bw)
      
      ## Old version 0.4.6
      #d <- ggplot()+geom_histogram(aes(x=obs, y=100*..count../sum(..count..), weight=weight, fill=Method), data=df.obs.threep3, binwidth=bw, position="dodge")+
      #  xlab("Percentage of subjects with DLTs")+ylab("Percent of trials")+ggtitle("DLTs")
      
      d <- ggplot()+geom_bar(aes(x=obs.rounded, y=100*..count.., weight=weight, fill=Method), data=df.obs.threep3, position="dodge")+
        xlab("Percentage of subjects with DLTs")+ylab("Percent of trials")+ggtitle("DLTs")
      
    }
    if(!is.null(file))
      pdf(paste(file, ".pdf", sep=""), ...)
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(2, 2)))
    vplayout <- function(x, y)  viewport(layout.pos.row=x, layout.pos.col=y)
    if(!is.null(a)) print(a, vp=vplayout(1, 1))  
    print(b, vp=vplayout(1, 2))
    print(c, vp=vplayout(2, 1))  
    print(d, vp=vplayout(2, 2))  
    if(!is.null(file))
      dev.off()
  }
}

#-----------------------------------------------------------------------
#    Plot function for an object of class threep3
# -----------------------------------


#' Plot the operating characteristics from a standard 3+3 trial
#' 
#' Plots of the operating characteristics obtained from a standard 3+3 trial, 
#' using \code{\link{threep3}}
#' 
#' This function plots the sample size distribution,  the experimentation
#' distribution,  the recommended dose distribution and the percentage of
#' subjects who experience the toxicity outcome (dose-limiting toxicity) for
#' the standard 3+3 trial.
#' 
#' @param x An object of class "threep3",  as returned by \code{\link{threep3}}.
#' @param file File name where the operating characteristic plot is stored,  in
#' a pdf format.
#' @param ... Further arguments passed to or from other methods
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
plot.threep3 <- function(x, file=NULL, ...){
  dose <- if(is.null(x$dose)) 1:length(x$truep) else x$dose
  dose.label <- if(is.null(x$dose)) "Dose levels" else "Dose"
  
  # sample size
  n.threep3 <- x$ssize
  df.n.threep3 <- data.frame(n=n.threep3, weight=x$prob)
  a <- ggplot()+stat_bin(aes(x=n, y=100*..density.., weight=weight), data=df.n.threep3, binwidth=1)+xlab("Sample size")+ylab("Percent")+ggtitle("Sample size")
  
  # experimentation
  exp.threep3 <- rep(dose, 10000*x$exp)
  df.exp.threep3 <- data.frame(exp=as.factor(exp.threep3))
  b <- ggplot()+geom_bar(aes(x=exp, y=..count../100), data=df.exp.threep3)+xlab(dose.label)+ylab("Percent")+ggtitle("Experimentation")
  
  # recommendation
  rec.threep3 <- dose[x$mtd]
  df.rec.threep3 <- data.frame(rec=factor(rec.threep3), weight=x$prob[x$mtd!=0])
  c <- ggplot()+geom_bar(aes(x=rec, y=100*..count.., weight=weight), data=df.rec.threep3)+xlab(dose.label)+ylab("Percent")+ggtitle("Recommendation")
  
  # observed DLTs
  obs.threep3 <- 100*x$dlt.no/x$ssize
  bw <- max(diff(range(obs.threep3[x$prob>0]))/30, 1)
  df.obs.threep3 <- data.frame(obs=obs.threep3, weight=x$prob, binwidth=bw)
  df.obs.threep3 <- subset(df.obs.threep3, df.obs.threep3$weight>0)
  df.obs.threep3$obs.rounded <- bw*round(obs.threep3/bw)
  d <- ggplot()+geom_bar(aes(x=obs.rounded, y=100*..count.., weight=weight), data=df.obs.threep3)+
    xlab("Percentage of subjects with DLTs")+ylab("Percent of trials")+ggtitle("DLTs")
  
  
  if(!is.null(file))
    pdf(paste(file, ".pdf", sep=""), ...)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2, 2)))
  vplayout <- function(x, y)  viewport(layout.pos.row=x, layout.pos.col=y)
  if(!is.null(a)) print(a, vp=vplayout(1, 1))  
  print(b, vp=vplayout(1, 2))
  print(c, vp=vplayout(2, 1))  
  print(d, vp=vplayout(2, 2))  
  if(!is.null(file))
    dev.off()
}


#-----------------------------------------------------------------------
#    Print function for an object of class bcrm
# -----------------------------------


#' Print information regarding a trial conducted using the Bayesian continuous
#' reassessment method
#' 
#' Print method for a trial or series of trials conducted using a
#' \code{\link{bcrm}} model.
#' 
#' If a single trial is conducted,  then the \code{\link{print}} function
#' currently produces summary information about the design used,  the data
#' observed,  current posterior estimates of toxicity,  and the next recommended
#' dose level.  If a simulation study is conducted,  then the following
#' operating characteristics are printed: \describe{ \item{Experimentation
#' proportion}{Proportion of patients recruited to each dose,  and to each true
#' region of toxicity,  across the simulated trials} \item{Recommendation
#' proportion}{Proportion of trials that recommend each of the dose levels as
#' the final maximum tolerated dose (i.e. with toxicity "closest" to the target
#' toxicity level),  and the associated regions of true toxicity for the
#' recommended MTDs} } If \code{trajectories = TRUE} then the dose level
#' administered and outcome observed are returned as matrices for every patient
#' (column) in every simulation (row).  If \code{threep3 = TRUE} then the
#' operating characteristics of the standard 3+3 design are displayed alongside
#' those of the \code{bcrm} design (see \code{\link{threep3}} for more
#' details).
#' 
#' @aliases print.bcrm print.bcrm.sim
#' @param x An object of class "bcrm" or "bcrm.sim" as returned by
#' \code{\link{bcrm}}
#' @param tox.cutpoints An optional argument passed to \code{print.bcrm.sim}
#' specifying the cutpoints of toxicity for which the operating characteristics
#' are to be categorised. Defaults to \code{seq(from=0, to=1, by=0.2)}
#' @param trajectories Should the individual simulation dose and outcome
#' trajectories be returned? Defaults to FALSE.
#' @param threep3 Should operating characteristics of a standard 3+3 rule-based
#' design be displayed alongside those from the \code{bcrm} design? Defaults to
#' \code{FALSE}.
#' @param ... Further arguments passed to or from other methods
#' @return The following two components are returned from
#' \code{print.bcrm.sim}: \item{exp}{A matrix with number of rows equal to the
#' number of doses,  and number of columns equal to the number of simulations.
#' Gives the experimentation proportions for each dose within each simulation.}
#' \item{rec}{A vector with length equal to the number of simulations,  giving
#' the recommended MTD for each simulation.}
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{bcrm}},  \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
print.bcrm <- function(x, ...){
  cat("\n Estimation method: ", x$method)
  
  cat("\n Target toxicity level: ", x$target.tox)
  
  ff.txt <- switch(x$ff
                   , ht="Hyperbolic Tangent"
                   , logit1="1-parameter logistic"
                   , power="1-parameter power"
                   , logit2="Two-parameter logistic")
  
  cat("\n Model: ", ff.txt)
  
  pa.txt <- switch(x$prior.alpha[[1]]
                   , "1"=paste("Gamma( Shape:", x$prior.alpha[[2]], ",  Scale:", x$prior.alpha[[3]], ")", sep="")
                   , "2"=paste("Uniform(", x$prior.alpha[[2]], ",  ", x$prior.alpha[[3]], ")", sep="")
                   , "3"=paste("Lognormal( Mean:", x$prior.alpha[[2]], ",  Variance:", x$prior.alpha[[3]], ")", sep="")
                   , "4"=paste("Log Multivariate Normal"))
  cat("\n Prior: ", pa.txt, "\n")
  if(x$prior.alpha[[1]]==4){
    message("Mean Vector: ")
    print(x$prior.alpha[[2]])
    cat("\nVariance-Covariance Matrix: \n")
    print(x$prior.alpha[[3]])
  }
  tab1 <- x$sdose
  names(tab1) <- x$dose
  cat("\n Standardised doses (skeleton): \n")
  print(tab1)
  
  if(x$constrain) { 
    if(is.null(x$dose)){
      cat("\n Modified (constrained) CRM used,  starting dose level: ", x$start, "\n") 
    } else {
      cat("\n Modified (constrained) CRM used,  starting dose: ", x$dose[x$start], "\n") 
    }
  } else { cat("\n Unmodified (unconstrained) CRM used \n") }
  
  if(!is.null(x$loss)){
    cat("\n Loss function given intervals of toxicity used to select next dose.")
    tab.lf <- x$loss
    names(tab.lf) <- levels(cut(0, breaks=c(0, x$tox.cutpoints, 1)))
    cat("\n Loss function: \n")
    print(tab.lf)
  } else if(x$pointest=="plugin"){
    cat("\n Plug-in estimate of probability of toxicity used to select next dose \n")
  } else if(x$pointest=="mean"){
    cat("\n Posterior mean estimate of probability of toxicity used to select next dose \n")
  } else { 
    cat("\n", 100*x$pointest, "percentile of (standardised) MTD distribution used to select next dose")
    cat("\n", 100*x$pointest, "percentile is:", x$ndose[[length(x$ndose)]]$target, "\n")
  }  
  
  tab <- rbind(x$tox+x$notox, x$tox)
  rownames(tab) <- c("n", "Toxicities")
  colnames(tab) <- x$dose
  names(dimnames(tab)) <- c("", "Doses")
  cat("\n Toxicities observed: \n")
  print(tab)
  
  if(x$method %in% c("exact", "exact.sim") & x$ff=="logit2"){
    tab2a <- signif(rbind(x$ndose[[length(x$ndose)]]$est), 3)
    rownames(tab2a) <- c("Point estimate")
  } else {
    tab2a <- signif(rbind(x$ndose[[length(x$ndose)]]$mean, x$ndose[[length(x$ndose)]]$sd, x$ndose[[length(x$ndose)]]$quantiles["50%", ]), 3)
    rownames(tab2a) <- c("Mean", "SD", "Median")
    tab2b <- signif(x$ndose[[length(x$ndose)]]$quantiles, 3)
    colnames(tab2b) <- x$dose
    names(dimnames(tab2b)) <- c("Quantiles", "Doses")
  }
  colnames(tab2a) <- x$dose
  names(dimnames(tab2a)) <- c("", "Doses")
  cat("\n Posterior estimates of toxicity: \n")
  print(tab2a)
  if(!(x$method %in% c("exact", "exact.sim") & x$ff=="logit2")){
    print(tab2b)
  }
  if(!is.null(x$loss)){
    tab3 <- rbind(x$ndose[[length(x$ndose)]]$est)
    colnames(tab3) <- x$dose
    cat("\n Posterior expected loss at each dose: \n")    
    print(tab3)
    tab4 <- x$ndose[[length(x$ndose)]]$probs
    colnames(tab4) <- x$dose
    rownames(tab4) <- levels(cut(0, breaks=c(0, x$tox.cutpoints, 1)))
    cat("\n Posterior probability of dose being in each toxicity interval")
    names(dimnames(tab4)) <- c("Toxicity intervals", "Doses")
    print(tab4)
  } else if(x$pointest=="plugin"){
    tab3 <- signif(rbind(x$ndose[[length(x$ndose)]]$est), 3)
    colnames(tab3) <- x$dose
    cat("\n Plug-in estimates of toxicity: \n")
    print(tab3)
  }
  if(is.null(x$dose)){
    cat("\n Next recommended dose level: ", x$ndose[[length(x$ndose)]]$ndose, "\n")
  } else {
    if(x$ndose[[length(x$ndose)]]$ndose!=0){
      cat("\n Next recommended dose: ", x$dose[x$ndose[[length(x$ndose)]]$ndose], "\n")
    } else {
      cat("\n No recommended dose levels are safe", "\n")
    }
  }
}

#-----------------------------------------------------------------------
#    Print function for an object of class bcrm.sim
# -----------------------------------
print.bcrm.sim <- function(x, tox.cutpoints=NULL, trajectories=FALSE,
                           threep3=FALSE, ...){
  if(trajectories){
    ## sample size
    n <- sapply(x, function(i){dim(i$data)[1]})
    dose.mat <- outcome.mat <- matrix(NA, nrow=length(x), ncol=max(n))
    for(i in 1:length(x)){
      dose.mat[i, 1:n[i]] <- x[[i]]$data$dose
      outcome.mat[i, 1:n[i]] <- x[[i]]$data$tox
    }
    return(list(doses=dose.mat, outcomes=outcome.mat))
  } else {
    # average sample size
    n.average <- mean(sapply(x, function(i){dim(i$data)[1]}))
    n.min <- min(sapply(x, function(i){dim(i$data)[1]}))
    n.max <- max(sapply(x, function(i){dim(i$data)[1]}))
    if(n.min==n.max){
      tab0 <- cbind(n.average)
      rownames(tab0) <- "Sample size"
      colnames(tab0) <- ""
    } else {
      tab0 <- cbind(n.average, n.min, n.max)
      rownames(tab0) <- "Sample size"
      colnames(tab0) <- c("Mean", "Minimum", "Maximum")
    }
    exp <- sapply(x, function(i){(i$tox+i$notox)/sum(i$tox+i$notox)})
    exp.tab <- c(NA, apply(exp, 1, mean))
    rec <- sapply(x, function(i){
      i$ndose[[length(i$ndose)]]$ndose
      })
    rec.tab <- prop.table(table(factor(rec, levels=0:length(x[[1]]$tox))))
    
    tab <- signif(rbind(exp.tab, rec.tab), 3)
    rownames(tab) <- c("Experimentation proportion", "Recommendation proportion")
    dose <- if(is.null(x[[1]]$dose)){1:length(x[[1]]$truep)} else {x[[1]]$dose}
    colnames(tab) <- c("No dose", dose)
    names(dimnames(tab)) <- c("", "Doses")
    if(is.null(tox.cutpoints)){
      tox.cutpoints <- seq(0, 1, by=0.2)
    } else {
      tox.cutpoints <- unique(c(0, tox.cutpoints, 1))
    }
    exp.tox <- prop.table(table(cut(unlist(sapply(x, function(i){rep(i$truep, (i$tox+i$notox))}, simplify=FALSE)), tox.cutpoints, include.lowest=T)))
    ## rec.tox updated on 04/07/17 so that trials that do not recommend a dose are classified as recommending a dose with 0% DLT rate
    rec.tox <- prop.table(table(cut(sapply(x, function(i){ifelse(i$ndose[[length(i$ndose)]]$ndose==0, 0, i$truep[i$ndose[[length(i$ndose)]]$ndose])}), tox.cutpoints, include.lowest=T)))
    tab2 <- signif(rbind(exp.tox, rec.tox), 3)
    rownames(tab2) <- c("Experimentation proportion", "Recommendation proportion")
    names(dimnames(tab2)) <- c("", "Probability of DLT")
    cat("Operating characteristics based on ", length(x), " simulations: \n \n")
    print(tab0)
    cat("\n")
    print(tab)
    cat("\n")
    print(tab2)
    if(threep3 & is.null(x[[1]]$threep3)){
      cat("\n Calculating 3+3 operating characteristics....\n")
      x[[1]]$threep3 <- threep3(x[[1]]$truep, x[[1]]$start)
    }
    if(threep3){
      cat("\n\n******************** 3+3 operating characteristics *****************\n")
      print.threep3(x[[1]]$threep3, tox.cutpoints=tox.cutpoints, dose=dose)
    }
  }
  invisible(list(rec=rec, exp=exp))
}

#-----------------------------------------------------------------------
#    Print function for an object of class threep3
# -----------------------------------


#' Print information regarding the operating characteristics of a standard 3+3
#' design
#' 
#' Print method for a 3+3 design specified using a \code{\link{threep3}}.
#' 
#' The following operating characteristics are printed for the standard 3+3
#' design: \describe{ \item{Sample size}{Mean,  minimum and maximum sample size
#' of the design} \item{Experimentation proportion}{Proportion of patients
#' recruited to each dose,  and to each true region of toxicity,  on average}
#' \item{Recommendation proportion}{Proportion of 3+3 trials that would
#' recommend each of the dose levels as the final maximum tolerated dose (see
#' \code{\link{threep3}} for definition of the MTD),  and the associated regions
#' of true toxicity for the recommended MTDs} \item{Average number of
#' patients}{The average number of patients dosed at each level} \item{Average
#' number of DLTs}{The average number of DLTs seen at each level} }
#' 
#' @param x An object of class "threep3" as returned by \code{\link{threep3}}
#' @param tox.cutpoints An optional argument passed to \code{print.threep3}
#' specifying the cutpoints of toxicity for which the operating characteristics
#' are to be categorised. Defaults to \code{seq(from=0, to=1, by=0.2)}
#' @param dose Optional vector of length \code{k} of actual doses for
#' presentation purposes
#' @param ... Further arguments passed to or from other methods
#' @author Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of
#' Cambridge,  UK)
#' @seealso \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
print.threep3 <- function(x, tox.cutpoints=NULL, dose=NULL, ...){
  if(is.null(dose)){
    dose <- 1:length(x$truep)
  }
  # average sample size
  n.average <- weighted.mean(x$ssize, x$prob)
  n.min <- min(x$ssize[x$prob>0])
  n.max <- max(x$ssize[x$prob>0])
  tab0 <- cbind(n.average, n.min, n.max)
  rownames(tab0) <- "Sample size"
  colnames(tab0) <- c("Mean", "Minimum", "Maximum")
  exp <- c(NA, x$exp)
  rec <- xtabs(x$prob~x$mtd)
  tab <- signif(rbind(exp, rec), 3)
  rownames(tab) <- c("Experimentation proportion", "Recommendation proportion")
  colnames(tab) <- c(paste("<", dose[1]), dose)
  names(dimnames(tab)) <- c("", "Doses")
  if(is.null(tox.cutpoints)){
    tox.cutpoints <- seq(0, 1, by=0.2)
  } else {
    tox.cutpoints <- unique(c(0, tox.cutpoints, 1))
  }
  exp.tox <- xtabs(x$exp~cut(x$truep, tox.cutpoints, include.lowest=T))
  rec.tox <- xtabs(rec[-1]~cut(x$truep, tox.cutpoints, include.lowest=T))
  tab2 <- signif(rbind(exp.tox, rec.tox), 3)
  rownames(tab2) <- c("Experimentation proportion", "Recommendation proportion*")
  names(dimnames(tab2)) <- c("", "Probability of DLT")
  
  tab3 <- signif(rbind(x$n.average, x$dlt.average), 3)
  rownames(tab3) <- c("Average number of patients", "Average number of DLTs")
  colnames(tab3) <- dose
  names(dimnames(tab3)) <- c("", "Doses")
  
  print(tab0)
  cat("\n")
  print(tab)
  cat("\n")
  print(tab2)
  cat("\n * Amongst those trials that recommend an MTD\n")
  cat("\n")
  print(tab3)
}

