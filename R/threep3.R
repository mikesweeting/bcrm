#####################################################################
#
# threep3 - Generates all possible 3+3 trial pathways,  probabilities
#            of occurrence,  MTD recommendation and sample size
# AUTHOR: Graham Wheeler
#
# ARGUMENTS
# truep - vector of probabilities for DLT at each dose level
#              to be investigated
# start - starting dose level for the design (defaults to 1)
#    dose    --> optional vector of dose labels (for printing and plotting purposes)
#
# VALUES
#       "prob" - probability of each trial occurring
#       "ssize" - sample size per trial
#       "mtd" - final MTD recommendation per trial
#     "exp" - experimentation proportions across all possible 3+3 trials
#       "truep" - true vector of probabilities
#
####################################################################
#####################################################################
#
# threep3 - Generates all possible 3+3 trial pathways,  probabilities
#            of occurrence,  MTD recommendation and sample size
# AUTHORS: Graham Wheeler,  Michael Sweeting
#
# ARGUMENTS
# truep - vector of probabilities for DLT at each dose level
#              to be investigated
# start - starting dose level for the design (defaults to 1)
#    dose    --> optional vector of dose labels (for printing and plotting purposes)
#
# VALUES
#       "prob" - probability of each trial occurring
#       "ssize" - sample size per trial
#       "mtd" - final MTD recommendation per trial
#     "exp" - experimentation proportions across all possible 3+3 trials
#       "truep" - true vector of probabilities
#
####################################################################



#' Calculate all possible trial pathways for the standard 3+3 design,  together
#' with their probability of occurring
#' 
#' All possible pathways of a standard 3+3 design are calculated and assigned a
#' probability of occurring. This facilitates the calculation of operating
#' characteristics,  using \code{\link{print.threep3}} and
#' \code{\link{plot.threep3}}.
#' 
#' The dose-escalation schema used here relates to that defined by Chang et al.
#' (2006),  which incorporates dose de-escalation also. Variations of this
#' design exist in the literature (see Storer 1989,  Reiner et al. 1999).
#' 
#' The first cohort of three patients are administered the starting dose
#' (usually the lowest dose). The trial then proceeds as follows: \itemize{
#' \item If none of the three patients experience a DLT,  then dose the next
#' three patients at the next highest dose level; \item If one of the three
#' patients last treated experiences a DLT,  then dose the next three patients
#' at the current dose level; \item If at least two patients in the first dose
#' level experience a DLT the trial is stopped for safety and no dose is
#' recommended; } Escalation / de-escalation rules to the next dose level for
#' subsequent cohorts proceed as follows: \itemize{ \item Escalate: If 0/3 or
#' at most 1/6 DLTs are observed in the current cohort AND the next highest
#' dose has not yet been tested; \item Stay at current dose level: If 1/3 DLTs
#' have been observed at this level. Dose a further three patients at the same
#' level; \item De-Escalate: If at least two out of three to six patients
#' experience DLTs at the current dose level AND fewer than six patients have
#' been dosed at the next lowest level }
#' 
#' If none of the rules above are satisfied then the trial stops. If the
#' current dose level has at most one DLT observed then this is claimed to be
#' the MTD,  otherwise the dose level below is deemed to be the MTD.
#' 
#' If dose-escalation extends to doses outside of that defined by \code{dose}, 
#' the MTD is determined to be the largest dose in \code{dose}.
#' 
#' @param truep A vector of length \code{k} (the number of doses being
#' considered in the trial),  with values equal to the true probabilities of
#' toxicity at the dose levels.
#' @param start Starting dose level. Defaults to 1,  i.e. the lowest dose level
#' @param dose Optional vector of length \code{k} of actual doses for
#' presentation purposes
#' @return \code{threep3} returns an object of class "threep3". The function
#' \code{\link{print}} (i.e. \code{\link{print.threep3}}) can be used to obtain
#' operating characteristics of the design used.
#' 
#' An object of class "threep3" is a list with the following components:
#' \item{prob}{A vector with the probabilities of each design occurring. As all
#' possible designs are calculated,  this vector sums to one} \item{ssize}{A
#' vector with the sample size of each design} \item{mtd}{A vector of dose
#' levels giving the recommended maximum tolerated dose (MTD) at the end of the
#' trial} \item{exp}{A vector of length \code{k} giving the average trial
#' experimentation proportions at each dose level} \item{dlt.no}{A vector with
#' the number of toxicities (DLTs) that occur in each trial} \item{truep}{The
#' true probabilities of toxicity at each dose level,  specified by the user}
#' \item{dose}{The actual doses as supplied in the function arguments}
#' \item{n.average}{The average number of patients dosed at each level}
#' \item{dlt.average}{The average number of DLTs experienced at each dose
#' level} \item{all.designs}{A matrix containing all possible 3+3 designs,  with
#' each row representing a different design. Columns labelled "d k" and "tox k"
#' represent the dose level and number of toxicities for the kth cohort, 
#' respectively.}
#' @author Graham Wheeler \email{graham.wheeler@@mrc-bsu.cam.ac.uk} (MRC
#' Biostatistics Unit,  Cambridge,  UK) and
#' 
#' Michael Sweeting \email{mjs212@@medschl.cam.ac.uk} (University of Cambridge, 
#' UK)
#' @seealso \code{\link{threep3}}
#' @references Sweeting M.,  Mander A.,  Sabin T. \pkg{bcrm}: Bayesian Continual
#' Reassessment Method Designs for Phase I Dose-Finding Trials. \emph{Journal
#' of Statistical Software} (2013) 54: 1--26.
#' \url{http://www.jstatsoft.org/article/view/v054i13}
#' 
#' Chang A.,  Ganz P.,  Hayes D.,  Kinsella T.,  Pass H.,  Schiller J.,  Stone R., 
#' Strecher V. \emph{Oncology: An Evidence-Based Approach}. Springer (2006).
#' 
#' Storer B. Design and Analysis of Phase I Clinical Trials. \emph{Biometrics}
#' (1989) 45: 925--937.
#' 
#' Reiner E.,  Paoletti X.,  O'Quigley J. Operating characteristics of the
#' standard phase I clinical trial design. \emph{Computational Statistics &
#' Data Analysis} (1999) 30: 303--315.
#' 
#' Neuenschwander B.,  Branson M.,  Gsponer T. Critical aspects of the Bayesian
#' approach to phase I cancer trials. \emph{Statistics in Medicine} (2008) 27:
#' 2420--2439.
#' @examples
#' 
#' ## What are the operating characteristics of a standard 3+3 design if we conside only the first 
#' ## 12 doses of the dose-escalation cancer trial example as described in Neuenschwander et al 2008.
#' ## Pre-defined doses
#' dose <- c(1, 2.5, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100)
#' ## Pre-specified probabilities of toxicity
#' p.tox0 <- c(0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.100, 0.170, 0.300, 0.400, 0.500)
#' 
#' \dontrun{
#' design.threep3 <- threep3(truep=p.tox0, start=1, dose=dose)
#' print(design.threep3)
#' plot(design.threep3)
#' }
#' 
#' @export threep3
threep3 <- function(truep, start=1, dose=NULL){
  # Check that dose is the same length as truep
  if(!is.null(dose) & length(dose)!=length(truep)) stop("Length of 'dose' must be the same as the length of 'truep'.")
  
  # Define number of doses and max. cohort number
  # 'mcplus1' used in computation
  doses <- length(truep)
  mcohort <- 2*doses
  mcplus1 <- mcohort+1
  
  # Check that starting dose is not top dose
  if(start==doses) stop("Starting dose cannot be uppermost dose")
  
  # Begin deriving pathways
  pmat <- as.data.frame(matrix(NA, nrow=1, ncol=2*mcplus1+1))
  colnames(pmat) <- c("stop", "desc", paste(c("d", "tox"), rep(1:mcohort, each=2)), paste("d", mcplus1))
  pmat[1, 1:3] <- c(0, 0, start)
  pmat <- pmat[rep(seq_len(nrow(pmat)),  rep(4, nrow(pmat))), ]
  pmat[, "tox 1"] <- c(0, 1, 2, 3)
  
  pmat[pmat[, "tox 1"]==0, "d 2"] <- start+1
  pmat[pmat[, "tox 1"]==1, "d 2"] <- start
  pmat[pmat[, "tox 1"]>1 & start>1, "d 2"] <- start-1
  pmat[pmat[, "tox 1"]>1 & start==1, "stop"] <- 1
  pmat[pmat[, "tox 1"]>1, "desc"] <- 1
  stopped.pmat <- pmat[pmat$stop==1, -2]
  all.designs <- stopped.pmat
  prob <- ssize <- mtd <- dlt.no <- NULL
  exp <- 0
  
  n.average <- dlt.average <- 0
  
  if(start==1){
    # Probabilties of stopped 3+3 trials occuring
    dose.mat <- stopped.pmat[, grep("d", names(stopped.pmat))]
    tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))]
    prob <- apply(matrix(dbinom(as.matrix(tox.mat), 3, truep[as.matrix(dose.mat)]), nrow=nrow(dose.mat)), 1, prod, na.rm=T)
    
    # Calculate sample size 
    ssize <- 3*apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]), 1, sum)
    
    # Determine MTD per stopped trial
    last.cohort <- apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]), 1, sum)
    last.drug.column <- paste("d", last.cohort)
    last.drug <- sapply(1:nrow(stopped.pmat), function(j){stopped.pmat[j, last.drug.column[j]]})
    previous.drug.column <- paste("d", last.cohort-1)
    previous.drug <- sapply(1:nrow(stopped.pmat), function(j){ifelse(previous.drug.column[j]=="d 0", 0, stopped.pmat[j, previous.drug.column[j]])})
    last.tox.column <- paste("tox", last.cohort)
    last.tox <- sapply(1:nrow(stopped.pmat), function(j){stopped.pmat[j, last.tox.column[j]]})
    mtd <- rep(NA, nrow(stopped.pmat))  
    mtd[last.tox==0] <- last.drug[last.tox==0]  
    mtd[last.tox==1 & previous.drug==last.drug] <- last.drug[last.tox==1 & previous.drug==last.drug]-1
    mtd[last.tox==1 & previous.drug!=last.drug] <- last.drug[last.tox==1 & previous.drug!=last.drug]
    mtd[last.tox>1] <- last.drug[last.tox>1]-1
    
    # Prob. that each dose is experimented on and trial occurs
    exp <- sapply(1:doses, function(j){sum(3*(stopped.pmat[, grep("d", names(stopped.pmat))]==j)*prob/ssize, na.rm=T)})
    # Average number of people dosed at each level
    n.average <- sapply(1:doses, function(j){sum(3*(stopped.pmat[, grep("d", names(stopped.pmat))]==j)*prob, na.rm=T)})
    
    # Number of subjects who have DLT per trial
    dlt.no <- apply(stopped.pmat[, grep("tox", names(stopped.pmat))], 1, sum, na.rm=T)
    # Average number of DLTs at each dose
    dlt.average <- sapply(1:doses, function(j){sum((stopped.pmat[, grep("tox", names(stopped.pmat))]*prob)[stopped.pmat[, grep("d", names(stopped.pmat))]==j], na.rm=T)})
    
    
  }
  
  for(i in 3:mcplus1){
    message(paste(round(100*i/mcplus1), "% complete", sep=""))
    dd <- as.character(paste("d", i))
    td <- as.character(paste("tox", i))
    dc <- as.character(paste("d", i-1))
    tc <- as.character(paste("tox", i-1))
    db <- as.character(paste("d", i-2))
    tb <- as.character(paste("tox", i-2))
    
    ## Creates new data.frame with 1,  2,  or 3 toxicities for every continued trial
    pmat <- pmat[rep(which(pmat[, "stop"]==0),  each=4), ]
    pmat[, tc] <- 0:3
    
    pmat[pmat[, tc]==0 & pmat[, "desc"]==0 & pmat[, dc]+1 <=doses, dd] <-  pmat[pmat[, tc]==0 & pmat[, "desc"]==0 & pmat[, dc]+1 <=doses, dc]+1
    pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==0, dd] <-  pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==0, dc]
    pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==1 & pmat[, dc]-1 >=1 , dd] <-  pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==1, dc]-1
    pmat[pmat[, tc]==0 & pmat[, "desc"]==1 , dd] <-  pmat[pmat[, tc]==0 & pmat[, "desc"]==1, dc]
    pmat[pmat[, tc]==1 & pmat[, "desc"]==1,  dd] <-  pmat[pmat[, tc]==1 & pmat[, "desc"]==1, dc]
    pmat[pmat[, tc]>1 & pmat[, dc]-1 >=1, "desc"] <-  1
    
    pmat[pmat[, tc]==1 & pmat[, "desc"]==0 & pmat[, tb]==1 & pmat[, dc]-1 >=1, "desc"] <-  1
    pmat[pmat[, tc]>1 & pmat[, dc]-1 >=1, dd] <-  pmat[pmat[, tc]>1 & pmat[, dc]-1>=1 , dc]-1
    
    
    excluding.dd <- names(pmat)[grepl("d ", names(pmat)) & names(pmat)!=dd]
    cnt <- apply(pmat[!is.na(pmat[, dd]), dd]==pmat[!is.na(pmat[, dd]), excluding.dd], 1, sum, na.rm=T)
    
    pmat[!is.na(pmat[, dd]), dd][cnt>1] <- NA
    pmat[is.na(pmat[, dd]), "stop"] <- 1
    stopped.pmat <- pmat[pmat$stop==1, -2]
    all.designs <- rbind(all.designs, stopped.pmat)
    if(dim(stopped.pmat)[1]>0){
      # Probabilties of stopped 3+3 trials occuring
      dose.mat <- stopped.pmat[, grep("d", names(stopped.pmat))[1:(i-1)]]
      tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))[1:(i-1)]]
      
      
      # Add these probabilities to prob
      prob.new <- apply(matrix(dbinom(as.matrix(tox.mat), 3, truep[as.matrix(dose.mat)]), nrow=nrow(dose.mat)), 1, prod, na.rm=T)
      prob <- c(prob, prob.new)
      
      # Calculate sample size and determine MTD per stopped trial
      # Add them to existing ssize and mtd vectors
      ssize.new <- rep(3*(i-1), nrow(stopped.pmat))
      ssize <- c(ssize, ssize.new)
      last.drug <- stopped.pmat[, dc]
      previous.drug <- stopped.pmat[, db]
      last.tox <- stopped.pmat[, tc]
      
      mtd.new <- rep(NA, nrow(stopped.pmat))
      mtd.new[last.tox==0] <- last.drug[last.tox==0]  
      mtd.new[last.tox==1 & previous.drug==last.drug] <- last.drug[last.tox==1 & previous.drug==last.drug]-1
      mtd.new[last.tox==1 & previous.drug!=last.drug] <- last.drug[last.tox==1 & previous.drug!=last.drug]
      mtd.new[last.tox>1] <- last.drug[last.tox>1]-1
      mtd <- c(mtd, mtd.new)
      
      # Prob. that each dose is experimented on and trial occurs (summing over all trials)
      exp <- exp+sapply(1:doses, function(j){sum(3*(stopped.pmat[, grep("d", names(stopped.pmat))]==j)*prob.new/ssize.new, na.rm=T)})
      # Average number of people dosed at each level
      n.average <- n.average+sapply(1:doses, function(j){sum(3*(stopped.pmat[, grep("d", names(stopped.pmat))]==j)*prob.new, na.rm=T)})
      
      
      # Number of subjects who have DLT per trial
      dlt.no <- c(dlt.no, apply(stopped.pmat[, grep("tox", names(stopped.pmat))], 1, sum, na.rm=T))
      # Average number of DLTs at each dose
      dlt.average <- dlt.average+sapply(1:doses, function(j){sum((stopped.pmat[, grep("tox", names(stopped.pmat))]*prob.new)[stopped.pmat[, grep("d", names(stopped.pmat))]==j], na.rm=T)})
      
      
    }
  }
  obj <- list(prob=prob, ssize=ssize, mtd=mtd, exp=exp, dlt.no=dlt.no, truep=truep, dose=dose, n.average=n.average, dlt.average=dlt.average, all.designs=all.designs)
  class(obj) <- "threep3"
  return(obj)
}
