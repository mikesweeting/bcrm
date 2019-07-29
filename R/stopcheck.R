# ----------------------------------------------------------------------
#     Check whether we have reached any stopping criteria
#
#     stop         --> list of stopping criteria (nmax,  nmtd,  precision,  nmin)
#     target.tox   --> Target Toxicity Level
#     ncurrent     --> current sample size
#     ndose     --> next dose information,  including posterior quantities
#     new.tox     --> No. toxicities to date
#     new.notox     --> No. non-toxicities to date
#     simulate     --> Is this a simulation study?
# ----------------------------------------------------------------------
stop.check <- function(stop, target.tox, ncurrent, ndose, new.tox, new.notox, simulate){
  answer1 <- answer2 <- answer3 <- answer5 <- FALSE
  answer4 <- TRUE
  if(!is.null(stop$nmin)){
    answer4 <- (ncurrent>=stop$nmin)
  } 
  if(!is.null(stop$nmax)){
    answer1 <- (ncurrent>=stop$nmax)
    if(answer1 & answer4 & !simulate)
      message("\n Stopping: Reached maximum sample size")
  }  
  if(!is.null(stop$nmtd)){
    answer2 <- (new.tox+new.notox)[ndose$ndose]>=stop$nmtd
    if(answer2 & answer4 & !simulate)
      message("\n Stopping: Reached maximum number at MTD estimate")
  }
  if(!is.null(stop$precision)){
    answer3 <-  stop$precision[1] <= ndose$quantiles["2.5%", ndose$ndose] & ndose$quantiles["97.5%", ndose$ndose]<= stop$precision[2]
    if(answer3 & answer4 & !simulate)
      message("\n Stopping: Reached required precision for MTD estimate")
  }
  if(!is.null(stop$safety)){
    ## If prob DLT of first dose being greater than TTL is > stop$safety then stop 
    answer5 <- ndose$quantiles[paste0(100*(1-stop$safety), "%"), 1]>=target.tox
    if(answer5 & answer4) {
      if(!simulate){
        message("\n Stopping: Safety criteria reached. No doses deemed safe")
      }
      ndose$ndose <- 0
    }
  }
  
  stop <- (answer1 | answer2 | answer3 | answer5) & answer4 
  return(list(stop=stop,  ndose=ndose))
}