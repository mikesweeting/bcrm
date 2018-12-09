
intFun <- function(X){
  attach(sys.frame(sys.parent(2)), 2)
  
  ncurrent <- sum(new.tox + new.notox)
  stopped <- stop.check(stop, ncurrent, ndose, new.tox, new.notox, simulate)

  new.tox <- as.numeric(xtabs(tox~factor(dose, levels=1:k), data=data))
  new.notox <- as.numeric(xtabs((1-tox)~factor(dose, levels=1:k), data=data))
  current <- as.numeric(data$dose[dim(data)[1]])
  alpha <- switch(method
                  , rjags=Posterior.rjags(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr)
                  , BRugs=Posterior.BRugs(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr)
                  , R2WinBUGS=Posterior.R2WinBUGS(new.tox, new.notox, sdose, ff, prior.alpha, burnin.itr, production.itr, bugs.directory)
                  , exact=Posterior.exact(new.tox, new.notox, sdose, ff, prior.alpha)
                  , exact.sim=Posterior.exact.sim(new.tox, new.notox, sdose, ff, prior.alpha, pointest)
  )
  
  ncurrent <- sum(new.tox + new.notox)
  newdata <- data
  
  results <- list(dose=dose, sdose=sdose, tox=new.tox, notox=new.notox, ndose=list(ndose), constrain=constrain, start=start, target.tox=target.tox, ff=ff, method=method, pointest=pointest, tox.cutpoints=tox.cutpoints, loss=loss, prior.alpha=prior.alpha, data=data)
  class(results) <- "bcrm"
  
  if(plot){
    plot(results, file)
  }
  
  interact <- crm.interactive(new.tox, new.notox, ncurrent, cohort, ndose, dose)
  if(interact$bk == TRUE){
    results$data <- newdata
    return(results)
  }
  y <- interact$y
  new.tox <- interact$tox
  new.notox <- interact$notox
  current <- interact$ans
  
  results$tox <- new.tox
  results$notox <- new.notox
  results$ndose[[length(results$ndose) + 1]] <- ndose
  results$data <- newdata
  
  if(plot){
    plot(results, file)
  }
  
}
