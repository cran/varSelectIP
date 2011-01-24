m0 <- function(Z){

  n0 <- sum(Z == 0)
  nn <- length(Z)
  n1 <- nn - n0

  #############################################################################
  Sigma.j <- diag(nn)
  lower <- c(rep(-Inf, n0), rep(0, n1))
  upper <- c(rep(0, n0), rep(Inf, n1))
    
  MLE.alpha <- qnorm(mean(Z))
  
  interval.length <- 0.3
  alphas <- seq(MLE.alpha - 6, MLE.alpha + 6, by=interval.length)
  
  output <- 0
  for (ii in 1:length(alphas)) {
      output <- output + pmvt(lower=lower, upper=upper, df=0, 
                corr=Sigma.j, delta=rep(alphas[ii], nn), abseps=1e-35, 
                algorithm=GenzBretz())[1]
  }
    
  return(output * interval.length)
  #############################################################################
}
