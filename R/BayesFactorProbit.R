BayesFactorProbit <- function(Z, covariates.retain=NULL, covariates.test, 
                              current.model){

  #############################################################################
  marginalZ.serial <- function() {
      Sigma.j <- diag(nn) + 2 * nn/k * H
      lower <- c(rep(-Inf, n0), rep(0, n1))
      upper <- c(rep(0, n0), rep(Inf, n1))
      
      del <- rep(0, times=nn)
      corrMat <- cov2cor(Sigma.j)
  
      interval.length <- 0.4
      MLE.alpha <- qnorm(mean(Z))
      alphas <- seq(MLE.alpha - 6, MLE.alpha + 6, by=interval.length)
  
      output <- 0
      for (ii in 1:length(alphas)) {
        lower.tmp <- (lower - alphas[ii])/sqrt(diag(Sigma.j))
        upper.tmp <- (upper - alphas[ii])/sqrt(diag(Sigma.j))
        output <- output + pmvt(lower=lower.tmp, upper=upper.tmp, df=0, 
                  corr=corrMat, delta=del, abseps=1e-35, 
                  algorithm=GenzBretz())[1]
      }
    
      return(output * interval.length)
  }
  #############################################################################
  #############################################################################

  if(is.null(covariates.retain) & (sum(current.model) == 0)) {
    return
  }

  n0 <- sum(Z == 0)
  nn <- length(Z)
  n1 <- nn - n0

  Xtmp <- cbind(rep(1, length(Z)), covariates.retain, 
                  covariates.test[ , which(current.model == 1)])
  XtX_inv <- solve(t(Xtmp) %*% Xtmp)

  H <- Xtmp %*% XtX_inv %*% t(Xtmp) # HAT matrix

  # total number of covariates (including intercept)
  k <- ncol(Xtmp)  
 
    output <- marginalZ.serial()
  
  return(output)
}
