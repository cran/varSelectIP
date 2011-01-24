# Y: A vector containing the (continuous) response.
# covariates.retain: A matrix of covariates that should always be in the 
#                     design matrix.
# covariates.test: A matrix of covariates to search through. It should be of 
#                   dimension (n \times p).
# current.model: A 0-1 vector of length p, whose Bayes Factor is to be found, 
#               versus the intercept-only model.

BayesFactorLinReg <- function(Y, covariates.retain=NULL, 
                              covariates.test, current.model){
  fn01 <- function(phi) {
    numerator <- (sin(phi))^(j-1)*(nn + (j+1)*(sin(phi))^2)^(0.5*(nn-j)) 
    denominator <- (nn*B + (j+1)*((sin(phi))^2))^(0.5*(nn-1))
    numerator/denominator
  }

  # If somehow we hit upon the intercept-only model, return 1.
  if(is.null(covariates.retain) & (sum(current.model) == 0)) {
    return(1)
  }
  else{  # compute the Bayes Factor for the current model agst intercept model.
    Xtmp <- cbind(rep(1, length(Y)), covariates.retain, 
                  covariates.test[ , which(current.model == 1)])
    XtX_inv <- ginv(t(Xtmp) %*% Xtmp)

    nn <- length(Y)		
    H <- Xtmp %*% XtX_inv%*% t(Xtmp)
    H0 <- matrix(1/nn, nrow=nn, ncol=nn)

    # total number of covariates (including intercept)
    j <- ncol(Xtmp)  

    # compute \mathscr{B}^n_{ij} from the paper (pp. 1212)
    # TODO: pre-compute the denominator here.
    B <- (Y %*% (diag(nn) - H) %*% Y)/(Y %*% (diag(nn) - H0) %*% Y)

    # compute expression (7) of paper
    output <- (2/pi)*(j+1)^((j-1)/2) * 
                integrate(fn01, lower=0, upper=pi/2)$value
  }
  return(output)
}
