#############################################################################
# Overall function that calls the rest
#############################################################################
varSelectIP <- function(response, covariates.retain=NULL, covariates.test, 
                     nsim, keep, q, a=0.2, model.type=c("probit", "reg"),
                     save.every=50, out.fname="models.csv",
                     parallel=FALSE, interactive=TRUE, nproc=2) {
 
  p <- ncol(covariates.test)
  n.obs <- nrow(covariates.test)

  if(is.null(covariates.retain))
    num.covariates.retain <- 0
  else {
    num.covariates.retain <- ifelse(class(covariates.retain) == "matrix",
      ncol(covariates.retain), 1) # Because covariates.retain could be a vec
  }

  # If q is missing, set it to be (total number of covariates to be tested - 1)
  # so that it is really a low-dimensional search.
  if (missing(q)) {
    q <- ncol(covariates.test) - 1
  }
  else {
    if (num.covariates.retain + q + 1 > n.obs) {
      stop("Number of observations < num.covariates.retain + q + 1 !\n")
    }
  }

  # Order the observations in the reponse column if we are testing a probit 
  # model.
  if (model.type == "probit") {
    tmp.df <- cbind(response, covariates.retain, covariates.test)
    tmp.df <- tmp.df[order(tmp.df[,1]), ] # decreasing = FALSE by default
    response <- tmp.df[,1]
    if(is.null(covariates.retain)) {
      covariates.test <- tmp.df[,-1]
    }
    else {
      covariates.retain <- tmp.df[,2:(1+num.covariates.retain)]
      covariates.test <- tmp.df[ ,(2+num.covariates.retain):ncol(tmp.df)]
    }
    m0.z <- m0(response)
  }
  else 
    m0.z <- NULL

  # Next line gives a matrix of "NA" values. The reason we have 2*nsim, is 
  # because at each step, we test TWO candidates, so both are stored.
  models.visited <- matrix(nrow=2*nsim, ncol=2, dimnames=list(NULL, 
    c("modelId", "BF")))

  # If parallel=TRUE, check if Rmpi is installed.
  if (parallel) {
    if (!("Rmpi" %in% installed.packages()[,1])) 
      stop("Rmpi is not installed - no parallel processing possible!\n")
    if (interactive) {
      # Check if Rmpi can be loaded properly
      rmpiLoaded <- require(Rmpi)
      if (rmpiLoaded) {
        mpi.spawn.Rslaves(nslaves=nproc)
      }
      else
        stop("Rmpi could not be loaded properly! Please load it before 
              calling varSelectIP().\n")
    }
    mpi.bcast.cmd(library(mvtnorm))
  }

  # Obtain starting point
  repeat{
    current.model <- sample(0:1, size = p, replace=TRUE)
    if (sum(current.model) <= q)
      break
  }
  if(model.type == "probit") {
    current.model.score <- BayesFactorProbit(response, covariates.retain,
      covariates.test, current.model, parallel)
  }
  else {
    current.model.score <- BayesFactorLinReg(response, covariates.retain,
      covariates.test, current.model)
  }
  current.gamma <- current.model
  current.active <- current.model
  if (sum(current.active) < q) {
    need.to.set <- q - sum(current.active)
    current.active[current.active!=1][1:need.to.set] <- 1
  }
  cat("\nStarting point:\n")
  cat("Initial Active Set: ", current.active, "\n", sep="")
  cat("Initial Gamma (latent): ", current.gamma, "\n", sep="")
  cat("Initial Model : ", current.model, "\n\n", sep="")
  
  # Start running the chain
  for (i in 1:nsim) {
    NewState <- NextModel(response, covariates.retain, covariates.test, 
      current.gamma, current.active, current.model, current.model.score, 
      a, model.type, parallel, models.visited)
    current.gamma <- NewState[[1]]
    current.active <- NewState[[2]]
    current.model <- NewState[[3]]
    current.model.score <- NewState[[4]]
    models.visited <- NewState[[5]]
  
    if (i %% save.every == 0)
      current.table <- writeToFile(models.visited, fname=out.fname, 
        keep=keep, p=p, m0.z)
    cat("Finished", i, "simulations", "\n", sep=" ")
  }

  final.table <- writeToFile(models.visited, fname=out.fname, keep=keep, 
                             p=p, m0.z)

  if(parallel)
    mpi.close.Rslaves(dellog=FALSE)
   
#  models.visited
  return(final.table)
}
