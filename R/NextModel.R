########################################################
# current.gamma: This is the 0-1 vector of length p, which specifies the latent 
#                model under consideration. It could have 
#                sum(current.gamma) > q.

NextModel <- function(response, covariates.retain=NULL, covariates.test, 
                      current.gamma, current.active, current.model, 
                      current.model.score, a=0.2, model.type, parallel, 
                      models.visited){

  p <- length(current.model)

  # Step 1 from page 12 of the probit paper.
  MetHastUpdateGamma <- function(current.gamma, current.active, current.model, 
                                 current.model.score, models.visited){
    # Generate proposal
    if (runif(1) < a) {
      proposed.gamma <- current.gamma
      proposed.gamma[which(current.active == 1)] <- 
        sample(0:1, size=sum(current.active), replace=TRUE)
    }
    else {
      proposed.gamma <- current.gamma
      index.to.shift <- sample(which(current.active == 1), size=1)
      proposed.gamma[index.to.shift] <- 1 - proposed.gamma[index.to.shift]
    }
  
    proposed.model <- rep(0, times = p)
    proposed.model[which((proposed.gamma == 1) & (current.active == 1))] <- 1

    # Compute acceptance probability
    if (identical(proposed.model, current.model)) {
      alpha01 <- 1
      proposed.model.score <- current.model.score
    }
    else {
      proposed.model.base10 <- convert2base10(proposed.model)
      tmp <- match(proposed.model.base10, models.visited[,"modelId"])
      if(!is.na(tmp)) {
        proposed.model.score <- models.visited[tmp, "BF"]
        # cat("Skipped!\n")
      }
      else {
        if (model.type == "probit") {
          proposed.model.score <- BayesFactorProbit(response, 
            covariates.retain, covariates.test, proposed.model, parallel)
        }
        else {
          proposed.model.score <- BayesFactorLinReg(response, 
            covariates.retain, covariates.test, proposed.model)
        }
      }
      models.visited[min(which(is.na(models.visited[,"modelId"]))), ] <- 
        c(convert2base10(proposed.model), proposed.model.score)

      alpha01 <- proposed.model.score/current.model.score
    }

    if (runif(1) < alpha01) {
      current.gamma <- proposed.gamma 
      current.model <- proposed.model
      current.model.score <- proposed.model.score
    }

    return(list(current.gamma, current.model, current.model.score, 
      models.visited))
  }

  # Step 2 from page 13 of the probit paper
  MetHastUpdateModel <- function(current.gamma, current.active, current.model, 
                                 current.model.score, models.visited){
    index.from.A <- sample(which(current.active == 1), size=1)
    index.from.NOT.A <- sample(which(current.active == 0), size=1)

    proposed.active <- current.active
    proposed.active[index.from.A] <- 0
    proposed.active[index.from.NOT.A] <- 1

    proposed.model <- rep(0, times = p)
    proposed.model[which((current.gamma == 1) & (proposed.active == 1))] <- 1

    # Compute acceptance probability
    if (identical(proposed.model, current.model)) {
      alpha02 <- 1
      proposed.model.score <- current.model.score
    }
    else {
      proposed.model.base10 <- convert2base10(proposed.model)
      tmp <- match(proposed.model.base10, models.visited[,"modelId"])
      if(!is.na(tmp)) { # If model has been visited before
        proposed.model.score <- models.visited[tmp, "BF"]
        # cat("Skipped!\n")
      }
      else {
        if (model.type == "probit") {
          proposed.model.score <- BayesFactorProbit(response, 
            covariates.retain, covariates.test, proposed.model, parallel)
        }
        else {
          proposed.model.score <- BayesFactorLinReg(response, 
            covariates.retain, covariates.test, proposed.model)
        }
      }
      models.visited[min(which(is.na(models.visited[,"modelId"]))), ] <- 
        c(convert2base10(proposed.model), proposed.model.score)

      alpha02 <- proposed.model.score/current.model.score
    }

    if (runif(1) < alpha02) {
      current.active <- proposed.active
      current.model <- proposed.model
      current.model.score <- proposed.model.score
    }

    return(list(current.active, current.model, current.model.score, 
      models.visited))
  }

  updated.gamma.model <- MetHastUpdateGamma(current.gamma, current.active, 
    current.model, current.model.score, models.visited)
  current.gamma <- updated.gamma.model[[1]]
  current.model <- updated.gamma.model[[2]]
  current.model.score <- updated.gamma.model[[3]]
  models.visited <- updated.gamma.model[[4]]

  updated.active.model <- MetHastUpdateModel(current.gamma, current.active, 
    current.model, current.model.score, models.visited)
  current.active <- updated.active.model[[1]]
  current.model <- updated.active.model[[2]]
  current.model.score <- updated.active.model[[3]]
  models.visited <- updated.active.model[[4]]

  return(list(current.gamma, current.active, current.model, 
    current.model.score, models.visited))
}
########################################################
