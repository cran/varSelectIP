\name{varSelectIP}
\alias{varSelectIP}
\title{
Objective Bayes Model Selection
}
\description{
This function will carry out a low-dimensional stochastic search in order 
to determine the ``best'' model, as measured by its posterior probability.
The types of model that this function can handle are probit and regression 
models. For full details on the model set-up and the stochastic search, 
please refer to the papers listed below.
}
\usage{
varSelectIP(response, covariates.retain = NULL, covariates.test, nsim, 
            keep, q, a = 0.2, model.type = c("probit", "reg"), 
            save.every = 50, out.fname = "models.csv")
}
\arguments{
  \item{response}{
    The vector of response values. If a probit model, this should be a binary 
    vector with the 0's coming before the 1's.
}
  \item{covariates.retain}{
    A matrix or a vector containing the covariates that should always be 
    retained when searching through all possible models.
}
  \item{covariates.test}{
    A matrix or a vector containing all the covariates that should  be taken 
    into consideration when searching through all possible models.
}
  \item{nsim}{
    The number of iterations of the stochastic search to run through.
}
  \item{keep}{
    The final number of models to report, along with their Bayes Factors.
}
  \item{q}{
    The maximum number of covariates to be included in each model considered. 
    These covariates will be chosen out of those in covariates.test above.
}
  \item{a}{
    The probability with which the entire set of active coefficients are 
    re-drawn. See page 12 of reference (2) for more details.
}
  \item{model.type}{
    This has to be either "probit" or "reg", specifying the type of model 
    to be fit.
}
  \item{save.every}{
    Specifies how often the models should be written out to a .csv file. This 
    allows a user to monitor progress of models found and to prevent loss of 
    effort in the case of power failure, etc.
}
  \item{out.fname}{
    The name of the .csv file to save models to.
}
}
\value{
A table containing the top models found is returned. Each row in the table 
represents a model. A 1 within a row indicates that that covariate, from 
within the covariates defined by covariates.test, is included in the model. A 
0 indicates that that particular covariate is left out of the model. The 
last column contains the Bayes Factor from comparing the specified model 
against the intercept-only model.
}
\references{
  Casella, G. and Giron, F.J. and Martinez, M.L. and Moreno, E. (2009) 
  Consistency of Bayesian Procedures for Variable Selection. _Annals of 
  Statistics_, *37*, 1207 - 1228.

  Leon-Novelo, L. and Moreno, E. and Casella, G. (2010) Objective Bayes 
  Model Selection in Probit Models.
  \url{http://www.stat.ufl.edu/~casella/Papers}
}
\author{
Gopal, V. and Novelo, L. L. and Casella, G.

Maintainer: Vik Gopal <viknesh.gopal@gmail.com>
}
\examples{
n <- 20 # number of observations
p <- 6 # total number of covariates
set.seed(0)

gene_expression <- matrix(runif(n*p)*4,nrow=n,ncol=p)
age <- sample(20:40,n,replace=TRUE)

truth_betavector <- c(-0.1, -.01, 1, -1, rep(0,p+2-4)) 
design <- cbind(1, age, gene_expression) # sets up the entire design matrix.

# Simulating the z-values and y-values and setting up the data-frame
y_tmp <- apply(design, 1, function(xi){rnorm(n=1, sum(xi * truth_betavector))})
y <- y_tmp[order(y_tmp)]
x <- design[order(y_tmp), -c(1:2)]
n0 <- sum(y<0)
n1 <- n-n0
z <- c(rep(0,n0),rep(1,n1))
mydata <- cbind(z, y, age, x)
colnames(mydata)<-c("z", "y", "age", paste("GE",1:p,sep=""))

# Linear regression function call:
varSelectIP(a=0.2, keep=32, covariates.retain=mydata[,3], model.type="reg", 
            q=5, covariates.test=mydata[,4:9], response=mydata[,2], nsim=25)

##NOT RUN
# Probit regression function call:
#varSelectIP(a=0.2, keep=2, covariates.retain=mydata[,3], model.type="probit", 
#           q=3, covariates.test=mydata[,4:7], response=mydata[,1], nsim=2)
}
\keyword{ models }
\keyword{ htest }
\keyword{ regression }
