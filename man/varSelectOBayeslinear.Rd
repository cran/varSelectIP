\name{varSelectOBayeslinear}
\alias{varSelectOBayeslinear}
\title{
Objective Bayes Model Selection and inference for Linear Models using mixtures of g-priors
}
\description{
This function carries out a random walk in order 
to determine the best model, as measured by its posterior probability.
The types of models that this function can handle are linear regression 
models with homoscedastic, independent errors. For full details on the model set-up and the random walk, 
please refer to the paper listed below.
}
\usage{
varSelectOBayeslinear(y,X,X0,
	type='IP',prior.par,EFF=TRUE,
	model.prior.type='Beta-Binomial',model.prior.par=c(1,1),
	RWflag=TRUE,N.draws=10^4,shuf=.2,start,RWpropflag=FALSE,
	inference='both',alpha=0.05,m.considered=100)
}
\arguments{
  \item{y}{
    The vector of response values.
}
  \item{X0}{
    A matrix containing the covariates that should always be 
    retained when searching through all possible models.  These constitute the base model.
}
  \item{X}{
    A matrix containing all the covariates that should  be taken 
    into consideration when searching through all possible models.  These constitute the test covariates.
}
  \item{type}{
	The kind of prior used to compute Bayes' Factors, taking the values 'NP','IP','ZS', or 'HG' for normal prior, intrinsic prior, Zellner-Siow prior, or hyper g-prior.  For the hyper g-prior, the prior on w=1/g is taken to be proportional to w^(-.5)*(b+w)^((a+1)/2).
}
  \item{prior.par}{
    The parameter of the g-prior, when appropriate.  The intrinsic and Zellner-Siow priors have no parameters.  The normal prior has parameter w=1/g, defaulting to 1.  The hyper g, has parameter (a,b), defaulting to (2,1).
}
  \item{EFF}{
    A 0-1 variable denoting whether the prior precision for the coefficients of a model should be scaled by the number of covariates.  Defaults to TRUE and a model with p covariates is scaled by (p+1)/n, if FALSE each model is scaled by 1/n.
}
  \item{model.prior.type}{
    The type of probability distribution for the model space, taking values 'Beta-Binomial', 'Binomial', or 'Uniform'.  The binomial and beta-binomial priors are formed by putting a Bernoulli(p) prior on inclusion, with the beta-binomial placing a Beta(a,b) prior on p.
}
  \item{model.prior.par}{
	Parameter for the model space prior.  For the binomial case it is 0<p<1, defaulting to 0.5.  For the beta-binomial case it is 0<a,b, defaulting to (1,1).
}
  \item{RWflag}{
    A 0-1 variable denoting whether to do a random walk (TRUE, the default) or enumerate the models (FALSE).
}
  \item{N.draws}{
	The number of simulations for the random walk, defaults to 10000.
}
  \item{shuf}{
    The percentage of times that the random walk kernel takes an independent draw from the prior, defaults to 0.2.
}
\item{start}{The starting model for the random walk.  Must be either a vector of 0s and 1s whose length is the number of test covariates, a collapsed string of such 0s and 1s, or a vector of column numbers.}
  \item{RWpropflag}{
    A 0-1 variable denoting whether to compute posterior probabilities by sample averages (TRUE) or by renormalization (FALSE, the default).
}
  \item{inference}{
  	A variable taking value of 'both' (the default), 'selected', 'averaged', or 'none' denoting the kind of posterior credible sets to compute for the regression coefficients.  All sets computed are quantile based.
}
\item{alpha}{The level for (1-alpha) credible sets.}
\item{m.considered}{The number of models to keep in order to for credible sets under model averaging.  If 0<m.considered<=1, this denotes keeping the top models that produce at least m.considered total posterior probability.}
}
\value{
\item{gamma}{Binary strings of the considered models ordered by their posterior probabilities.
}
\item{p}{The number of covariates (including base covariates) in the considered models.
}
\item{logPrior}{The log of the prior probabilities of the considered models.
}
\item{Rsquared}{The coefficients of determination for the considered models.
}
\item{logB}{The log of the Bayes factors of the considered models to the base model.
}
\item{prob}{The posterior probabilities of the considered models.}
\item{selected.model}{The binary string of the selected model.}
\item{selected.model.prob}{The posterior probability of the selected model.}
\item{prob0}{The posterior probabilities that each regression coefficient is 0.}
\item{means}{A matrix with two columns.  The first is the model averaged posterior means of the regression coefficients.  The second is the posterior mean conditioned on the selected model.}
\item{credsets.selected}{Quantile based credible sets for the selected model.  This is a matrix with rows for the regression coefficients and columns for upper and lower bounds.}
\item{credsets.averaged}{Quantile based credible sets under model averaging.  This is a list with an element for each regression coefficient.  The elements of the list are themselves lists containing the intervals of the credible set and the posterior probabilities of the intervals.  The credible sets comprise an open interval and possibly a point mass at 0.}
\item{resides}{A matrix with two columns.  The first column is the residuals from the model averaged mean.  The second column is the residuals from the mean conditioned on the selected model.}
}
\references{
  Womack, A. J. and Leon-Novelo, L. L. and Casella, G. (2013) Inference from Intrinsic Bayes\' Procedures
Under Model Selection and Uncertainty
  \url{http://www.stat.ufl.edu/~ajwomack/WNC-IP.pdf}
}
\author{
Andrew Womack

Maintainer: Andrew Womack <ajwomack@ufl.edu>
}
\examples{
n <- 20 # number of observations
p <- 6 # total number of covariates

gene_expression <- matrix(runif(n*p)*4,nrow=n,ncol=p)
age <- sample(20:40,n,replace=TRUE)

truth_betavector <- c(-0.1, -.01, 1, -1, rep(0,p+2-4)) 
design <- cbind(c(1), age, gene_expression) # sets up the entire design matrix.

# Simulating the y-values and setting up the data-frame
y <- apply(design, 1, function(xi){rnorm(n=1, sum(xi * truth_betavector))})
X <- design[, -c(1:2)]
X0<-design[,1:2]

# Linear regression function call:
varSelectOBayeslinear(y,X,X0,RWflag=FALSE,m.considered=.99)
}
\keyword{ linear models }
\keyword{ regression }
\keyword{ g prior }
