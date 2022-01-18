\name{Bayesiantreg}
\alias{Bayesiantreg}
\title{Bayesiantreg}
\description{Function to do Bayesian t Regression: joint mean and variance modeling and estimation of the degrees of freedom}
\usage{
Bayesiantreg(y, x, z, nsim, bini, bpri, Bpri, gini, gpri,Gpri, glini, glpri, 
           type, apriori, propuesta, Maxi=NULL,
           lambda = NULL, p = NULL, burn, jump, graph1 = TRUE, graph2 = TRUE,
           graph3 = TRUE)
}
\arguments{
  \item{y}{object of class matrix, with the dependent variable.}
  \item{x}{object of class matrix, with the variables for modelling the mean.}
  \item{z}{object of class matrix, with the variables for modelling the precision.}
  \item{nsim}{a number that indicate the number of iterations.}
  \item{bini}{a vector with the initial values of beta.}
  \item{bpri}{a vector with the values of the mean of the prior of beta.}
  \item{Bpri}{a matrix with the values of the variance of the prior of beta.}
  \item{gini}{a vector with the initial values of gamma.}
  \item{gpri}{a vector with the values of the mean of the prior of gamma.}
  \item{Gpri}{a matrix with the values of the variance of the prior of gamma.}
  \item{glini}{a vector with the initial value of the degrees of freedom.}
  \item{glpri}{a vector with the value of the the prior of the degrees of freedom.}
  \item{type}{a vector that can take the value "D" if the prior for the degrees of freedom considered 
	      as discrete or "C" if it is continuous.}
  \item{apriori}{when type is "D", it is a vector that can take the values of "poi" for a Poisson prior or 
                 "unif" for a uniform prior. When type is "C", it is a vector that can take the values of 
		 "exp" for the exponential prior, "unif" for the uniform prior or "J2" for the Jeffrey's 
		 prior.}
  \item{propuesta}{when type is "D", it is a vector that can take the values of "poi" for a Poisson proposal, 
                 "unif" for a uniform proposal or by default the proposal made by Marin and Cepeda (_). 
                 When type is "C", it is a vector that can take the values of "exp" for the exponential 
		 proposal, "unif" for the uniform proposal, "J2" for the Jeffrey's proposal or by default 
		 the proposal made by Marin and Cepeda (_).}
  \item{Maxi}{a number indicating the maximum value for the uniform prior an the uniforme proposal.}
  \item{lambda}{a number indicating the mean parameter value for the Poisson prior an the Poisson proposal.}
  \item{p}{a number indicating the parameter value for the Jeffrey's prior an the Jeffrey's proposal.}
  \item{burn}{a proportion that indicate the number of iterations to be burn at the beginning of the chain.}
  \item{jump}{a number that indicate the distance between samples of the autocorrelated the chain, 
	      to be excluded from the final chain.}
  \item{graph1}{if it is TRUE present the graph of the chains without jump and burn.}
  \item{graph2}{if it is TRUE present the graph of the chains with jump and burn.}
  \item{graph3}{if it is TRUE present the graph of the standardized residuals, the the standardized residuals against 
		the lineal predictor, the pseudo deviance residuals and the pseudo deviance residuals against the 
		lineal predictor.}
 }
\details{
The bayesian t regression allows the joint modelling of mean and variance and the estimation of the degrees of freedom 
of a t distributed variable, as is proposed in Marin and Cepeda (_), with identical link for the mean and logarithmic for the variance,
and differents discrete and continuous aproach for the degrees of freedom. 
}
\value{
 object of class bayesbetareg with:
 \item{coefficients}{object of class matrix with the estimated coefficients of beta, gamma and degrees of freedom.}
 \item{interv}{object of class matrix with the estimated confidence intervals of beta, gamma and the degrees of freedom.}
 \item{fitted.values}{object of class matrix with the fitted values of y.}
 \item{residuals}{object of class matrix with the residuals of the regression.}
 \item{residualsstd}{object of class matrix with the standardized residuals of the regression.}
 \item{residualsdev}{object of class matrix with the pseudo deviance residuals of the regression.}
 \item{variance}{object of class matrix with the variance terms of the regression.}
 \item{beta.mcmc}{object of class matrix with the complete chains for beta.}
 \item{gamma.mcmc}{object of class matrix with the complete chains for gamma.}
 \item{gl.mcmc}{object of class matrix with the complete chains for the degrees of freedom.}
 \item{beta.mcmc.burn}{object of class matrix with the chains for beta after the burned process.}
 \item{gamma.mcmc.burn}{object of class matrix with the chains for gamma after the burned process.}
 \item{gl.mcmc.burn}{object of class matrix with the chains for the degreees of freedom after the burned process.}
 \item{loglik}{the logaritmic of the liklihood of the model.}
 \item{AIC}{AIC of the model.}
 \item{BIC}{BIC of the model.}
 \item{DIC}{BIC of the model.}
 \item{PseudoDeviance}{Pseudo deviance criteria of the model as is proposed by Marin and Cepeda (_).}
 \item{arb}{acceptance percentage for beta.}
 \item{arg}{acceptance percentage for gamma.}
 \item{argl}{acceptance percentage for the degrees of freedom.}
 \item{call}{Call.} 
}
\references{
1. Marin and Cepeda-Cuervo (_). A Bayesian regression model for the non-standardized t distribution with location, 
scale and degrees of freedom parameters. Unpublished}
\author{
 Margarita Marin \email{mmarinj@unal.edu.co},
 Edilberto Cepeda-Cuervo \email{ecepedac@unal.edu.co}  

2. Cepeda-Cuervo E. (2001). Modelagem da variabilidade em modelos lineares generalizados. Unpublished
Ph.D. tesis. Instituto de Matematicas. Universidade Federal do Rio do Janeiro.

3. Cepeda C., E. and Gamerman D. (2001). Bayesian Modeling of Variance Heterogeneity in Normal
Regression Models. Brazilian Journal of Probability and Statistics. 14, 207-221
      
}

\examples{


n <- 10
X1 <- runif(n,0,10)
X2 <- runif(n,5,10)
X3 <- runif(n,10,15)


y1 <- c(0.09, 1.68, -2.43, 0.23, 2.94, 1.50, 3.40, 2.22, 0.28, -0.17)

betas <- c(0,0,0,0)
gammas <- c(0,0,0)
gl <- 3

x <- cbind(rep(1,n),X1,X2,X3)
z <- cbind(rep(1,n),X2,X3)
y <- y1

Bpri <- diag(rep(100,4))
bpri <- rep(0,4)

Gpri <- diag(rep(10,3))
gpri <- rep(0,3)

glpri <- 7

propuesta <- "unif2"
apriori <- "unif"
tipo <- "D"

Maxi <- 100
nsim <- 50

bini=bpri
gini=gpri
glini=glpri

reg1 <- Bayesiantreg(y, x, z, nsim=nsim, bini, bpri, 
                     Bpri, gini, 
                     gpri,Gpri, glini, glpri, 
                     type=tipo, apriori=apriori, 
                     propuesta=propuesta, 
                     Maxi=Maxi,burn=0.3, jump=3, 
                    graph1 = TRUE, graph2 = TRUE, graph3 = TRUE)

}
\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Meancovariance modelling}
\keyword{t}
\keyword{Degrees of Freedom}


