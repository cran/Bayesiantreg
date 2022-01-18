\name{glpost}
\alias{glpost}
\title{Posterior value of the degrees of freedom}
\description{Calculate a value for posterior density of the degrees of freedom parameter}
\usage{
glpost(y, x, z, betas.ini, gammas.ini, gl.ini, Maxi, lambda, p, prior, type)
}
\arguments{
  \item{y}{object of class matrix, with the dependent variables.}
  \item{x}{object of class matrix, with the variables for modelling the mean.}
  \item{z}{object of class matrix, with the variables for modelling the variance.}
  \item{betas.ini}{a vector with the proposal beta parameters.}
  \item{gammas.ini}{a vector with the proposal gamma parameters.}
  \item{gl.ini}{a vector with the proposal degrees of freedom parameter.}
  \item{Maxi}{a number indicating the maximum value for the uniform prior an the uniforme proposal}
  \item{lambda}{a number indicating the mean parameter value for the Poisson prior an the Poisson proposal}
  \item{p}{a number indicating the parameter value for the Jeffrey's prior an the Jeffrey's proposal}
  \item{type}{a vector that can take the value "D" if the prior for the degrees of freedom considered 
	      as discrete or "C" if it is continuous.}
  \item{prior}{when type is "D", it is a vector that can take the values of "poi" for a Poisson prior or 
                 "unif" for a uniform prior. When type is "C", it is a vector that can take the values of 
		 "exp" for the exponential prior, "unif" for the uniform prior or "J2" for the Jeffrey's 
		 prior.}

}
\details{
Generate the posterior density for the degrees of freedom proposed by Marin and Cepeda (_).
}
\value{
\item{value}{a value with the posterior denity for the degrees of freedom}
}
\references{
1. Marin and Cepeda-Cuervo (_). A Bayesian regression model for the non-standardized t distribution with location, 
scale and degrees of freedom parameters. Unpublished

2. Cepeda-Cuervo E. (2001). Modelagem da variabilidade em modelos lineares generalizados. Unpublished
Ph.D. tesis. Instituto de Matematicas. Universidade Federal do Rio do Janeiro.

3. Cepeda C., E. and Gamerman D. (2001). Bayesian Modeling of Variance Heterogeneity in Normal
Regression Models. Brazilian Journal of Probability and Statistics. 14, 207-221
}

\author{
 Margarita Marin \email{mmarinj@unal.edu.co},
 Edilberto Cepeda-Cuervo \email{ecepedac@unal.edu.co}        
}
\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Meancovariance modelling}
\keyword{t}
\keyword{Degrees of Freedom}