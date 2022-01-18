\name{gammakernel}
\alias{gammakernel}
\title{the probability of a gamma parameter from the probability density funcion defined by old parameters.}
\description{evaluate the probability of a gamma parameter from the probability density function defined by old parameters.}
\usage{
gammakernel(y, x, z, betas.ini, gammas.now, gammas.old, gl.ini, gpri, Gpri)
}
\arguments{
  \item{y}{object of class matrix, with the dependent variable}
  \item{x}{object of class matrix, with the variables for modelling the mean}
  \item{z}{object of class matrix, with the variables for modelling the variance}
  \item{betas.ini}{a vector with the beta parameters that define the old p.d.f}
  \item{gammas.now}{a vector with the gamma parameters - new parameters - to evaluate in the old p.d.f}
  \item{gammas.old}{a vector with the gamma parameters that define the old p.d.f}
  \item{gl.ini}{a vector with the degrees of freedom parameters that define the old p.d.f}
  \item{gpri}{a vector with the initial values of gamma}
  \item{Gpri}{a matrix with the initial values of the variance of gamma}
}
\details{
Evaluate the probability of a gamma parameter from the probability density function defined by old parameters, according with the model proposed by Marin and Cepeda-Cuervo (_). 
}
\value{ 
 \item{value}{a vector with the probability for the gamma parameter from the probability density function defined by old parameters.}
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
