\name{dpostg}
\alias{dpostg}
\title{Posterior density of gamma}
\description{Propose a value for posterior density of the gamma parameter}
\usage{
dpostg(y, x, z, betas, gammas, gl, gpri, Gpri)
}
\arguments{
  \item{y}{object of class matrix, with the dependent variables.}
  \item{x}{object of class matrix, with the variables for modelling the mean.}
  \item{z}{object of class matrix, with the variables for modelling the variance.}
  \item{betas}{a vector with the proposal beta parameters.}
  \item{gammas}{a vector with the proposal gamma parameters.}
  \item{gl}{a vector with the proposal degrees of freedom parameter.}
  \item{gpri}{a vector with the values of the mean of the prior of gamma.}
  \item{Gpri}{a matrix with the values of the variance of the prior of gamma.}
}
\details{
Generate the posterior density for the gamma proposed by Marin and Cepeda (_).
}
\value{
\item{value}{a value with the posterior denity for gamma}
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
