\name{summary.Bayesiantreg}
\alias{summary.Bayesiantreg}
\title{summary of the Bayesian t regression}
\description{Summarized the Bayesian Bayesian t Regression: joint mean and variance modeling and estimation of the degrees of freedom}
\usage{
\method{summary}{Bayesiantreg}(object, \ldots)
}
\arguments{
  \item{object}{an object of class Bayesiantreg}
   \item{\ldots}{not used.}
}
\value{
\item{call}{Call}
\item{coefficients}{Coefficients.}
 \item{AIC}{AIC of the model.}
 \item{BIC}{BIC of the model.}
 \item{DIC}{BIC of the model.}
 \item{PseudoDeviance}{Pseudo deviance criteria of the model as is proposed by Marin and Cepeda (_).}
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
