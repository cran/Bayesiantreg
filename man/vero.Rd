\name{vero}
\alias{vero}
\title{Loglikelihood for the t model}
\description{Calculate the loglikelihood for the t model}
\usage{
vero(y, mu, sigma2, grados)
}
\arguments{
  \item{y}{object of class matrix, with the dependent variables.}
  \item{mu}{object of class matrix, with the mean of the model.}
  \item{sigma2}{object of class matrix, with the variace of the model.}
  \item{grados}{a vector with the degrees of freedom of the model.}

}
\details{
Calculate the loglikelihood for the t model as proposed by Marin and Cepeda (_).
}
\value{
\item{l}{a value with the loglikelihood for the t model}
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