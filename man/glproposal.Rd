\name{glproposal}
\alias{glproposal}
\title{A proposal for degrees of freedom parameter}
\description{Propose a value for the degrees of freedom parameter}
\usage{
glproposal(gl.ini, lambda, p, Maxi, matriz, propuesta, type)
}
\arguments{
  \item{gl.ini}{a vector with the previous proposal degrees of freedom parameter}
  \item{lambda}{a number indicating the mean parameter value for the Poisson prior an the Poisson proposal.}
  \item{p}{a number indicating the parameter value for the Jeffrey's prior an the Jeffrey's proposal.}
  \item{Maxi}{a number indicating the maximum value for the uniform prior an the uniforme proposal.}
  \item{matriz}{a matrix generate by the function tabla of the bayesiantreg package.}
  \item{propuesta}{when type is "D", it is a vector that can take the values of "poi" for a Poisson proposal, 
                 "unif" for a uniform proposal or by default the proposal made by Marin and Cepeda (_). 
                 When type is "C", it is a vector that can take the values of "exp" for the exponential 
		 proposal, "unif" for the uniform proposal, "J2" for the Jeffrey's proposal or by default 
		 the proposal made by Marin and Cepeda (_).}
  \item{type}{a vector that can take the value "D" if the prior for the degrees of freedom considered 
	      as discrete or "C" if it is continuous.}

}
\details{Generate a proposal for the gamma parameter according to the model proposed by Marin and Cepeda-Cuervo (_).}
\value{
 \item{gl.pro}{a number with the proposal for the degrees of freedom parameter.}
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