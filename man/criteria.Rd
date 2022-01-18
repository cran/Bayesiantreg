\name{criteria}
\alias{criteria}
\title{criteria for comparison the bayesian t regression models}
\description{Performs the comparison criterias for the Bayesian t Regression}
\usage{
criteria(object,...)
}
\arguments{

 \item{object}{object of class "Bayesiantreg"}
 \item{\ldots}{not used.}
}
\details{
This function allows to extract the information criteria from the model AIC, BIC, DIC and pseudo-deviance. 
}
\value{
 \item{loglik}{the logaritmic of the liklihood of the model}
 \item{AIC}{the AiC criteria}
 \item{BIC}{the BIC criteria}
 \item{DIC}{the DIC criteria}
 \item{PseudoDeviance}{the pseudo deviance criteria of the model as is proposed in Marin and Cepeda (_).}
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
 Edilberto Cepeda-Cuervo \email{ecepedac@unal.edu.co},     
}

\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Meancovariance modelling}
\keyword{t}
\keyword{Degrees of Freedom}
