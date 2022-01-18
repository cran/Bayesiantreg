dpostg <- function(y, x, z, betas, gammas, gl, gpri, Gpri) {
  
  if(is.null(x)|is.null(z)){
    stop("There is no data")
  }
  
  if(nrow(x)!=nrow(z)){
    stop("The mean and precision data have a different size")
  }
  
  if(!is.vector(gpri)){
    stop("the initial parameters for beta must be a vector")
  }
  
  if(ncol(Gpri)!=nrow(Gpri)){
    stop("The initial covariance matrix for beta is not square")
  }
  
  
  mu <- x%*%betas
  sigma2 <- as.vector(exp(z%*%gammas))
  gl <- gl
  n <- nrow(y)
  vero <- (gamma((gl+1)/2)/(gamma(1/2)*gamma(gl/2)*(gl*sigma2)^0.5))*(1+(y-mu)^2/(sigma2*gl))^(-(gl+1)/2)  
  L <- prod(vero)
  P <- dmvnorm(t(gammas), gpri, Gpri) # priori
  value=L*P
  value
}