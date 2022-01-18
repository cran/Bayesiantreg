gammakernel <- function(y, x, z, betas.ini, gammas.now, gammas.old, gl.ini, gpri, Gpri) {
  
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
  
  
  mu <- x%*%betas.ini
  sigma2 <- exp(z%*%gammas.old)
  gl <- gl.ini
  y.mono1 <- z%*%gammas.old - 1 + ((gl-2)/gl)*(y-mu)^2/sigma2
  ifelse(gl>4,  vary.now1 <- 2*(gl-1)/(gl-4), vary.now1 <- 8)
  Gpos <- ginv(ginv(Gpri)+ t(z)%*%ginv(diag(as.vector(vary.now1),nrow(z),nrow(z)))%*%z)
  Gpos <- as.matrix(forceSymmetric(as.matrix(Gpos)))
  gpos <-Gpos%*%(ginv(Gpri)%*%gpri + t(z)%*%ginv(diag(as.vector(vary.now1),nrow(z),nrow(z)))%*%y.mono1)
  dmvnorm(t(gammas.now),gpos,Gpos) #These functions provide the density function for the multivariate normal
  # distribution with mean equal to mean and covariance matrix sigma.
}