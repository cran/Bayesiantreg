mukernel <- function(y, x, z, betas.now, betas.old, gammas.ini, gl.ini, bpri, Bpri) {
  
  if(is.null(x)|is.null(z)){
    stop("There is no data")
  }
  
  if(nrow(x)!=nrow(z)){
    stop("The mean and precision data have a different size")
  }
  
  if(!is.vector(bpri)){
    stop("the initial parameters for beta must be a vector")
  }
  
  if(ncol(Bpri)!=nrow(Bpri)){
    stop("The initial covariance matrix for beta is not square")
  }
  
  
  sigma2 = exp(z%*%gammas.ini)
  gl = gl.ini
  Tao <- gl*sigma2/(gl-2)
  Bpos <- ginv(ginv(Bpri)+ t(x)%*%ginv(diag(as.vector(Tao)))%*%x)
  Bpos <- as.matrix(forceSymmetric(as.matrix(Bpos)))
  bpos <- Bpos%*%(ginv(Bpri)%*%bpri + t(x)%*%ginv(diag(as.vector(Tao)))%*%y)
  dmvnorm(t(betas.now),bpos,Bpos) #These functions provide the density function for the multivariate normal 
  #distribution with mean equal to mean and covariance matrix sigma.
}