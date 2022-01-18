glpost <- function(y, x, z, betas.ini, gammas.ini, gl.ini, Maxi, lambda, p, prior, type) {
  
  if(is.null(x)|is.null(z)){
    stop("There is no data")
  }
  
  if(nrow(x)!=nrow(z)){
    stop("The mean and precision data have a different size")
  }
  
  mu <- x%*%betas.ini
  sigma2 <- exp(z%*%gammas.ini)
  grados <- gl.ini
  
  if(type=="D"){
    
    if(prior=="pois"){
      grados <- grados
      P <- dpois(grados,lambda)
      grados <- grados
    } else {
      P <- dunif(grados, 3, Maxi) # priori hay q poner las otras opciones
    } 
    
  } else{
    
    if(prior=="exp"){
      P <- dexp(grados,lambda) # priori hay q poner las otras opciones
    } else if(prior=="J2"){
      P <- dJ2(grados,p)
    } else {
      P <- dunif(grados, 3, Maxi) # priori hay q poner las otras opciones 
    } 
  }
  
  vero <- (gamma((grados+1)/2)/(gamma(1/2)*gamma(grados/2)*(grados*sigma2)^0.5))*(1+(y-mu)^2/(sigma2*grados))^(-(grados+1)/2)
  L <- prod(vero)
  value=L*P
  value
} 
