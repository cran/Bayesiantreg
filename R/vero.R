vero <- function(y, mu, sigma2, grados){
  l <- prod((gamma((grados+1)/2)/(gamma(1/2)*gamma(grados/2)*(grados*sigma2)^0.5))*(1+(y-mu)^2/(sigma2*grados))^(-(grados+1)/2))
  l <- log(l)
  return(l)
}