pJ2 <- function(gl.ini, p){
  a <-p
  dJ3 <- function(gl.ini){dJ2(gl.ini,a)}
  J1I <- integrate(dJ3,0,gl.ini)$value/integrate(dJ3,0,Inf)$value
  return(J1I)
}