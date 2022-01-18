dJ2 <- function(gl.ini, p){
  a <- (gl.ini/(gl.ini+3))^0.5
  b <- trigamma(gl.ini/2)-trigamma((gl.ini+1)/2) - (2*(gl.ini+3)/(gl.ini*(gl.ini+1)^2))
  c <- (gl.ini+1)/(gl.ini+3)
  J1 <- a*(b^0.5)*(c^(p/2))
  return(J1)
}