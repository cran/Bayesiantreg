tabla <- function(min, max, p){    
  secue <- seq(min,max,by=min)
  matriz <- matrix(0,nrow=length(secue),1)
  
  for(i in 1:length(secue)){
    matriz[i,] <- pJ2(secue[i],p)
  }
  return(matriz)
}