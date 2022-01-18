rJ2 <- function(n, matriz, min, max){
  secue <- seq(min,max,by=min)
  quantil <-runif(n)
  grados <- matrix(0,nrow=length(quantil),ncol=1)
  
  for(i in 1:length(quantil)){
    if (quantil[i]<min(matriz)){
      grados[i,] <- secue[1]
    } else if(quantil[i]>max(matriz)){
      grados[i,] <- secue[length(secue)]        
    } else {
      grados[i,] <- secue[(which(matriz==max(matriz[matriz<quantil[i]]))+1)]
    }
  }
  grados <- as.numeric(grados)
  return(grados)
}
