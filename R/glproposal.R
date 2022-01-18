glproposal<-function(gl.ini, lambda, p, Maxi, matriz, propuesta, type) {
  
  if(type=="D"){
    if(propuesta=="poi"){
      gl.pro <- rpois(1,lambda)+3
    } else if(propuesta=="unif"){
      gl.pro <- round(runif(1,3,Maxi))
    } else{
      if(gl.ini==3){
        gl.pro <- round(runif(1,3,5))
      } else{
        gl.pro <- round(runif(1,(gl.ini-1),(gl.ini+1)))    
      }
    } 
    
  } else{
    
    if(propuesta=="exp"){
      gl.pro <- rexp(1,lambda)+2.0001
    } else if(propuesta=="J2"){
      gl.pro <- rJ2(1,matriz,0.001,100)+2.0001
      #a <- 3
    } else if(propuesta=="unif"){
      gl.pro <- runif(1,2.0001,Maxi)
    } else{
      if (gl.ini<4.0001){gl.ini <- runif(1,2.0001,6)}else{
        gl.pro <- runif(1,(gl.ini-2),(gl.ini+2))
      }
    }
  }
  
}
