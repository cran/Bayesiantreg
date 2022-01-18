BayesiantregEst <-
  function(y, x, z, nsim, bini, bpri, Bpri, gini, gpri, Gpri, glini, glpri, 
           type, apriori, propuesta, Maxi=NULL,
           lambda = NULL, p = NULL, burn, jump, graph1 = TRUE, graph2 = TRUE,
           graph3 = TRUE){
    
    if(is.null(x)|is.null(z)|is.null(y)){
      stop("There is no data")
    }

    if(burn> 1 | burn < 0){
      stop("Burn must be a proportion between 0 and 1")
    }

    if(nsim <= 0){
      stop("the number of simulations must be greater than 0")
    }

    if(jump < 0|jump > nsim){
      stop("Jumper must be a positive number lesser than nsim")
    }
  
    #cambiamos la dimensión de y, x, x
      
    y <- as.matrix(y)
    
    #cambioy <- nchar(max(round(y)))-2
    if(nchar(max(round(y)))==1){
      y <- as.vector(y)
      cambioy <- 0
    }else{
      y <- y/10^(nchar(max(round(y)))-2)
      y <- as.vector(y)
      cambioy <- nchar(max(round(y)))-2
    }

    cambiox <- matrix(0,1,ncol(x))
    for (i in 1:ncol(x)){
      if(nchar(max(round(x[,i])))==1){
        x[,i] <- x[,i]
        cambiox[,i] <- 0 
      } else {
        cambiox[,i] <- nchar(max(round(x[,i])))-2    
        x[,i] <- x[,i]/10^(nchar(max(round(x[,i])))-2)
      }
    }
    
    cambioz <- matrix(0,1,ncol(z))
    for (i in 1:ncol(z)){
      
      if(nchar(max(round(z[,i])))==1){
        z[,i] <- z[,i]
        cambioz[,i] <- 0 
      } else {
        cambioz[,i] <- nchar(max(round(z[,i])))-2   
        z[,i] <- z[,i]/10^(nchar(max(round(z[,i])))-2)
      }
    }
 
    #aqui empieza el algoritmo
       
    ind1<-rep(0,nsim)
    ind2<-rep(0,nsim)
    ind3<-rep(0,nsim)
    
    if (is.null(bini)){
      betas.ind <- matrix(bpri,nrow=ncol(x))
    }else{
      betas.ind <- matrix(bini,nrow=ncol(x))
    }
    
    if (is.null(gini)){
      gammas.ind <- matrix(gpri,nrow=ncol(z))
    }else{
      gammas.ind <- matrix(gini,nrow=ncol(z))
    }
    
    if (is.null(glini)){
      #gl.ind <- matrix(glpri,nrow=1)
      gl.ind <- glpri
    }else{
      #gl.ind <- matrix(glini,nrow=1)
      gl.ind <- glini
    }
    
    beta.mcmc<-matrix(NA,nrow=nsim,ncol=ncol(x))
    gamma.mcmc<-matrix(NA,nrow=nsim,ncol=ncol(z))  
    gl.mcmc<-matrix(NA,nrow=nsim,ncol=1)
    
    prior <- apriori
    
    if(prior=="J2"){
      tablaa <- tabla(0.001,100,p)
    }
    
    for(i in 1:nsim) {
      
      #Betas
      
      betas.sim <- matrix(muproposal(y,x,z,betas.ind,gammas.ind,gl.ind,bpri,Bpri),nrow=ncol(x))
      gammas.sim <- matrix(gammaproposal(y,x,z,betas.ind,gammas.ind,gl.ind,gpri,Gpri),nrow=ncol(z))
      gl.sim <- glproposal(gl.ind,lambda,p,Maxi,tablaa, propuesta, type=type)
      
      q1.mu <- mukernel(y,x,z,betas.sim,betas.ind,gammas.sim,gl.sim,bpri,Bpri)
      q2.mu <- mukernel(y,x,z,betas.ind,betas.sim,gammas.ind,gl.ind,bpri,Bpri)
      p1.mu<-dpostb(y,x,z,betas.sim,gammas.ind,gl.ind,bpri,Bpri)
      p2.mu<-dpostb(y,x,z,betas.ind,gammas.ind,gl.ind,bpri,Bpri)
      
      alfa1<-((p1.mu/p2.mu)*(q1.mu/q2.mu))
      
      if(is.na(alfa1)==TRUE |alfa1==Inf){
        alfa1=10
      }
      
      Mu.val<-min(1,alfa1)
      u<-runif(1)
      if (u <=Mu.val) {
        betas.ind <- betas.sim
        ind1[i] = 1
      }
      
      beta.mcmc[i,]<-betas.ind
      beta.mcmc <- as.ts(beta.mcmc)
      
      q1.gamma <- gammakernel(y,x,z,betas.sim,gammas.sim,gammas.ind,gl.sim,gpri,Gpri)
      q2.gamma <- gammakernel(y,x,z,betas.ind,gammas.ind,gammas.sim,gl.ind,gpri,Gpri)
      p1.gamma<-dpostg(y,x,z,betas.ind,gammas.sim,gl.ind,gpri,Gpri)
      p2.gamma<-dpostg(y,x,z,betas.ind,gammas.ind,gl.ind,gpri,Gpri)  
      
      
      alfa2<-((p1.gamma/p2.gamma)*(q1.gamma/q2.gamma))
      
      if(is.na(alfa2)==TRUE |alfa2==Inf){
        alfa2=10
      }
      
      Gamma.val<-min(1,alfa2)
      u<-runif(1)
      if (u <=Gamma.val) {
        gammas.ind <- gammas.sim
        ind2[i] = 1
      }
      gamma.mcmc[i,]<-gammas.ind
      gamma.mcmc <- as.ts(gamma.mcmc)
      
      
      p1.gl<-glpost(y, x, z,betas.ind,gammas.ind,gl.sim,Maxi,lambda,p,apriori,type=type)
      p2.gl<-glpost(y, x, z,betas.ind,gammas.ind,gl.ind,Maxi,lambda,p,apriori,type=type)    
      
      alfa3<-(p1.gl/p2.gl)
      
      if(is.na(alfa3)==TRUE |alfa3==Inf){
        alfa3=10
      }
      
      Gl.val<-min(1,alfa3)
      u<-runif(1)
      if (u <=Gl.val) {
        gl.ind <- gl.sim
        ind3[i] = 1
      }
      gl.mcmc[i,]<-gl.ind
      gl.mcmc <- as.ts(gl.mcmc)
      
      
      if (i%%100 == 0)
        cat("Burn-in iteration : ", i, "\n")
    }
    
    #volvemos las variables a su dimensión original
    
    for (i in 1:ncol(beta.mcmc)){
      beta.mcmc[,i] <- beta.mcmc[,i]*(10^cambioy)/(10^cambiox[,i])
    }
    
    gamma.mcmc[,1] <- log(10^(2*cambioy))+gamma.mcmc[,1]
    
    for (i in 2:ncol(gamma.mcmc)){
      gamma.mcmc[,i] <- gamma.mcmc[,i]/(10^cambioz[,i])
    }   
   
    #Beta y Gamma estimations
    Bestimado <- colMeans(beta.mcmc)
    Gammaest <- colMeans(gamma.mcmc)
    #Glest <- mfv(gl.mcmc)
    Glestm1 <- mean(gl.mcmc)
    Glestm2 <- median(gl.mcmc)
    

    
    #Credibility intervals for beta
    B1 <- matrix(0, ncol(x),1)
    B2 <- matrix(0, ncol(x),1)
    
    
    for(i in 1:ncol(x)){
      B1[i,]<-quantile(beta.mcmc[,i],0.025)
      B2[i,]<-quantile(beta.mcmc[,i],0.975)
      B <- cbind(B1,B2)
    }
    
    # Credibility intervals for gamma
    
    G1 <- matrix(0, ncol(z),1)
    G2 <- matrix(0, ncol(z),1)
    
    for(i in 1:ncol(z)){
      G1[i,]<-quantile(gamma.mcmc[,i],0.025)
      G2[i,]<-quantile(gamma.mcmc[,i],0.975)
      G <- cbind(G1,G2)
    }
    
    U1 <- matrix(0, 1,1)
    U2 <- matrix(0, 1,1)
    
    U1<-quantile(gl.mcmc,0.025)
    U2<-quantile(gl.mcmc,0.975)
    U <- cbind(U1,U2)
    
    burn <- burn
    beta.mcmc.burn <- as.ts(beta.mcmc[-c(1:(nsim*burn)),])
    gamma.mcmc.burn <- as.ts(gamma.mcmc[-c(1:(nsim*burn)),])
    gl.mcmc.burn <- as.ts(gl.mcmc[-c(1:(nsim*burn))])
    
    if(jump>0){
      salta <- seq(1, nrow(beta.mcmc.burn), by=jump)
      beta.mcmc.burn <- as.ts(beta.mcmc.burn[salta,])
      gamma.mcmc.burn <- as.ts(gamma.mcmc.burn[salta,])
      gl.mcmc.burn <- as.ts(gl.mcmc[salta])
      
    }
    
    colMeans(beta.mcmc.burn)
    colMeans(gamma.mcmc.burn)
    median(gl.mcmc.burn)
    mean(gl.mcmc.burn)
    
    if (graph1==TRUE) {
      
      for(i in 1:ncol(x)){
        dev.new()
        ts.plot(beta.mcmc[,i], main=paste("Complete chain for beta",i), xlab="number of iterations", ylab=paste("parameter beta",i))
      }
      
      for(i in 1:ncol(z)){
        dev.new()
        ts.plot(gamma.mcmc[,i], main=paste("Complete chain for gamma",i), xlab="number of iterations", ylab=paste("parameter gamma",i))
        
      }
      
      dev.new()
      ts.plot(gl.mcmc, main="Complete chain for degrees of freedom", xlab="number of iterations", ylab="parameter for the degrees of freedom")
      
    } else{
    }
    
    
    if (graph2==TRUE) {
      
      for(i in 1:ncol(x)){
        dev.new()
        ts.plot(beta.mcmc.burn[,i], main=paste("Burn chain for beta",i), xlab="number of iterations", ylab=paste("parameter beta",i))
        
      }
      
      for(i in 1:ncol(z)){
        dev.new()
        ts.plot(gamma.mcmc.burn[,i], main=paste("Burn chain for gamma",i), xlab="number of iterations", ylab=paste("parameter gamma",i))
        
      }
      
      ts.plot(gl.mcmc.burn, main="Burn chain for degrees of freedom", xlab="number of iterations", ylab="parameter for the degrees of freedom")
      
    } else{
    }
    
    #Credibility intervals for beta
    B1 <- matrix(0, ncol(x),1)
    B2 <- matrix(0, ncol(x),1)
    
    
    for(i in 1:ncol(x)){
      B1[i,]<-quantile(beta.mcmc.burn[,i],0.025)
      B2[i,]<-quantile(beta.mcmc.burn[,i],0.975)
      Bburn <- cbind(B1,B2)
    }
   
    # Credibility intervals for gamma
    
    G1 <- matrix(0, ncol(z),1)
    G2 <- matrix(0, ncol(z),1)
    
    for(i in 1:ncol(z)){
      G1[i,]<-quantile(gamma.mcmc.burn[,i],0.025)
      G2[i,]<-quantile(gamma.mcmc.burn[,i],0.975)
      Gburn <- cbind(G1,G2)
    }
    
    
    U1 <- matrix(0, 1,1)
    U2 <- matrix(0, 1,1)
    
    U1<-quantile(gl.mcmc.burn,0.025)
    U2<-quantile(gl.mcmc.burn,0.975)
    Uburn <- cbind(U1,U2)
    
    
    # Acceptance Rates 
    arb <- sum(ind1)/nsim
    arg <- sum(ind2)/nsim
    argl <- sum(ind3)/nsim
    
    #Estimaciones adicionales
    mu <- x%*%Bestimado
    sigma2 <- exp(z%*%Gammaest)
    gl1 <- Glestm1
    gl2 <- Glestm2
    varianza <- (gl1/(gl1-2))*sigma2
    loglik <- vero(y,mu,sigma2,as.numeric(Glestm1))
    AIC <- 2*(ncol(x)+ncol(z)+1)-2*loglik
    BIC <- log(length(y))*(ncol(x)+ncol(z))-2*loglik
    residuos <- y-mu
    residuosstd <- residuos/varianza
    
    deviance <- -2*(devero(y,mu,sigma2,as.numeric(gl1))-devero(y,y,sigma2,as.numeric(gl1)))
    deviance <- sign(y-mu)*sqrt(deviance)
    
    if(graph3==TRUE){
      dev.new()
      plot(residuosstd, main="Standardized Residuals")
      dev.new()
      plot(residuosstd, mu, main="Standardized Residuals vs lineal predictor", ylab="XB")
      dev.new()
      plot(deviance, main="Pseudo deviance Residuals")
      dev.new()
      plot(deviance, mu, main="Pseudo deviance Residuals vs lineal predictor", ylab="XB")
      
    }else{}
    
    Deviance <- -2*(vero(y,mu,sigma2,as.numeric(gl1))-vero(y,y,sigma2,as.numeric(gl1)))
    Deviance
    
    expD <- matrix(0,nrow(beta.mcmc),1)
    for(i in 1:nrow(expD)){
      expD[i,] <- -2*vero(y,x%*%beta.mcmc[i,],exp(z%*%gamma.mcmc[i,]),as.numeric(gl.mcmc[i,]))
    }
    
    if(any(is.infinite(expD)==TRUE)==TRUE){
      Dhat <- mean(expD[-which(is.infinite(expD)==TRUE),], na.rm = TRUE)
    }else{
      Dhat <- mean(expD, na.rm = TRUE)
    }
    
    pd <- Dhat - (-2*loglik)
    DIC <- pd + Dhat
    DIC2 <- (-2*loglik) + 2*pd
  
    list(y=y, x=x, z=z, Bestimado=Bestimado, Gammaest=Gammaest, 
         Glestm1=Glestm1,Glestm2=Glestm2,
         B=B, G=G, U=U, Bburn=Bburn, Gburn=Gburn,
         Uburn=Uburn, arb=arb, arg=arg, argl=argl,
         yestimado=mu, variance=varianza, residuals=residuos,
         residualsstd=residuosstd,residualsPseudodev=deviance,
         beta.mcmc=beta.mcmc, gamma.mcmc=gamma.mcmc, gl.mcmc=gl.mcmc,
         beta.mcmc.burn=beta.mcmc.burn, gamma.mcmc.burn=gamma.mcmc.burn, 
         gl.mcmc.burn=gl.mcmc.burn,loglik=loglik, AIC=AIC, 
         BIC=BIC, DIC=DIC, PseudoDeviance=Deviance)
  }
