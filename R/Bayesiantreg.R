Bayesiantreg <-
  function(y, x, z, nsim, bini, bpri, Bpri, gini, gpri, Gpri, glini, glpri, 
           type, apriori, propuesta, Maxi=NULL,
           lambda = NULL, p = NULL, burn, jump, graph1 = TRUE, graph2 = TRUE,
           graph3 = TRUE){
    
    est1 <- BayesiantregEst(y, x, z, nsim, bini, bpri, Bpri, gini, gpri,Gpri, glini, glpri, 
                           type, apriori, propuesta, Maxi=Maxi,
                           lambda = lambda, p = p, burn, jump, graph1, 
                           graph2, graph3)
    est <- list()
    
    est$coefficients <- matrix(c(est1$Bestimado,est1$Gammaest,est1$Glestm1,est1$Glestm2))
    names <- c(colnames(est1$x),colnames(est1$z),"mean", "median")
    par <-rep(c("beta.","gamma.","gl."),c(ncol(est1$x),ncol(est1$z),2))
    rownames(est$coefficients) <-paste(par,names,sep="")
    #est$desv <-  matrix(c(est$DesvBeta, est$DesvGamma))
    est$interv<- rbind(est1$B,est1$G,est1$U,c("",""))
    est$fitted.values <- est1$yestimado
    est$residuals <- est1$residuos
    est$residualsstd <- est1$residuosstd
    est$residualsdev <- est1$residuosdev
    est$variance <- est1$variance
    est$beta.mcmc<-est1$beta.mcmc
    est$gamma.mcmc<-est1$gamma.mcmc
    est$gl.mcmc<-est1$gl.mcmc
    est$beta.mcmc.burn<-est1$beta.mcmc.burn
    est$gamma.mcmc.burn<-est1$gamma.mcmc.burn
    est$gl.mcmc.burn<-est1$gl.mcmc.burn
    est$loglik <- est1$loglik
    est$AIC <- est1$AIC
    est$BIC <- est1$BIC
    est$DIC <- est1$DIC
    est$PseudoDeviance <- est1$PseudoDeviance
    arb <- est1$arb 
    arg <- est1$arg 
    argl <- est1$argl
    est$y<-est1$y
    est$x<-est1$x
    est$z<-est1$z
    
    
    est$call <- match.call()
    
    class(est) <- "Bayesiantreg"
    
    est 
    
  }

