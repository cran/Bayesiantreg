criteria <-
  function(object, ...){

    loglik <- object$loglik
    AIC <- object$AIC
    BIC <- object$BIC
    DIC <- object$DIC
    PseudoDeviance <-object$PseudoDeviance 
    
    criterias <- as.matrix(c(loglik, AIC, BIC, DIC, PseudoDeviance))
    rownames(criterias) <- c("loglik", "AIC", "BIC", "DIC", "PseudoDeviance")
    return(criterias)
  }