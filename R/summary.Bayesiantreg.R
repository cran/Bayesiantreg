summary.Bayesiantreg <-
function(object, ...){
  
   TAB <- cbind(Coefficient = object$coefficients,
                L.CredIntv = as.numeric(object$interv[,1]),
                U.CredIntv = as.numeric(object$interv[,2])
  )
  
  colnames(TAB) <- c("Estimate", "L.CredIntv",  "U.CredIntv")
  
  
  AIC <- object$AIC
  BIC <- object$BIC
  DIC <- object$DIC
  PseudoDeviance <-object$PseudoDeviance 
  
  Info <- as.matrix(c(AIC, BIC, DIC, PseudoDeviance))
  rownames(Info) <- c("AIC", "BIC", "DIC", "PseudoDeviance")
  
  res <- list(call=object$call, coefficients=TAB, Criteria=Info)  
  
  class(res) <- "summary.Bayesiantreg"
  res  
}
