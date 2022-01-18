print.summary.Bayesiantreg <-
function(x, ...){
  
  cat (" \n ################################################################
            ###                    Bayesian t Regression                  ###
            ################################################################ \n")
  
  cat("\n Call: \n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, digits=4)
  
  cat("\n AIC: \n")
  print(x$Criteria[1])
  
  cat("\n BIC: \n")
  print(x$Criteria[2])
  
  cat("\n DIC: \n")
  print(x$Criteria[3])
  
  cat("\n Pseudo Deviance: \n")
  print(x$Criteria[4])
  
  }
