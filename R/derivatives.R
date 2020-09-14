#'@export
dCdt <- function(t, state, parameters){

  with(as.list(c(state,parameters)), {

    # Create an interpolation function: for any given tau,
    # the function estimates channel concentration at t minus tau
    Cch <- numeric(100)
    CchTauVals <- seq(tauMin, tauMax, length.out = 100)
    for(i in 1:100){
      if(CchTauVals[i]>=t){
        Cch[i] <- 0
      }else{
        Cch[i] <- lagvalue(t-CchTauVals[i])[1]
      }
    }
    CchInterpFunc <- approxfun(
      x = CchTauVals,
      y = Cch
    )
    # Calculate derivative of average channel concentration with respect to time
    dCCH <- (q/(sCHsHZ*sHZ)) * (
      integrate(
        integrand_firstOrderExit,
        lower = tauMin,
        upper = tauMax,
        tauMin = tauMin,
        tauMax = tauMax,
        alpha = alpha,
        t = t,
        CchInterpFunc = CchInterpFunc,
        k = k
      )$value - CCH)

    # Calculate derivative of average hyporheic zone concentration with respect to time
    dCHZ <- (q/sHZ)*(CCH -
      integrate(
        integrand_firstOrderExit,
        lower = tauMin,
        upper = tauMax,
        tauMin = tauMin,
        tauMax = tauMax,
        alpha = alpha,
        t = t,
        CchInterpFunc = CchInterpFunc,
        k = k
      )$value) -
      k * (integrate(
        integrand_firstOrderUptake,
        lower = tauMin,
        upper = tauMax,
        tauMin = tauMin,
        tauMax = tauMax,
        alpha = alpha,
        t = t,
        CchInterpFunc = CchInterpFunc,
        k = k
      )$value)

    # Return derivatives in a list as required by dede
    list(c(dCCH, dCHZ))
  })
}
