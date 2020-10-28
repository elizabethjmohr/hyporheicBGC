#'@export
dCdt <- function(t, y, parms){

  with(as.list(c(y,parms)),{

    # Create an interpolation function: for any given tau,
    # the function estimates channel concentration at t minus tau
    laggedCCH <- numeric(100)
    tausForInterp <- seq(0, tauMax, length.out = 100)
    for(i in 1:100){
      if(t-tausForInterp[i]<=0){
        laggedCCH[i] <- CHZ0
      }else{
        laggedCCH[i] <- lagvalue(t-tausForInterp[i],1)
      }
    }
    CchInterpFunc <- approxfun(
      x = tausForInterp,
      y = laggedCCH
    )
    # Calculate derivative of average channel concentration with respect to time
    dCCH <- (q/(sCHsHZ*sHZ)) * (
      integrate(
        RTDintegrand,
        lower = tauMin,
        upper = tauMax,
        tauMin = tauMin,
        tauMax = tauMax,
        alpha = alpha,
        t = t,
        CchInterpFunc = CchInterpFunc,
        k = k
      )$value - CCH)

    # Return derivatives in a list as required by dede solver
    list(dCCH)
  })
}

