#'@export
#'@param sCHsHZ Ratio of channel volume to hyporheic zone volume
simulateSAM <- function(tauMin,
                        tauMax,
                        alpha,
                        k,
                        sHZ,
                        sCHsHZ,
                        CCH0,
                        CHZ0,
                        times = seq(0, tauMax, length.out = 1000)){

  # Specify model parameters
  parameters <- c(
    tauMin = tauMin,
    tauMax = tauMax,
    alpha = alpha,
    k = k,
    sHZ = sHZ,
    sCHsHZ = sCHsHZ,
    q = q(tauMin = tauMin, tauMax = tauMax, alpha = alpha, sHZ = sHZ),
    CHZ0 = CHZ0
  )

  # Specify initial conditions in the channel
  state <- c(
    CCH = CCH0
  )

  # Solve the channel model using a solver for delay differential equations
  out <- dede(y = state,
              times = times,
              func = dCdt,
              parms = parameters)

  # Define a function that uses interpolation to estimate channel concentration
  # at any time, based on model output
  CchInterpFunc <- approxfun(
    x = out[,1],
    y = out[,2]
  )

  # Expand domain of interpolation function to include times before
  # the beginning of the model run
  interpTimes <- seq(-tauMax, max(out[,1]), length.out = 1000)

  laggedCCH <- numeric(1000)
  for(i in 1:1000){
    if(interpTimes[i]< 0){
      laggedCCH[i] <- CHZ0
    }else{
      laggedCCH[i] <- CchInterpFunc(interpTimes[i])
    }
  }

  CchInterpFunc <- approxfun(
    x = interpTimes,
    y = laggedCCH
  )

  output <- data.frame(out) %>%
    # Calculate hyporheic zone concentrations at each time
    mutate(CHZ = sapply(times,
                        ChzFunc,
                        tauMin = parameters["tauMin"],
                        tauMax = parameters["tauMax"],
                        alpha = parameters["alpha"],
                        k = parameters["k"],
                        CchInterpFunc = CchInterpFunc)) %>%
    # Bind rows for ancillary data
    bind_rows(
      data.frame(
        time = seq(-tauMax/4, -0.001, length.out = 10),
        CCH = rep(CHZ0, times = 10),
        CHZ = rep(CHZ0, times = 10)
      )
    )
  return(output)
}
