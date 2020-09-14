#'@export
E <- function(tau, tauMin, tauMax, alpha){
  (tau^(-alpha) - tauMax^(-alpha))*(-alpha+1)/
    (tauMax^(-alpha+1)-tauMin^(-alpha+1)-tauMax^(-alpha)*(tauMax - tauMin)*(1-alpha))
}

#'@export
W <- function(tau, tauMin, tauMax, alpha){
  (tauMax^(-alpha+1) - tau^(-alpha+1)- (tauMax^(-alpha)*(tauMax-tau)*(-alpha+1)))/
    (tauMax^(-alpha+1) - tauMin^(-alpha+1)- (tauMax^(-alpha)*(tauMax-tauMin)*(-alpha+1)))
}

#'@export
W_int <- function(tauMin, tauMax, alpha){
  result <- integrate(W, lower = tauMin, upper = tauMax, tauMin = tauMin, tauMax = tauMax, alpha = alpha)
  result$value
}

#'@export
I <- function(tau, tauMin, tauMax, alpha){
  W(tau, tauMin, tauMax, alpha)/W_int(tauMin, tauMax, alpha)
}

#'@export
q <- function(tauMin, tauMax, alpha, sHZ = 1){
  sHZ/W_int(tauMin, tauMax, alpha)
}

#'@export
integrand_firstOrderUptake <- function(tau, tauMin, tauMax, alpha, t, CchInterpFunc, k){
  I(tau, tauMin, tauMax, alpha) * CchInterpFunc(tau) *exp(-k*tau)
}

#'@export
integrand_firstOrderExit <- function(tau, tauMin, tauMax, alpha, t, CchInterpFunc, k){
  E(tau, tauMin, tauMax, alpha) * CchInterpFunc(tau) *exp(-k*tau)
}


