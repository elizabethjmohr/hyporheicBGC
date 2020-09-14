rm(list = ls())

library(deSolve)
library(hyporheicBGC)
library(ggplot2)

# Specify model parameters
# Note: sCHsHZ is the ratio of channel volume to HZ volume
parameters <- c(
  tauMin = 0.01,
  tauMax = 2,
  alpha = 1.6,
  k = 0,
  sHZ = 1,
  sCHsHZ = 1,
  q = q(tauMin = 0.01, tauMax = 2, alpha = 1.6)
)

# Specify the initial conditions
state <- c(
  CCH = 3,
  CHZ = 0
)

# Generate a vector of times at which model output is desired
times <- seq(0, parameters["tauMax"]*2, length.out = 200)

# Solve the model with numerical integration
out <- dede(y = state, times = times, func = dCdt, parms = parameters)

# Plot the results
ggplot(data.frame(out),aes(x = time, y = CHZ)) +
  geom_line(size = 2, aes(color = "black")) +
  geom_line(size = 2, aes(x = time, y = CCH,  color = "blue"))+
  ylab("Concentration")+
  xlab("Time")+
  theme_bw()+
  scale_color_identity(guide = guide_legend(override.aes = list(color = c("black", "blue"))),
                       name = NULL,
                       breaks = c("black", "blue"),
                       labels = c("Hyporheic Zone", "Channel"))

