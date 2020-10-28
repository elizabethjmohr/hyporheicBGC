rm(list = ls())
library(hyporheicBGC)
library(deSolve)
library(tidyverse)
library(ggplot2)

kZeroLowAlpha <- simulateSAM(
  tauMin = 0.01,
  tauMax = 2,
  alpha = 1.3,
  k = 0,
  sHZ = 1,
  sCHsHZ = 1,
  CCH0 = 3,
  CHZ0 = 0
) %>%
  mutate(simulation = "kZeroLowAlpha")

kZeroHighAlpha <- simulateSAM(
  tauMin = 0.01,
  tauMax = 2,
  alpha = 1.9,
  k = 0,
  sHZ = 1,
  sCHsHZ = 1,
  CCH0 = 3,
  CHZ0 = 0
) %>%
  mutate(simulation = "kZeroHighAlpha")

kNonZeroLowAlpha <- simulateSAM(
  tauMin = 0.01,
  tauMax = 2,
  alpha = 1.3,
  k = 0.1,
  sHZ = 1,
  sCHsHZ = 1,
  CCH0 = 3,
  CHZ0 = 0
) %>%
  mutate(simulation = "kNonZeroLowAlpha")

kNonZeroHighAlpha <- simulateSAM(
  tauMin = 0.01,
  tauMax = 2,
  alpha = 1.9,
  k = 0.1,
  sHZ = 1,
  sCHsHZ = 1,
  CCH0 = 3,
  CHZ0 = 0
)%>%
  mutate(simulation = "kNonZeroHighAlpha")

allOutput <- bind_rows(
  kZeroLowAlpha,
  kZeroHighAlpha,
  kNonZeroLowAlpha,
  kNonZeroHighAlpha
) %>%
  mutate(k = case_when(simulation %in% c("kZeroLowAlpha", "kZeroHighAlpha") ~ "0",
                       TRUE ~ "0.1")) %>%
  mutate(alpha = case_when(simulation %in% c("kZeroLowAlpha", "kNonZeroLowAlpha") ~ "1.3",
                       TRUE ~ "1.9"))

# Plot the results
p1<- ggplot(kNonZeroLowAlpha,aes(x = time, y = CHZ)) +
  geom_line(aes(color = "black")) +
  geom_line(aes(x = time, y = CCH,  color = "dodgerblue"))+
  ylab("Concentration")+
  xlab("Time")+
  ggtitle(expression(paste(alpha, "= 1.3, k = 0.1")))+
  theme_bw()+
  scale_color_identity(guide = guide_legend(override.aes = list(color = c("black", "dodgerblue"))),
                       name = NULL,
                       breaks = c("black", "dodgerblue"),
                       labels = c("Hyporheic Zone", "Channel"))

quartz(width = 6, height = 4)
p1

filteredOutput <- allOutput %>%
  filter(time >= 0)
pWrap <- ggplot(filteredOutput,aes(x = time, y = CCH)) +
  geom_line() +
  ylab("Concentration")+
  xlab("Time")+
  theme_bw() +
  facet_wrap(vars(simulation))
quartz(width = 6, height = 4)
pWrap

pLines <- ggplot(filteredOutput,aes(x = time, y = CCH)) +
  geom_line(aes(color = alpha, linetype = k)) +
  ylab("Concentration")+
  xlab("Time")+
  theme_bw() +
  scale_color_manual(name = expression(paste(alpha)),
                     values = c("dodgerblue", "midnightblue"))


quartz(width = 6, height = 4)
pLines

