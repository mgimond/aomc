library(tidyverse)

dat <- read_table("wave.out")

# B = 0.019
# Note that the wavelength number here is simply the
# scattering coefficient number

dat$bB <- dat$wavelength * 0.019

M <- lm(R ~ bB, dat)
coef(M)[2]

ggplot(dat, aes(bB, R)) + geom_point() +
  stat_smooth(method = "lm", se=FALSE)
