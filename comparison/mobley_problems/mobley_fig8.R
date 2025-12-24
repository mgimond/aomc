library(tidyverse)
source("make_long_rad.R")
# Figure 8 comparison
#aomc <- read_csv("prob2_aomc_rad.csv")  # This is the rad formatted df from make_long_rad
aomc <- make_long_rad("../mobley2/")
mobley <- read_csv("mobley_fig8.csv")

# Tau = 0 (Depth = 0)
# Wavelength # 2 is the case where omega = 0.9
fig8_aomc_tau0 <- aomc %>% filter(depth == 0, 
                                  Wavelength == 2, 
                                  phi %in% c(0,180)) %>% 
  mutate(alpha = ifelse(phi == 180, alpha +180, -(alpha -180)))

fig8_mobley_tau0 <- mobley %>% filter(tau == 0)

ggplot(fig8_aomc_tau0, aes(x=alpha, y = radiance, lty = "AOMC")) + 
  geom_line() +
  geom_line(data = fig8_mobley_tau0, aes(x=alpha, y = radiance, lty="Mobley"),
            col = "grey30") +
  scale_y_log10(
    breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100),
    limits = c(0.0001, 100),
    labels = scales::label_number(drop0trailing = TRUE),
    minor_breaks = FALSE) +
  guides(y = guide_axis_logticks())+
  scale_x_continuous(breaks = seq(0, 360, 45))

# Tau = 5 (Depth = 0.5)
# Wavelength # 2 is the case where omega = 0.9
fig8_aomc_tau5 <- aomc %>% filter(depth == 0.5, 
                                  Wavelength == 2, 
                                  phi %in% c(0,180)) %>% 
  mutate(alpha = ifelse(phi == 180, alpha +180, -(alpha -180)))

fig8_mobley_tau5 <- mobley %>% filter(tau == 5)

ggplot(fig8_aomc_tau5, aes(x=alpha, y = radiance, lty = "AOMC")) + geom_line() +
  geom_line(data = fig8_mobley_tau5, aes(x=alpha, y = radiance, lty="Mobley"),
            col = "grey30") +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1),
    limits = c(0.001, 1),
    labels = scales::label_number(drop0trailing = TRUE),
    minor_breaks = FALSE) +
  guides(y = guide_axis_logticks())+
  scale_x_continuous(breaks = seq(0, 360, 45))


# Tau = 20 (Depth = 2.0)
# Wavelength # 2 is the case where omega = 0.9
fig8_aomc_tau20 <- aomc %>% filter(depth == 2.0, 
                                  Wavelength == 2, 
                                  phi %in% c(0,180)) %>% 
  mutate(alpha = ifelse(phi == 180, alpha +180, -(alpha -180)))

fig8_mobley_tau20 <- mobley %>% filter(tau == 20)

ggplot(fig8_aomc_tau20, aes(x=alpha, y = radiance, lty = "AOMC")) + geom_line() +
  geom_line(data = fig8_mobley_tau20, aes(x=alpha, y = radiance, lty="Mobley"),
            col = "grey30") +
  scale_y_log10(
    breaks = c(0.0001, 0.001, 0.01),
    limits = c(0.0001, 0.01),
    labels = scales::label_number(drop0trailing = TRUE),
    minor_breaks = FALSE) +
  guides(y = guide_axis_logticks())+
  scale_x_continuous(breaks = seq(0, 360, 45))

# Combined plots

ggplot(fig8_aomc_tau0, aes(x=alpha, y = radiance, lty = "AOMC")) + geom_line() +
  geom_line(data = fig8_mobley_tau0, aes(x=alpha, y = radiance, lty="Mobley"), 
            col = "grey30") +
  geom_line(data = fig8_aomc_tau5, aes(x=alpha, y = radiance, lty = "AOMC"),
            col = "grey30") +
  geom_line(data = fig8_mobley_tau5, aes(x=alpha, y = radiance, lty="Mobley"), 
            col = "grey30") +
  geom_line(data = fig8_aomc_tau20, aes(x=alpha, y = radiance, lty = "AOMC"), 
            col = "grey30") +
  geom_line(data = fig8_mobley_tau20, aes(x=alpha, y = radiance, lty="Mobley"), 
            col = "grey30") +
  scale_y_log10(
    breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100),
    limits = c(0.0001, 100),
    labels = scales::label_number(drop0trailing = TRUE),
    minor_breaks = FALSE) +
  guides(y = guide_axis_logticks())+
  scale_x_continuous(breaks = seq(0, 360, 45))

