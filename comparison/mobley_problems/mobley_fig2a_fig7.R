library(tidyverse)
source("make_long_aop.R")

mob <- readxl::read_excel("Fig2_fig7.xlsx")

# Figure 2a
aomc2a <- make_long_aop("../mobley1/.") %>% 
  filter(aop %in% c("Ed", "Eou", "Lu")) %>% 
  filter(wavelength == 2) %>% 
  filter(depth <= 4,
         value > 0) %>%  #for log plot 
  arrange(aop, depth) %>% 
  mutate(depth = ifelse(depth == -1, -0.2, depth),
         dataset = "AOMC")

mob2a <- mob %>% 
  filter(figure == "fig2a")  %>% 
  arrange(aop, depth)%>% 
  mutate(depth = ifelse(depth == -1, -0.2, depth),
         dataset = "Mobley")

ggplot(aomc2a, aes(value, depth, col = aop, lty = dataset)) + 
  geom_path(lwd = 2, alpha = 0.1) + # Preserves order of variables, otherwise, it orders by x
  geom_path(data =mob2a, aes(value,depth, col = aop, lty=dataset)) +
   scale_y_reverse(limits = c(1, -0.2))+
  scale_x_log10(
    breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100),
    limits = c(0.0001, 1),
    labels = scales::label_number(drop0trailing = TRUE),
    minor_breaks = FALSE) +
  guides(x = guide_axis_logticks())


# Figure 7
aomc2a <- make_long_aop("../mobley6/.") %>% 
  filter(aop %in% c("Ed", "Eou", "Lu")) %>% 
  filter(wavelength == 1) %>% 
  filter(depth <= 4,
         value > 0) %>%  #for log plot 
  arrange(aop, depth) %>% 
  mutate(depth = ifelse(depth == -1, -0.4, depth),
         dataset = "AOMC")


mob2a <- mob %>% 
  filter(figure == "fig7")  %>% 
  arrange(aop, depth) %>% 
  mutate(depth = ifelse(depth == -1, -0.4, depth),
         dataset = "Mobley")

ggplot(aomc2a, aes(value, depth, col = aop, lty = dataset)) + 
  geom_path(lwd = 2, alpha = 0.1) + # Preserves order of variables, otherwise, it orders by x
  geom_path(data =mob2a, aes(value,depth, col = aop, lty=dataset)) +
  geom_point(data =mob2a, aes(value,depth, col = aop)) +
  scale_y_reverse()+
  scale_x_log10(
    #   breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100),
    #   limits = c(0.0001, 100),
    labels = scales::label_number(drop0trailing = TRUE),
    minor_breaks = FALSE) +
  guides(x = guide_axis_logticks())
