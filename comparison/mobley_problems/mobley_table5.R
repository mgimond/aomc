library(tidyverse)
library(scales)
library(kableExtra)

source("make_long_aop.R")

mob <- readxl::read_excel("Mobley_comparisons.xlsx", sheet = "long_form")


pad_percent <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%06.", digits, "f%%"), x * 100))
}


# Problem 1
aop <- make_long_aop("../mobley1/")
w0.2 <- aop %>% 
  filter(aop %in% c("Ed", "Eou", "Lu")) %>% 
  filter(wavelength %in% c(1)) %>% 
  mutate(depth = round(depth,1)) %>% 
  filter(depth%in% c(0.8, 4, 8)) 
w0.9 <- aop %>% 
  filter(aop %in% c("Ed", "Eou", "Lu")) %>% 
  filter(wavelength %in% c(2)) %>% 
  mutate(depth = round(depth,1)) %>% 
  filter(depth %in% c(0.1, 0.5, 1)) 

prob1 <- rbind(w0.2, w0.9) %>% 
  mutate(omega = ifelse(wavelength == 1, 0.2, 0.9)) %>% 
  select(-wavelength)


table <- mob %>% 
  filter(problem == 1) %>%
  mutate(omega = ifelse(wavelength == 1, 0.2, 0.9)) %>% 
  select(-problem, -wavelength) %>% 
  left_join(prob1, by = c("omega", "depth", "aop"),
            suffix = c("_mob", "_aomc")) %>% 
  mutate(rd = (value_mob - value_aomc) / value_mob ,
         rd = pad_percent(rd, digits =2),
         CV = pad_percent(CV, digits =2)) %>% 
  select(aop, depth, omega, value_mob, value_aomc, CV, rd)

knitr::kable(table)
knitr::kable(table %>% filter(omega==0.2) %>% select(-omega)) 
knitr::kable(table %>% filter(omega==0.9) %>% select(-omega))   

# Problem 2
aop <- make_long_aop("../mobley2/")
w0.2 <- aop %>% 
  filter(aop %in% c("Ed", "Eou", "Lu")) %>% 
  filter(wavelength %in% c(1)) %>% 
  mutate(depth = round(depth,1)) %>% 
  filter(depth%in% c(0.8, 4, 8)) 
w0.9 <- aop %>% 
  filter(aop %in% c("Ed", "Eou", "Lu")) %>% 
  filter(wavelength == 2) %>% 
  mutate(depth = round(depth,1)) %>% 
  filter(depth %in% c(0.1, 0.5, 1)) 

prob1 <- rbind(w0.2, w0.9) %>% 
  mutate(omega = ifelse(wavelength == 1, 0.2, 0.9)) %>% 
  select(-wavelength)

table <- mob %>% 
  filter(problem == 2) %>%
  mutate(omega = ifelse(wavelength == 1, 0.2, 0.9)) %>% 
  select(-problem, -wavelength) %>% 
  left_join(prob1, by = c("omega", "depth", "aop"),
            suffix = c("_mob", "_aomc")) %>% 
  mutate(rd = (value_mob - value_aomc) / value_mob ,
         rd = pad_percent(rd, digits =2),
         CV = pad_percent(CV, digits =2)) %>% 
  select(aop, depth, omega, value_mob, value_aomc, CV, rd)

knitr::kable(table)
knitr::kable(table %>% filter(omega==0.2) %>% select(-omega)) 
knitr::kable(table %>% filter(omega==0.9) %>% select(-omega))   


# Problem 6
aop <- make_long_aop("../mobley6/.")
w0.2 <- aop %>% 
  filter(aop %in% c("Ed", "Eou", "Lu")) %>% 
  filter(wavelength %in% c(1)) %>% 
  mutate(depth = round(depth,1)) %>% 
  filter(depth%in% c(0.8, 4)) 

prob1 <- w0.2 %>% 
  mutate(omega = ifelse(wavelength == 1, 0.2, 0.9)) %>% 
  select(-wavelength)

table <- mob %>% 
  filter(problem == 6) %>%
  mutate(omega = ifelse(wavelength == 1, 0.2, 0.9)) %>% 
  select(-problem, -wavelength) %>% 
  left_join(prob1, by = c("omega", "depth", "aop"),
            suffix = c("_mob", "_aomc")) %>% 
  mutate(rd = (value_mob - value_aomc) / value_mob ,
         rd = pad_percent(rd, digits =2),
         CV = pad_percent(CV, digits =2)) %>% 
  select(aop, depth, omega, value_mob, value_aomc, CV, rd)

knitr::kable(table)
knitr::kable(table %>% filter(omega==0.2) %>% select(-omega)) 

