# Script to estimate GPP from O2 measurements
library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(streamMetabolizer)
library(brms)
library(mgcv)
library(rstan)
library(purrr)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)
source("./code/metabolism_functions.R")
'%ni%' <- Negate('%in%')
# library(rnoaa)
# library(rwunderground)
theme_set(theme_minimal())
# import oxygen data 
o2_data = read.csv(file = "./data/o2_temp.csv", TRUE) %>%
  dplyr::mutate(date_time = as.POSIXct(date_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "America/Chicago"))

temp_data = read.csv(file = "./data/temperature.csv", TRUE) %>%
  dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y"),
    date_time = as.POSIXct(paste0(date," ",paste0(time,":00")), format = "%Y-%m-%d %H:%M:%S", tz = "America/Chicago")) %>%
  dplyr::select(tank, date_time, temp_deg_f) %>%
  dplyr::mutate(temp_c = F_to_C(fahr = temp_deg_f)) 

# split into list and merge temps to fill in holes
o2_data_tanks = o2_data %>%
  named_group_split(tank)
temp_data_tanks = temp_data %>%
  named_group_split(tank)

o2_data = purrr::map2(o2_data_tanks, temp_data_tanks, ~.x %>% dplyr::mutate(tempFill = unname(merge_temps(o2Df = .x, tempDf = .y)))) %>% bind_rows(.id = 'tank') %>%
  dplyr::mutate(temp_deg_c = tempFill) %>% dplyr::select(-tempFill) %>%
  dplyr::mutate(tank = as.integer(tank),
                # calculate the DO percent saturation
                o2_do_sat = calc_DO_sat(temp.water = temp_deg_c,
                                                           pressure.air = (28.96/0.029530)),
                o2_do_pct_sat = (o2_do_mg_l/o2_do_sat)*100)

# separate unique metabolism "runs"
grp1_dates = as.Date(c("2022-06-21","2022-06-22"), tz = attr(o2_data$date_time, 'tzone'))
grp2_dates = as.Date(c("2022-06-27","2022-06-28"), tz = attr(o2_data$date_time, 'tzone'))

o2_data = o2_data %>%
  dplyr::mutate(date = as.Date(date_time, tz = attr(date_time, "tzone")),
                run = case_when(date %in% grp1_dates ~ "run1",
                                date %in% grp2_dates ~ "run2",
                                TRUE ~ NA_character_),
                TOD = case_when(lubridate::hour(date_time) < 12 ~ "morning",
                                TRUE ~ "evening")) %>%
  group_by(tank, run) %>%
  dplyr::mutate(day = as.integer(date - min(date)+1)) %>%
  group_by(tank, run, day, TOD) %>%
  arrange(date_time) %>%
  dplyr::mutate(time_point = 1:n()) %>%
  ungroup %>% group_by(tank, run) %>%
  dplyr::mutate(run_hr = as.numeric(date_time - min(date_time))/3600)
# import experimental design metadata 

meta_data = read.csv(file = "./data/treatments.csv", TRUE)

# 
exp_data = o2_data %>%
  left_join(meta_data %>%
              dplyr::select(tank, temp_treat, nutrient_treat), by = "tank")

exp_metadata = exp_data %>%
  group_by(tank) %>%
  dplyr::summarise(across(c(temp_treat, nutrient_treat), ~unique(.x)))

rm(list = ls()[ls() %ni% c("exp_data", "exp_metadata", "site_coords", ls.str(mode = 'function'))])
