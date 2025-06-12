# Load o2 data and combine with metadata

source("./code/06_clean_metabolism_data.R")
#### 06_clean_metabolism_data.R loads: #####
# `exp_data` = the o2 and temperature data and metadata
# `grp1_dates` = the dates in the first sampling
# `grp2_dates`= the dates in the second sampling
# functions for gpp and er estimation
##
#### Data plotting -----
exp_data %>%
  ggplot(aes(x = date_time, y = o2_do_mg_l))+
  geom_point(aes(color = temp_treat), size = 1.1)+
  geom_line(aes(linetype = nutrient_treat))+
  scale_x_datetime(date_breaks = "1 day", date_labels = "%m/%y")+
  scale_color_manual(values = c("blue","red"))+
  scale_linetype_manual(values = c("solid","dotted"))+
  facet_wrap(~tank)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

### Estimate GPP from full interpolated o2 signal ----
# fit a model to each of the diel curves to interpolate o2
tank1_1 = exp_data %>% named_group_split(tank) %>% .[[1]] %>% filter(run == 'run1') %>% data.frame %>% ungroup %>% droplevels()

tank1_2 = exp_data %>% named_group_split(tank) %>% .[[1]] %>% filter(run == 'run2') %>% data.frame %>% ungroup

tank1_1 %>%
  ggplot()+
  geom_point(aes(x = date_time, y = o2_do_mg_l))+
  geom_smooth(aes(x =date_time, y = o2_do_mg_l), method = 'gam',  se = TRUE)

tank1_2 %>% 
  ggplot()+
  geom_point(aes(x = date_time, y = o2_do_mg_l))+
  geom_smooth(aes(x = date_time, y = o2_do_mg_l), method = 'loess', se = TRUE)

### Run the gams for temp and o2
# functions are found in metabolism-functions.R script
# tankList %>% purrr::walk(~.x %>% purrr::walk(~fit_o2_gam(.x)))
# tankList %>% purrr::walk(~.x %>% purrr::walk(~fit_temp_gam(.x)))


# estimate metabolism from bayesian model with LakeMetabolizer
tankNameList = exp_data %>% dplyr::select(tank, run) %>% 
  dplyr::mutate(names = paste0("t",tank,"_",run)) %>% 
  ungroup %>% dplyr::select(names)%>% 
  distinct %>% 
  unlist %>% 
  sapply(.,function(a) gsub("run","",a)) %>% unname

# Estimate metabolism parameters from modal o2 and temp estimates ----
# purrr::walk(tankNameList, ~estimateContinuous(tankID = .x))

# Get list of files from all metabolism models ----
metabFiles = list.files("./data/models",".*bayes.*.rds", full.names = TRUE)

## Extract the summaries and relevant plots ----
expMetab = metabFiles %>% purrr::map(~extract_metab(metabModel = .x)) %>%
  bind_rows %>%
  dplyr::mutate(tankMod = as.numeric(gsub("t(\\d{1,2})_\\d{1}","\\1", tank)),
                run = gsub("t\\d{1,2}_(\\d{1})","\\1", tank)) %>%
  dplyr::select(-tank) %>%
  dplyr::rename(tank = 'tankMod') %>%
  left_join(exp_metadata, by ='tank')

expMetab %>%
  ggplot()+
  geom_boxplot(aes(x = temp_treat, y = GPP, fill = nutrient_treat))+
  geom_hline(aes(yintercept = 0))+
  theme()+
  facet_wrap(~run)

expMetab %>%
  ggplot()+
  geom_boxplot(aes(x = temp_treat, y = R, fill = nutrient_treat))+
  geom_hline(aes(yintercept = 0))+
  facet_wrap(~run)

expMetab %>% 
  summarise(GPP = mean(GPP), .by = c('temp_treat','nutrient_treat','run'))
# convert all the gpp coefficient estimates for each tank by multiplying GPP coefficient by par and summing across all time ----
# debugonce(convert_metab)
expMetabFull = metabFiles %>% purrr::map(~convert_metab(metabModel = .x)) %>% setNames(.,nm = tankNameList)

expMetabFull %>% 
  bind_rows() %>% 
  mutate(run = gsub('t\\d{1,2}_(\\d{1})','\\1', summary$tank)) %>% 
  mutate(tank = gsub("t(\\d{1,2})_\\d{1}", "\\1", summary$tank)) %>% 
  merge(exp_data %>%
          ungroup %>%
          dplyr::select(tank, temp_treat, nutrient_treat) %>%
          dplyr::mutate(tank = as.character(tank)) %>% 
          distinct, by = 'tank') %>% 
  summarise(GPP = mean(gpp.out),
            GPP_l = quantile(gpp.out, 0.025),
            GPP_u = quantile(gpp.out, 0.975),
            R = mean(r.out),
            R_l = quantile(r.out, 0.025),
            R_u = quantile(r.out, 0.975),
            NEP = mean(nep.out),
            NEP_l = quantile(nep.out, 0.025),
            NEP_u = quantile(nep.out, 0.975),
            .by = c('temp_treat','nutrient_treat','run'))

expMetab %>%
  ggplot()+
  geom_boxplot(aes(x = temp_treat, y = NEP, fill = nutrient_treat))+
  geom_hline(aes(yintercept = 0))+
  facet_wrap(~run)

saveRDS(expMetabFull, "./data/gppPosteriors.rds")

expMetabSumm = expMetabFull %>% 
  map(~.x %>% pluck('summary')) %>% 
  bind_rows() %>% 
  mutate(run = gsub('t\\d{1,2}_(\\d{1})','\\1', tank)) %>% 
  mutate(tank = gsub("t(\\d{1,2})_\\d{1}", "\\1", tank)) %>% 
  merge(exp_data %>%
          ungroup %>%
          dplyr::select(tank, temp_treat, nutrient_treat) %>%
          dplyr::mutate(tank = as.character(tank)) %>% 
          distinct, by = 'tank')

expMetabSumm %>% 
  mutate(tank = as.numeric(tank)) %>% 
  ggplot()+
  # geom_errorbar(aes(x = run, ymin = (GPP-GPPsd), ymax = (GPP+GPPsd)), width = 0.1)+
  # geom_errorbar(aes(x = run, ymin = (R+Rsd), ymax = (R-Rsd)), width = 0.1)+
  geom_col(aes(x = run, y = GPP), fill = 'darkgreen', width = 0.5, alpha = 0.5)+
  geom_col(aes(x = run, y = R), fill = 'red', width = 0.5, alpha = 0.5)+
  geom_hline(aes(yintercept = 0))+
  facet_wrap(~tank)
  
#### Estimate metabolism from Dawn-Dusk-Dawn -----
# estimate patterns of GP, NP, and R
exp_dataList = exp_data %>%
  ungroup %>%
  named_group_split(tank)

# debug(estimateDuskDawn)
dusk_dawnMetEstimates = exp_dataList %>%
  purrr::map(~.x %>% named_group_split(run) %>%
               purrr::map(~estimateDuskDawn(.x)) %>%
               bind_rows(.id = 'run')) %>%
  bind_rows(.id = 'tank') %>%
  merge(exp_data %>%
          ungroup %>%
          dplyr::select(tank, temp_treat, nutrient_treat) %>%
          dplyr::mutate(tank = as.character(tank)) %>% 
          distinct, by = 'tank')

# create boxplot of NP
dusk_dawnMetEstimates %>%
  ggplot()+
  geom_boxplot(aes(x = temp_treat, y = NP_mg_o2_l_hr, fill = nutrient_treat)) +
  scale_x_discrete(name = "Temperature treatment")+
  scale_y_continuous(name = expression("Net production ( mg"~O[2]~L^-1~hr^-1~")"),
                     expand = c(0.01,0.01))+
  theme(legend.position = 'inside',
        legend.position.inside = c(1,1),
        legend.justification = c(1,1))+
  facet_wrap(~run)

# create boxplot of R
dusk_dawnMetEstimates %>%
  ggplot()+
  geom_boxplot(aes(x = temp_treat, y = Rnight_mg_o2_l_hr, fill = nutrient_treat))+
  scale_x_discrete(name = "Temperature treatment")+
  scale_y_continuous(name = expression("Respiration ( mg"~O[2]~L^-1~hr^-1~")"),
                     limits = c(NA, 0), expand = c(0.01,0.01))+
  theme(legend.position = "inside",
        legend.position.inside = c(1,1),
        legend.justification = c(1,1))+
  facet_wrap(~run)

# create boxplot of gpp
dusk_dawnMetEstimates %>%
  ggplot()+
  geom_boxplot(aes(x = temp_treat, y = GP_mg_o2_l_hr, fill = nutrient_treat))+
  scale_x_discrete(name = "Temperature treatment")+
  scale_y_continuous(name = expression("Gross production ( mg"~O[2]~L^-1~hr^-1~")"),
                     limits = c(0, NA), expand = c(0.01,0.01))+
  theme(legend.position = 'inside',
        legend.position.inside = c(1,1),
        legend.justification = c(1,1))+
  facet_wrap(~run)

dusk_dawnMetEstimates %>%
  ggplot()+
  geom_boxplot(aes(x = temp_treat, y = (GP_mg_o2_l_hr/abs(Rnight_mg_o2_l_hr)), fill = nutrient_treat))+
  scale_x_discrete(name = "Temperature treatment")+
  scale_y_continuous(name = expression("GPP:ER"),
                     limits = c(0.5,NA),expand = c(0.01,0.01))+
  geom_hline(aes(yintercept = 1))+
  theme(legend.position = 'inside',
        legend.position.inside = c(1,1),
        legend.justification = c(1,1))+
  facet_wrap(~run)

# scatter plot of GP & R
dusk_dawnMetEstimates %>%
  ggplot()+
  geom_point(aes(x = GP_mg_o2_l_hr, y = abs(Rnight_mg_o2_l_hr), color = temp_treat, fill = nutrient_treat), shape = 21, size =3, stroke = 1.3)+
  geom_abline()+
  geom_smooth(aes(x = GP_mg_o2_l_hr, y = abs(Rnight_mg_o2_l_hr)), method = 'lm', se = FALSE)+
  scale_y_continuous(name = expression("Respiration ( -mg"~O[2]~L^-1~hr^-1~")"),
                     limits = c(0,NA), expand = c(0.01,0.01))+
  scale_x_continuous(name = expression("Gross production ( mg"~O[2]~L^-1~hr^-1~")"),
                     limits = c(0,NA), expand = c(0.01,0.01))+
  scale_color_manual(values = c("blue","red"))+
  scale_fill_manual(values = c("blue","red"))+
  facet_wrap(~run)

dusk_dawnMetEstimates %>% 
  summarise(GPP = mean(GP_mg_o2_l_hr*24), .by = c('temp_treat','nutrient_treat','run'))

# compare the two methods to estimating
expMetabSumm %>% 
  select(tank,run, GPP, R) %>%
  mutate(method = 'continuous') %>% 
  group_by(tank, run, method) %>% 
  tidyr::pivot_longer(cols = GPP:R, names_to = 'variable', values_to = 'value') %>% 
  bind_rows(dusk_dawnMetEstimates %>%
              mutate(method = 'duskdawn',
                     run = gsub('run(\\d{1})','\\1', run),
                     GPP = GP_mg_o2_l_hr*24,
                     R = Rnight_mg_o2_l_hr*24) %>%
              select(tank, run, method, GPP, R) %>% 
              group_by(tank, run, method) %>% 
              tidyr::pivot_longer(cols = GPP:R, names_to = 'variable', values_to = 'value')) %>% 
  ungroup %>%
  tidyr::unite(col = 'measure', method, variable ) %>%
  tidyr::pivot_wider(id_cols = c('tank', 'run'), names_from = measure, values_from = value) %>% 
  ggplot()+
  geom_point(aes(x = continuous_GPP, y = duskdawn_GPP))+
  geom_abline()+
  facet_wrap(~run)

expMetabSumm %>% 
  mutate(tank = as.numeric(tank)) %>% 
  ggplot()+
  geom_point(data = dusk_dawnMetEstimates %>% mutate(run = gsub("run(\\d{1})","\\1", run),
                                                     tank = as.numeric(tank)), aes(x = run, y = (GP_mg_o2_l_hr*24)), color = 'darkgreen')+
  geom_point(data = dusk_dawnMetEstimates %>% mutate(run = gsub("run(\\d{1})","\\1", run),
                                                     tank = as.numeric(tank)), aes(x = run, y = (Rnight_mg_o2_l_hr*24)), color = 'red')+ 
  geom_col(aes(x = run, y = GPP), fill = 'darkgreen', width = 0.5, alpha = 0.5)+
  geom_col(aes(x = run, y = R), fill = 'red', width = 0.5, alpha = 0.5)+
  geom_hline(aes(yintercept = 0))+
  facet_wrap(~tank)

