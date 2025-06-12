#'
#'
#'
F_to_C = function(fahr){
  (fahr - 32) * 5/9
}


#'
#'
named_group_split = function (.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), 
                                              sep = " / ")))
  grouped %>% group_split() %>% rlang::set_names(names)
}

#'
#'
#'
merge_temps = function(o2Df = NULL, tempDf = NULL, tempCol = 'temp_c', maxDiff = 2){
  
  o2_timeCol = do.call("c",o2Df[,sapply(o2Df, function(x) inherits(x, 'POSIXt'))])
  temp_timeCol = do.call("c",tempDf[,sapply(tempDf, function(x) inherits(x, 'POSIXt'))])
  
  hourDiff = sapply(o2_timeCol, function(x) temp_timeCol[which.min(abs(temp_timeCol-x))]-x)
  hourDiff_mod = sapply(hourDiff, function(x) ifelse(x > maxDiff,NA,x))
  timeDiff = sapply(o2_timeCol, function(x) which.min(abs(temp_timeCol-x)))
  tempFill = unlist(tempDf[timeDiff,tempCol])
  
  return(tempFill)
}


#'
#'
#'
estimateDuskDawn = function(df = NULL,...){
  require(plyr)
  require(dplyr)
  require(magrittr)
  date_timeCol = names(df[sapply(df, function(x) inherits(x, "POSIXt"))])
  # dayCol = names(df[sapply(df, function(x) inherits(x, "difftime"))])
  day = 'day'
  # set the variables for day of experiment, time of day and 
  df = df %>% ungroup #%>%
  #   dplyr::mutate(date = as.Date(date_time),
  #                 !!day := (julian(as.Date(!!as.symbol(date_timeCol), tz=attr(df[[date_timeCol]],"tzone"))) - min(julian(as.Date(!!as.symbol(date_timeCol), tz=attr(df[[date_timeCol]],"tzone"))))+1),
  #                 TOD = case_when(lubridate::hour(!!as.symbol(date_timeCol)) < 12 ~ "morning",
  #                                 TRUE ~ "evening")) %>%
  #   group_by(day, TOD) %>%
  #   dplyr::mutate(!!as.symbol(date_timeCol) := sort(!!as.symbol(date_timeCol)),
  #                 time_point = 1:n())
  
  # test that >1 day of data exist for mornings
  if(length(unique(df[which(grepl("morning",df[["TOD"]])),df[[day]]])) <= 1) stop("GPP calculations require >1 day of oxygen measurements")
  
  # Calculate NP
  dayVec = unique(df[[day]])
  morningDf = df[which(df[["time_point"]] == 1 & df[["TOD"]] == "morning"),]
  NP_mg_o2_l_hr = unlist(sapply(dayVec[-1], function(x) morningDf[which(morningDf[[day]] == x),"o2_do_mg_l"] - morningDf[which(morningDf[[day]]== (x-1)), "o2_do_mg_l"]))/unlist(sapply(dayVec[-1], function(x) morningDf[which(morningDf[[day]] == x), date_timeCol]-morningDf[which(morningDf[[day]] == (x-1)), date_timeCol]))
  # calculate R @ night
  nightDf = df[which(df[['TOD']] == 'evening'),]
  nightDf = nightDf[which(nightDf[['time_point']] == max(nightDf[['time_point']])),]
  Rnight_mg_o2_l_hr = unlist(sapply(dayVec[-length(dayVec)], function(x) nightDf[which(nightDf[[day]] == x), "o2_do_mg_l"] - morningDf[which(morningDf[[day]] == (x+1)),"o2_do_mg_l"]))/unlist(sapply(dayVec[-length(dayVec)], function(x) nightDf[which(nightDf[[day]] == x), date_timeCol] - morningDf[which(morningDf[[day]] == (x+1)), date_timeCol]))
  
  # estimate GP as the sum of NP & R
  
GP_mg_o2_l_hr = NP_mg_o2_l_hr + (abs(Rnight_mg_o2_l_hr))
  
  # return this list of GP, NP, and R estimates
  return(list(NP_mg_o2_l_hr = NP_mg_o2_l_hr,
              Rnight_mg_o2_l_hr = Rnight_mg_o2_l_hr,
              GP_mg_o2_l_hr = GP_mg_o2_l_hr))
}


#'
#'
#'
estimateContinuous = function(tankID = NULL,...){
require(LakeMetabolizer)
require(streamMetabolizer)
# set variables for metadata
site_coords = c(42.80049, -96.92677)
runID = substr(tankID, nchar(tankID), nchar(tankID))
# id the relevant file names
o2File = list.files(path = "./data/models/", pattern = paste0("o2_GAM_",tankID,".rds"), full.names = TRUE)
tempFile = list.files(path = "./data/models/", pattern = paste0("temp_GAM_",tankID,".rds"), full.names = TRUE)

# read the files in
o2Model = readRDS(o2File)
tempModel = readRDS(tempFile)

# make prediction DF

modelTimes = with(o2Model$data,
                  expand.grid(run_hr = seq(min(run_hr), max(run_hr), length = 100))) 

modelDf = modelTimes %>%
  bind_cols(predict(o2Model,
                    modelTimes,
                    se.fit = FALSE,
                    type = 'response')) %>%
  dplyr::select(run_hr, do_obs = 'Estimate') %>%
  bind_cols(predict(tempModel,
                    modelTimes,
                    se.fit = FALSE,
                    type = 'response')) %>%
  dplyr::select(run_hr, do_obs, temp_obs = 'Estimate') %>%
  dplyr::mutate(date_time = case_when(runID == 1 ~ as.POSIXct(run_hr*3600, origin = as.POSIXct("2022-06-21 05:30:00")),
                                      runID == 2 ~ as.POSIXct(run_hr*3600, origin = as.POSIXct("2022-06-27 05:30:00"))),
                solar_time = streamMetabolizer::calc_solar_time(date_time, longitude = site_coords[2]),
                light = streamMetabolizer::calc_light(solar_time,
                                                      latitude = site_coords[1],
                                                      longitude = site_coords[2]),
                do_sat = calc_DO_sat(temp_obs,
                                     pressure.air = (28.96/0.029530)),
                K600 = k.cole(data.frame(datetime = solar_time,
                                         wnd = wind.scale.base(1,1)))[['k600']],
                k_gas = k600.2.kGAS.base(k600 = K600, temp_obs, gas = "O2"), 
                z_mix = 0.3) %>%
  dplyr::mutate(across(matches('do_obs'), ~ifelse(.x <=0,0.01,.x))) %>%
  dplyr::select(date_time, do_obs, temp_obs, do_sat, light, k_gas, z_mix )

bayes.res = metab.bayesian(do.obs = modelDf[['do_obs']],
                           do.sat = modelDf[['do_sat']],
                           k.gas = modelDf[['k_gas']],
                           z.mix = modelDf[['z_mix']],
                           irr = modelDf[['light']],
                           wtr = modelDf[['temp_obs']],
                           datetime = modelDf[['date_time']])
saveRDS(bayes.res, file = paste0("./data/models/metab_bayes_",tankID,".rds"))
return(list(tankID = bayes.res,
            modelDf = modelDf))

}

#'
#'
# covariance matrix of outcomes
cor_matrix <- function(x, r, v = rep(1, length(x)), na.rm = FALSE) {
  mat <- diag(v)
  se <- sqrt(v)
  se[is.na(se)] <- 0
  if (length(x) > 1L) {
    for (i in 2:nrow(mat)) {
      for (j in 1:(i-1)) {
        if (x[i] == x[j]) {
          mat[i, j] <- mat[j, i] <- r * se[i] * se[j] 
        }
      }
    }
  }
  dimnames(mat) <- list(1:nrow(mat), 1:ncol(mat))
  if (na.rm) {
    keep <- !is.na(diag(mat))
    mat <- mat[keep, keep]
  }
  mat
}


fit_o2_gam = function(df = NULL,...){
  require(brms)
  require(mgcv)
  # if(!is.formula(form) | is.na(as.formula(form))) stop("`form` must be formula object or coercible")
  
  tankID =paste0("t",unique(df$tank))
  runID = gsub("run","",unique(df$run))
  file_name = paste0("./data/models/o2_GAM_",tankID,"_",runID,".rds")
  
  formDefault = bf(o2_do_mg_l ~ s(run_hr, bs = 'tp'))
  formOtherwise = bf(o2_do_mg_l ~ s(run_hr, bs = 'tp', k = 9))
  
  if(nrow(df) > 10){
    model = brm(formDefault,
                data = df,
                family = gaussian(),
                cores = 5, seed = 42,
                iter = 1.5e4, warmup = 1e4, thin = 10,
                file = file_name,
                control = list(adapt_delta = 0.99),
                save_pars = save_pars(all = TRUE),
                sample_prior = TRUE,
                file_refit = 'on_change',
                backend = 'cmdstanr')
  } else{
    model = brm(formOtherwise,
                data = df,
                family = gaussian(),
                cores = 5, seed = 42,
                iter = 1.5e4, warmup = 1e4, thin = 10,
                file = file_name,
                control = list(adapt_delta = 0.99),
                save_pars = save_pars(all = TRUE),
                sample_prior = TRUE,
                file_refit = 'on_change',
                backend = 'cmdstanr')
  }
  # 
  #     ,
  #                   otherwise = function(...){
  #                     brm(formOtherwise,
  #                         data = df,
  #                         family = gaussian(),
  #                         cores = 5, seed = 42,
  #                         iter = 1.2e4, warmup = 1e4, thin = 10,
  #                         file = file_name,
  #                         sample_prior = TRUE,
  #                         control = list(adapt_delta = 0.99),
  #                         save_pars = save_pars(all = TRUE),
  #                         file_refit = 'on_change',
  #                         backend = 'cmdstanr')
  #                     },
  #                   quiet = TRUE
  #                   )
}

#'
#'
#'
fit_temp_gam = function(df = NULL,...){
  require(brms)
  require(mgcv)
  # if(!is.formula(form) | is.na(as.formula(form))) stop("`form` must be formula object or coercible")
  
  tankID =paste0("t",unique(df$tank))
  runID = gsub("run","",unique(df$run))
  file_name = paste0("./data/models/temp_GAM_",tankID,"_",runID,".rds")
  
  formDefault = bf(temp_deg_c ~ s(run_hr, bs = 'tp'))
  formOtherwise = bf(temp_deg_c ~ s(run_hr, bs = 'tp', k = 9))
  
  if(nrow(df) > 10){
    model = brm(formDefault,
                data = df,
                family = gaussian(),
                cores = 6, seed = 42, chains = 3,
                threads = threading(2),
                iter = 1.5e4, warmup = 1e4, thin = 10,
                file = file_name,
                control = list(adapt_delta = 0.99),
                save_pars = save_pars(all = TRUE),
                sample_prior = TRUE,
                file_refit = 'on_change',
                backend = 'cmdstanr')
  } else{
    model = brm(formOtherwise,
                data = df,
                family = gaussian(),
                cores = 6, seed = 42, chains = 3,
                threads = threading(2),
                iter = 1.5e4, warmup = 1e4, thin = 10,
                file = file_name,
                control = list(adapt_delta = 0.99),
                save_pars = save_pars(all = TRUE),
                sample_prior = TRUE,
                file_refit = 'on_change',
                backend = 'cmdstanr')
  }
  # 
  #     ,
  #                   otherwise = function(...){
  #                     brm(formOtherwise,
  #                         data = df,
  #                         family = gaussian(),
  #                         cores = 5, seed = 42,
  #                         iter = 1.2e4, warmup = 1e4, thin = 10,
  #                         file = file_name,
  #                         sample_prior = TRUE,
  #                         control = list(adapt_delta = 0.99),
  #                         save_pars = save_pars(all = TRUE),
  #                         file_refit = 'on_change',
  #                         backend = 'cmdstanr')
  #                     },
  #                   quiet = TRUE
  #                   )
}

#'
#'
#'
extract_metab = function(metabModel = NULL,...){
tankID = gsub(".*_(\\w{1}\\d{1,2}_\\d{1}).rds","\\1",sapply(strsplit(metabModel, "/"),"[",length(unlist(strsplit(metabModel,"/")))))
x = readRDS(metabModel)
return(data.frame(tank = tankID,
                  x$metab,
                  x$metab.sd))
}

#'
#'
#'
convert_metab = function(metabModel = NULL,...){
  tankID = gsub(".*_(\\w{1}\\d{1,2}_\\d{1}).rds","\\1",sapply(strsplit(metabModel, "/"),"[",length(unlist(strsplit(metabModel,"/")))))
  x = readRDS(metabModel)
  n_obs = environment(x[["model"]][["model"]][["data"]])[["data"]][["N"]]
  gpp.coef = x$model$BUGSoutput$sims.list$C[,1]
  r.coef = x$model$BUGSoutput$sims.list$C[,2]
  par = environment(x[["model"]][["model"]][["data"]])[["data"]][["U"]][,1]
  lntemp = environment(x[["model"]][["model"]][["data"]])[["data"]][["U"]][,2]
  gpp.out = sapply(gpp.coef, function(a) mean(a*par)*n_obs)
  r.out = sapply(r.coef, function(a) mean(a * lntemp)*n_obs)
  nep.out = gpp.out + r.out
  return(list(summary = data.frame(tank = tankID,
                                   x$metab,
                                   x$metab.sd),
              gpp.out = gpp.out,
              r.out = r.out,
              nep.out = nep.out))
}
