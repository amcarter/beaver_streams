###################################################################
# Script to run streamMetabolizer on prepared datasets
##################################################################

#Reinstall unitted and streamMetabolizer if needed
# remotes::install_github('appling/unitted', force = TRUE)
# remotes::install_github("GLEON/LakeMetabolizer", force = TRUE)

#load all packages
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(streamMetabolizer)
library(lubridate)
library(tidyverse)
library(dygraphs)

files <- list.files('data/prepared/')
# load datasets: ####
dat <- read_csv('data/prepared/EBU.csv')

# Visualize the data #####
head(dat) ; tail(dat)

dat %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)',
            light='PAR\n(umol m^-2 s^-1)')
dat %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

#sM model specs--------------- ####
bayes_name <- mm_name(type='bayes', pool_K600='normal',
                      err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_specs <- specs(bayes_name, burnin_steps = 1000, saved_steps = 1000,
                     K600_daily_meanlog_sdlog = 0.05)
#Run models
for(i in 1:length(files)){

    ff <- files[i]
    dat <- read_csv(paste0('data/prepared/', ff))
    fit <- streamMetabolizer::metab(bayes_specs, data = dat)
    saveRDS(fit, paste0('data/metab_fits/',
                        substr(ff, 1, nchar(ff)-3), 'rds'))
}




####
# Rerun models for ELK sites with tighter parameters on K600 sigma sigma:
#sM model specs--------------- ####
bayes_name <- mm_name(type='bayes', pool_K600='normal',
                      err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_specs <- specs(bayes_name, burnin_steps = 1000, saved_steps = 1000,
                     K600_daily_meanlog_sdlog = 0.01)
#Run models
ELK <- files[c(5,6)]
for(i in 1:length(ELK)){

    ff <- ELK[i]
    dat <- read_csv(paste0('data/prepared/', ff))
    fit <- streamMetabolizer::metab(bayes_specs, data = dat)
    saveRDS(fit, paste0('data/metab_fits/',
                        substr(ff, 1, nchar(ff)-4), 'K600sigsig_01.rds'))
}

comp_met <- data.frame()

for(i in 1:length(ELK)){
    ff <- ELK[i]
    site <- substr(ff, 1, nchar(ff)-5)
    location <- str_match( ff, '([UD])\\.csv$')[,2]
    fit <- readRDS(paste0('data/metab_fits/', substr(ff, 1, nchar(ff)-4),
                          'K600sigsig_01.rds'))
    met <- fit@fit$daily %>%
        select(date, K600 = K600_daily_mean, GPP_Rhat, ER_Rhat,
               K600_Rhat = K600_daily_Rhat) %>%
        mutate(site = site,
               location = location)

    met <- predict_metab(fit) %>%
        left_join(met, by = 'date') %>%
        filter(!is.na(GPP)) %>%
        relocate(site) %>%
        relocate(c(msgs.fit, warnings, errors), .after = K600_Rhat)
    p <- plot_DO_preds(fit) + ggtitle(paste(site, location, sep = " "))
    print(p)
    comp_met <- bind_rows(comp_met, met)
}

compiled_met <- comp_met %>%
    mutate(distance = case_when(location == 'U' ~ -100,
                                location == 'D' ~ 5),
           location = case_when(location == 'U' ~ 'upstream',
                                location == 'D' ~ 'downstream'),
           location = factor(location, levels = c('upstream', 'downstream')))

write_csv(compiled_met, 'data/metab_fits/compiled_metabolism_K600sigsig_01.csv')


# Rerun models for EBU and WBXD sites with looser priors on K600 sigma sigma:
#sM model specs--------------- ####
bayes_name <- mm_name(type='bayes', pool_K600='normal',
                      err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_specs <- specs(bayes_name, burnin_steps = 1000, saved_steps = 1000,
                     K600_daily_meanlog_sdlog = 0.5)
#Run models
EB_WBX <- files[c(4,11)]
for(i in 1:length(EB_WBX)){

    ff <- EB_WBX[i]
    dat <- read_csv(paste0('data/prepared/', ff))
    fit <- streamMetabolizer::metab(bayes_specs, data = dat)
    saveRDS(fit, paste0('data/metab_fits/',
                        substr(ff, 1, nchar(ff)-4), 'K600sigsig_5.rds'))
}

comp_met <- data.frame()

for(i in 1:length(EB_WBX)){
    ff <- EB_WBX[i]
    site <- substr(ff, 1, nchar(ff)-5)
    location <- str_match( ff, '([UD])\\.csv$')[,2]
    fit <- readRDS(paste0('data/metab_fits/', substr(ff, 1, nchar(ff)-4),
                          'K600sigsig_1.rds'))
    met <- fit@fit$daily %>%
        select(date, K600 = K600_daily_mean, GPP_Rhat, ER_Rhat,
               K600_Rhat = K600_daily_Rhat) %>%
        mutate(site = site,
               location = location)

    met <- predict_metab(fit) %>%
        left_join(met, by = 'date') %>%
        filter(!is.na(GPP)) %>%
        relocate(site) %>%
        relocate(c(msgs.fit, warnings, errors), .after = K600_Rhat)
    p <- plot_DO_preds(fit) + ggtitle(paste(site, location, sep = " "))
    print(p)
    comp_met <- bind_rows(comp_met, met)
}

compiled_met <- comp_met %>%
    mutate(distance = case_when(location == 'U' ~ -100,
                                location == 'D' ~ 5),
           location = case_when(location == 'U' ~ 'upstream',
                                location == 'D' ~ 'downstream'),
           location = factor(location, levels = c('upstream', 'downstream')))

write_csv(compiled_met, 'data/metab_fits/compiled_metabolism_K600sigsig_5.csv')

