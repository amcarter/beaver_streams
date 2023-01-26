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

# load datasets: ####
dat <- read_csv('data/prepared/EB100_1_7.7.19.csv')

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
bayes_specs <- specs(bayes_name, burnin_steps = 1000, saved_steps = 1000)
#Run models
dat <- read_csv('data/prepared/EB100_1_7.7.19.csv')
fit <- streamMetabolizer::metab(bayes_specs, data = dat)
saveRDS(fit, 'data/metab_fits/EB100_1_7.7.19.rds')

dat <- read_csv('data/prepared/EBD5_1_7.7.19.csv')
fit <- streamMetabolizer::metab(bayes_specs, data = dat)
saveRDS(fit, 'data/metab_fits/EBD5_1_7.7.19.rds')

dat <- read_csv('data/prepared/ELKD5_1_7.3.19.csv')
fit <- streamMetabolizer::metab(bayes_specs, data = dat)
saveRDS(fit, 'data/metab_fits/ELKD5_1_7.3.19.rds')

dat <- read_csv('data/prepared/ELKU100_1_7.3.19.csv')
fit <- streamMetabolizer::metab(bayes_specs, data = dat)
saveRDS(fit, 'data/metab_fits/ELKU100_1_7.3.19.rds')

dat <- read_csv('data/prepared/WBXD5_1_7.9.19.csv')
fit <- streamMetabolizer::metab(bayes_specs, data = dat)
saveRDS(fit, 'data/metab_fits/WBXD5_1_7.9.19.rds')

dat <- read_csv('data/prepared/WBXU100_1_7.9.19.csv')
fit <- streamMetabolizer::metab(bayes_specs, data = dat)
saveRDS(fit, 'data/metab_fits/WBXU100_1_7.9.19.rds')


