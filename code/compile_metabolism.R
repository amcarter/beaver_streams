##########################################################
# Extract and compile metabolism estimates from model fits
##########################################################

library(streamMetabolizer)
library(tidyverse)

# get list of model fits
files <- list.files('data/metab_fits/')
files <- grep('compiled', files, value = TRUE, invert = TRUE)
files <- grep('K600', files, value = TRUE, invert = TRUE)
comp_met <- data.frame()

for(i in 1:length(files)){
    site <- substr(files[i], 1, nchar(files[i])-5)
    location <- str_match( files[i], '([UD])\\.rds$')[,2]
    fit <- readRDS(paste0('data/metab_fits/', files[i]))
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
    p <- plot_DO_preds(fit)
    print(p)
    comp_met <- bind_rows(comp_met, met)
}

compiled_met <- comp_met %>%
    mutate(distance = case_when(location == 'U' ~ -100,
                                location == 'D' ~ 5),
           location = case_when(location == 'U' ~ 'upstream',
                                location == 'D' ~ 'downstream'),
           location = factor(location, levels = c('upstream', 'downstream')))

write_csv(compiled_met, 'data/metab_fits/compiled_metabolism_K600sigsig_0.05.csv')


compiled_met05<- read_csv( 'data/metab_fits/compiled_metabolism_K600sigsig_0.05.csv') %>%
    mutate(K600sigsig = 0.05)
compiled_met1<- read_csv( 'data/metab_fits/compiled_metabolism_K600sigsig_0.1.csv') %>%
    mutate(K600sigsig = 0.1)
compiled_met<- read_csv( 'data/metab_fits/compiled_metabolism.csv')%>%
    mutate(K600sigsig = 0.7)

compiled_met %>%
    mutate(site = paste0(site, distance)) %>%
    ggplot(aes(K600, ER, col = site)) +
    geom_point() +
    geom_smooth(method = 'lm', se = FALSE) +
    theme_classic()

compiled_met %>%
    mutate(K600_Rhat = case_when(K600_Rhat > 1.1 ~ 'Rhat > 1.1',
                                 K600_Rhat <= 1.1 ~ 'Rhat < 1.1'))%>%
    ggplot(aes(date, K600, col = K600_Rhat, group = location)) +
    geom_line(col = 'grey') +
    geom_point() +
    facet_wrap(site ~ ., scales = 'free_x') +
    ylab('metabolism (g O2/m2/d)') +
    theme_bw()

compiled_met %>%
    ggplot(aes(date, GPP)) +
    geom_point(col = 'forestgreen') +
    geom_line(col = 'forestgreen') +
    geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), col = 'forestgreen')+
    geom_point(aes(y = ER), col = 'brown3') +
    geom_line(aes(y = ER), col = 'brown3') +
    geom_errorbar(aes(ymin = ER.lower, ymax = ER.upper), col = 'brown3')+
    facet_grid(site ~ location, scales = 'free_x') +
    ylab('metabolism (g O2/m2/d)') +
    theme_bw()

bind_rows(compiled_met, compiled_met1, compiled_met05)%>%
    ggplot(aes(date, ER, col = factor(K600sigsig))) +
    geom_line() + geom_point() +
    facet_grid(site~location, scales = 'free_x')
