###################################################
# Extract and plot fits from stream metabolizer
###################################################

library(streamMetabolizer)

    fit <- readRDS('data/metab_fits/EBD5_1_7.7.19.rds')
plot_DO_preds(fit)
fit@fit$daily %>% select(date, K600 = K600_daily_mean, ends_with('Rhat'))
predict_metab(fit)
plot_metab_preds(fit)
