################################################################################
# This script prepares raw O2 data from yellowstone beaver streams by
# downloading air pressure data, calculating DO saturation, and
# adjusting to solar time with modeled light.
################################################################################

#Reinstall unitted and streamMetabolizer if needed
# remotes::install_github('appling/unitted', force = TRUE)
# remotes::install_github("USGS-R/streamMetabolizer", force = TRUE)
# remotes::install_github("GLEON/LakeMetabolizer", force = TRUE)
# remotes::install_github("streampulse/StreamPULSE")

#load all packages
library(tidyverse)
library(streamMetabolizer)
library(LakeMetabolizer)
library(StreamPULSE)
library(lubridate)
library(zoo)

# This function downloads data the nearest noaa station:
# lat <- 44.94998
# lon <- -109.4337
# air_pres <- StreamPULSE:::FindandCollect_airpres(lat, lon, '2019-01-01', '2019-12-02')
# air_pres <- data.frame(DateTime_UTC = seq(air_pres$DateTime_UTC[1],
#                                           air_pres$DateTime_UTC[nrow(air_pres)],
#                                           by = 'min')) %>%
#     left_join(air_pres) %>%
#     mutate(air_kPa = na.approx(air_kPa, na.rm = F)) %>% select(-air_temp)
# write_csv(air_pres, 'data/NOAA_air_pressure.csv')
air_pres <- read_csv('data/NOAA_air_pressure.csv')

# Read in data file and format:
files <- list.files(path = 'data/raw/')
for(i in 1:length(files)){

    dat <- read_csv(paste0('data/raw/', files[i]))

    # reformat data columns
    dat <- dat %>% slice(-1)%>%
        select(-DO.sat, -date.time_Q_Z)%>%
        mutate(MST.time = as.POSIXct(MST.time, tz = 'MST',
                                     tryFormats = c('%m/%d/%Y %H:%M:%S',
                                                    '%m/%d/%Y %H:%M')),
               DateTime_MST = round_date(MST.time, unit = 'minute'),
               across(c(DO.obs, depth, temp.water, light, discharge, elev,
                        dist.d100, obsDO.satprct),
                          as.numeric)) %>%
        tidyr::extract(lat, into = c('degs', 'mins'), regex = '([0-9]+). ([0-9]+\\.[0-9]+)') %>%
        mutate(lat = as.numeric(degs) + as.numeric(mins)/60) %>%
        tidyr::extract(long, into = c('degs', 'mins'), regex = '([0-9]+). ([0-9]+\\.[0-9]+)') %>%
        mutate(long = -as.numeric(degs) + as.numeric(mins)/60) %>%
        select(-MST.time, -degs, -mins)


    # round times to nearest minute and convert to solar time
    dat <- dat %>%
        mutate(DateTime_UTC = with_tz(DateTime_MST, tz = 'UTC'),
               solar.time = convert_UTC_to_solartime(
                   DateTime_UTC, dat$long[1],
                   time.type = c("mean solar")))

    # # Generate light
    dat$light<-calc_light(dat$solar.time,
                          latitude = dat$lat[i],
                          longitude = dat$long[i])

    # Add NOAA air pressure data
    dat <- left_join(dat, air_pres, by = 'DateTime_UTC') %>%
        mutate(air_pres_mbar = air_kPa *10)

    # calculate DO saturation based on temp and pressure
    dat$DO.sat <- o2.at.sat.base(temp = dat$temp.water,
                            # baro = dat$air_pres_mbar,
                            altitude = dat$elev[1])

    #Compile necessary parts for sM model
    mod_dat <- select(dat,
                      solar.time, DO.obs, DO.sat,
                      depth, temp.water, light)

    # Make the time steps 15 minutes

    write_csv(mod_dat, paste0('data/prepared/', files[i]))

}


