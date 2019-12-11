#Download and install required packages
install.packages("biwavelet")
install.packages("raster")
library(tidyverse)
library(biwavelet)
library(lubridate)
library(padr)
library(gridExtra)
library(xts)
library(sf)
library(raster)

##############################################################################

#Project research questions:
# *What was the maximum salinity of the Hurricane Florence and Dorian storm surges and how long did it take for salinity levels to return to average site conditions following each event? 
#   *During the flash drought of 2019 what patterns do we see in our salinity data?
#   *Can wavelet coherance diagrams show us both hurricanes and the flash drought in time frequency space? 



################################################################################

#import and tidy the site data-------
#step one. Import data for AR0
read_csv("data/AR0_Full_111819.csv") -> ar0
#change datetime column to datetime format in R
ar0$record_date_time <- mdy_hm(ar0$record_date_time)
#fill data gaps so there are none. observations of filled times will be NA
ar0_no_date_gap <- ar0%>% pad %>% fill_by_value(NA)
#mutate in salinity column, units = ppt
ar0 %>%
  mutate(spc_ms_cm = `specific_conductance_microScm-1`/1000) %>%
  mutate(salinity = 0.008+(-0.1692*(sqrt(spc_ms_cm/53.0647665)))+(0.47837957*spc_ms_cm)+(14.091*((spc_ms_cm/53.0647665)^1.5))+(-7.0261*((spc_ms_cm/53.0647665)^2))+(2.7081*((spc_ms_cm/53.0647665)^2.5))+(8.60585198*(0.0005+(-0.0056*(sqrt(spc_ms_cm/53.0647665)))+(-0.0066*(spc_ms_cm/53.0647665))+(-0.0375*((spc_ms_cm/53.0647665)^1.5))+(0.0636*((spc_ms_cm/53.0647665)^2))+(-0.0144*((spc_ms_cm/53.0647665)^2.5)))))%>%
  mutate(site = "ar0") -> ar0

#create a subset of only 2018 and 2019 data
ar0_2019 <- subset(ar0, format(as.Date(record_date_time),"%Y")==2019)
ar0_2018 <- subset(ar0, format(as.Date(record_date_time),"%Y")==2018)
ar0_2018_2019 <- full_join(ar0_2018,ar0_2019)
#aggregate the data to an hourly time scale for faster plotting and processing
ar0_2018_2019_hourly <- xts(ar0_2018_2019, order.by=ar0_2018_2019$record_date_time)
ends <- endpoints(ar0_2018_2019_hourly,'hours',1) 
period.apply(ar0_2018_2019_hourly,ends ,mean)->ar0_2018_2019_hourly
data.frame(date_time=index(ar0_2018_2019_hourly), coredata(ar0_2018_2019_hourly))->ar0_2018_2019_hourly
#remove last row from dataframe else the timestep will not be equal for wavelet analysis
head(ar0_2018_2019_hourly,-1)->ar0_2018_2019_hourly
#exploratory data analysis- plot time series-----

##add more options to this plot
##change labels, change theme, at title, look into changing the x axis label format -> FIGURE 1
ggplot(ar0_2018_2019_hourly, mapping = aes(x=date_time, y=salinity))+
  geom_point()
#
#
#create site map -> FIGURE 2-----
#
#read in spatial data 
read_sf("data/AGU_sites.shp") -> sites
read_sf("data/AR01.shp")->ar01_point
read_sf("data/app_canal_clp.shp")->canals
raster("data/elev_APP3m.tif")->app_dem
#check projection codes
st_crs(sites)
st_crs(ar01_point)
st_crs(canals)
#Assign epsg codes
st_transform(sites, 4326)->sites
st_transform(ar01_point,4326->ar01_point)
st_transform(canals,4326)->canals
projectRaster(app_dem,crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#what are average salinity conditions for the site? these are the conditions we will use as a benchmark for determining how long salinity stayed elevated. 

mean(ar0_2018_2019_hourly$salinity,na.rm = TRUE) -> site_avg_salt


#hurricane analyis-----
#Hurricane Florence made landfall in NC on 2018/09/14 as a category 1 hurricane. 
#Hurricane Dorian made landfall in NC on 209/09/06 as a cat 1 hurricane.
#step 1. take a month long subset surrounding each hurricane

florence <- subset(ar0_2018_2019_hourly, date_time> ymd_hms('2018-09-01 00:45:00') & date_time < ymd_hms('2019-01-01 00:45:00'))
dorian <-subset(ar0_2018_2019_hourly, date_time> ymd_hms('2019-09-06 00:45:00') & date_time < ymd_hms('2019-10-06 00:45:00'))

#exploratory plots of salinity data from hurricanes
ggplot(florence, mapping = aes(x=date_time, y=salinity))+
  geom_point()

ggplot(dorian, mapping = aes(x=date_time, y=salinity))+
  geom_point()
#What was the maximum salinity for each event? 

  
subset(florence,date_time>ymd_hms('2018-09-14 00:45:00') & date_time > ymd_hms('2018-10-01 0:45:00')) ->flor_period
max(flor_period$salinity, na.rm=TRUE) -> max_florence_salinity
subset(dorian,date_time>ymd_hms('2019-09-06 00:45:00') & date_time > ymd_hms('2019-09-20 0:45:00')) ->dor_period
max(dor_period$salinity) -> max_dorian_salinity

#determine how long it took for salinity to return to average site salinity
#create seperate object to drop the NAs - interpolation for small gaps? 
florence_salt <- na.omit(flor_period$salinity)
elevated_salinity_dur <-  rep(NA, length(florence_salt))
dur = 0
for (record in 1:length(elevated_salinity_dur)){ 
  if (florence_salt[record] < 5){
    if (dur>0){
      elevated_salinity_dur[record] <- dur 
      dur = 0
    }


  }
  else{
    dur=dur+1
  }
  
}

#find a way to work in a duration below the threshold idea





# florence %>%
#   filter(salinity >= 5)->elevated_salt
# 
# datenums <- as.numeric(elevated_salt$date_time)
# 
# delta_t <- diff(datenums)/3600
# append(delta_t,0, after = length(dt))->delta_t
#     
# elevated_salt %>%
#   mutate(time_diff = delta_t)->elevated_salt
# 
# elevated_salinity_dur = 0
# for (record in elevated_salt){
#   while (delta_t == 1){
#     elevated_salinity_dur = elevated_salinity_dur+1
#   }
# }
# 
# for (record in florence){
#   if salintity < threshold, 
# }
#if_else not written correctly. write the pseudo-code first and then write the code. code goal ^ from the start of the record, determine when the salinity first becomes elevated and start to count the duration there. have the for loop break when the salinity threshold is crossed again. would be simplest way to do it for isolated events, but will need to be modified when looking at the entire timeseries. 

#wtc test -> use only targeted analysis periods for the wavelet analysis for this. computer is not powerful enough to run the code -> FIGURE 3 -------

data.frame(ar0_2018_2019_hourly$date_time, ar0_2018_2019_hourly$salinity)->salinity_ts
data.frame(ar0_2018_2019_hourly$date_time, ar0_2018_2019_hourly$water_level_m)->lvl_ts
salinity_lvl_wtc <- wtc(salinity_ts, lvl_ts)

