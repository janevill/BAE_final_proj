#Download and install required packages
install.packages("biwavelet")
install.packages("raster")
install.packages("rgdal")
library(tidyverse)
library(biwavelet)
library(lubridate)
library(padr)
library(gridExtra)
library(xts)
library(sf)
library(raster)
library(rgdal)
library(ggspatial)

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
  geom_point(color ="black", fill = "grey", pch = 21)+
  guides(fill=FALSE)+
  labs(x = "", y= "Salinity (ppt)", caption = "fig2a. Time series data for 2018 and 2019")+
  theme_classic() -> fig2a
  
#
#
#create site map -> FIGURE 2-----
library(RColorBrewer)
#read in spatial data 
read_sf("data/spatial_data/AR01.shp")->ar01_point
read_sf("data/spatial_data/panhand_canals.shp")->canals
raster("data/spatial_data/panhandle_3m.tif")->panhandle_dem
#check projection codes
st_crs(ar01_point)
st_crs(canals)
projection(panhandle_dem)
#Assign epsg codes
st_transform(ar01_point,CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))->ar01_point
st_transform(canals,CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))->canals
projectRaster(panhandle_dem,crs="+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
#uncomment the above line if the raster has no projection. this line takes a very long time to run
#
#turn raster layer into datafrane
dem_points<-data.frame(rasterToPoints(panhandle_dem))
colnames(dem_points) <- c("X","Y","Elevation")

#create sitemap
n=3
breaks = seq(min(dem_points$Elevation),max(dem_points$Elevation), length.out = n)
ggplot()+
  geom_raster(dem_points, mapping = aes(x= X, y = Y, fill = Elevation))+
  scale_fill_gradientn(limits = c(0,5),
  na.value = "black",
  colors = c("black", "white"),
  name= "Elevation (m)")+
  labs(x = "", y = "", caption = "Pan Handle Portion of the Albemarle - Pamlico Peninsula \n located in northeastern North Carolina")+
  theme(axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.title.y=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank())+
  geom_sf(data=canals,color = "blue")+
  geom_sf(data=ar01_point, fill = "red", color = "black",pch=21, size = 4, stroke = 1)+
  scale_x_continuous(breaks = breaks)+
  annotation_scale(location = "br")+
  annotation_north_arrow(location="tr",pad_x=unit(0.25, "cm"), style=north_arrow_fancy_orienteering())+
  theme_classic()
#summarize Data-------
data_summary <- data.frame(Variable= c("Specific Conductance","Water Level", "Salinity" ), Units=c("microS/cm", "m","ppt"))
kable(data_summary)  

#what are average salinity conditions for the site? these are the conditions we will use as a benchmark for determining how long salinity stayed elevated. ------

mean(ar0_2018_2019_hourly$salinity,na.rm = TRUE) -> site_avg_salt


#hurricane analyis-----
#Hurricane Florence made landfall in NC on 2018/09/14 as a category 1 hurricane. 
#Hurricane Dorian made landfall in NC on 209/09/06 as a cat 1 hurricane.
#step 1. take a month long subset surrounding each hurricane

florence <- subset(ar0_2018_2019_hourly, date_time> ymd_hms('2018-09-01 00:45:00') & date_time < ymd_hms('2019-01-01 00:45:00'))
dorian <-subset(ar0_2018_2019_hourly, date_time> ymd_hms('2019-09-06 00:45:00') & date_time < ymd_hms('2019-10-06 00:45:00'))

#exploratory plots of salinity data from hurricanes
ggplot(florence, mapping = aes(x=date_time, y=salinity))+
  geom_point(aes(fill = florence$water_level_m),color ="black", pch = 21)+
  labs(x = "", y= "Salinity (ppt)", caption = "fig3. Hurricane Florence (2018)", fill="Water Level (m)")+
  scale_fill_continuous(type = "viridis")+
  theme_classic() +
  theme(plot.caption = element_text(hjust = 0, face = "italic"))#->fig2b


ggplot(dorian, mapping = aes(x=date_time, y=salinity))+
  geom_point(aes(fill = dorian$water_level_m),color ="black", pch = 21)+
  labs(x = "", y= "Salinity (ppt)", caption = "fig2c. Hurricane Dorian (2019)", fill = "Water Level (m)")+
  scale_fill_continuous(type = "viridis")+
  #scale_x_date(labels=date_format("%b-%y")) +
  theme_classic() +
  theme(plot.caption = element_text(hjust = 0, face = "italic"))->fig2c

grid.arrange(fig2a,fig2b,fig2c, nrow =3)

#What was the maximum salinity the period of analysis


  
filter(florence,date_time>=ymd_hms('2018-09-14 00:45:00') & date_time <ymd_hms('2018-10-15 0:45:00')) ->flor_period
max(flor_period$salinity, na.rm=TRUE) -> max_florence_salinity
filter(dorian,date_time>=ymd_hms('2019-09-06 00:45:00') & date_time < ymd_hms('2019-09-20 0:45:00')) ->dor_period
max(dor_period$salinity) -> max_dorian_salinity

#create tables with nominal stats for florence and dorain
flor_period_table_data <- data.frame(flor_period$date_time, flor_period$water_level_m, flor_period$water_temperature_degC, flor_period$salinity, flor_period$mean_wind_speed_ms.1)
florence_summary <- summary(flor_period_table_data)
kable(florence_summary)


dor_period_table_data <- data.frame(dor_period$date_time, dor_period$water_level_m, dor_period$water_temperature_degC, dor_period$salinity, dor_period$mean_wind_speed_ms.1)
dorian_summary <- summary(dor_period_table_data)
kable(dorian_summary)

ggplot(flor_period, mapping=aes(x=date_time, y=salinity))+
  geom_point()
#
ggplot(dor_period, mapping=aes(x=date_time, y=salinity))+
  geom_point()
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
elevated_salinity_dur <- na.omit(elevated_salinity_dur)
#the increased salinity from hurricane florence lingered for about 23 days
#repeat for hurricane Dorian
#
dor_salt <- na.omit(dor_period$salinity)
elev_salt_dur_dor <-  rep(NA, length(dor_salt))
dur = 0
for (record in 1:length(elev_salt_dur_dor)){ 
  if (dor_salt[record] < 5){
    if (dur>0){
      elev_salt_dur_dor[record] <- dur 
      dur = 0
    }
  }
  else{
    dur=dur+1
  }
}
elev_salt_dur_dor <- na.omit(elev_salt_dur_dor)
#elevated salt from hurricane dorian only lingered 19 hours 
#
#find a way to work in a duration below the threshold idea
#
#
##   *During the flash drought of 2019 what patterns do we see in our salinity data?
##  The 2019 Flash Drought occurred from 2019-09-20 00:45:00 to 2019-10-07 23:45:00
##   plot the flash drought
##   turn into a dual plot with precipitation on the y-axis
filter(ar0_2018_2019_hourly,date_time>=ymd_hms("2019-09-20 00:45:00") & date_time <ymd_hms('2019-10-07 23:45:00')) ->flash_d_period

ggplot(flash_d_period)+
  geom_point(mapping = aes(x=date_time, y= salinity),size = 2, fill = "red",pch =21,rcolor = "black")+
  labs(x= "",y= "Salinity (ppt)", caption = "fig 3. Salinity during the flash drought period of 2019")+
  guides(fill=FALSE)+
  theme_classic()
#to understand the pattern here, lets calculate the duration of the drought induced saltwater intrusion events

flash_d_salt <- na.omit(flash_d_period$salinity)
flash_d_int <-  rep(NA, length(flash_d_salt))
int = 0
for (record in 1:length(flash_d_int)){ 
  if (flash_d_salt[record] >1){
    if (int>0){
      flash_d_int[record] <- int 
      int = 0
    }
  }
  else{
    int=int+1
  }
}
flash_d_int <- na.omit(flash_d_int)
mean(flash_d_int)

#compute continuous wavelet transform to find at what timescale we see significant variability in the flash drought record
filter(ar0_2018_2019_hourly,date_time>=ymd_hms("2019-09-28 00:45:00") & date_time <ymd_hms('2019-10-07 23:45:00')) ->cwt_fd_period
data.frame(cwt_fd_period$date_time, cwt_fd_period$salinity)->fd_salinity_ts
fd_salinity_ts$cwt_fd_period.salinity[is.na(fd_salinity_ts$cwt_fd_period.salinity)] <- -999
cwt_fd <- wt(fd_salinity_ts)
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(cwt_fd, plot.cb=TRUE,xaxt="n" )
axis(1,at = c(min(flash_d_period$date_time),median(flash_d_period$date_time),max(flash_d_period$date_time)),labels = c(min(flash_d_period$date_time),median(flash_d_period$date_time),max(flash_d_period$date_time)))


#average duration of events is close to 12 hours leading us to believe this is  a lunar tide phenomena. We can use wavelets to test this hypothesis. If we see a clear band across the 
#continuous wavelet transform, it tells us there is significant variability in our data across this time-scale. 
#
##   *Can wavelet coherance diagrams show us both hurricanes and the flash drought in time frequency space? 

#wtc test -> use only targeted analysis periods for the wavelet analysis for this. computer is not powerful enough to run the code -> FIGURE 3 -------

data.frame(flor_period$date_time, flor_period$salinity)->salinity_ts
data.frame(flor_period$date_time, flor_period$water_level_m)->lvl_ts
#change all NAs to -99999
salinity_lvl_wtc <- wtc(salinity_ts, lvl_ts)
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(salinity_lvl_wtc, plot.cb = TRUE,plot.phase = TRUE, main = "Salinity vs. Depth Wavelet Coherence \n Hurricane Florence",xaxt = 'n')
axis(1,at = c(min(flor_period$date_time),median(flor_period$date_time),max(flor_period$date_time)),labels = c(min(flor_period$date_time),median(flor_period$date_time),max(flor_period$date_time)))


data.frame(dor_period$date_time, dor_period$salinity)->dor_salinity_ts
data.frame(dor_period$date_time, dor_period$water_level_m)->dor_lvl_ts
#change all NAs to -99999
dor_salinity_lvl_wtc <- wtc(dor_salinity_ts, dor_lvl_ts)
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(dor_salinity_lvl_wtc, plot.cb = TRUE,plot.phase = TRUE, main = "Salinity vs. Depth Wavelet Coherence \n Hurricane Dorian",xaxt = 'n')
axis(1,at = c(min(dor_period$date_time),median(dor_period$date_time),max(dor_period$date_time)),labels = c(min(dor_period$date_time),median(dor_period$date_time),max(dor_period$date_time)))
