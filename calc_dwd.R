library(tidyverse)
library(stringr)
library(dplyr)
library(rgdal)
library(raster) 
library(rgeos)
library(sf)
library(mlogit)

library(fGarch) # skewed normal dist
library(abind) # combine arrays
# Statistical mode
library(DescTools) 

library(rdwd)
library(pbapply)

## Load data monthly gridded data from DWD database

# Load string vector with all vailable directories (CDC)
data("gridIndex")

# Filter growing season subset, only directories with monthly gridded data from March till 
# September for the 10 years prior to 2020 are relevant (see Stetter,Sauer (2022))
years=as.character(2010:2019)

base_dwd <- as.data.frame(gridIndex) %>%
 mutate("month"=as.numeric(str_sub(gridIndex,-09,-8)),
        "year"=as.numeric(str_sub(gridIndex,-13,-10))) %>% 
 filter(year %in% years, 
        month %in% 3:10)


# growing season precipitation string index
precip_index <- base_dwd  %>% 
 filter(str_detect(gridIndex, "monthly_precipitation"))

# growing season mean temperature string index
temp_index <- base_dwd  %>% 
 filter(str_detect(gridIndex, "monthly_air_temp_mean"))

# Download monthly data
raster_rain <- stack(dataDWD(precip_index[[1]], base=gridbase, joinbf=TRUE, dir=locdir(), dividebyten=FALSE))
raster_temp <- stack(dataDWD(temp_index[[1]], base=gridbase, joinbf=TRUE, dir=locdir()))

# Assign prjection 
# (https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_mean/DESCRIPTION_gridsgermany_monthly_air_temperature_mean_en.pdf)
crs(raster_rain) <- crs(raster_temp) <- "EPSG:31467"

## add choice experiment data to allocate observations in space
dat <- readRDS("../dat/choice_data/data_calc.rds")
shapeBY <- st_read("../dat/spatial_data/bayern_ex.shp")[,"SCH"]

## Match projections
shapeBY <- st_transform(shapeBY, st_crs(raster_temp))

## Delete unnecessary files to save memory
rm(gridIndex, precip_index, temp_index)



## 1.5 Shrink raster dataset
raster_temp <- mask(raster_temp, shapeBY)
raster_rain <- mask(raster_rain, shapeBY)


## 2.1 Yearly average temperature growing season
### Annual GS data at raster level
list_temp <- pblapply(years, function(i) 
 raster_temp[[tidyselect::vars_select(
  names(raster_temp), 
  contains(i))]])

list_temp <- lapply(list_temp, function(x) 
 stackApply(x, indices = 1, fun='mean'))

#--- convert to an sf object ---#
list_temp <- pblapply(list_temp, function(i)
 st_as_sf(
  as(i, 'SpatialPolygonsDataFrame')
 )
)


# Make dataframe
df_temp <- list_temp %>% 
 purrr::map(dplyr::select, 1) %>% 
 purrr::reduce(st_join, join=st_equals)
colnames(df_temp) <- c(paste0("temp_", years), "geometry") 


saveRDS(df_temp, '../dat/dwd_data/df_temp_dwd.rds')

# remove redundant files to save memory space
rm(list_temp, raster_temp)

## 2.2 Yearly precip sum - growing season
### Annual GS data at raster level
list_rain <- lapply( years, function(i) 
 raster_rain[[tidyselect::vars_select(
  names(raster_rain), 
  contains(i))]])

list_rain <- pblapply(list_rain, function(x) 
 stackApply(x, indices = 1, fun='sum'))

#--- convert to an sf object ---#
list_rain <- pblapply(list_rain, function(i)
 st_as_sf(
  as(i, 'SpatialPolygonsDataFrame')
 )
)

# Make dataframe
df_rain <- list_rain1 %>% 
 purrr::map(dplyr::select, 1) %>% 
 purrr::reduce(st_join, join=st_equals)
colnames(df_rain) <- c(paste0("rain_", years), "geometry") 


saveRDS(df_temp, '../dat/dwd_data/df_rain_dwd.rds')



#### Other indicators: daily data ####
growing_season <- paste0(rep(years, each=8), ".", c(rep(0,7), 1) , c(3:9, 0), ".")
shapeBY <- st_transform(shapeBY, st_crs(raster_temp))


#### Load daily data ####
files_rain_daily <- list.files('../dat/dwd_data/rain',
                               pattern='*.nc', 
                               full.names=TRUE)


list_rain_d <- lapply(files_rain_daily, stack)


files_Tmax_daily <- list.files('../dat/dwd_data/tmax',
                               pattern='*.nc', 
                               full.names=TRUE)


list_Tmax_d <- lapply(files_Tmax_daily, stack)


shapeBY <- st_transform(shapeBY, st_crs(list_rain_d[[1]]))



#### R20 ####
list_rain_gs <- lapply(list_rain_d, function(i) i[[tidyselect::vars_select(
 names(i),
 contains(growing_season))]])


list_rain_gs <- pblapply(list_rain_gs, function(i)  crop(i, shapeBY))


list_r20 <- pblapply(list_rain_gs, function(x) 
 calc(x, function(i) r20(i, 20)))

list_r20 <- pblapply(list_r20, function(i)  mask(i, shapeBY))

#--- convert to an sf object ---#
list_r20 <- pblapply(list_r20, function(i)
 st_as_sf(
  as(i, 'SpatialPolygonsDataFrame')
 )
)

# Make dataframe
df_r20 <- b %>% 
 purrr::map(dplyr::select, 1) %>% 
 purrr::reduce(st_join, join=st_equals)


#### Dry days ####
dd <- function(x,t){
 z <- sum(x < t)
 return(z)
}


list_dd <- pblapply(list_rain_gs, function(x) 
 calc(x, function(i) dd(i, 1)))


list_dd <- pblapply(list_dd, function(i)  mask(i, shapeBY))

#--- convert to an sf object ---#
list_dd <- pblapply(list_dd, function(i)
 st_as_sf(
  as(i, 'SpatialPolygonsDataFrame')
 )
)

# Make dataframe
df_dd <- list_dd %>% 
 purrr::map(dplyr::select, 1) %>% 
 purrr::reduce(st_join, join=st_equals)
colnames(df_dd) <- c(paste0("X", 2000:2019), "geometry") 



## 2.3 Yearly maximum temperature growing season
### Yearly max temp data at raster level
list_Tmax_gs <- lapply(list_Tmax_d, function(i) i[[tidyselect::vars_select(
 names(i),
 contains(growing_season))]])


list_Tmax_gs <- pblapply(list_Tmax_gs, function(i)  crop(i, shapeBY))


## formula count hot days
hd <- function(x,t){
 z <- sum(x > t)
 return(z)
}

## calculate hot days in ag season
list_hd <- pblapply(list_Tmax_gs, function(x) 
 calc(x, function(i) hd(i, 29)))


list_hd <- pblapply(list_hd, function(i)  mask(i, shapeBY))

#--- convert to an sf object ---#
list_hd <- pblapply(list_hd, function(i)
 st_as_sf(
  as(i, 'SpatialPolygonsDataFrame')
 )
)

# Make dataframe
df_hd <- list_hd %>% 
 purrr::map(dplyr::select, 1) %>% 
 purrr::reduce(st_join, join=st_equals)
colnames(df_hd) <- c(paste0("X", years), "geometry") 




ggplot() +
 geom_sf(data = df_hd[1:2,],aes(), col="red") +
 geom_sf(data = df_dd[1:20,],aes(), col="blue")+
 geom_sf(data = a,aes(), col="green")

a <- st_join(df_dd[1:20,], df_hd[1:5,], join = st_intersects, largest=T)

 
 
 