library(tidyverse)
library(stringr)
library(dplyr)
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
library(parallel)
library(terra)
library(gstat)
library(stars)

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
raster_rain <- stack(dataDWD(precip_index[[1]], base=gridbase, joinbf=TRUE, dir=locdir()))
raster_temp <- stack(dataDWD(temp_index[[1]], base=gridbase, joinbf=TRUE, dir=locdir()))

# Assign prjection 
# (https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_mean/DESCRIPTION_gridsgermany_monthly_air_temperature_mean_en.pdf)
crs(raster_rain) <- crs(raster_temp) <- "EPSG:31467"

## add choice experiment data to allocate observations in space
dat <- readRDS("../dat/choice_data/data_calc.rds")
shapeBY <- st_read("../dat/spatial_data/bayern_ex.shp")[,"SCH"]

## Match projections
shapeBY <- st_transform(shapeBY, st_crs(raster_rain))

## Delete unnecessary files to save memory
rm(gridIndex, precip_index, temp_index)



## 1.5 Shrink raster dataset
raster_temp <- crop(raster_temp, shapeBY)
raster_rain <- crop(raster_rain, shapeBY)

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
# Make parallel because slow
cl <- makeCluster(5)
clusterEvalQ(cl, {
 library(sf)
 library(pbapply)
 library(raster)
 })

clusterExport(cl, c("list_rain"))

list_rain <- pblapply(list_rain, function(i)
 st_as_sf(
  as(i, 'SpatialPolygonsDataFrame')
 ), cl = cl
)

stopCluster(cl)


# Make dataframe
df_rain <- list_rain %>% 
 purrr::map(dplyr::select, 1) %>% 
 purrr::reduce(st_join, join=st_equals)

colnames(df_rain) <- c(paste0("rain_", years), "geometry") 

# correct for wrong unit in rdwd package
df_rain[,1:10] <- st_drop_geometry(df_rain)*10

# remove all observations outside Bavaria (==NA)
df_rain <- filter(df_rain, rain_2010 !=0)

saveRDS(df_rain, '../dat/dwd_data/df_rain_dwd.rds')

rm(list_rain, raster_rain)




#### Other indicators: daily data ####
growing_season <- paste0(rep(years, each=8), ".", c(rep(0,7), 1) , c(3:9, 0), ".")

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
## formula count heavy rain days
r20 <- function(x,t){
 z <- sum(x > t)
 return(z)
}


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
df_r20 <- list_r20 %>% 
 purrr::map(dplyr::select, 1) %>% 
 purrr::reduce(st_join, join=st_equals)

colnames(df_r20) <- c(paste0("r20_", years), "geometry") 

saveRDS(df_r20, '../dat/dwd_data/df_r20_dwd.rds')

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
colnames(df_dd) <- c(paste0("dd_", years), "geometry") 

saveRDS(df_dd, '../dat/dwd_data/df_dd_dwd.rds')

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
colnames(df_hd) <- c(paste0("hd_", years), "geometry") 

saveRDS(df_hd, '../dat/dwd_data/df_hd_dwd.rds')


#### Combine data ####
df_rain <- st_transform(df_rain, st_crs(df_dd))
df_temp <- st_transform(df_temp, st_crs(df_dd))


ggplot() +
 geom_sf(data = df_hd[1:2,],aes(), col="red") +
 geom_sf(data = df_dd[1:35,],aes(), col="blue") +
 geom_sf(data = df_r20[1:35,],aes(), col="yellow") +
 geom_sf(data = df_temp[1:35,],aes(), col="green")



cl <- makeCluster(5)

clusterEvalQ(cl, {
 library(sf)
 library(pbapply)
})

clusterExport(cl, c("df_dd", "df_temp", "df_rain"))

df_temp <- pblapply(1:10, function(i) 
 st_interpolate_aw(df_temp[,i], df_dd, extensive=F), cl=cl) %>% 
 purrr::reduce(st_join, join=st_equals)

df_rain <- pblapply(1:10, function(i) 
 st_interpolate_aw(df_rain[,i], df_dd, extensive=F), cl=cl) %>% 
 purrr::reduce(st_join, join=st_equals)

stopCluster(cl)

ggplot() +
 geom_sf(data = df_hd[1:2,],aes(), col="red") +
 geom_sf(data = df_dd[1:35,],aes(), col="blue") +
 geom_sf(data = df_r20[1:35,],aes(), col="yellow") +
 geom_sf(data = df_temp[1:35,],aes(), col="green")



## Build dataframe
df_dwd <- st_join(df_dd, df_hd, join = st_intersects, largest=T)
df_dwd <- st_join(df_dwd, df_r20, join = st_equals)
df_dwd <- st_join(df_dwd, df_temp, join = st_equals)
df_dwd <- st_join(df_dwd, df_rain, join = st_equals)

## Filter out NAs that come from merging
filter_na <-with(df_dwd, which(is.na(rain_2010) |
                                is.na(dd_2010) |
                                is.na(hd_2010) |
                                is.na(r20_2010) |
                                is.na(temp_2010) ))

df_dwd <- df_dwd[-filter_na,]

saveRDS(df_dwd, '../dat/dwd_data/df_dwd_cleaned.rds')




ggplot() +
 geom_sf(data = df_dwd,aes(), col="red") 


#### Build relevant weather history dataframe ####

df_dwd <- df_dwd %>% 
 mutate(
  ## Precip sum
  rain_1_3 = rowMeans(select(st_drop_geometry(df_dwd), rain_2017:rain_2019)),
  rain_4_10 = rowMeans(select(st_drop_geometry(df_dwd), rain_2010:rain_2016)),
  ## Mean temp
  temp_1_3 = rowMeans( select(st_drop_geometry(df_dwd), temp_2017:temp_2019)),
  temp_4_10 = rowMeans(select(st_drop_geometry(df_dwd), temp_2010:temp_2016)),
  ## Dry days sum
  dd_1_3 = rowMeans( select(st_drop_geometry(df_dwd), dd_2017:dd_2019)),
  dd_4_10 = rowMeans( select(st_drop_geometry(df_dwd), dd_2010:dd_2016)),
  ## Heavy rain
  r20_1_3 = rowMeans( select(st_drop_geometry(df_dwd), r20_2017:r20_2019)),
  r20_4_10 = rowMeans( select(st_drop_geometry(df_dwd), r20_2010:r20_2016)),
  ## Hot days sum
  hd_1_3 = rowMeans( select(st_drop_geometry(df_dwd), hd_2017:hd_2019)),
  hd_4_10 = rowMeans( select(st_drop_geometry(df_dwd), hd_2010:hd_2016))
 )

#### Mean-center data ####
# The originally used data in Stetter, Sauer (2022) was mean-centered. 
# We use same means as in the original paper to make datasets compatible.

# Load original processed weather data based on zip codes
weather_stetter_sauer <-  readRDS("../dat/choice_data/fd_weather_plz.rds") 

colnames(weather_stetter_sauer) <- str_replace(
 colnames(weather_stetter_sauer), "ann_", "") 

old_names <- c("rain", "temp", "dd", "r20", "hd")

center_weather <- cbind(
 ## 1-3 lags
 weather_stetter_sauer %>% filter(year>=2017 & year<=2019) %>%
  aggregate(., by = list(.$plz), FUN = mean) %>% 
  dplyr::select(Group.1, rain:hd) %>%
  rename_at(vars(colnames(.)), ~ c("zip", paste0(old_names, "_1_3" ))),
 
 ## 4-10 lag
 weather_stetter_sauer %>% filter(year>=2010 & year<=2016) %>%
  aggregate(., by = list(.$plz), FUN = mean) %>% 
  dplyr::select(rain:hd) %>%
  rename_at(vars(colnames(.)), ~ c(paste0(old_names, "_4_10" )))
) %>% 
 dplyr::select(-zip) %>% 
 colMeans()

# re-order variables to align with created raster weather dataset
center_weather <- center_weather[c("rain_1_3", "rain_4_10",
                                   "temp_1_3",  "temp_4_10",
                                   "dd_1_3", "dd_4_10",
                                   "r20_1_3", "r20_4_10",
                                   "hd_1_3", "hd_4_10")] 

df_dwd_centered <- scale(st_drop_geometry(df_dwd) %>% select(rain_1_3:hd_4_10), 
                         center_weather, 
                         scale = FALSE)
df_dwd_centered <- as.data.frame(df_dwd_centered)

colnames(df_dwd_centered) <- paste0(colnames(df_dwd_centered), "_c")

# Save data, which will be test data for MAB
saveRDS(df_dwd_centered, '../dat/dwd_data/df_dwd_centered.rds')

saveRDS(cbind(df_dwd, df_dwd_centered), '../dat/dwd_data/df_dwd_final.rds')




#### Land cover ####
## Load dwd (test) dataset (from above)
df_grid <- readRDS('../dat/dwd_data/df_dwd_final.rds')


# Load land cover 10m grid
files_terra <- list.files('../dat/terrascope_10m',
                          pattern='*.tif', 
                          full.names=TRUE)

# Load
ic <- lapply(files_terra, rast)

# Align prjections
shape <- st_transform(df_grid, crs(ic[[1]]))

# Build single raster of the 6 raster files within Bavaria
ic <- sprc(lapply(ic, function(i) crop(i, shape)))
r <- mosaic(ic)

# Convert polygon to terra format
shape <- vect(shape)


# Make function to 1. see if 1x1km contains any cropland and 2. calculate share of cropland within grid 1x1km grid
calc_crop_share <- function(cells){ 
 x <- crop(r, shape[cells,])
 y <- ifel(x == 40, 1, NA)
 res1 <- zonal(cellSize(y,unit="ha"), y, sum)[[2]]
 res2 <- expanse(shape[cells,], unit="ha")
 sharec <- res1/res2
 sharec <- ifelse(sharec>1, 1, sharec)
 sharec <- ifelse(is_empty(sharec), 0, sharec)
 
 y <- ifel(x == 10 | x == 20, 1, NA)
 res1 <- zonal(cellSize(y,unit="ha"), y, sum)[[2]]
 res2 <- expanse(shape[cells,], unit="ha")
 sharef <- res1/res2
 sharef <- ifelse(sharef>1, 1, sharef)
 sharef <- ifelse(is_empty(sharef), 0, sharef)
 
 res <- data.frame(share_crop = sharec,
                   share_forest=sharef)
 return(res)
}


crop_tree_share <- 
 as.data.frame(t(
  pbsapply(1:nrow(df_grid), function(i)
   calc_crop_share( cells=i))
 ))

df_grid$share_crop <-  unlist(crop_tree_share$share_crop)
df_grid$share_forest <-  unlist(crop_tree_share$share_forest)



ggplot() +
geom_sf(data = df_grid, aes(fill=share_forest), col=NA)

saveRDS(df_grid, '../dat/dwd_data/df_dwd_final_lc.rds')



# 5 Calculate expected utility

## 5.1 Load choice model from Stetter, Sauer (2023)
# Load original processed weather data based on zip codes
res_stetter_sauer <-  readRDS("../dat/choice_data/res_cor_13_410.rds") 


## 5.2 Predicted utility when constant coefficients are assumed
# Extract model matrix
mm <- model.matrix(res_stetter_sauer)

# Extract mean coefficents
coef_stetter_sauer <- coef(res_stetter_sauer)[1:27]

# Return margin follows log distr. -> exponent
coef_stetter_sauer["DB"] <- exp(coef_stetter_sauer["DB"])

coefsubset <- intersect(names(coef_stetter_sauer), colnames(mm))

mm <- mm[, coefsubset]
utiltiy_homo <- 
 crossprod(
  t(mm[, coefsubset]), 
  coef_stetter_sauer[coefsubset]
 )


## 5.3 Predicted utility when farm-level coefficients are assumed
# Prepare indiv. coefficients for row-wise multiplication
ind_par <- as.matrix(t(purrr::map_dfr(seq_len(36), ~res_stetter_sauer$indpar)))

# Calculate  w/ heterogeneous coefficients
utility_hetero <- sapply(1:nrow(mm), function(i)
 mm[i,] %*%
  ind_par[-1,i]
)

# Normalize
utility_hetero_norm <- scales::rescale(utility_hetero, to=c(0,1))



#### Spatial WTA ####
#### Load shape data 
shape <- sf::st_read("../dat/plz_spatial/plz-3stellig.shp")

## subset only relevant zip codes
shape <- shape[startsWith(shape$plz, "63") | 
                startsWith(shape$plz, "8") | 
                startsWith(shape$plz, "9"),]

shape5 <- sf::st_read("../dat/plz_spatial/plz-gebiete.shp")
shape5 <- shape5[ which(shape5$plz %in% 
                         dat$location[which(nchar(dat$location)==5)]),]

shape$note <- shape5$note <- NULL
shapeAll <- rbind(shape, shape5)

shapeAll <- shapeAll %>% 
 group_by(plz) %>%
 summarise(geometry = sf::st_union(geometry)) %>%
 ungroup() %>% 
 rename(zip=plz)

# centroid
a <- st_centroid(shapeAll)


shapeBY <- st_read("../dat/spatial_data/bayern_ex.shp")[,"SCH"]
shapeBY <- st_transform(shapeBY, st_crs(a))


foo <- left_join(dat, a, by="zip") %>% 
 st_as_sf %>% 
 select(ID) %>% 
 rename(id=ID)

res_stetter_sauer <-  readRDS("../dat/choice_data/res_cor_13_410.rds") 


foo1 <- left_join(foo, res_stetter_sauer$indpar, by="id")

dd <- as.data.frame(res_stetter_sauer$model) %>% 
 select(rain_1_3_c:hd_4_10_c)
dd$id <- res_stetter_sauer$model$idx$ID

foo2 <- left_join(dd, foo1, by="id") %>% 
 st_as_sf


indpar_src <- foo2 %>% 
 select(ends_with("c.1")) %>% 
 st_drop_geometry

indpar_ac <- foo2 %>% 
 select(ends_with("c.2")) %>% 
 st_drop_geometry

var_weather <- foo2 %>% 
 select(rain_1_3_c:hd_4_10_c) %>% 
 st_drop_geometry

foo2$wta_src <- (foo2$X.Intercept..1 +
                  sapply(1:7128, function(i)
                   as.matrix(var_weather[i,]) %*% t(indpar_src[i,]))  ) /
 foo2$DB 

foo2$wta_ac <- (foo2$X.Intercept..2 +
                 sapply(1:7128, function(i)
                  as.matrix(var_weather[i,]) %*% t(indpar_src[i,])) ) /
 foo2$DB


a <- st_transform(a, st_crs(df_grid))
foo2 <- st_transform(foo2, st_crs(df_grid))
test <- stars::st_as_stars(df_grid$geometry)


dataCV <- foo2 %>% 
 group_by(id) %>% 
 summarise_at(vars(wta_src, wta_ac), list(mean))

wta_ac = st_as_sf(idw(wta_ac ~ 1, dataCV, test, idp = 2, nmax=198))
wta_src = st_as_sf(idw(wta_src ~ 1, dataCV, test, idp = 2, nmax=198))

wta_grid <- cbind(wta_ac[,1], wta_src[,1]) %>% 
 select(1:3) %>% 
 rename("wta_ac"=1, 
        "wta_src"=2)


#split data frame into n equal-sized data frames
list_small_grid <- split(df_grid, factor(sort(rank(row.names(df_grid))%%6000)))


list_final_wta <- pblapply(list_small_grid, function(i)
 st_interpolate_aw(st_crop(wta_grid, st_bbox(i)), i, extensive=F)
)

df_wta <- do.call(rbind, list_final_wta)

saveRDS(df_wta, '../dat/data_processed/df_wta.rds')


library("metR")
ggplot(pivot_longer(df_wta, 1:2)) +
 geom_sf(aes(fill=value), col=NA) +
 facet_grid(~name) +
 scale_fill_binned(type="viridis",
                   breaks=my_breaks(c(min(df_wta$wta_src), max(df_wta$wta_src)))) +
 theme_void()


#### 6 Build training dataframe ####
df_train <- res_stetter_sauer$model  %>%  
 as.data.frame %>% 
 dplyr::select(DB:hd_4_10_c) %>% 
 dplyr::mutate(reward_homo=utiltiy_homo,
               reward_hetero=utility_hetero,
               reward_hetero_norm=utility_hetero_norm,
               decision_arm=rep(c("SRC", "AC", "SQ"),nrow(mm)/3),
               wta_ac = foo2$wta_ac,
               wta_src = foo2$wta_ac)



## Landscape weighted norm utility
# Crop by (inverse of) crop share or 1-share
# SRC and AC by (inverse of) forest share
# df_train$reward_hetero_norm_inverse <- 
#  with(df_train,ifelse(decision_arm=="SQ",
#                       reward_hetero_norm*1/ls_shares_zip$share_crop,
#                       reward_hetero_norm*1/ls_shares_zip$share_forest))
# 
# 
# df_train$reward_hetero_norm_share <- 
#  with(df_train,ifelse(decision_arm=="SQ",
#                       reward_hetero_norm*(1-ls_shares_zip$share_crop),
#                       reward_hetero_norm*(1-ls_shares_zip$share_forest)))

write.csv(df_train, file = "../dat/data_processed/df_train", row.names = F)


# 7 Build test dataframe
df_grid <- readRDS('../dat/dwd_data/df_dwd_final_lc.rds')
df_wta <- readRDS('../dat/data_processed/df_wta.rds')


df_test <- cbind(df_grid, df_wta) %>% 
 filter(df_grid$share_crop!=0) %>% 
 dplyr::select(contains("_c"), -share_crop, wta_ac) %>% 
 st_drop_geometry() 

write.csv(df_test, file = "../dat/data_processed/df_test", row.names = F)
