library(spdep)

coef_mean <- coef(res_stetter_sauer)[1:27]

coef_src <- subset(coef_mean, str_detect(names(coef_mean), "_c:1"))
coef_ac <- subset(coef_mean, str_detect(names(coef_mean), "_c:2"))

dnearneigh(df_grid[1:10,], d1=0, d2=0.75*max_1nn, row.names=IDs)


map <- df_grid[1:50,]
map$ID <- 1:50
point <- st_centroid(df_grid$geometry[1:2])
circle <- st_buffer(point, dist = 5000)

ggplot() +
 geom_sf(data = map$ge) +
 geom_sf(data = point) +
 geom_sf(data = circle, fill=NA) +
 geom_sf_text(data= map, aes(label = ID),
              nudge_x = -300, nudge_y = 300) +
 theme_void()

# Check out units of measurement of projection
metacapa::projection_units(df_grid) # m

# distance neighbors
neighbor_list <- 
 include.self(
  dnearneigh(
   st_centroid(df_grid$geometry), 
   d1=0, 
   d2=15500
  )
 )

# Show cell 1 neighbors
map$nn1 <- as.character(ifelse(map$ID %in% neighbor_list[[1]], 1, NA))
map$nn2 <- as.character(ifelse(map$ID %in% neighbor_list[[2]], 2, NA))
map$nn3 <- as.character(ifelse(map$ID %in% neighbor_list[[3]], 3, NA))
map$nn4 <- as.character(ifelse(map$ID %in% neighbor_list[[4]], 4, NA))
map$nn5 <- as.character(ifelse(map$ID %in% neighbor_list[[5]], 5, NA))

map_long <- pivot_longer(map, c(nn1:nn5)) %>% filter(!is.na(value))

ggplot() +
 geom_sf(data = map_long, aes(fill=value)) +
 geom_sf(data = point) +
 geom_sf(data = circle, fill=NA) +
 geom_sf_text(data= map, aes(label = ID),
              nudge_x = -300, nudge_y = 300) +
 facet_wrap(~value, ncol=2) +
 theme_bw()

a <- pbsapply(1:nrow(df_test), function(i)
df_test[neighbor_list[[i]],] %>% colMeans(na.rm = T))

df_test_buffer15 <- as.data.frame(t(a))



df_grid <- readRDS('../dat/dwd_data/df_dwd_final_lc.rds')

lll <- df_grid %>% 
 mutate(wta_src=as.matrix(df_test) %*% coef_src / exp(coef_mean["DB"])) %>% 
 ggplot() +
 geom_sf(aes(fill=pred_land_use), col = 0) +
 theme_void(base_size = 30) +
 scale_fill_viridis_b(option="magma", breaks = (c(0, 200)))



ggplot() +
 geom_sf(data = weather_stetter_sauer$) 
