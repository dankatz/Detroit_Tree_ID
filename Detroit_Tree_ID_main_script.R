#Detroit tree ID: main script for manuscript

### data set descriptions, preprocessing, etc ######################################
## 2017 Lidar ======================================================================
# 2017 LiDAR was provided by Everett Root, slightly ahead of release schedule
# It is now publicly available by request, details available here:
# https://semcog.org/Portals/0/ElevationQuickFacts.pdf
# The files are large, so I have it stored on an external hard drive 
# E:\lidar\2017 LiDAR
# preprocessing was originally done in this script: CHM_script180628.R 
# lots of additional things are stored there

# Covert .LAS to .LAZ to speed up computation --------------------------------------
library(lidR)
library(rlas)
library(raster)

working_directory <- "D:/lidar/2017 LiDAR/Classified_LAS"
setwd("D:/lidar/2017 LiDAR/Classified_LAS")
las_list = list.files(getwd(), pattern=".las$", full.names=FALSE)
for(i in 1:707){writelax(las_list[i])} #write lax files (index of las files) to speed computation

# Create a CHM----------------------------------------------------------------------
th = seq(from = 0, to = 100, by = 5)
kernel = matrix(1,3,3)

myAnalyze <- function(las){
  chm = grid_tincanopy(las, 2, thresholds = th, max_edge = c(0,5), subcircle = 1.0)
  chm_rast = as.raster(chm)
  chm_rast_s = raster::focal(chm_rast, w = kernel, fun = median)
  chm_rast_s = raster::focal(chm_rast_s, w = kernel, fun = median)
  
  extent_orig <- extent(chm_rast_s)
  extent_new <- c((extent_orig@xmin + 10),
                  (extent_orig@xmax - 10),
                  (extent_orig@ymin + 10),
                  (extent_orig@ymax - 10))
  chm_rast_sc <- crop(chm_rast_s, extent_new)
  
  
  writeRaster(chm_rast_sc, filename = paste("SDM_", xmin(chm_rast_sc), "_", 
                                            ymin(chm_rast_sc), ".tif", sep =""), overwrite = TRUE)
}

ctg = catalog("D:/lidar/2017 LiDAR/Classified_LAS")
plot(ctg)
buffer(ctg) <- 10
cores(ctg) <- 4L
by_file(ctg) <- TRUE
#tiling_size(project) <- 2500

catalog_apply(ctg, myAnalyze, select = "xyzr") #takes ~55 hours

# create a DEV ---------------------------------------------------------------------
setwd("D:/lidar/2017 LiDAR/CHM_test")

myAnalyze <- function(las){
  tile_DEV <- grid_terrain(las, res = 2, method = "delaunay")
  DEV <- as.raster(tile_DEV)
  writeRaster(DEV, filename = paste("DEV_", xmin(tile_DEV), "_", 
                                    ymin(tile_DEV), ".tif", sep =""), overwrite = TRUE)
}

ctg = catalog("D:/lidar/2017 LiDAR/CHM_test")
plot(ctg)
buffer(ctg) <- 10
cores(ctg) <- 4L
by_file(ctg) <- TRUE

#ctg2 <- ctg[1:2]

catalog_apply(ctg, myAnalyze) #takes ~X hours

# create derived variables for tree ID ---------------------------------------------

## Street tree database ===========================================================
# The street tree database was provided by MI DNR, Davey, and Detroit
# The MOU and downloaded data are stored here:
# C:\Users\dsk856\Box\MIpostdoc\Detroit spatial data\Detroit Street Tree Inventory feb 2016
# E:\other Detroit data\street trees

#script used to filter out  trees from street tree dataset that are unsuitable for inclusion

setwd("D:/tree_classification")
rm(list = ls()) 
library(httr)
library(raster)
library(jpeg)
library(dplyr)
library(sf) #install.packages("sf")

#read in street trees
t <- st_read("streettrees_to_filter.shp")

head(t)
unique(t$genus)

f <- filter(t, trees09_17 == 1, DBH > 4)
f$taxa <- as.character(f$genus)
f$taxa[f$taxa == "Thuja"] <- "Cupressaceae"
f$taxa[f$taxa == "Juniperus"] <- "Cupressaceae"
f$taxa[f$SPP == "Acer platanoides"] <- "Acer_platanoides"

f <- filter(f, taxa != "stump")
f <- filter(f, nSDM_2017_ > 16.4 & ht_m_09 > 5)
#unique(f$taxa)

f$dbhsub20_heightover40 <- NA
f$dbhsub20_heightover40[f$nSDM_2017_ > 40 & f$DBH < 20] <- "overtopped"
f$dbhsub20_heightover40[f$ht_m_09 > 12 & f$DBH < 20] <- "overtopped"

f2 <- filter(f, is.na(dbhsub20_heightover40))
# ggplot(f2, aes(x = DBH, y = nSDM_2017_)) + theme_bw() + facet_wrap(~taxa) + geom_point(alpha = 0.1)+
#   #geom_quantile() + 
#   geom_density_2d() + 
#   geom_smooth(col = "blue", method = "lm")

#group tree taxa with sparse data together as "other"
trees_per_taxa <- group_by(f2, taxa) %>% summarize(n_trees = n()) 
trees_per_taxa$geometry <- NULL
f2 <- left_join(f2, trees_per_taxa)

#f <-  filter(f, n_trees > 100 | taxa == "Cupressaceae" | taxa == "Magnolia")
f2$taxa[f2$n_trees < 100] <- "other"
f2$taxa500 <- f2$taxa
f2$taxa500[f2$n_trees < 500] <- "other" #group_by(f2, taxa500) %>% summarize(n_trees = n()) 


#unique(f2$taxa)
write_sf(f2, "street_trees_filtered180326.shp")


#unique(f2$SPP[f2$taxa == "Acer"])

### Table 1: street tree database in Detroit ######################################################
library(readr)
library(dplyr)
library(forcats)
st_trees <- read_csv("C:/Users/dsk856/Box/MIpostdoc/Detroit spatial data/Detroit Street Tree Inventory feb 2016/street_trees_with_coordinates_only.csv")
# hist(st_trees$DBH)
# names(st_trees)

st_trees <- st_trees %>% 
              mutate(genus = gsub( " .*$", "", SPP),
                     group = genus,
                     group = replace(x = group, SPP == "Acer platanoides", "Acer platanoides")) %>%
              filter(genus != "vacant",
                     genus != "stump") %>%
              #case_when(SPP == "Acer platanoides" ~ group = "Acer platanoides") %>%
              mutate(BA = (DBH ^2) * 0.005454154) #%>% #convert DBH in inches to BA in square ft
              
total_ba <- sum(st_trees$BA)
table_1 <- st_trees %>% group_by(group) %>%
              summarize(BA_sum = sum(BA),
                        BA_rel = BA_sum / total_ba,
                        n_trees = n()) %>%
              mutate(group = replace(x = group, BA_rel < 0.002, "other")) %>%
              mutate(group = replace(x = group, group == "unknown", "other")) %>%
              mutate(group = replace(x = group, group == "Malus", "other")) %>% #looks like it was just on the cutoff point; 
              #and wasn't included in the analysis. If I re-do the analysis drop the infrequent genera including:
              #Ailanthus, Celtis, Ginkgo, Morus
              group_by(group) %>%
              summarize(BA_sum = sum(BA_sum),
                        BA_rel = sum(BA_rel),
                        n_trees = sum(n_trees)) 
table_1 %>% summarize(BA_rel_total = sum(BA_rel), #totals for table
                      n_trees_total = sum(n_trees))

#what are the common species in a genus?
st_trees %>% filter(genus == "Ginkgo") %>% group_by(SPP) %>% summarize(n = n())



## World View 2 image ==============================================================
# Three WV2 images were acquired from Digital Globe via a data grant from 
# UM's Clark Library. Contact people were Mara Blake and, later, Nicole Scholtz
# we received a few folders:
# DigitalGlobeImagery_v1 contains imagery at level 3D which IS orthorectified  
# DigitalGlobeImagery_v2 contains imagery at level 2A which IS NOT orthorectified
# Shannon was in charge of atmospheric correction; to do so he used Erdas Imagine
# I manually made a cloud mask, it's here: "C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/WV2_cloud_cover.shp"

library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(RStoolbox)
library(randomForest)
library(sp)
library(stringr)

#RMSE compared to LiDAR
# lidar_wv2_rmse <- read.csv("C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/imagery_rmse/rmse_lidar_wv2_200426.txt",
#                            stringsAsFactors = FALSE, sep = '\t')
# dist_1 <- pointDistance(cbind(lidar_wv2_rmse$X1, lidar_wv2_rmse$Y1), cbind(lidar_wv2_rmse$X2, lidar_wv2_rmse$Y2), lonlat=FALSE)
# (sqrt(mean(dist_1^2))) #rmse

# WV2 spectral index creation -----------------------------------------------------
#originally from this script: wv2_derived_spectral_indices_creation190206.R
setwd("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics")
#wp <- stack("mosaic11jun13164022pan.img")
wm <- stack("mosaic11jun13164023mul_reflec.img")

setwd("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics/derived_spectral_indices")

wm_derived_ndvi75 <- (wm$mosaic11jun13164023mul_reflec.7 - wm$mosaic11jun13164023mul_reflec.5)/
  (wm$mosaic11jun13164023mul_reflec.7 + wm$mosaic11jun13164023mul_reflec.5)
writeRaster(wm_derived_ndvi75, "wm_derived_ndvi75.tif")

wm_derived_ndvi76 <- (wm$mosaic11jun13164023mul_reflec.7 - wm$mosaic11jun13164023mul_reflec.6)/
  (wm$mosaic11jun13164023mul_reflec.7 + wm$mosaic11jun13164023mul_reflec.6)
writeRaster(wm_derived_ndvi76, "wm_derived_ndvi76.tif")

wm_derived_ndvi86 <- (wm$mosaic11jun13164023mul_reflec.8 - wm$mosaic11jun13164023mul_reflec.6)/
  (wm$mosaic11jun13164023mul_reflec.8 + wm$mosaic11jun13164023mul_reflec.6)
writeRaster(wm_derived_ndvi86, "wm_derived_ndvi86.tif")

wm_derived_ndvi85 <- (wm$mosaic11jun13164023mul_reflec.8 - wm$mosaic11jun13164023mul_reflec.5)/
  (wm$mosaic11jun13164023mul_reflec.8 + wm$mosaic11jun13164023mul_reflec.5)
writeRaster(wm_derived_ndvi85, "wm_derived_ndvi85.tif")

wm_derived_ndvi84 <- (wm$mosaic11jun13164023mul_reflec.8 - wm$mosaic11jun13164023mul_reflec.4)/
  (wm$mosaic11jun13164023mul_reflec.8 + wm$mosaic11jun13164023mul_reflec.4)
writeRaster(wm_derived_ndvi84, "wm_derived_ndvi84.tif")

wm_derived_ndvi61 <- (wm$mosaic11jun13164023mul_reflec.6 - wm$mosaic11jun13164023mul_reflec.1)/
  (wm$mosaic11jun13164023mul_reflec.6 + wm$mosaic11jun13164023mul_reflec.1)
writeRaster(wm_derived_ndvi61, "wm_derived_ndvi61.tif")

wm_derived_ndvi65 <- (wm$mosaic11jun13164023mul_reflec.6 - wm$mosaic11jun13164023mul_reflec.5)/
  (wm$mosaic11jun13164023mul_reflec.6 + wm$mosaic11jun13164023mul_reflec.5)
writeRaster(wm_derived_ndvi65, "wm_derived_ndvi65.tif")

wm_derived_845 <- wm$mosaic11jun13164023mul_reflec.8 /
  (wm$mosaic11jun13164023mul_reflec.4 + wm$mosaic11jun13164023mul_reflec.5)
writeRaster(wm_derived_845, "wm_derived_ndvi845.tif")

wm_derived_cari <- (wm$mosaic11jun13164023mul_reflec.6 - wm$mosaic11jun13164023mul_reflec.5) -
  ((0.2 * (wm$mosaic11jun13164023mul_reflec.6 - wm$mosaic11jun13164023mul_reflec.3)))
writeRaster(wm_derived_cari, "wm_derived_cari.tif")



# extract multispectral imagery and the derived spectral indices --------------------------------
library(velox)
## load WV2 imagery (will clip trees to middle swath) ===================================
setwd("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics")
wm <- stack("mosaic11jun13164023mul_reflec.img")

## load in tree canopy objects (that are considered street trees) =======================
setwd("D:/tree_classification/segmentation/")
t <- st_read("tree_can_obj_190326_d.shp") #polygons created from tree segmentation in R
t2 <- st_crop(t, wm) #only select polygons within the datastack area
t2 <- dplyr::filter(t2, area < 500 & area > 4)

### extracting values from existing rasters of spectral indices###
setwd("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics/derived_spectral_indices")
wv_list <- grep(".tif", str_sub(dir(), start = -4)) 
wv_list <- dir()[wv_list]
wv_list_length <- length(wv_list)

for(r in 1: (wv_list_length + 1)){
  
  #derived indices
  if(r < wv_list_length + 1){
    wv_rast <- stack(wv_list[r])  #wv_rast <- stack(wv_list[1])
    #wv_rast <- crop(wv_rast, wm)
    wv_rast_vx <- velox(wv_rast) #takes a couple min to run, maybe 6 GB memory
    wv_rast_names <- names(wv_rast)
  }
  
  #raw multispectral bands
  if(r == wv_list_length + 1){
    wv_rast <- wm
    wv_rast_vx <- velox(wv_rast) #takes a couple min to run, maybe 6 GB memory
    wv_rast_names <- names(wv_rast)
  }
  
  start_time <- Sys.time()
  
  for(i in 1:nrow(t2)){ #loop through all of the objects
    t2_small <- t2[i,]
    obj_data_vx <- wv_rast_vx$copy() #extent(obj_data_vx)
    obj_data_vx$crop(t2_small) #extent(obj)
    
    bands_d_mean <- obj_data_vx$extract(sp = t2_small, fun=mean)
    bands_d_mean_df <- as.data.frame(bands_d_mean)
    names(bands_d_mean_df) <- paste0(wv_rast_names, "_d_mn")
    bands_d_sd <- obj_data_vx$extract(t2_small, fun=sd)
    bands_d_sd_df <- as.data.frame(bands_d_sd)
    names(bands_d_sd_df) <- paste0(wv_rast_names, "_d_sd")
    bands_d_median <- obj_data_vx$extract(t2_small, fun=median)
    bands_d_median_df <- as.data.frame(bands_d_median)
    names(bands_d_median_df) <- paste0(wv_rast_names, "_d_md")
    dp_focal_poly <- cbind(bands_d_sd_df, bands_d_mean_df, bands_d_median_df)
    if( i == 1){dp <- dp_focal_poly}
    if( i > 1){dp <- bind_rows(dp, dp_focal_poly)}
  } #end band extraction loop
  
  print(Sys.time() - start_time)
  
  if( r == 1){dp_results <- bind_cols(t2, dp)}
  if( r > 1){dp_results <- bind_cols(dp_results, dp)}
  
}#end raster loop

#save results
#dp_results_backup <- dp_results #names(dp_results) #dp_results <- dp_results_backup

#saving the removal of NAs for the next stage
#dp_results[, 34:84][is.na(dp_results[, 34:84])] <- 0 #throwing an error message, but still seems to work
#dp_results$wm_derived_ndvi61_d_md[is.na(dp_results$wm_derived_ndvi61_d_md)]
dp_results2 <- dp_results[,c(5, 15, 26, 34:84)]
dp_results2$geometry <- NULL
write.csv(dp_results2, "wv_spect_ind_190327.csv")



## NearMap imagery =================================================================
# NearMap data was downloaded (Clark Library had a subscription)
# data and original scripts are stored here:
# C:\Users\dsk856\Box\MIpostdoc\Detroit spatial data\Nearmap
# rmse for a nearmap image compared to the 2017 LiDAR (with 50 control points)
# lidar_nearmap_rmse <- read.csv("C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/imagery_rmse/rmse_lidar_nearmap181025.txt",
#                            stringsAsFactors = FALSE, sep = '\t')
# dist_1 <- pointDistance(cbind(lidar_nearmap_rmse$X1, lidar_nearmap_rmse$Y1), 
#                         cbind(lidar_nearmap_rmse$X2, lidar_nearmap_rmse$Y2), lonlat=FALSE)
# (sqrt(mean(dist_1^2))) #rmse: 1.04 m

## figuring out which tiles are in Detroit ------------------------------------------
#the purpose of this script was to determine which Nearmap tiles are within Detroit
#this is the 4th version of this script, it now uses a 200 m buffer around Detroit
#for use with the script nearmap_api_downloadscript.R
#rm(list = ls()) 


#which tiles are in the general region for zoom level 17, written 
#regional x and y coordinates, gotten from nearmap
# lower left: x = 35208  y =48540, #top right: x= 35370, y = 48425
x_min <- 35208
x_max <- 35370
y_min <- 48425 
y_max <- 48540


#creating a matrix of all tile x coordinates within the region
D_region_x_tile <- matrix(nrow = x_max - x_min, ncol =  y_max - y_min)
for(x in x_min:x_max){
  for(y in y_min:y_max){
    D_region_x_tile[(x - x_min), (y - y_min)] <- x
  }
}

#creating a matrix of all tile y coordinates
D_region_y_tile <- matrix(nrow = x_max - x_min, ncol =  y_max - y_min)
for(x in x_min:x_max){
  for(y in y_min:y_max){
    D_region_y_tile[(x - x_min), (y - y_min)] <- y
  }
}

#creating a dataframe of all tile coordinates
D_region_points <- data.frame(as.vector(D_region_x_tile), as.vector(D_region_y_tile))
names(D_region_points) <- c("x_tile", "y_tile")


#loading and converting Detroit boundary in to tile coordinate system
library('ggplot2')
zoom <- 17
setwd("C:/Users/dsk586/Box/MIpostdoc/Detroit spatial data/Nearmap")
D_outline_points <- read.csv('Detroit_boundary_points_sorted_200mbuffer.csv')  #previously 'Detroit_boundary_points_sorted.csv'
plot(D_outline_points$long, D_outline_points$lat)

#using equations found at the url below to convert between long/lat and tile coordinates: 
#https://groups.google.com/forum/#!topic/Google-Maps-API/NICY9wcl_JY

#convert Detroit outline longitude in to google tile coordinates
D_outline_points$x_temp = round (256*(2^(zoom-1)))+(D_outline_points$lon*((256*(2^zoom))/360))
D_outline_points$x_tile <- floor(D_outline_points$x_temp/256)
# Absolute X to tile X (where x = absolute x pixel coord; x2 = relative
#                       x pixel (on tile); x3 = the tile's x coord):

#convert Detroit outline latitude in to google tile coordinates
D_outline_points$y_temp <- sin((D_outline_points$lat*pi)/180)
D_outline_points$y_temp = round (256*(2^(zoom-1)))+
  ((.5*log((1+D_outline_points$y_temp)/(1-D_outline_points$y_temp)))*((-256*(2^zoom))/(2*pi)))
D_outline_points$y_tile <- floor(D_outline_points$y_temp / 256)

D_outline_points_longlat <- data.frame(D_outline_points$x_tile, D_outline_points$y_tile)
names(D_outline_points_longlat) <- c("x_tile", "y_tile")

plot(D_outline_points_longlat$x_tile, D_outline_points_longlat$y_tile)

#checking which points in the region are within Detroit
library('mgcv')
D_region_points$arepointsin <- in.out(bnd = as.matrix(D_outline_points_longlat), x = as.matrix(D_region_points)) #test whether points are in or not

ggplot(D_region_points, aes(x= x_tile, y = y_tile, color = arepointsin)) + geom_point()

D_tiles_within_boundary <- subset(D_region_points, arepointsin == TRUE)
D_tiles_within_boundary <- D_tiles_within_boundary[,1:2] 

plot(D_tiles_within_boundary)

#write final results to a csv for use with the nearmap download API
write.csv(D_tiles_within_boundary, "coords_tiles_within_Detroit_z17.csv")


#which tiles are in the general region for zoom level 18
#regional x and y coordinates, gotten from nearmap
# lower left: x = 70375  y =97102, #top right: x= 70724, y = 96855 #lat long == -82.88, lat == 42.46

#creating a matrix of all tile x coordinates within the region
D_region_x_tile <- matrix(nrow = 70724 - 70375, ncol =  97102 - 96855)
for(x in 70375:70724){
  for(y in 96855:97102){
    D_region_x_tile[(x - 70375), (y -96855)] <- x
  }
}

#creating a matrix of all tile y coordinates
D_region_y_tile <- matrix(nrow = 70724 - 70375, ncol =  97102 - 96855)
for(x in 70375:70724){
  for(y in 96855:97102){
    D_region_y_tile[(x - 70375), (y -96855)] <- y
  }
}

#creating a dataframe of all tile coordinates
D_region_points <- data.frame(as.vector(D_region_x_tile), as.vector(D_region_y_tile))
names(D_region_points) <- c("x_tile", "y_tile")

#loading and converting Detroit boundary in to tile coordinate system
zoom <- 18
D_outline_points <- read.csv('Detroit_boundary_points_sorted_200mbuffer.csv')  #previously 'Detroit_boundary_points_sorted.csv'
#plot(D_outline_points$long, D_outline_points$lat)

#using equations found at the url below to convert between long/lat and tile coordinates: 
#https://groups.google.com/forum/#!topic/Google-Maps-API/NICY9wcl_JY

#convert Detroit outline longitude in to google tile coordinates
D_outline_points$x_temp = round (256*(2^(zoom-1)))+(D_outline_points$lon*((256*(2^zoom))/360))
D_outline_points$x_tile <- floor(D_outline_points$x_temp/256)
# Absolute X to tile X (where x = absolute x pixel coord; x2 = relative
#                       x pixel (on tile); x3 = the tile's x coord):

#convert Detroit outline latitude in to google tile coordinates
D_outline_points$y_temp <- sin((D_outline_points$lat*pi)/180)
D_outline_points$y_temp = round (256*(2^(zoom-1)))+
  ((.5*log((1+D_outline_points$y_temp)/(1-D_outline_points$y_temp)))*((-256*(2^zoom))/(2*pi)))
D_outline_points$y_tile <- floor(D_outline_points$y_temp / 256)

D_outline_points_longlat <- data.frame(D_outline_points$x_tile, D_outline_points$y_tile)
names(D_outline_points_longlat) <- c("x_tile", "y_tile")

plot(D_outline_points_longlat$x_tile, D_outline_points_longlat$y_tile)

#checking which points in the region are within Detroit
library('mgcv')
D_region_points$arepointsin <- in.out(bnd = as.matrix(D_outline_points_longlat), x = as.matrix(D_region_points)) #test whether points are in or not
ggplot(D_region_points, aes(x= x_tile, y = y_tile, color = arepointsin)) + geom_point()

D_tiles_within_boundary <- subset(D_region_points, arepointsin == TRUE)
D_tiles_within_boundary <- D_tiles_within_boundary[,1:2] 
plot(D_tiles_within_boundary)

#write final results to a csv for use with the nearmap download API
write.csv(D_tiles_within_boundary, "coords_tiles_within_Detroit.csv")



## downloading tiles ----------------------------------------------------------------
# this script is available here: nearmap_api_downloadscript.R
# the basic steps are: download each tile, convert .JPEG to .TIF with the proper extent, stich the image together

## creating derived spectral indices for Nearmap ------------------------------------

rm(list = ls()) 
gc()

library(sf)
library(raster)
library(dplyr)
library(RStoolbox)
library(randomForest)
library(sp)
library(stringr)
library(velox)


## load WV2 imagery (will clip trees to middle swath)
setwd("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics")
wm <- stack("mosaic11jun13164023mul_reflec.img")

# load in tree canopy objects (that are considered street trees) 
setwd("D:/tree_classification/segmentation/")
t <- st_read("tree_can_obj_190326_d.shp") #polygons created from tree segmentation in R
t2 <- st_crop(t, wm) #only select polygons within the datastack area
t2 <- dplyr::filter(t2, area < 500 & area > 4)
t2_paste <- dplyr::select(t2, UnqID, SPP, DBH, geometry)
#head(t2)

# extracting values from existing rasters of spectral indices
setwd("D:/aerial imagery/nearmap/derived_spectral_indices")
nearmap_list <- grep(".tif", str_sub(dir(), start = -4)) 
nearmap_list <- dir()[nearmap_list]

for(r in 1: 56){
  near_rast <- stack(nearmap_list[r])  #near_rast <- stack(nearmap_list[1]) #plot(near_rast)
  near_rast <- crop(near_rast, wm)
  near_rast_vx <- velox(near_rast) #takes a couple min to run, maybe 6 GB memory
  near_rast_names <- names(near_rast)
  
  start_time <- Sys.time()
  
  for(i in 1:nrow(t2)){ #loop through all of the objects
    t2_small <- t2[i,]
    obj_data_vx <- near_rast_vx$copy() #extent(obj_data_vx)
    obj_data_vx$crop(t2_small) #extent(obj)
    
    bands_d_mean <- obj_data_vx$extract(sp = t2_small, fun=mean)
    bands_d_mean_df <- as.data.frame(bands_d_mean)
    names(bands_d_mean_df) <- paste0(near_rast_names, "_d_mn")
    bands_d_sd <- obj_data_vx$extract(t2_small, fun=sd)
    bands_d_sd_df <- as.data.frame(bands_d_sd)
    names(bands_d_sd_df) <- paste0(near_rast_names, "_d_sd")
    bands_d_median <- obj_data_vx$extract(t2_small, fun=median)
    bands_d_median_df <- as.data.frame(bands_d_median)
    names(bands_d_median_df) <- paste0(near_rast_names, "_d_md")
    dp_focal_poly <- cbind(bands_d_sd_df, bands_d_mean_df, bands_d_median_df)
    if( i == 1){dp <- dp_focal_poly}
    if( i > 1){dp <- bind_rows(dp, dp_focal_poly)}
  } #end band extraction loop
  
  print(Sys.time() - start_time)
  
  if( r == 1){dp_results <- bind_cols(t2_paste, dp)}
  if( r > 1){dp_results <- bind_cols(dp_results, dp)}
  
}#end raster loop
dp_results$geometry <- NULL
#save results #getwd()
setwd("D:/aerial imagery/nearmap/derived_spectral_indices") 
write.csv(dp_results, "nearmap_spect_ind_190327.csv", row.names = FALSE)



### segmentation of tree canopies ##################################################
library(lidR)
library(rlas)
library(rgdal)

## segmenting full raster (SDM tiles) for Detroit ==================================
library(stringr)
tile_rast <- raster("D:/lidar/2017 LiDAR/nSDM_2017_UTM17.tif")
#split the raster in arcgis (it was easier than in R)
setwd("D:/tree_classification/segmentation/nSDM2017_UTM17_tiles")
all_files_subfiles <- dir() #all files, including subfiles
all_files <- all_files_subfiles[str_sub(all_files_subfiles, -4, -1) == ".TIF"] #only select .TIF

#loop to load file
current_max_treeID <- 0
for(i in 1:length(all_files)){
  #for(i in 1318:1324){ #tile 1319 isn't working
  focal_rast <- raster(all_files[i])
  
  #dealing with tiles that don't have data
  if(maxValue(focal_rast) > 5){ #there seem to be some empty tiles that have values greater than 0  #plot(focal_rast)
    focal_rast <- focal_rast + 0.01 #makes sure that the file is read as double (required by tree_detection)
    f = function(x) { x * 0.02 + 4 } #function for variable moving window size
    ttops <- tree_detection(focal_rast, algorithm = lmf(ws = f, hmin = 10, shape = "circular"))
    crowns = dalponte2016(focal_rast, ttops, th_tree = 10, th_seed = 0.45, th_cr = 0.5, max_cr = 40)()
    crowns <- crowns + current_max_treeID #making sure that treeID is unique across tiles
    
    tile_number <- as.numeric(gsub("([0-9]+).*$", "\\1", str_sub(all_files[i], start = 17, end = 20)))
    file_name <- paste("treeID_tile_", tile_number,".tif", sep = "")
    file_dir_name <- paste("../nSDM2017_UTM17_seg_rast/",file_name, sep = "")
    writeRaster(crowns, file_dir_name, overwrite = TRUE)
    current_max_treeID <- maxValue(crowns)
  }
}#results in ~2,500,000 polygons 

#rasters combined in arcgis using mosaic to new raster D_treeID_190113.shp
#this raster was too big for arcgis to do raster to polygon, so I broke it up into 9 pieces
#D:\tree_classification\segmentation\nSDM2017_UTM17_seg_rast_chunks
#these chunks were then converted to shapefile using raster to polygon
#output: D:\tree_classification\segmentation\nSDM2017_UTM17_seg_rast_chunks\rast_to_polygons

#merge doesn't seem to be working in arcgis (too big?) so I'm going to try selecting out the relevant polygons in R


##intersect polygons with filtered street tree data for each segmentation chunk ========================
library(sf)
library(stringr)

#load tree locations
setwd("D:/tree_classification")
t <- st_read("street_trees_filtered180326.shp") #filtering criteria are in script: "filtering_Davey_tree_points_to_create_training.R"
#plot(t)

#load in a chunk of segmentation for Detroit
setwd("D:/tree_classification/segmentation/nSDM2017_UTM17_seg_rast_chunks/rast_to_polygons")
dir()
all_files_subfiles <- dir() #all files, including subfiles
all_files <- all_files_subfiles[str_sub(all_files_subfiles, -4, -1) == ".shp"] #only select .TIF

#loop through each chunk
st_crs(t) <- st_crs(seg) #not quite sure why this had to be done, but it seems to work okay

for(i in 1:4){ #this seems to choke because of bad geometry in chunk 4
  #for(i in 5:9){  #the work-around is to run it first for chunks 1 -3, manually fix it with the section below
  #and then run it for chunks 5-9
  seg <- st_read(all_files[i])
  t_chunk <- st_crop(t, seg)
  
  #use this to fix the bad geometry in chunk 4
  #plot(t_chunk)
  #bad features in chunk 4:
  # seg_validity <- st_is_valid(seg)  #chunk 4
  # seg$seg_validity <- seg_validity
  # seg <- filter(seg, seg_validity == TRUE)
  # seg$seg_validity <- NULL
  
  tree_poly_points_all <- st_join(y = t_chunk, x = seg) # matches for those circles covering a point
  tree_poly_points <- tree_poly_points_all[!is.na(tree_poly_points_all$DBH),] #getting rid of polygons that didn't have a tree inside
  #plot(st_geometry(test2))
  tree_poly_points$area <- as.numeric(st_area(tree_poly_points))
  tree_poly_points <- dplyr::filter(tree_poly_points, area < 1500) #filter out non-trees (trees won't be bigger than ~1000m2)
  
  #append to the previous chunks
  if(i == 1){tree_poly_points_chunks <- tree_poly_points}else{
    tree_poly_points_chunks <- rbind(tree_poly_points_chunks, tree_poly_points)}
}

plot(st_geometry(tree_poly_points_chunks))

#remove polygons that have the same geometry (i.e., 2 points within same polygon)
unique_geometry <- tree_poly_points_chunks %>% dplyr::select(geometry) %>% duplicated()
unique_geometry_rev <- tree_poly_points_chunks %>% dplyr::select(geometry) %>% duplicated(fromLast = TRUE)
tree_poly_points_chunks$unique_geometry_forward <- unique_geometry
tree_poly_points_chunks$unique_geometry_backward <- unique_geometry_rev
tree_poly_points_chunks_unique <- dplyr::filter(tree_poly_points_chunks, 
                                                unique_geometry_forward == FALSE & #
                                                  unique_geometry_backward == FALSE) #need to go both directions

setwd("D:/tree_classification/segmentation/")
st_write(tree_poly_points_chunks_unique, "tree_can_obj_190326_d.shp")

### assembling all predictor variables #############################################
#Object based identification of trees in Detroit
#setwd("//ITS-AD-DFS02.adsroot.itcs.umich.edu/MWS/dept/seas/Ibanez Lab/Dan Katz/POLLEN POSTDOC/Detroit spatial data/Nearmap")

library(sf)
library(raster)
library(dplyr)
#library(RStoolbox)
library(randomForest)
library(sp)

rm(list = ls()) 
#gc()

### assemble object dataframe ################################################################
#load in the middle WV2 swath for analysis
setwd("E:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics")
wm <- stack("mosaic11jun13164023mul_reflec.img")

#load in canopy polygons 
setwd("E:/tree_classification/segmentation/")
t <- st_read("tree_can_obj_190326_d.shp") #polygons created from tree segmentation in R
t2 <- st_crop(t, wm) #only select polygons within the datastack area
t2 <- dplyr::filter(t2, area < 500 & area > 4)

#load in the various derived variables
#nearmap
#created in: "nearmap_object_vars_190320.R"
setwd("D:/aerial imagery/nearmap/derived_spectral_indices") #getwd()
nearmap_spect <- read.csv("nearmap_spect_ind_190327.csv", stringsAsFactors = FALSE) #31678 * 2

## wv2
#created in: "wv2_derived_spectral_indices_creation.R"
setwd("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics/derived_spectral_indices")
wv2_spect <- read.csv("wv_spect_ind_190327.csv", stringsAsFactors = FALSE) #includes lots of NAs

#created in "wv2_texture_indices_object.R"
setwd("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics/derived_textural_indices")
wv2_text <- read.csv("wv_text_ind_190329.csv", stringsAsFactors = FALSE)


## lidar 2017
setwd("D:/lidar/2017 LiDAR")
lid17 <- read.csv("lid17_ind_190412.csv") #indices for raw lidar, intensity, nSDM - DTM
lid17$X <- NULL
lid17_tex <- read.csv("lid17_tex_190415.csv") #lidar texture
lid17_tex$X <- NULL
names(lid17_tex) <- sub("wp_", "lid17_", names(lid17_tex)) #quick fix for column names that weren't renamed correctly in texture script
#lid17_inten <- read.csv("") #lidar intensity
lid17_tex_inten <- read.csv("lid17_tex_inten_190412.csv") #lidar intensity texture
lid17_tex_inten$X <- NULL
names(lid17_tex_inten) <- sub("wp_", "lid17_inten_", names(lid17_tex_inten)) #quick fix for column names that weren't renamed correctly in texture script


## join the variables into one large df
t3 <- left_join(t2, wv2_spect)
t3 <- left_join(t3, wv2_text)
t3 <- left_join(t3, nearmap_spect)
t3 <- left_join(t3, lid17)
t3 <- left_join(t3, lid17_tex)
t3 <- left_join(t3, lid17_tex_inten)

## load in variables based on the independent variable csvs that were created for the predictions -----------
# #work around for re-doing a section of the analysis
# library(readr)
# library(purrr)
# 
# setwd("E:/tree_classification/predictions/csv_for_pred_vars")
# csvs_for_pred <- dir()#[1:3]
# 
# compile_csvs <- function(file_name){
#   read_csv(file_name, col_types = cols(Id = col_character(), grdcd = col_character())) %>%
#   select(-c("X1", "grond", "tree", "build", "NA_px", "X", "nSDM_2017_", "nSDM_201_1", "nSDM_201_2", "ply__"))
# }
# 
# all_csvs <- map_dfc(csvs_for_pred, compile_csvs) #this is a little messy since the Id and grdcd columns are duplicated.
# # double check that all of the grdcd columns are the same (no offsets)
# # select(all_csvs, contains("grdcd")) %>% sample_n(10000) %>% t(.) %>% as.data.frame %>% distinct() 
# # this check works by making sure that there aren't differences in a sample of the grdcd values
# 
# #remove the duplicate grdcd and Id columns
# cols_to_remove_grdcd <- names(all_csvs) %>% grep(pattern = "grdcd", x = ., value = TRUE)
# cols_to_remove_grdcd <- cols_to_remove_grdcd[2:length(cols_to_remove_grdcd)]
# cols_to_remove_Id <- names(all_csvs) %>% grep(pattern = "Id", x = ., value = TRUE)
# cols_to_remove_Id <- cols_to_remove_Id[2:length(cols_to_remove_Id)]
# cols_to_remove <- c(cols_to_remove_grdcd, cols_to_remove_Id)
# all_csvs2 <- all_csvs %>% select(-cols_to_remove) %>%
#             mutate(grdcd = as.numeric(grdcd),
#                    Id = as.numeric(Id))


unique_polys <- unique(t2$grdcd)
all_pred_vars_for_analysis <- filter(all_pred_vars, grdcd %in% unique_polys) 
# write_csv(all_pred_vars_for_analysis, 
#   "C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/all_data_for_analysis200502.csv")

# join the independent variables with the analysis dataset
t3 <- left_join(t2, all_pred_vars_for_analysis)


#t3 <- left_join(t2, all_csvs2)
#write_csv(t3, "C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/polygons_with_imagery_for_analysis.csv")

#remove the polygons that are under the cloud mask
clouds <- read_sf("C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/WV2_cloud_cover.shp") 
clouds <- st_transform(clouds, crs(t3)) #the crs *seems* to be exactly the same between them, but somewhere it isn't
cloud2 <- clouds %>% mutate(cloud = "cloud") %>% dplyr::select(cloud)

t3_cloud <- st_join(t3, cloud2)
t3_cloud <- filter(t3_cloud, is.na(cloud)) %>%
              dplyr::select(-cloud)

library(forcats)
t3_nogeo <- t3_cloud
st_geometry(t3_nogeo) <- NULL #speeds up processing
t3_nogeo <- t3_nogeo %>% mutate(taxa = fct_recode(taxa, conifer = "Picea", conifer = "Pinus")) %>%
            mutate(taxa2 = fct_lump_min(taxa, min = 100, other_level = "other")) #only including taxa that have more than 100 polygons in focal area
         
# count_taxa <- t3_nogeo %>% group_by(taxa2) %>% summarize(count_taxa = n())  #number of polygons from each taxon
# print(count_taxa, n = nrow(count_taxa))
#t3 <- left_join(t3, count_taxa, by = "taxa")
#t4 <- dplyr::filter(t3, count_taxa > 100) #only including taxa that have more than 100 polygons in focal area

#t3$taxa2 <- t3$taxa
#t3$taxa2[t3$count_taxa < 100] <- "other"

#str(t3$taxa2)
t4 <- t3_nogeo
#t5 <- sample_n(t4, 10000)

### classification #################################################################
set.seed(100)
dp <- t4
dp$geometry <- NULL
dp$count_taxa <- NULL
dp2 <- dp[!is.na(dp$nearmap_181025_gei_d_sd),]

#names(which(sapply(dp2, anyNA)))
sapply(dp2, function(y) sum(length(which(is.na(y))))) #count number of NAs in each column

dp2 <- dp2 %>% replace(., is.na(.), 0)
### TRY THIS OUT
#try missing forest algorithm in order to impute missing values

#dp2 <- dp2[Reduce('&', lapply(dp2, is.finite)),] #remove any NAs/ INF values

dp2$taxa <- as.character(dp2$taxa)
# dp2$taxa[dp2$taxa == "Aesculus" | dp2$taxa == "Ailanthus" |dp2$taxa == "Catalpa"| 
#            dp2$taxa == "Pyrus"] <- "other"

#making a completely balanced dataset
# dp3 <- dp2
# #dp3$taxa <- as.character(dp3$taxa)
# dp3$taxa <- as.factor(as.character(dp3$taxa)) 
# group_by(dp3, taxa) %>% summarize(count = n())


# # #original single taxon version
# group_by(dp2, taxa) %>% summarize(count = n())
# focal_taxon <- "Quercus"
#dp3 <- dp2
#dp3$taxa <- as.character(dp3$taxa)
# dp3$taxa[dp3$taxa != focal_taxon] <- "other"
# dp3$taxa <- as.factor(as.character(dp3$taxa))
#dp3 <- group_by(dp3, taxa) %>% sample_n(n_focal)
#dp3$taxa <- as.factor(as.character(dp3$taxa))
dp3 <- dp2 #names(dp2)
dp3$train0_test1 <- rbinom(n = nrow(dp3) , size = 1, prob = 0.3) #?rbinom
dp3_for_modeling <- dp3[, c(33:215, 218:336)] #names(dp3)
names(dp3_for_modeling)
dp3_for_modeling$taxa <- as.factor(as.character(dp3_for_modeling$taxa2)) ##TrainSet2 <- cbind(TrainSet2, TrainSet$taxa)
dp3_for_modeling$taxa2 <- NULL
#names(dp3_for_modeling)

#TrainSet2 <- TrainSet[, c(34:374)] #MAKE SURE ITS ONLY REMOTE SENSING VARS LEFT HERE

TrainSet <- filter(dp3_for_modeling,  train0_test1== 0)
TrainSet %>% group_by(taxa) %>% summarize(n_taxa = n())  %>% print(n=40)

ValidSet <- filter(dp3_for_modeling, train0_test1 == 1)
ValidSet %>% group_by(taxa) %>% summarize(n_taxa = n()) %>% print(n=40)




#n_focal <- filter(TrainSet, taxa == focal_taxon) %>% NROW()
#TrainSet %>% group_by(taxa) %>% summarize(n_taxa = n())
#TrainSet <- group_by(TrainSet, taxa) %>% sample_n(n_focal)

#ValidSet %>% group_by(taxa) %>% summarize(n_taxa = n())
#ValidSet <- group_by(ValidSet, taxa) %>% sample_n(200, replace = TRUE)

#train <- sample(nrow(dp3), 0.7*nrow(dp3), replace = FALSE)
#TrainSet <- dp3[train,] #length(unique(TrainSet$taxa)) #levels(TrainSet$taxa)
#TrainSet %>% group_by(taxa) %>% summarize(n_taxa = n())
#ValidSet <- dp3[-train,] #length(unique(ValidSet$taxa))

#  TrainSet <- group_by(TrainSet, taxa) %>% sample_n(800, replace = TRUE)




# Split into Train and Validation sets
# Training Set : Validation Set = 70 : 30 (random)
# set.seed(100)
# train <- sample(nrow(dp3), 0.7*nrow(dp3), replace = FALSE)
# TrainSet <- dp3[train,] #length(unique(TrainSet$taxa)) #levels(TrainSet$taxa)

#ValidSet <- dp3[-train,] #length(unique(ValidSet$taxa))

#TrainSet[is.infinite(TrainSet)]

#str(TrainSet)
#names(TrainSet)
#TrainSet <- dp2[train,] #length(unique(TrainSet$taxa)) #levels(TrainSet$taxa)
# TrainSet2 <- TrainSet[, c(34:374)] #MAKE SURE ITS ONLY REMOTE SENSING VARS LEFT HERE
#TrainSet2 <- dplyr::select(TrainSet, contains("nearmap_"))

# TrainSet2$taxa <- as.factor(as.character(TrainSet$taxa2)) ##TrainSet2 <- cbind(TrainSet2, TrainSet$taxa)
# TrainSet2$taxa2 <- NULL
# names(TrainSet2)

#adding in Boruta to select only the useful variables: spoiler- all variables were listed as important
# library("Boruta")
# boruta_train <- Boruta(taxa~., data = TrainSet, doTrace = 2) #takes 5 hours
# print(boruta_train)
# attributes_to_include <- getSelectedAttributes(boruta_train)
# plot(boruta_train)
# boruta_df <- attStats(boruta_train)
# setwd("//ITS-AD-DFS02.adsroot.itcs.umich.edu/MWS/dept/seas/Ibanez Lab/Dan Katz/POLLEN POSTDOC/trees/tree_identificaiton")
# write.csv(boruta_df, "boruta_df_190418.csv") #getwd()
# print(bank_df)
# plotImpHistory(boruta_train)

start_time <- Sys.time() #takes about an hour 
#TrainSet2 <- TrainSet2[,c(1:318,340:341)]#names(TrainSet2) #not including certain variable classifications in model
model1 <- randomForest(taxa ~ ., data = TrainSet, importance = TRUE)
model1
print(Sys.time() - start_time)
#plot(model1)

# model2 <- randomForest(taxa ~ ., data = TrainSet, ntree = 500, mtry = 6, importance = TRUE)
# model2

# Predicting on train set
#predTrain <- predict(model1, TrainSet, type = "class")
#str(predTrain)
# Checking classification accuracy
#table(predTrain, TrainSet$taxa)  



### model assessment ###############################################################
# Predicting on Validation set
#ValidSet$taxa2 <- as.character(ValidSet$taxa) #ValidSet$taxa <- NULL #ValidSet$taxa <- ValidSet$taxa2
predValid <- predict(model1, ValidSet, type = "class")
#?predict
# Checking classification accuracy
#mean(predValid == ValidSet$taxa)                    
# test <- (as.matrix(table(predValid,ValidSet$taxa)))
# test

mean(predValid == ValidSet$taxa) 
confusion_df <- as.data.frame.matrix(table(predValid,ValidSet$taxa)) #rows = predicted, col = actual
confusion_df$row_sum <- rowSums(confusion_df)
confusion_df[ (nrow(confusion_df) +1) ,] <- colSums(confusion_df)

#user accuracy
confusion_matrix <- as.data.frame.matrix(table(predValid,ValidSet$taxa))
confusion_df$user_accuracy <- c(unlist(diag(as.matrix(confusion_matrix)))/confusion_df$row_sum[1:18], NA)
confusion_df$user_accuracy <- round(confusion_df$user_accuracy, 2)

#producer accuracy
producer_accuracy <- c(unlist(diag(as.matrix(confusion_matrix)))/confusion_df[19,])
producer_accuracy <- round(unlist(producer_accuracy), 2)
confusion_df <- rbind(confusion_df, producer_accuracy)
#View(confusion_df)

model_importance_df <- as.data.frame(importance(model1))        
View(model_importance_df) #names(model_importance_df)
varImpPlot(model1) 
#hist(model_importance_df$MeanDecreaseGini)
#str(model_importance_df)

save(model1, file = "C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/d_rf200503.RData")




### check the contributions of each dataset ###################################################
setwd("C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/dataset_selection")

#list of all datasets to do separately
trainset_names <- sort(names(TrainSet))
names_wv2 <- grep(pattern = "wm_derived|mosaic11|taxa|train0", x = trainset_names, value = TRUE)
names_wv2_raw <- grep(pattern = "mosaic11|taxa|train0", x = trainset_names, value = TRUE)
names_wv2_indices <- grep(pattern = "wm_derived|taxa|train0", x = trainset_names, value = TRUE)
names_nearmap <- grep(pattern = "nearmap|taxa|train0", x = trainset_names, value = TRUE)
names_lidar <- grep(pattern = "SDM_DTM|nSDM_2017|lidar2017|taxa|train0", x = trainset_names, value = TRUE)
names_lidar_int <- grep(pattern = "lidar2017|taxa|train0", x = trainset_names, value = TRUE)
names_lidar_nearmap <- grep(pattern = "nearmap|SDM_DTM|nSDM_2017|lidar2017|taxa|train0", x = trainset_names, value = TRUE)
names_lidar_wv2 <-  grep(pattern = "mosaic11|wm_derived|SDM_DTM|nSDM_2017|lidar2017|taxa|train0", x = trainset_names, value = TRUE)
names_all<- names(TrainSet)
# names_wv2_raw <- names(TrainSet)[c(28:51, 341)]
# names_wv2_indices <- names(TrainSet)[c(1:27, 341)]
# names_wv2_tex <- names(TrainSet)[c(52:72, 341)]
# names_nearmap <- names(TrainSet)[c(73:288, 341)]
# names_lidar <- names(TrainSet)[c(289:339, 341)]
# names_lidar_int <- names(TrainSet)[c(292:294, 341)]
# names_lidar_tex <- names(TrainSet)[c(298:318, 341)]
# names_lidar_int_tex <- names(TrainSet)[c(319:339, 341)]
# names_lidar_nearmap <- names(TrainSet)[c(73:339, 341)]
# names_lidar_wv2 <- names(TrainSet)[c(1:72, 289:339, 341)]
# names_all <- names(TrainSet)[c(1:339, 341)]
# names_all_no_tex <- names(TrainSet)[c(1:51, 73:297, 341)]

list_models <- list("names_wv2" = names_wv2, 
                    "names_wv2_raw" = names_wv2_raw, 
                    "names_wv2_indices" = names_wv2_indices,
                    #"names_wv2_tex" = names_wv2_tex,
                    "names_nearmap" = names_nearmap,
                    "names_lidar" = names_lidar,
                    "names_lidar_int" = names_lidar_int,
                    #"names_lidar_tex" = names_lidar_tex,
                    #"names_lidar_int_tex" = names_lidar_int_tex,
                    "names_lidar_nearmap" = names_lidar_nearmap,
                    "names_lidar_wv2" = names_lidar_wv2,
                    "names_all" = names_all
                    #"names_all_no_tex" = names_all_no_tex 
)

for(i in 1:9){#start loop#for(i in 1:12){#start loop  #getwd()
  start_time <- Sys.time() #takes about an hour 
  focal_model <- list_models[i]
  focal_model_name <- names(focal_model)
  print(focal_model_name)
  
  focal_model_data <- dplyr::select(TrainSet, unlist(focal_model))
  names(focal_model_data) <- unlist(focal_model) #names_wv2 #fixing the naming issue by going back to the original
  #names(focal_model_data)
  
  #run model
  model1 <- randomForest(taxa ~ ., data = focal_model_data, importance = TRUE)
  predValid <- predict(model1, ValidSet, type = "class")
  print(mean(predValid == ValidSet$taxa))
  
  #confusion matrix
  confusion_df <- as.data.frame.matrix(table(predValid,ValidSet$taxa)) #rows = predicted, col = actual
  confusion_df$row_sum <- rowSums(confusion_df) #View(confusion_df)
  confusion_df[ (nrow(confusion_df) +1) ,] <- colSums(confusion_df)
  confusion_df_file_name <- paste0("confus_df_", focal_model_name, ".csv")
  write.csv(confusion_df, file = confusion_df_file_name, row.names = TRUE) #?write.csv
  
  #model importance df
  model_importance_df <- as.data.frame(randomForest::importance(model1))   #View(model_importance_df) 
  confusion_df_file_name <- paste0("model_import_", focal_model_name, ".csv")
  write.csv(model_importance_df, confusion_df_file_name)
  print(Sys.time() - start_time)
}#end loop

## create a table from the results ===============================
library("psych")
setwd("C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/dataset_selection")
confus_file_list <- dir()[grep("confus_df", dir())]
n_models <- length(confus_file_list)

#create dataframe to hold results

models_results <- data.frame(model_name = rep(NA, n_models), 
                             model_accuracy = rep(NA, n_models),
                             model_user_acc = rep(NA, n_models),
                             model_prod_acc = rep(NA, n_models),
                             kappa_stat = rep(NA, n_models))

#loop through all files
for(i in 1:length(confus_file_list)){
  #i <- 1
  #load confusion matrix
  confusion_df <- read.csv(confus_file_list[i])
  confusion_df_small <- confusion_df[1:(nrow(confusion_df ) - 1), 1:(ncol(confusion_df) - 1)]
  confusion_df_small$X <- NULL
  confusion_df_small <- as.matrix(confusion_df_small)
  
  #model_name
  models_results$model_name[i] <- sub("confus_df_names_", "", confus_file_list[i])
  models_results$model_name[i] <- sub(".csv", "", models_results$model_name[i])
  
  #total accuracy
  models_results$model_accuracy[i] <- sum(diag(confusion_df_small))/sum(confusion_df_small)
  
  #users accuracy
  #diag(confusion_df_small) / rowSums(confusion_df_small)
  model_user_acc <- mean( diag(confusion_df_small) / rowSums(confusion_df_small), na.rm = TRUE)
  models_results$model_user_acc[i] <- model_user_acc
  
  #producers accuracy
  #diag(confusion_df_small) / colSums(confusion_df_small)
  model_prod_acc <- mean( diag(confusion_df_small) / colSums(confusion_df_small), na.rm = TRUE)
  models_results$model_prod_acc[i] <- model_prod_acc
  
  #kappa
  kappa_stat <- cohen.kappa(confusion_df_small) #use the estimate for unweighted kappa
  models_results$kappa_stat[i] <- kappa_stat$kappa 
  
}

#save model diagnostic table
models_results
library(stargazer)
stargazer(models_results, type = "html", summary = FALSE, rownames = FALSE, out = "summary_table.doc")

#save confusion matrix table for all 
confus_matrix_all <- read.csv("confus_df_names_all.csv")
stargazer(confus_matrix_all, type = "html", summary = FALSE, rownames = FALSE, out = "confusion_matrix_table.doc")

#save importance values table
importance_vals_all <- read.csv("model_import_names_all.csv")
importance_vals_all <- arrange(importance_vals_all, by = -MeanDecreaseAccuracy)
importance_vals_all <- rename(importance_vals_all, "variable" = "X")
importance_vals_all_trunc <- importance_vals_all[1:50, ]
stargazer(importance_vals_all_trunc, type = "html", summary = FALSE, rownames = FALSE, out = "SI_variable_importance_table.doc")

### prediction #####################################################################
#predict tree identity across all of Detroit 
#the model this is based on was run in "tree_classification_obj_D_190419.R"
# this was originally in the script: "tree_classification_predictions_D190807.R"

#set up work environment
library(sf)
library(raster)
library(tidyr)
library(dplyr)
library(velox)
library(stringr)
rm(list = ls()) 


## assemble dataset for predictions ============================================
# load in ALL objects from segmentation 
#load in a chunk of segmentation for Detroit: created in "tree_segmentation190326.R"
setwd("D:/tree_classification/segmentation/nSDM2017_UTM17_seg_rast_chunks/rast_to_polygons")
#dir()
all_files_subfiles <- dir() #all files, including subfiles
all_files <- all_files_subfiles[str_sub(all_files_subfiles, -4, -1) == ".shp"] 

# lidar 
setwd("D:/lidar/2017 LiDAR")
lid17 <- raster("nSDM_2017_UTM17.tif") #lid17 <- NULL
tbg <- raster("trees_ground_build_v2_UTM17_v2_within_D.tif")

# load WV2 imagery (will clip trees to middle swath) 
setwd("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics")
wm <- stack("mosaic11jun13164023mul_reflec.img") # wm <- NULL


# load in tree canopy objects (that are considered street trees) 
setwd("D:/tree_classification/segmentation/")
t <- st_read("tree_can_obj_190326_d.shp") #polygons created from tree segmentation in R

## start loop to go through all of the objects in chunks 

## extract lidar nSDM for all objects in that chunk of the shapefile ==========================
setwd("D:/tree_classification/segmentation/nSDM2017_UTM17_seg_rast_chunks/rast_to_polygons")
for(j in 1:1){ #probly less than 12 hours to run for all chunks
  start_time <- Sys.time()
  
  seg <- st_read(all_files[j])
  #names(seg)
  #seg[1,1]
  
  # seg2 <- st_crop(seg, wm) #only select polygons within the datastack area
  # t2 <- dplyr::filter(t2, area < 500 & area > 4)
  
  ### extracting values from existing rasters of spectral indices###
  #derived indices
  lid_rast <- lid17  #plug in each of the lidar rasters here #wv_rast <- stack(wv_list[1])
  lid_rast <- crop(lid17, wm)
  lid_rast_vx <- velox(lid_rast) #takes a couple min to run, maybe 6 GB memory
  lid_rast_names <- names(lid_rast)
  
  #extract lidar data for each object
  for(i in 1:nrow(seg)){ #loop through all of the objects#nrow(seg)){ #loop through all of the objects
    seg_small <- seg[i,]
    obj_data_vx <- lid_rast_vx$copy() #extent(obj_data_vx)
    obj_data_vx$crop(seg_small) #extent(obj)
    
    bands_d_mean <- obj_data_vx$extract(sp = seg_small, fun=mean)
    bands_d_mean_df <- as.data.frame(bands_d_mean)
    names(bands_d_mean_df) <- paste0(lid_rast_names, "_d_mn")
    bands_d_sd <- obj_data_vx$extract(seg_small, fun=sd)
    bands_d_sd_df <- as.data.frame(bands_d_sd)
    names(bands_d_sd_df) <- paste0(lid_rast_names, "_d_sd")
    bands_d_median <- obj_data_vx$extract(seg_small, fun=median)
    bands_d_median_df <- as.data.frame(bands_d_median)
    names(bands_d_median_df) <- paste0(lid_rast_names, "_d_md")
    dp_focal_poly <- cbind(bands_d_sd_df, bands_d_mean_df, bands_d_median_df)
    if( i == 1){dp <- dp_focal_poly}
    if( i > 1){dp <- bind_rows(dp, dp_focal_poly)}
  } #end band extraction loop #40 min for 158000 obs.  
  
  #save results from object chunk
  # if(j == 1){dp_results <- bind_cols(seg, dp)}
  # if(j > 1){
  dp_results_chunkj <- bind_cols(seg, dp)
  #dp_results <- bind_rows(dp_results, dp_results_chunkj) 
  #}
  
  #save results
  object_chunk_lidar_filename <- paste0("object_lidar_data_chunk_",j,".csv")
  dp_results_save <- dp_results_chunkj
  dp_results_save$geometry <- NULL
  write.csv(dp_results_save, object_chunk_lidar_filename)
  print(Sys.time() - start_time)
} # end object chunk loop



# extract whether each object is classified as tree, building, ground 
#for all objects in that chunk of the shapefile 

setwd("D:/tree_classification/segmentation/nSDM2017_UTM17_seg_rast_chunks/rast_to_polygons")
for(j in 1:9){ #~5 hours to run for all chunks
  start_time <- Sys.time()
  
  seg <- st_read(all_files[j]) #seg <- st_read(all_files[5])
  #names(seg)
  #seg[1,1]
  
  # seg2 <- st_crop(seg, wm) #only select polygons within the datastack area
  # t2 <- dplyr::filter(t2, area < 500 & area > 4)
  
  ### extracting values from existing rasters of spectral indices###
  #derived indices
  tbg_rast <- tbg  #plug in each of the tbg rasters here #wv_rast <- stack(wv_list[1])
  tbg_rast <- crop(tbg, wm)
  tbg_rast_vx <- velox(tbg_rast) #takes a couple min to run, maybe 6 GB memory
  tbg_rast_names <- names(tbg_rast)
  
  #extract tbg data for each object
  for(i in 1:nrow(seg)){ #loop through all of the objects#nrow(seg)){ #loop through all of the objects
    seg_small <- seg[i,] #seg_small <- seg[13987, ] #413829
    obj_data_vx <- tbg_rast_vx$copy() #extent(obj_data_vx)
    obj_data_vx$crop(seg_small) #extent(seg_small) #plot(seg_small)
    #plot(tbg_rast, xlim =c(327249, 327400), ylim =c(4700800, 4700930))
    #plot(seg_small, add = TRUE)
    
    tbg_vals <- obj_data_vx$extract(sp = seg_small, fun= NULL)
    tbg_vals_unlist <- unlist(tbg_vals)
    
    NA_pixels <- as.data.frame(0)
    if(length(na.omit(tbg_vals_unlist)) == 0){
      NA_pixels <- as.data.frame(length(is.na(tbg_vals_unlist)))
    }
    names(NA_pixels) <- "NA_pixels"
    ground_pixels <- as.data.frame(length(tbg_vals_unlist[tbg_vals_unlist == 1]))
    names(ground_pixels) <- "ground"
    tree_pixels <- as.data.frame(length(tbg_vals_unlist[tbg_vals_unlist == 2]))
    names(tree_pixels) <- "tree"
    build_pixels <- as.data.frame(length(tbg_vals_unlist[tbg_vals_unlist == 3]))
    names(build_pixels) <- "build"
    dp_focal_poly <- cbind(ground_pixels, tree_pixels, build_pixels, NA_pixels)
    if( i == 1){dp <- dp_focal_poly}
    if( i > 1){dp <- bind_rows(dp, dp_focal_poly)}
  } #end band extraction loop #40 min for 158000 obs.  
  
  #save results from object chunk
  # if(j == 1){dp_results <- bind_cols(seg, dp)}
  # if(j > 1){
  dp_results_chunkj <- bind_cols(seg, dp)
  #dp_results <- bind_rows(dp_results, dp_results_chunkj) 
  #}
  
  #save results
  object_chunk_tbg_filename <- paste0("object_tbg_data_chunk_",j,".csv")
  dp_results_save <- dp_results_chunkj
  dp_results_save$geometry <- NULL
  write.csv(dp_results_save, object_chunk_tbg_filename)
  print(Sys.time() - start_time)
} # end object chunk loop

# # #load in a completed chunk and check it out manually in arcgis
# seg <- st_read(all_files[5])
# test <- read.csv("object_tbg_data_chunk_5.csv")
# test2 <- left_join(seg, test)
# test2$tree2 <- 0
# test2$tree2[test2$tree > test2$build & 
#               test2$tree > test2$ground & 
#               test2$tree > test2$NA_pixels ] <- "1"
# write_sf(test2, "test6.shp")



### load all of the processed csv chunks, combine, and filter for height and tree ################
setwd("D:/tree_classification/segmentation/nSDM2017_UTM17_seg_rast_chunks/rast_to_polygons")
#dir()
all_files_subfiles <- dir() #all files, including subfiles
all_files <- all_files_subfiles[str_sub(all_files_subfiles, -4, -1) == ".csv"] 
tbg_files <- all_files[grep("tbg", all_files)]
lidar_files <- all_files[grep("lidar", all_files)]
seg_files <- all_files_subfiles[str_sub(all_files_subfiles, -4, -1) == ".shp"] 

all_t_seg <- NULL
for(i in 1:9){
  tbg_chunk <- read.csv(tbg_files[i]) #names(test)
  tbg_chunk$X <- NULL
  lidar_chunk <- read.csv(lidar_files[i])
  tbglid_chunk <- left_join(tbg_chunk, lidar_chunk,  by = c("Id", "gridcode"))
  tbglid_chunk_f <- filter(tbglid_chunk, tree > build & 
                             tree > ground &
                             tree > NA_pixels)
  tbglid_chunk_f <- filter(tbglid_chunk_f, nSDM_2017_UTM17_d_md > 15 & nSDM_2017_UTM17_d_md < 150)
  
  if(nrow(tbglid_chunk_f) > 1){ #protect against the chunks that don't overlap with this wv2 area
    tbglid_chunk_f$poly_is_tree <- "tree"
    seg_chunk <- st_read(seg_files[i])
    seg_chunk_rast <- left_join(seg_chunk, tbglid_chunk_f)
    seg_chunk_rast_f <- filter(seg_chunk_rast, poly_is_tree == "tree")
    
    if(exists("all_t_seg") == FALSE){all_t_seg <- seg_chunk_rast_f}
    if(exists("all_t_seg") == TRUE){all_t_seg <- rbind(all_t_seg, seg_chunk_rast_f)}
  }
} #hist(all_t_seg$nSDM_2017_UTM17_d_md, breaks = 100)
setwd("D:/tree_classification/segmentation")
write_sf(all_t_seg, "all_tree_polygons_inD190508.shp") #save results





# load in rasters for all predictor vars 
#clear out workspace
rm(list = ls()) 
gc() #velox is pretty memory hungry, need to manually clear space for it
setwd("D:/tree_classification/predictions/csv_for_pred_vars")

wm <- stack("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics/mosaic11jun13164023mul_reflec.img") # wm <- NULL
all_t_seg <- st_read("D:/tree_classification/segmentation/all_tree_polygons_inD190508.shp") 

#load in the required files 
# lidar
lid_list <- c("D:/lidar/2017 LiDAR/nSDM_2017_UTM17.tif",
              "D:/lidar/2017 LiDAR/lidar2017_intensity_utm.tif",
              "D:/lidar/2017 LiDAR/SDM_DTM_2017_UTM17.tif")

# wv2 
wv2_list <- dir("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics/derived_spectral_indices",
                full.names = TRUE)
wv2_list <- wv2_list[str_sub(wv2_list, -4, -1) == ".tif"] 
wv2_list_raw <- dir("D:/satellite imagery/DigitalGlobe_2011_2012/DigitalGlobe_mosaics/raw_bands", 
                    full.names = TRUE)
wv2_list_raw <- wv2_list_raw[str_sub(wv2_list_raw, -4, -1) == ".tif"] 
wv2_list_full <- c(wv2_list, wv2_list_raw)

# nearmap 
near_list <- dir("D:/aerial imagery/nearmap/derived_spectral_indices",
                 full.names = TRUE)
near_list <- near_list[str_sub(near_list, -4, -1) == ".tif"] 
#remove hsv rasters (they're just the first raster from the stack of 3)
near_list_nohsv <- near_list[grep(pattern = "hsv", x = near_list, invert = TRUE)]

#add hsv individual rasters
near_hsv_list <- dir("D:/aerial imagery/nearmap/derived_spectral_indices/nearmap_hsv/", 
                     full.names = TRUE)
near_list_full <- c(near_list_nohsv, near_hsv_list)

#all rasters #note that some are as bands within a raster stack #nearmap hsv
all_rasters <- c(lid_list, wv2_list_full, near_list_full)
#grep(vars_to_extract_median,all_rasters)
#all_rasters_md <- all_rasters[2]
# 

#THIS IS IN A SECTION BELOW, NEED TO RUN IT FIRST
vars_to_extract_median2 <- str_sub(vars_to_extract_median, 1, -8)
# grep(c("gei","grvi"), all_rasters)

#double checking that all rasters were extracted
existing_files <- dir("D:/tree_classification/predictions/csv_for_pred_vars/")
existing_files <- sub("poly_", "", existing_files)
existing_files <- sub(".csv", "", existing_files)
vars_to_extract_median2 <- vars_to_extract_median2[!vars_to_extract_median2 %in% existing_files]


#getting all rasters that were messed up in SJBs extraction and pairing that with the vector
all_rasters2 <- unique (grep(paste(vars_to_extract_median2,collapse="|"), all_rasters, value=TRUE))
all_rasters <- all_rasters2



### extracting values from existing rasters of spectral indices###
total_start_time <- Sys.time()
for(r in 1:length(all_rasters)){#  1){    #length(all_rasters)){ #  # #start raster loop
  start_time <- Sys.time()
  
  focal_rast <- raster(all_rasters[r]) #r <- 1
  focal_rast <- crop(focal_rast, wm)
  focal_rast_vx <- velox(focal_rast) #focal_rast_vx <- NULL #gc()
  focal_rast_names <- names(focal_rast)
  
  
  
  #for(i in 1:1000){#
  for(i in 1:nrow(all_t_seg)){ #loop through all of the objects #for(i in 1:1000){ #
    poly_to_extract_small <- all_t_seg[i,]
    obj_data_vx <- focal_rast_vx$copy() #extent(obj_data_vx)
    obj_data_vx$crop(poly_to_extract_small) #extent(obj)
    
    bands_d_mean <- obj_data_vx$extract(sp = poly_to_extract_small, fun=mean)
    bands_d_mean_df <- as.data.frame(bands_d_mean)
    names(bands_d_mean_df) <- paste0(focal_rast_names, "_d_mn")
    bands_d_sd <- obj_data_vx$extract(poly_to_extract_small, fun=sd)
    bands_d_sd_df <- as.data.frame(bands_d_sd)
    names(bands_d_sd_df) <- paste0(focal_rast_names, "_d_sd")
    bands_d_median <- obj_data_vx$extract(poly_to_extract_small, fun=median)
    bands_d_median_df <- as.data.frame(bands_d_median)
    names(bands_d_median_df) <- paste0(focal_rast_names, "_d_md")
    dp_focal_poly <- cbind(bands_d_sd_df, bands_d_mean_df, bands_d_median_df)
    if( i == 1){dp <- dp_focal_poly}
    if( i > 1){dp <- bind_rows(dp, dp_focal_poly)}
  } #end band extraction loop
  
  #dp_results <- bind_cols(all_t_seg[1:1000, ], dp) 
  dp_results <- bind_cols(all_t_seg[, ], dp) 
  dp_results$geometry <- NULL
  csv_save_name <- paste0("D:/tree_classification/predictions/csv_for_pred_vars/","poly_",focal_rast_names,".csv")
  write.csv(dp_results, csv_save_name)
  
  #for memory usage
  focal_rast_vx <- NULL #
  gc()
  print(paste("finished:", focal_rast_names))
  print(Sys.time() - start_time)
}#end raster loop

Sys.time() - total_start_time 



###extract data from csvs and make a large df for RF #######################################
setwd("D:/tree_classification/predictions/csv_for_pred_vars")
csvs_for_pred <- dir()
#test <- read.csv(csvs_for_pred[1], nrows = 1000)
myfiles = do.call(cbind, lapply(csvs_for_pred, 
                                function(x) read.csv(x, stringsAsFactors = FALSE, nrows = 1000)))

myfiles2 <- myfiles[, !duplicated(colnames(myfiles), fromLast = TRUE)] 
#rename the hsv bands to match t4 ("_" vs ".")
names(myfiles2) <- gsub(pattern = "hsv_", replacement = "hsv.", x = names(myfiles2))

#compare names against dataframe from predictions DF
#LOAD IN t4 from the tree_class_obj_D_190419.R, takes about a min
#names(t4) # nearmap_181025_hsv.3_d_mn
#names(myfiles2) #nearmap_181025_hsv_d_mn
names(myfiles2)[names(myfiles2) %in% names(t4) == FALSE]
names(t4)[names(t4) %in% names(myfiles2) == FALSE]





## predict object identity using pre-run RF model ============================================
library(randomForest)
all_t_seg <- st_read("E:/tree_classification/segmentation/all_tree_polygons_inD190508.shp")

#load rf model
load("C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/d_rf200503.RData")
print(model1)

#load in extracted data that was extracted by SJB on ~6/5/19
library(readr)
library(dplyr)
library(purrr)
zonal_means <- read_csv("E:/tree_classification/predictions/sjb_rast_extraction/katz_zonal_mean.txt") #, n_max = 100000
zonal_medians <- read_csv("E:/tree_classification/predictions/sjb_rast_extraction/katz_zonal_median.txt") 
zonal_sd <- read_csv("E:/tree_classification/predictions/sjb_rast_extraction/katz_zonal_std.txt")

names(zonal_means) <-  c(names(zonal_means[1:15]), paste0(names(zonal_means[16:115]), "_mn")) #head(zonal_means)
names(zonal_medians) <- c(names(zonal_medians[1:5]), paste0(names(zonal_medians[6:105]), "_md"))
names(zonal_sd) <- c(names(zonal_sd[1:5]), paste0(names(zonal_sd[6:105]), "_sd"))

zonal_extract <- left_join(zonal_means, zonal_medians)
zonal_extract <- left_join(zonal_extract, zonal_sd)

#names of columns that SJB didn't get effectively
test <- zonal_extract %>% 
  sample_n(10000) %>% 
  select(-c(ply__)) #remove the column that was character
  
bad_col_fn <- function(dataset_col){
  n_na <- sum(is.na(dataset_col))
  n_small_n <- length(dataset_col[abs(dataset_col) < 10])
  quantiles_return <- quantile(dataset_col, na.rm = TRUE)
  return(c(n_na, n_small_n, quantiles_return, names(dataset_col)))
}

bad_cols_from_sjb <- map_dfr(test, bad_col_fn) %>% 
         t() %>% 
         as.data.frame() %>%
         mutate(var_name = as.character(row.names(test3))) %>% 
         filter(abs(V4) < 0.0001) %>% #if more than 25% of a column was zeros it was considered bad
         select(var_name) %>%
         filter(var_name != "grond" & var_name != "build" & var_name != "NA_px") %>%
         unlist()

zonal_extract <- zonal_extract %>% select( -one_of(bad_cols_from_sjb)) #removing the bad columns

zonal_extract$idtxt <- NULL
zonal_extract$OBJECTID <- NULL
zonal_extract$Shape_Area <- NULL
zonal_extract$Shape_Length <- NULL
zonal_extract$train0_test1 <- 0
zonal_extract$X.1 <- 0

# #create a lookup table to translate between zonal_extract names and original names 
# test2 <- sort(names(myfiles2))
# test3 <- sort(names(zonal_extract))
# writeClipboard(test3)

#rename column names to match what the original names were
name_vectors <- read.csv("E:/tree_classification/predictions/sjb_rast_extraction/name_match_table190610.csv",
                         stringsAsFactors = FALSE)
names(zonal_extract) = name_vectors$names_orig[match(names(zonal_extract), name_vectors$names_sjb)]

head(zonal_extract)
##make sure that the names line up between zonal_extract and the originally used dataset
#intersect(test2, test3)
#names(myfiles2)[names(myfiles2) %in% names(test) == FALSE]
#names(t4)[names(t4) %in% names(myfiles2) == FALSE]


# adding in the variables that SJB didn't get succesfully
library(readr)
library(purrr)

setwd("E:/tree_classification/predictions/csv_for_pred_vars")
csvs_for_pred <- dir()#[1:3]

compile_csvs <- function(file_name){
  read_csv(file_name, col_types = cols(Id = col_character(), grdcd = col_character())) %>%
    select(-c("X1", "grond", "tree", "build", "NA_px", "X", "nSDM_2017_", "nSDM_201_1", "nSDM_201_2", "ply__"))
}

all_csvs <- map_dfc(csvs_for_pred, compile_csvs) #this is a little messy since the Id and grdcd columns are duplicated.
# double check that all of the grdcd columns are the same (no offsets)
# select(all_csvs, contains("grdcd")) %>% sample_n(10000) %>% t(.) %>% as.data.frame %>% distinct() 
# this check works by making sure that there aren't differences in a sample of the grdcd values

#remove the duplicate grdcd and Id columns
cols_to_remove_grdcd <- names(all_csvs) %>% grep(pattern = "grdcd", x = ., value = TRUE)
cols_to_remove_grdcd <- cols_to_remove_grdcd[2:length(cols_to_remove_grdcd)]
cols_to_remove_Id <- names(all_csvs) %>% grep(pattern = "Id", x = ., value = TRUE)
cols_to_remove_Id <- cols_to_remove_Id[2:length(cols_to_remove_Id)]
cols_to_remove <- c(cols_to_remove_grdcd, cols_to_remove_Id)
all_csvs2 <- all_csvs %>% select(-cols_to_remove) %>%
  mutate(grdcd = as.numeric(grdcd),
         Id = as.numeric(Id))

#combine the good columns from SJB with the good columns from me
cols_from_csvs2 <- colnames(all_csvs2)
good_cols_from_sjb <- colnames(zonal_extract)
cols_not_needed <- good_cols_from_sjb[good_cols_from_sjb %in% cols_from_csvs2]
cols_not_needed <- cols_not_needed[3:141]
zonal_extract_join <- zonal_extract %>% select(-cols_not_needed) %>% select(-c("grond", "tree", "build","NA_px", "X", 
                                                                               "nSDM_2017_", "nSDM_201_1", "nSDM_201_2",
                                                                               "ply__", "X.1", "train0_test1"))
all_pred_vars <- bind_cols(all_csvs2, zonal_extract_join)

#write_csv(all_pred_vars, "E:/tree_classification/predictions/all_predictor_variables200502.csv")
#all_pred_vars <- read_csv("E:/tree_classification/predictions/all_predictor_variables200502.csv")
all_pred_vars2 <- all_pred_vars %>% replace(., is.na(.), 0)
all_pred_vars2$train0_test1 <- 0

#predict myfiles2 tree identity
#myfiles2$train0_test1 <- 0
#predValid <- predict(model1, ValidSet, type = "class")
preds <- as.data.frame(predict(model1, all_pred_vars2, type = "class"))
#test <- as.data.frame(predict(model1, myfiles2, type = "class"))
names(preds) <- "predicted_taxon"
preds2 <- bind_cols(all_pred_vars2, preds)
preds3 <- dplyr::select(preds2, grdcd, Id, predicted_taxon)
preds4 <- left_join(preds3, all_t_seg) #the join converts it to a tibble
preds5 <- st_as_sf(preds4) #explicitly turning it back into an sf

#setwd("")
write_sf(preds5, "E:/tree_classification/predictions/pred200504.shp")
write_sf(preds5, "C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/predictions/pred190715.shp")

summary(preds5$predicted_taxon)
#
preds_sub <- dplyr::filter(preds5, predicted_taxon == "Acer_platanoides") %>% dplyr::select(predicted_taxon)
#preds_sub <- preds5[1:10000,3]
preds_sub_small <- st_crop(preds5, c(xmin= 321234, xmax= 321534, ymin= 4691234, ymax= 4691534)) %>% 
  dplyr::select(predicted_taxon)

plot(preds_sub_small, border = "gray")


#merged the predictions with street tree data, double check the accuracy now
#read in shapefile
tree_pred <- st_read("D:/tree_classification/predictions/pred190715_joinedstreettrees.shp")
names(tree_pred)
tree_pred$prdctd_[1:100]
unique(tree_pred$genus[1:100])
#compare accuracy

# ## tree pollen production equations plug in
# test2 <- tree_pred %>% dplyr::select(prdctd_, SPP, genus, DBH, nSDM_2017_) %>% filter(genus != "<NA>")
# test2$area <- st_area(test2)
# test2$area[test2$prdctd == "Acer_platanoides" & !is.na(test2$prdctd)]
# test2$area <- as.numeric(test2$area)
# test2$predpollen <- NA
# 
# test2$predpollen[test2$prdctd == "Quercus" & !is.na(test2$prdctd)] <- 
#   test2$area[test2$prdctd == "Quercus" & !is.na(test2$prdctd)] * 0.22 + 3.93  
# 
# test2$predpollen[test2$prdctd == "Gleditsia" & !is.na(test2$prdctd)] <- 
#   test2$area[test2$prdctd == "Gleditsia" & !is.na(test2$prdctd)] * 0.46 -3.18  
# 
# test2$predpollen[test2$prdctd == "Platanus" & !is.na(test2$prdctd)] <- 
#   test2$area[test2$prdctd == "Platanus" & !is.na(test2$prdctd)] * 0.07 - 1.04
# 
# p_Acpl <- filter(tree_pred, prdctd_ == "Acer_platanoides" & !is.na(prdctd_)) %>% 
#   dplyr::select(area, predpollen)
# p_Acpl$predpollen <- p_Acpl$area * 0.05 + 0.61  
# 
# 
# d_rast <- raster(ncol = 50, nrow = 50)
# extent(d_rast) <- extent(test)
# #a2_small_extent <- extent(c(13280000, 13280000 + 30000, 265242.3,  300536.3))
# #a2_rast <- crop(a2_rast, a2_small_extent)
# test <- rasterize(p_Acpl, d_rast, field = "predpollen", fun = sum)
# pixel_area <- (res(test)[1] * res(test)[2]) #get pixel area in m2
# test <- (test/ pixel_area) * 1000000 #convert to pol/m2 (was originally in millions)
# test[test < 0] <- 0 #make sure that none of the trees have values below 0
# #plot(test)
# 
# 
# write_sf(p_Acpl, "D:/tree_classification/predictions/p_acpl190806.shp")
# sf_write
# plot(p_Acpl)
# 
# 
# names(test2)
# hist(test2$predpollen)
# 
# #predicted pollen surface 
# 
# 
# test2$geometry <- NULL
# test2$genus_st <- as.character(test2$genus)
# test2$genus_st[test2$SPP == "Acer platanoides"] <- "Acer_platanoides"
# test3 <- dplyr::select(test2, prdctd_, genus_st)
# test4 <- as.data.frame(table(test3))

#mean(predValid == ValidSet$taxa) 
# confusion_df <- as.data.frame.matrix(table(test3$prdctd_, test3$genus_st)) #rows = predicted, col = actual
# confusion_df$row_sum <- rowSums(confusion_df)
# confusion_df[ (nrow(confusion_df) +1) ,] <- colSums(confusion_df)
# 
# ggplot(test2, aes(x = prdctd_, y = genus_st)) + geom_point()

#names(ValidSet)
#names(myfiles2)
#names(myfiles2)[names(myfiles2) %in% names(ValidSet) == FALSE]



### stats for predictions across Detroit ########################################
library(sf)
# read in predictions 
preds5 <- st_read("E:/tree_classification/predictions/pred200504.shp")
#preds5 <- st_read("C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/predictions/pred190715.shp")
preds5$area <- st_area(preds5)
preds5_nogeo <- preds5
preds5_nogeo$geometry <- NULL

names(preds5)
preds5 %>% filter(is.na(predicted_taxon)) #%>% plot()


Table5 <- preds5_nogeo %>% group_by(predicted_taxon) %>%
  summarize(n_trees = n(),
            area_sum = round(sum(area), 0)) %>% 
  arrange(-n_trees)

tree_pred_total_n_trees <- sum(Table5$n_trees)
tree_pred_total_area <- sum(Table5$area_sum)

names(preds5)
test <- unlist(as.numeric(t2$area))
summary(test[test >4 & test < 400])

head(preds5)
