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


## World View 2 image ==============================================================
# Three WV2 images were acquired from Digital Globe via a data grant from 
# UM's Clark Library. Contact people were Mara Blake and, later, Nicole Scholtz
# we received a few folders:
# DigitalGlobeImagery_v1 contains imagery at level 3D which IS orthorectified  
# DigitalGlobeImagery_v2 contains imagery at level 2A which IS NOT orthorectified
# Shannon was in charge of atmospheric correction; to do so he used Erdas Imagine


library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(RStoolbox)
library(randomForest)
library(sp)
library(stringr)


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

### classification #################################################################

### model assessment ###############################################################

### prediction #####################################################################