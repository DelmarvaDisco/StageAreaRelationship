#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Wetland stage-area relationship demo
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 6/21/2022
#Purpose: Demonstrate wetland stage-area relationship calculations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup workspace------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Download packages of interest 
library(readxl)
library(stars)
library(fasterize)
library(mapview)
library(whitebox)
#Was having trouble with installing whitebox. Added this step.
whitebox::install_whitebox()
library(sf)
library(raster)
library(purrr)
library(igraph)
library(tidyverse)

#Define master projection (UTM 17)
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Define scratch workspace
scratch_dir<-"data/II_temp/"

#download DEM site locations
dem<-raster("data/III_output/dem_all.tif")

sites<-st_read('data/III_output/sites_all.shp') %>% 
  filter(str_detect(Site_ID, "SW"))

#Plot dem for funzies
mapview(dem) + mapview(sites)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Prep DEM for Analysis -----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Export DEM to workspace
writeRaster(dem, paste0(scratch_dir,"dem.tif"), overwrite=T)

#Fill single cell pits
wbt_fill_single_cell_pits(
  dem="dem.tif", 
  output = "dem_fill.tif", 
  wd = scratch_dir)

#Use gaussian filter to smooth wetlands
wbt_gaussian_filter(
  input = "dem_fill.tif", 
  output = "dem_filter.tif", 
  sigma = 2, 
  wd = scratch_dir)

#Read raster back into workspace
dem_filter<-raster(paste0(scratch_dir,"dem_filter.tif"))
crs(dem_filter)<-p

# mapview(dem) + mapview(dem_filter)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Define "Root" GIWs --------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Use the Stochastic Deprresion Tool to identify deprresions-----------------
#Export fitlered DEM to workspace
# writeRaster(dem_filter, file = paste0(scratch_dir,"dem_filter.tif"), overwrite=T)

#Apply stochastic depression analysis tool
RNGkind(sample.kind = "Rounding")
set.seed(100)
# !!! Note this step takes a very long time to run. Output is on google drive inside the
# !!! II_temp folder file name = giws.tif. It will save some time. 
# !!! Took ~1.5 hrs on 16 GB of RAM. Be prepared to do other things while this model 
# !!! Runs in the background. 
wbt_stochastic_depression_analysis(
  dem = "dem_filter.tif", 
  output = "giws.tif", 
  wd = scratch_dir,
  rmse = 0.18, 
  range = 10, 
  iterations = 1000,
  verbose_mode = TRUE)

#3.2 Define depression based on threshold of occurrence in stochastic procedure----
#Reclass raster (any depression that was delineated less than 80% of the time is out!)
wbt_reclass(
  input = "giws.tif", 
  output = "reclass_80.tif", 
  reclass_vals = "'0;0;0.8'", 
  wd = scratch_dir)
wbt_reclass(
  input = "reclass_80.tif", 
  output = "reclass_80.tif", 
  reclass_vals = "'1;0.8;1'", 
  wd = scratch_dir)

#Reclass with 95% threshold
wbt_reclass(
  input = "giws.tif", 
  output = "reclass_95.tif", 
  reclass_vals = "'0;0;0.95'", 
  wd = scratch_dir)
wbt_reclass(
  input = "reclass_95.tif", 
  output = "reclass_95.tif", 
  reclass_vals = "'1;0.95;1'", 
  wd = scratch_dir)

#Reclass with a 50% threshold
wbt_reclass(
  input = "giws.tif", 
  output = "reclass_50.tif", 
  reclass_vals = "'0;0;0.5'", 
  wd = scratch_dir)
wbt_reclass(
  input = "reclass_50.tif", 
  output = "reclass_50.tif", 
  reclass_vals = "'1;0.5;1'", 
  wd = scratch_dir)

#Reclass with 97% threshold
wbt_reclass(
  input = "giws.tif",
  output = "reclass_97.tif",
  reclass_vals = "'0;0;0.97'",
  wd = scratch_dir)
wbt_reclass(
  input = "reclass_97.tif",
  output = "reclass_97.tif",
  reclass_vals = "'1;0.97;1'",
  wd = scratch_dir)

#Reclass with 70% threshold
wbt_reclass(
  input = "giws.tif",
  output = "reclass_70.tif",
  reclass_vals = "'0;0;0.7'",
  wd = scratch_dir)
wbt_reclass(
  input = "reclass_97.tif",
  output = "reclass_97.tif",
  reclass_vals = "'1;0.7;1'",
  wd = scratch_dir)

#Reclass with 25% threshold
wbt_reclass(
  input = "giws.tif",
  output = "reclass_25.tif",
  reclass_vals = "'0;0;0.25'",
  wd = scratch_dir)
wbt_reclass(
  input = "reclass_25.tif",
  output = "reclass_25.tif",
  reclass_vals = "'1;0.25;1'",
  wd = scratch_dir)



#Convert 0 to NA
giws_80<-raster(paste0(scratch_dir,"reclass_80.tif")) 
giws_80<-raster::clump(giws_80)
giws_95<-raster(paste0(scratch_dir,"reclass_95.tif")) 
giws_95<-raster::clump(giws_95)
giws_50<-raster(paste0(scratch_dir,"reclass_50.tif")) 
giws_50<-raster::clump(giws_50)
giws_97<-raster(paste0(scratch_dir,"reclass_97.tif"))
giws_97<-raster::clump(giws_97)
giws_70<-raster(paste0(scratch_dir,"reclass_70.tif"))
giws_70<-raster::clump(giws_70)
giws_25<-raster(paste0(scratch_dir,"reclass_25.tif"))
giws_25<-raster::clump(giws_25)


#Export as polygon
giws_80[giws_80==1]<-NA
giws_80<- giws_80 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
giws_95[giws_95==1]<-NA
giws_95<- giws_95 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
giws_50[giws_50==1]<-NA
giws_50<- giws_50 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
giws_97[giws_97==1]<-NA
giws_97<- giws_97 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
giws_70[giws_70==1]<-NA
giws_70<- giws_70 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
giws_25[giws_25==1]<-NA
giws_25<- giws_25 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)


#Write polygon shapes to workspace
st_write(giws_80, paste0(scratch_dir, "giws_80.shp"), delete_layer=TRUE)
st_write(giws_95, paste0(scratch_dir, "giws_95.shp"), delete_layer=TRUE)
st_write(giws_50, paste0(scratch_dir, "giws_50.shp"), delete_layer=TRUE)
st_write(giws_97, paste0(scratch_dir, "giws_97.shp"), delete_layer=TRUE)
st_write(giws_70, paste0(scratch_dir, "giws_70.shp"), delete_layer=TRUE)
st_write(giws_25, paste0(scratch_dir, "giws_25.shp"), delete_layer=TRUE)

#3.3 Filter depressions---------------------------------------------------------

#Filter by area and P:A Ratio for 80% stochastic threshold
wbt_polygon_perimeter(
  input="giws_80.shp", 
  wd = scratch_dir)
wbt_polygon_area(
  input="giws_80.shp", 
  wd = scratch_dir)
giws_80<-st_read(paste0(scratch_dir, "giws_80.shp"))
giws_80<-giws_80 %>%
  #Remove small depressions
  filter(AREA>250) %>%
  #Remove oddly large depressions
  filter(AREA<1e5) %>%
  #Remove ditched depressions
  mutate(p_a_ratio = AREA/PERIMETER) %>%
  filter(p_a_ratio>2)

#Filter by area and P:A Ratio for the 95% stochastic threshold
  wbt_polygon_perimeter(
    input="giws_95.shp", 
    wd = scratch_dir)
  wbt_polygon_area(
    input="giws_95.shp", 
    wd = scratch_dir)
  giws_95<-st_read(paste0(scratch_dir, "giws_95.shp"))
  giws_95<-giws_95 %>%
    #Remove oddly large depressions
    filter(AREA<1e5) %>%
  #Remove ditched depressions
  mutate(p_a_ratio = AREA/PERIMETER) %>%
    filter(p_a_ratio>2)
  
#Filter by attributes for the 50% stochastic threshold
  wbt_polygon_perimeter(
    input="giws_50.shp", 
    wd = scratch_dir)
  wbt_polygon_area(
    input="giws_50.shp", 
    wd = scratch_dir)
  giws_50<-st_read(paste0(scratch_dir, "giws_50.shp"))
  giws_50<-giws_50 %>%
    #Remove small depressions
    filter(AREA>250) %>%
    #Remove oddly large depressions
    filter(AREA<1e5) %>%
    #Remove ditched depressions
    mutate(p_a_ratio = AREA/PERIMETER)
    # filter(p_a_ratio>2)
  
#Filter by attributes for the 70% stochastic threshold
  wbt_polygon_perimeter(
    input="giws_70.shp", 
    wd = scratch_dir)
  wbt_polygon_area(
    input="giws_70.shp", 
    wd = scratch_dir)
  giws_70<-st_read(paste0(scratch_dir, "giws_70.shp"))
  giws_70<-giws_70 %>%
    #Remove small depressions
    filter(AREA>250) %>%
    #Remove oddly large depressions
    filter(AREA<1e5) %>%
    #Remove ditched depressions
    mutate(p_a_ratio = AREA/PERIMETER)
  # filter(p_a_ratio>2)
  
  
#Filter by attributes for the 97% stochastic threshold
  wbt_polygon_perimeter(
    input="giws_97.shp", 
    wd = scratch_dir)
  wbt_polygon_area(
    input="giws_97.shp", 
    wd = scratch_dir)
  giws_97<-st_read(paste0(scratch_dir, "giws_97.shp")) %>% 
    mutate(p_a_ratio = AREA/PERIMETER)
  
#Filter by attributes for the 25% stochastic threshold
  wbt_polygon_perimeter(
    input="giws_25.shp", 
    wd = scratch_dir)
  wbt_polygon_area(
    input="giws_25.shp", 
    wd = scratch_dir)
  giws_25<-st_read(paste0(scratch_dir, "giws_25.shp"))
  giws_25<-giws_25 %>%
    #Remove small depressions
    filter(AREA>250) %>%
    #Remove oddly large depressions
    filter(AREA<1e5) %>%
    #Remove ditched depressions
    mutate(p_a_ratio = AREA/PERIMETER)
  # filter(p_a_ratio>2)
  
  
#No need to filter depressions from giws_99. Even significant wetlands can be large. 
p <- "+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
   
st_crs(giws_97) <- p
st_crs(giws_25) <- p
st_crs(giws_50) <- p
st_crs(giws_80) <- p
st_crs(giws_95) <- p
st_crs(giws_70) <- p
crs(dem) <- p

mapview(giws_50, alpha = 0.1) + mapview(giws_80, alpha = 0.4) +
mapview(giws_95, alpha = 0.6) + mapview(giws_70, alpha = 0.25) + mapview(giws_25, alpha = 0.05) + mapview(dem)

#3.4 Connect depressions to research sites -------------------------------------

#Give ID based on row number
giws_50$WetID<-seq(1, nrow(giws_50))
giws_80$WetID<-seq(1, nrow(giws_80))
giws_95$WetID<-seq(1, nrow(giws_95))
giws_97$WetID<-seq(1, nrow(giws_97))
giws_70$WetID<-seq(1, nrow(giws_70))


#Create function to find giw closest to each point
site_fun<-function(n){

  #Select site
  site<-sites[n,]

  #Select giws within 50 m of site
  giw_x<-giws_x[st_buffer(site, 10),]

  #Define distances between site and giws
  giw_x<-giw_x %>%
    mutate(dist = st_distance(site, giw_x, by_element = T))

  #Filter to closest giw and export
  export<-giw_x %>%
    st_drop_geometry() %>%
    filter(dist == min(dist, na.rm=T)) %>%
    dplyr::select(WetID) %>%
    mutate(Site_ID = site$Site_ID)

  #export results
  export
}

#Apply site matching function to each stochastic threshold

#50% threshold
giws_x <- giws_50
site_giw<-lapply(FUN = site_fun, 
                    X = seq(1, nrow(sites))) %>%
  bind_rows()
giws_50_matched <-left_join(giws_50, site_giw) %>% drop_na()
#Clean up environment
rm(giws_50, giws_x, site_giw)

#80% threshold
giws_x <- giws_80 
site_giw<-site_giw<-lapply(FUN = site_fun, 
                           X = seq(1, nrow(sites))) %>%
  bind_rows()
giws_80_matched <- left_join(giws_80, site_giw) %>% drop_na()
#Clean up environment
rm(giws_80, giws_x, site_giw)

#95% threshold
giws_x <- giws_95
site_giw<-site_giw<-lapply(FUN = site_fun, 
                           X = seq(1, nrow(sites))) %>%
  bind_rows()
giws_95_matched <- left_join(giws_95, site_giw) %>% drop_na()
#Clean up environment
rm(giws_95, giws_x, site_giw)

#97% threshold
giws_x <- giws_97
site_giw<-site_giw<-lapply(FUN = site_fun, 
                           X = seq(1, nrow(sites))) %>%
  bind_rows()
giws_97_matched <- left_join(giws_97, site_giw) %>% drop_na()
#Clean up environment
rm(giws_97, giws_x, site_giw)

#70% threshold
giws_x <- giws_70
site_giw<-site_giw<-lapply(FUN = site_fun, 
                           X = seq(1, nrow(sites))) %>%
  bind_rows()
giws_70_matched <- left_join(giws_70, site_giw) %>% drop_na()
#Clean up environment
rm(giws_70, giws_x, site_giw)

# Look at mapview
mapview(giws_50_matched) + mapview(giws_80_matched) + 
mapview(giws_95_matched) + mapview(giws_97_matched) +
mapview(giws_70_matched) + mapview(dem)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Create stage-area relationships -------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.1 Create function to develop stage area relationship -------------------
inundation_fun<-function(n){

  #Select GIW of interest
  giw<-giws_x %>% slice(n)
  
  #gather information about raster
  r<-res(dem)[1]
  ext<-raster(extent(giw), res=r)
  dem_giw<-rasterize(giw, ext, field=1)
  dem_giw<-dem_giw*dem
  
  #Create Minimum Raster
  dem_giw_min<-dem_giw*0+minValue(dem_giw)
  
  #Create function to return conditional raster
  Con<-function(condition, trueValue, falseValue){
    return(condition * trueValue + (!condition)*falseValue)
  }
  
  #Create function to calculate inundation area/volume
  inundate<-function(z){
    #Create area and volume rasters on rise in water level (z)
    area   <- Con(dem_giw>(dem_giw_min+z),0,1)
    volume <- (((z+dem_giw)-dem_giw_min)*area)*(r^2)
    
    #Calculate area and volume based on rasters
    output <- tibble(
      Site_ID = giw$Site_ID,
      z, #Depth (m)
      area_m = cellStats(area, 'sum')*(r^2), #area (m^2)
      volume_m3 =cellStats(volume, 'sum')) #volume (m^3)
    
    #Export Ouptut
    output
  }  
  
  #Apply inundation for different stage values
  df<-lapply(seq(z_low, z_high, 0.01), inundate) %>% bind_rows()
  
  #export
  df
}

# Inundation function for 80% threshold (low water; NO merged sites)
z_low <- 0
z_high <- 1.25

giws_x <- giws_80_matched

inundate_80 <- lapply(FUN = inundation_fun,
                      X = seq(1, nrow(giws_x))) %>% 
  bind_rows()
#Clean up environment
rm(z_low, z_high, giws_x)

# Inundation function for 80% stochastic threshold (typical water levels; some merges)
# z_low <- .51
# z_high <- 1.1
# giws_x <- giws_80_matched
# 
# inundate_80 <- lapply(FUN = inundation_fun,
#                       X = seq(1, nrow(giws_x))) %>%
#   bind_rows()
# 
# rm(z_low, z_high, giws_x)
# 
# # Inundation function for 70% Stochastic threshold (high water levels: tons of merging)
# 
# z_low <- 1.11
# z_high <- 1.5
# giws_x <- giws_70_matched
# inundate_70 <- lapply(FUN = inundation_fun,
#                       X = seq(1, nrow(giws_x))) %>%
#   bind_rows()
# 
# rm(z_low, z_high, giws_x)

# Combine different water levles
# inundate_full <- bind_rows(inundate_70, inundate_80, inundate_97)

# 4.3 Quick check ---------------------------------------------------------

hipso_plot <- ggplot(data = inundate_80 %>% 
                       filter(Site_ID %in% c("TS-SW", "BD-SW", "DK-SW", "ND-SW",
                                             "FN-SW")),
                   mapping = aes(x = z,
                                 y = area_m,
                                 color = Site_ID)) +
geom_line(size = 1.5) +
scale_y_continuous() +
# scale_x_continuous(limits = c(0,1.1)) +
theme_bw()

(hipso_plot)

#6.0 Export --------------------------------------------------------------------
write_csv(inundate_80, "docs/stage_area_relationships.csv")
