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
library(tidyverse)

#Define master projection (UTM 17)
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Define scratch workspace
scratch_dir<-"data/II_temp/"

#download DEM site locations
dem<-raster("data/III_output/dem_jr.tif")
sites<-st_read('data/III_output/jr_sites.shp')

#Plot dem for funzies
# mapview(dem)

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

mapview(dem) + mapview(dem_filter)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Define "Root" GIWs --------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Use the Stochastic Deprresion Tool to identify deprresions-----------------
#Export fitlered DEM to workspace
# writeRaster(dem_filter, file = paste0(scratch_dir,"dem_filter.tif"), overwrite=T)

#Apply stochastic depression analysis tool
RNGkind(sample.kind = "Rounding")
set.seed(100)
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


#Convert 0 to NA
giws_80<-raster(paste0(scratch_dir,"reclass_80.tif")) 
giws_80<-raster::clump(giws_80)
giws_95<-raster(paste0(scratch_dir,"reclass_95.tif")) 
giws_95<-raster::clump(giws_95)
giws_50<-raster(paste0(scratch_dir,"reclass_50.tif")) 
giws_50<-raster::clump(giws_50)
giws_97<-raster(paste0(scratch_dir,"reclass_97.tif"))
giws_97<-raster::clump(giws_97)


#Export as polygon
giws_80[giws_80==1]<-NA
giws_80<- giws_80 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
giws_95[giws_95==1]<-NA
giws_95<- giws_95 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
giws_50[giws_50==1]<-NA
giws_50<- giws_50 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
giws_97[giws_97==1]<-NA
giws_97<- giws_97 %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)

#Write polygon shapes to workspace
st_write(giws_80, paste0(scratch_dir, "giws_80.shp"), delete_layer=TRUE)
st_write(giws_95, paste0(scratch_dir, "giws_95.shp"), delete_layer=TRUE)
st_write(giws_50, paste0(scratch_dir, "giws_50.shp"), delete_layer=TRUE)
st_write(giws_97, paste0(scratch_dir, "giws_97.shp"), delete_layer=TRUE)

#Check out the different layers
mapview::mapview(giws_80) + mapview::mapview(giws_95) +
mapview::mapview(giws_50) + mapview::mapview(giws_97)

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
  
#Filter by attributes for the 99% stochastic threshold
  wbt_polygon_perimeter(
    input="giws_97.shp", 
    wd = scratch_dir)
  wbt_polygon_area(
    input="giws_97.shp", 
    wd = scratch_dir)
  giws_97<-st_read(paste0(scratch_dir, "giws_97.shp"))
#No need to filter depressions from giws_99. Even significant wetlands can be large. 
  
mapview(giws_97) + mapview(giws_50) + mapview(giws_80) + mapview(giws_95)


#3.4 Connect depressions to research sites -------------------------------------
#Give ID based on row number
# giws_97$WetID<-seq(1, nrow(giws_97))
# 
# #Create function to find giw closest to each point
# site_fun<-function(sites, giws){
#   
#   #Row numbers for each site
#   n <- seq(1, nrow(sites))
# 
#   # #Select site
#   for(i in sites[n,]){
#   
#   #Select giws within 50 m of site
#   giws<-giws[st_buffer(site, 50),]
#   
#   #Define distances between site and giws
#   giws<-giws %>% 
#     mutate(dist = st_distance(site, giws, by_element = T))
#   
#   #Filter to closest giw and export
#   export<-giws %>% 
#     st_drop_geometry() %>% 
#     filter(dist == min(dist, na.rm=T)) %>% 
#     select(WetID) %>% 
#     mutate(Site_ID = site$Site_ID)
#   
#   #export results
#   export
#   }
#   
#   rbind(export)
#   
# }
# 
# 
# #apply function
# site_giws_97 <-map_dfr(.f = site_fun, .x = sites, .y = giws_97) 
# 
# site <- sites[n, ]

#join to giws
giws_97<-left_join(giws_97, site_giw) %>% drop_na()

mapview(giws_97) + mapview(giws_95) + mapview(giws_80) + mapview(giws_50) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Create stage-area relationships -------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.1 Create function to develop stage area relationship -------------------
inundation_fun<-function(n){

  #Select GIW of interest
  giw<-giws %>% slice(n)
  
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
  df<-lapply(seq(0,2, 0.01), inundate) %>% bind_rows()
  
  #export
  df

}

#4.2 Apply function to individual GIWs -----------------------------------------
df<-lapply(
    X = seq(1,nrow(giws_97)), 
    FUN = inundation_fun) %>% 
  bind_rows()

hipso_plot <- ggplot(data = df, 
                     mapping = aes(x = z,
                                   y = area_m,
                                   color = Site_ID)) +
  geom_line(size = 1.5) +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,1.1)) +
  theme_bw()

(hipso_plot)

#4.3 Export --------------------------------------------------------------------
write_csv(df, "docs/jr_stage_area_relationships.csv")
