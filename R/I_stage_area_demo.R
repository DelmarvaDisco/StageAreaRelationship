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
library(sf)
library(raster)
library(tidyverse)

#Define master projection (UTM 17)
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Define scratch workspace
scratch_dir<-"data/II_temp/"

#download DEM site locations
dem<-raster("data/III_output/dem_jr.tif")
sites<-st_read('data/III_output/jr_sites.shp')

#Plot dem for funzies
mapview(dem)

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

#plot Updated dem
mapview(dem_filter)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Define "Root" GIWs --------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Use the Stochastic Deprresion Tool to identify deprresions-----------------
#Export fitlered DEM to workspace
writeRaster(dem_filter, paste0(scratch_dir,"dem_filter.tif"), overwrite=T)

#Apply stochastic depression analysis tool
set.seed(100)
wbt_stochastic_depression_analysis(
  dem = "dem_filter.tif", 
  output = "giws.tif", 
  wd = scratch_dir,
  rmse = 0.18, 
  range = 10, 
  iterations = 100)

#3.2 Define depression based on threshold of occurence in stochastic procedure----
#Reclass raster (any depression that was delineated less than 80% of the time is out!)
wbt_reclass(
  input = "giws.tif", 
  output = "reclass.tif", 
  reclass_vals = "'0;0;0.8'", 
  wd = scratch_dir)
wbt_reclass(
  input = "reclass.tif", 
  output = "reclass.tif", 
  reclass_vals = "'1;0.8;1'", 
  wd = scratch_dir)

#Convert 0 to NA
giws<-raster(paste0(scratch_dir,"reclass.tif")) 
giws<-raster::clump(giws)

#Export as polygon
giws[giws==1]<-NA
giws<- giws %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)

#Write polygon shapes to workspace
st_write(giws, paste0(scratch_dir, "giws.shp"), delete_layer=TRUE)

#3.3 Filter depressions---------------------------------------------------------
#Filter by area and P:A Ratio
wbt_polygon_perimeter(
  input="giws.shp", 
  wd = scratch_dir)
wbt_polygon_area(
  input="giws.shp", 
  wd = scratch_dir)
giws<-st_read(paste0(scratch_dir, "giws.shp"))
giws<-giws %>%
  #Remove small depressions
  filter(AREA>250) %>%
  #Remove oddly large depressions
  filter(AREA<1e5) %>%
  #Remove ditched depressions
  mutate(p_a_ratio = AREA/PERIMETER) %>%
  filter(p_a_ratio>2)

#3.4 Connect depressions to research sites -------------------------------------
#Give ID
giws$WetID<-seq(1, nrow(giws))

#Create function to find giw clostes to each point
site_fun<-function(n){

  #Select site
  site<-sites[n,]
  
  #Select giws within 50 m of site
  giw<-giws[st_buffer(site, 50),]
  
  #Define distances between site and giws
  giw<-giw %>% 
    mutate(dist = st_distance(site, giw, by_element = T))
  
  #Filter to closest giw and export
  export<-giw %>% 
    st_drop_geometry() %>% 
    filter(dist == min(dist, na.rm=T)) %>% 
    select(WetID) %>% 
    mutate(Site_ID = site$Site_ID)
  
  #export results
  export
}

#apply function
site_giw<-lapply(
    FUN = site_fun, 
    X = seq(1, nrow(sites))) %>% 
  bind_rows()

#join to giws
giws<-left_join(giws, site_giw) %>% drop_na()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Create stage-area relationships -------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.1 Create function to develop stage discharge relationship -------------------
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
    X = seq(1,nrow(giws)), 
    FUN = inundation_fun) %>% 
  bind_rows()

#4.3 Export --------------------------------------------------------------------
write_csv(df, "docs/jr_stage_area_relationships.csv")
