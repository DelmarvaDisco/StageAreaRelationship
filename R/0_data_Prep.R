#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Data Prep
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 6/21/2022
#Purpose: Prep DEM and site for analysis at Jones Road
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup workspace -----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Download packages of interest 
library(readxl)
library(stars)
library(fasterize)
library(whitebox)
library(sf)
library(mapview)
library(raster)
library(tidyverse)

#Define master projection (UTM 17)
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Download relevant data 
dem<-raster("data/I_raw/2007_1m_DEM.tif")
jr_clip<-st_read("data/I_raw/JR_boundary.shp")

#Download site attributes
site_data_path <- paste0("data/I_raw/Site_Directory_Core.xlsx")

sheet_names <- excel_sheets(path = site_data_path)  

sheet_names <- sheet_names[1:3] %>% 
  as.list()

#Read in site data
jr_sites <- lapply(sheet_names, 
                    function(x) read_excel(path = site_data_path, 
                                           sheet = x)) %>% 
  reduce(rbind)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Clip DEM to workspace -----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Process DEM
dem[dem<0]<-NA
dem<-projectRaster(dem, crs=p)  

#reproject JR boundary file
jr_proj<-st_transform(jr_clip, p)

#Create buffer around JR boundary
jr_buffer<-st_buffer(jr_proj, dist = 500)

#Crop DEM
dem_crop<-crop(dem, jr_buffer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Convert site data to sf object --------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jr_sites<- jr_sites %>%
  filter(str_detect(`Site Name`,"Wetland")) %>% 
  filter(str_detect(`Catchment`, "Baltimore Corner")) %>% 
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude = as.numeric(Latitude)) %>%
  st_as_sf( 
    coords = c("Longitude", "Latitude"), 
    crs = '+proj=longlat +datum=WGS84 +no_defs') %>% 
  st_transform(., p)

mapview(dem_crop) + mapview(jr_sites)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Export data---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
writeRaster(dem_crop, 'data/III_output/dem_jr.tif', overwrite = TRUE)
st_write(jr_sites, 'data/III_output/jr_sites.shp', overwrite = TRUE, append = FALSE)
