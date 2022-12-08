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
library(gdalUtilities)

#Define master projection (UTM 17)
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Download relevant data 
dem<-raster("data/I_raw/2007_1m_DEM.tif")

jr_clip <- st_read("data/I_raw/JR_boundary.shp")
jl_clip <- st_read("data/I_raw/JL_boundary.shp")

#Download site attributes
site_data_path <- paste0("data/I_raw/Site_Directory_Core.xlsx")

sheet_names <- excel_sheets(path = site_data_path)  

sheet_names <- sheet_names[1:3] %>% 
  as.list()

#Read in site data
sites <- lapply(sheet_names, 
                function(x) read_excel(path = site_data_path, 
                                       sheet = x)) %>% 
  reduce(rbind)

rm(sheet_names, site_data_path) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Clip DEM to workspace -----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Process DEM
dem[dem<0]<-NA

#Had some trouble running this function, received and error message. Works on 
# the 2nd try.
dem<-projectRaster(dem, crs=p, 
                   #Runs really slow needed to specify resolution 
                   #went as granular as I could. 
                   res = 1)

#reproject boundary files
jr_proj <- st_transform(jr_clip, p)
jl_proj <- st_transform(jl_clip, p)

#Create buffer around JR & JL boundary
jr_buffer <- st_buffer(jr_proj, dist = 500)
jl_buffer <- st_buffer(jl_proj, dist = 500)

#Crop DEM
dem_crop_jr <- crop(dem, jr_buffer)
dem_crop_jl <- crop(dem, jl_buffer)

#Combine the crops
dem_crop_all <- merge(dem_crop_jl, dem_crop_jr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Convert site data to sf object --------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sites_all <- sites %>% 
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude = as.numeric(Latitude)) %>%
  st_as_sf( 
    coords = c("Longitude", "Latitude"), 
    crs = '+proj=longlat +datum=WGS84 +no_defs') %>% 
  st_transform(., p) %>% 
  # Remove the river and the Ag ditches. 
  filter(!Site_ID %in% c("CR-SW", "AG-SW", "TR-SW"))
  

jr_sites<- sites_all %>%
  filter(str_detect(`Site Name`,"Wetland")) %>% 
  filter(str_detect(`Catchment`, "Baltimore Corner")) 

jl_sites<- sites %>%
  filter(str_detect(`Site Name`,"Wetland")) %>% 
  filter(Catchment %in% c("Jackson Lane", "Beetree Rd")) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Export data---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

writeRaster(dem_crop_jr, 'data/III_output/dem_jr.tif', overwrite = TRUE)

writeRaster(dem_crop_jl, "data/III_output/dem_jl.tif", overwrite = TRUE)

writeRaster(dem_crop_all, "data/III_output/dem_all.tif", overwrite = TRUE)

#Write well points as shapefiles
st_write(sites_all, "data/III_output/sites_all.shp",
         overite = TRUE, append = FALSE)

st_write(jr_sites, 'data/III_output/jr_sites.shp', 
         overwrite = TRUE, append = FALSE)

st_write(jl_sites, "data/III_output/jl_sites.shp",
         overwrite = TRUE, append = FALSE)

