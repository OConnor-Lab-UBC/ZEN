# emmett's note
######################################################
#'
#' @title  Code to extract annual temperature data from HadSST
#'
#' @author Jarrett Byrnes
#' @author jarrett.byrnes at umb.edu
#'
#' @log
#' 
#' 5/16/2016 - MO adapting for ZEN
#'  2/25/2016 - fixed up script to run on a multicore server
#'  2/25/2016 - purrr to do means of cells
#'  2/25/2016 - Parallel extraction
#'  2/24/2016 - First version, still having extract problems
#'
######################################################

###### 0) Libraries and the like
#The Hadleyverse
library(dplyr)
library(ncdf4)
library(purrr)
library(readr)
library(tidyr)
library(rgdal)
library(lubridate)

#others
library(raster)
library(doParallel)
library(rgeos)
library(data.table)

#Set raster options to extract more easily
rasterOptions(maxmemory=5e+08,chunksize=5e+07) 

###### 1) Load the unique lat/longs into zen.coords
zen.coords <- read.csv("../site data/ZEN_2014_Lat&Long.csv")

#filter out marine sites, as they are the only ones with SST
unique_lat_long <- zen.coords %>%
  filter(!is.na(zen.coords$Latitude)) %>%
  filter(!is.na(zen.coords$Longitude)) #%>%
#filter(zen.coords) %>%
#filter(BioType == "marine")

###### 2) Load and rasterize HADSST data set from 1950 - 2015
hadsst <- raster::brick("../../HadISST_sst.nc")
raster::NAvalue(hadsst) <- -1000 #make sure we don't have any super small NAs

###### 3a) Test plot to make sure things line up!
ull_points <- SpatialPoints(cbind(unique_lat_long$Longitude, unique_lat_long$Latitude),
                            proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

ull_points <- spTransform(ull_points, projection(hadsst))
plot(hadsst[[1]], ylim = c(20,100))
plot(ull_points, pch = 19, add=T)

###### 4) Extract data at only those lat/longs used
area_raster <- raster::area(hadsst) #gets area in sq km
radii <- sqrt(raster::extract(area_raster, ull_points))/2*1000

safe_mean <- function(x){
  if(sum(is.na(x))==length(x)) return(NA)
  if(sum(x == -1000, na.rm=T) == length(x)) return(NA)
  return(max(x, na.rm=T))
  
}

safe_colMeans <- function(x){
  if(class(x)=="numeric") return(x)
  return(colMeans(x, na.rm=T))
  
}

######Parallelize the extractions
nchunks <- 10
cl <- makeCluster(nchunks)
### CLUSTER ### Create the copies of R running in parallel.
### CLUSTER ### Outfile is for where error is stored.
registerDoParallel(cl)
### CLUSTER ### To use the foreach function

breaks <- round(seq(from = 1, to = nrow(unique_lat_long), length.out=nchunks+1))

hadsst_extracted <- foreach(i = 1:nchunks+1, .combine=rbind) %dopar% {
  ### CLUSTER ### .combine: what does once the analysis in paralel is finished
  ### CLUSTER ### %dopar% "do in paralell"
  
  #which points
  idx <- seq(breaks[i-1], breaks[i], by=1)
  
  #subset of points
  ull <- as.data.frame(ull_points)[idx, ]
  
  #correct for overlap between pieces
  if(i>2) ull <- ull[-1,]
  
  #where are these points - with a 3 degree buffer
  ull_ext <- raster::buffer(ull_points[idx,], width = 3)
  
  hadsst_cropped <- crop(hadsst, ull_ext)
  
  hadsst_vals <- raster::extract(x = hadsst_cropped, 
                                 y = cbind(ull[,1], ull[,2]),
                                 buffer=radii[idx]*1.5,  na.rm=FALSE)
  
  hadsst_df <- map_df(hadsst_vals, function(x) data.frame(t(safe_colMeans(x))))
  
  #look at it
  #  hadsst_df[,1:3]
  
  #fix NaN
  hadsst_df[is.na(hadsst_df)] <- NA
  
  #return the df
  hadsst_df
}

stopImplicitCluster()
stopCluster(cl)
### CLUSTER ### to stop the clustering process.


hadsst_frame <- cbind(unique_lat_long, hadsst_extracted) %>%
  tidyr::gather(DateName, tempC, -Latitude, -Longitude, -len) %>%
  mutate(DateName = gsub("X", "", as.character(DateName))) %>%
  mutate(Year = year(parse_date_time(DateName, orders="ymd"))) %>%
  mutate(tempC=ifelse(tempC == -1000, NA, tempC))


###### 5) Write out temp kelp data as an intermediate step
write.csv(hadsst_frame, "../../workspace/derived_data/hadsst_at_latlongs.csv", row.names=F)