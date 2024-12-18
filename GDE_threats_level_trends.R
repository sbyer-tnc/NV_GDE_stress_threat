# Groundwater level data
# Consolidate data from NWIS and NDWR sites to claculate a trend in groundwater level for each site
# Calculate trends using modified mann-kendall package
#
#---------------------------------------------------------------------
# Load packages

#library(rgdal)
#library(sp)
#library(raster)
#library(readOGR)
library(sf)
library(terra)
library(mblm)
library(modifiedmk)


#library(rgeos)
#library(tidyverse)
#library(dplyr)
#library(ggplot2)

# USGS-derived packages to retrieve data in R
# https://owi.usgs.gov/R/dataRetrieval.html#1 (tutorial)
#install.packages("dataRetrieval")
library(dataRetrieval)
#GRAN_pkg <- available.packages(contrib.url("https://owi.usgs.gov/R"))
#names(GRAN_pkg[,1]) # Geological Survey R Archive Network 

#---------------------------------------------------------------------
#---------------------------------------------------------------------

# custom function for filtering USGS data for sites that meet data criteria (see below) and calculate trend over time
# Function has a couple of checks to ensure valid data are included
# Returns a dataframe for the groundwater site with trend statistics (p-value, sens-slope, min year, max yea, etc.) 


# Site criteria:
# - groundwater level data from 5 or more unique years between 1984 - 2021
# e.g. a site will be included if it has one or more observations from 1999, 2000, 2001, 2013, 2018, 2020
# if a site meets the criteria, the function will calculate a trend and statistics

annWL <- function(site_number){
  x <- s[s$site_no==site_number,]
  x$lev_va <- -1 * x$lev_va # Reverse groundwater values so sens-slope trend is negative where depth to gw is increasing
  # If there are zero (less than one) observations, do not include site in analysis
  if (nrow(x) < 1){
    print(cat(site_number,"does not have enough annual readings to calculate trend"))
    df <- data.frame(SITENO = site_number, New_Pval = NA, Old_Pval = NA,
                     n_effective = NA, n_raw = NA,
                     Sens_Slope = NA, MinYear = NA, MaxYear = NA,
                     Notes = "Not enough annual measurements")
  } 
  # If there are only no-data (null) observations, do not include site in analysis
  else if (is.na(sum(x$lev_va) == TRUE)){
    print(cat(site_number,"only has NULL readings"))
    df <- data.frame(SITENO = site_number, New_Pval = NA, Old_Pval = NA,
                     n_effective = NA, n_raw = NA,
                     Sens_Slope = NA, MinYear = NA, MaxYear = NA,
                     Notes = "Only NULL measurement values")
  }
  # Continue analysis if sites has valid observations
  # Calculates annual average if there are more than one observations per year
  else {
    x <- x %>% mutate(DATE = as.Date(lev_dt, format = "%Y-%m-%d"))
    annavg <- aggregate(lev_va ~ cut(DATE, "1 year"), x, mean)
    annavg$YEAR <- as.numeric(substr(annavg$`cut(DATE, "1 year")`, 1, 4))
    # Calculate trend (sens-slope) if there are 5 or more years of observations
    if (var(annavg$lev_va) > 0 & nrow(annavg) >= 5){
      modMK <- data.frame(t(mmkh(annavg$lev_va)))
      df <- data.frame(SITENO = site_number, New_Pval = modMK$new.P.value, Old_Pval = modMK$old.P.value,
                       n_effective = modMK$N.N., n_raw = nrow(annavg[which(!is.na(annavg$lev_va)),]),
                       Sens_Slope = modMK$Sen.s.slope, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), Notes = NA)
    } 
    # Do not calculate trend there is just one water level observation
    else if (var(annavg$lev_va) == 0 & nrow(annavg) >= 5){
      print(cat(site_number,"Only has one water level value - cannot run modified mann-kendall"))
      df <- data.frame(SITENO = site_number, New_Pval = NA, Old_Pval = NA,
                       n_effective = NA, n_raw = nrow(annavg[which(!is.na(annavg$lev_va)),]),
                       Sens_Slope = 0, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), 
                       Notes = "Only one water level value; trend == 0")
    } 
    # Do not calculate trend if there are not enough observations
    else {
      print(cat(site_number,"does not have enough annual readings to calculate trend"))
      df <- data.frame(SITENO = site_number, New_Pval = NA, Old_Pval = NA,
                       n_effective = NA, n_raw = NA,
                       Sens_Slope = NA, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR),
                       Notes = "Not enough annual measurements")
    }
    
    return(df)}
}

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Get USGS data from NWIS database

# Use Hydrographic Areas (HAs) to spatially-filter groundwater sites
# Could use other polygons (i.e. WBDs aka HUCs)
gdb <- "K:\\GIS3\\Projects\\GDE\\Maps\\GDE_Threats\\GDE_Threats.gdb" # Location (geodatabase in this case) for hydrographics areas (HAs)
fc_list <- st_layers(gdb)
ha <- sf::st_read(gdb, layer = "hydrographic_basin_boundaries")
ha <- sf::st_transform(ha, crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs ")
plot(ha['HYD_AREA'])

# You can read in polygons from any format (gdb feature class, shapefile, etc.)
# Just needs to be written to an sf object

# State of Nevada is too large to grab all water level data at once using dataRetieval package
# Loop through all basins to:
# - run readNWISdata() function using each HA to select groundwater sites
# - download groundwater level observations for the sites
# - run annWL() function to calculate trend, or note whether sites don't meet criteria


# Calculate trends and format data
# Create empty data frame to populate
annWL_df <- data.frame(SITENO = as.character(), New_Pval = as.numeric(), Old_Pval = as.numeric(),
                       n_effective = as.numeric(), n_raw = as.numeric(), Sens_Slope = as.numeric(),
                       MinYear = as.numeric(), MaxYear = as.numeric(), Notes = as.character())

# Use for-loop and annWL() to populate empty data frame with trend stats
for (j in seq(1, nrow(ha))){
  # Get sites for one basin
  # Create bounding box for the HA to get sites
  poly <- ha[j,]
  av_box <- poly@bbox
  mybox <- round(c(av_box[1,1], av_box[2,1], av_box[1,2], av_box[2,2]), 4) # bounding box lat/long coordiates
  av_gw <- readNWISdata(bBox = mybox, service="gwlevels", parameterCd="72019") # Sites with water level below surface measurements
  av_sites <- as.numeric(av_gw$site_no)
  if (is.na(av_gw$site_no)){
    print(cat(poly$HYD_AREA_N, "Has no valid sites"))
    # Careful, this might skip sites in hyd areas
    # Might be better to remove problematic areas from the loop
  } else{
    av_sites_loc <- readNWISsite(av_sites) # Get location data for sites
    av_points <- SpatialPointsDataFrame(coords = data.frame(av_sites_loc$dec_long_va, av_sites_loc$dec_lat_va),
                                        data = av_sites_loc)
    # Clip sites to Area of Interest
    crs(av_points) <- crs(poly) # set coord system of site spoints to be the same as selecitng polygon
    #crs(av_points) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") # Set CRS of spatialPointsDataFrame
    nwis <- av_points[poly, ] # Grab points within basin polygon
    if (nrow(nwis) == 0){
      print(cat(poly$HYD_AREA, "Does not have any NWIS sites"))
      df <- data.frame(SITENO = poly$HYD_AREA, New_Pval = NA, Old_Pval = NA,
                       n_effective = NA, n_raw = NA,
                       Sens_Slope = NA, MinYear = NA, MaxYear = NA,
                       Notes = "No NWIS sites in basin")
      annWL_df <- rbind(annWL_df, df)
    }
    else {
      nwis_sites <- nwis@data$site_no # Site numbers of sites in polygon
      # Get groundwater level data for final list of sites in given time period
      start.date <- "1984-01-01"
      end.date <- "2021-12-31"
      s <- readNWISgwl(nwis_sites, startDate = start.date, endDate = end.date)
      # Use for-loop and annWL() to populate empty data frame with trend stats
      for (i in seq(1, nrow(nwis))){
        siteno <- nwis[i,]$site_no
        bad_sites <- c("393233115594801", "404808116220801")
        if (siteno %in% bad_sites){
          print("YIKES don't use this one")
        } else{
          site_df <- annWL(site_number = siteno)
          annWL_df <- rbind(annWL_df, site_df) 
        }
      }
    }
  }
}
annWL_df_NWIS <- annWL_df


# Remove sites where siteid is BAD
# ex. nchar(as.character(annWL_df_NWIS$SITENO[5611]))
rem_index <- vector()
for(i in seq(1, nrow(annWL_df_NWIS))){
  if(nchar(as.character(annWL_df_NWIS$SITENO[i])) < 7){
    print("Bad site ID")
    print(as.vector(annWL_df_NWIS$SITENO[i]))
    rem_index <- c(rem_index, i)
  }
}
annWL_df_NWIS$SITENO[rem_index] # Print bad site numbers
annWL_df_NWIS <- annWL_df_NWIS[-c(rem_index),] # Remove bad site numbers from df
dim(annWL_df)
dim(annWL_df_NWIS)

# Remove a couple additional sites
rem_sites <- c(364640114050301, 364727114045601) # These sites have 2 entries each - remove one at least
annWL_df_NWIS[which(annWL_df_NWIS$SITENO %in% rem_sites),]
annWL_df_NWIS[c(12173, 12181),]
annWL_df_NWIS <- annWL_df_NWIS[-c(12173, 12181),]
dim(annWL_df_NWIS)


#----------------------------------
# Append CNRWA data to NWIS sites
# Central NV Regional Water Authority has additional data not-yet integrated in NWIS
cnrwa <- read.csv("E:\\RCF Data Recovery\\Recovered Files\\P00\\GDE_Threats\\Hydrology\\cnrwa_add_wells.csv")
head(cnrwa)
cnrwa <- cnrwa %>% dplyr::mutate(DATE = as.Date(Date, format = "%m/%d/%Y"))
cnrwa$YEAR <- as.numeric(substr(cnrwa$DATE, 1, 4)) 
head(cnrwa)

# See which sites have additional data from the 2021 CNRWA report
# Sites that already exist in the trend calculations
# Calculate new trend with additional data and REPLACE
check_nwis <- as.vector(unique(cnrwa[which(cnrwa$site_no %in% annWL_df_NWIS$SITENO),]$site_no))
check_nwis


# Function to find and replace NWIS data
new_trend <- function(site_number){
  og_site_trend <- annWL_df_NWIS[annWL_df_NWIS$SITENO==site_number,]
  cnrwa_site_data <- cnrwa[cnrwa$site_no==site_number,]
  #print(og_site_trend$MaxYear)
  #print(max(cnrwa_site_data$YEAR))
  if ((og_site_trend$MaxYear < max(cnrwa_site_data$YEAR)) | is.na(og_site_trend$MaxYear)){
    print(paste("Re-calculate and replace site trend", site_number, sep=": "))
    # Get OG NWIS gw level readings
    start.date <- "1984-01-01"
    end.date <- "2021-12-31"
    nwis_lvls <- readNWISgwl(site_number, startDate = start.date, endDate = end.date)
    nwis_lvls$lev_va <- -1 * nwis_lvls$lev_va
    # Append CNRWA gw levels - okay if it repeats a reading
    # i.e. 2010 readings from both - will be averaged to annual 2010 reading
    append_lvls <- full_join(nwis_lvls, cnrwa_site_data, by=c("site_no" = "site_no", "lev_va" = "lvl_below_surf", 
                                                              "lev_dt" = "DATE"))
    annavg <- aggregate(lev_va ~ cut(lev_dt, "1 year"), append_lvls, mean)
    annavg$YEAR <- as.numeric(substr(annavg$`cut(lev_dt, "1 year")`, 1, 4))
    annavg <- annavg %>% dplyr::filter(YEAR >= 1984)
    if (var(annavg$lev_va) > 0 & nrow(annavg) >= 5){
      modMK <- data.frame(t(mmkh(annavg$lev_va)))
      df <- data.frame(SITENO = site_number, New_Pval = modMK$new.P.value, Old_Pval = modMK$old.P.value,
                       n_effective = modMK$N.N., n_raw = nrow(annavg[which(!is.na(annavg$lev_va)),]),
                       Sens_Slope = modMK$Sen.s.slope, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), 
                       Notes = "USGS groundwater levels supplemented by 2021 CNRWA report data")
    } else if (var(annavg$lev_va) == 0 & nrow(annavg) >= 5){
      print(cat(site_number,"Only has one water level value - cannot run modified mann-kendall"))
      df <- data.frame(SITENO = site_number, New_Pval = NA, Old_Pval = NA,
                       n_effective = NA, n_raw = nrow(annavg[which(!is.na(annavg$lev_va)),]),
                       Sens_Slope = 0, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), 
                       Notes = "Only one water level value; trend == 0")
    } else {
      print(cat(site_number,"does not have enough annual readings to calculate trend"))
      df <- data.frame(SITENO = site_number, New_Pval = NA, Old_Pval = NA,
                       n_effective = NA, n_raw = NA,
                       Sens_Slope = NA, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR),
                       Notes = "Not enough annual measurements")
    }
    return(df)
  }
}

# Recalculate
nwis_replace <- data.frame(SITENO = as.character(), New_Pval = as.numeric(), Old_Pval = as.numeric(),
                           n_effective = as.numeric(), n_raw = as.numeric(), Sens_Slope = as.numeric(),
                           MinYear = as.numeric(), MaxYear = as.numeric(), Notes = as.character())
for(i in seq(1, length(check_nwis))){
  s <- check_nwis[i]
  x <- new_trend(s)
  nwis_replace <- rbind(nwis_replace, x)
}
nwis_replace

# Replace
for(i in seq(1, nrow(nwis_replace))){
  x <- as.vector(nwis_replace$SITENO[i])
  print(x)
  nwis_index <- which(annWL_df_NWIS$SITENO==x)
  print(nwis_index)
  annWL_df_NWIS[nwis_index,] <- nwis_replace[i,]
}
annWL_df_NWIS[which(annWL_df_NWIS$SITENO==403515114571701),]
annWL_df[which(annWL_df$SITENO==403515114571701),]
#nwis_replace[which(nwis_replace$SITENO==as.vector(annWL_df_NWIS[659,]$SITENO)),]
#which(nwis_replace$SITENO==annWL_df_NWIS[659]$SITENO)


#----------------------------------
# See which NEW sites are in the 2021 CNRWA report
# Sites that I dont already have in the trend calculations
`%notin%` <- Negate(`%in%`)
add_nwis <- unique(cnrwa[c(which(cnrwa$site_no %notin% annWL_df_NWIS$SITENO)),]$site_no)
add_nwis <- add_nwis[-c(1, 3, 4, 10)]
add_nwis

nwis_new <- data.frame(SITENO = as.character(), New_Pval = as.numeric(), Old_Pval = as.numeric(),
                           n_effective = as.numeric(), n_raw = as.numeric(), Sens_Slope = as.numeric(),
                           MinYear = as.numeric(), MaxYear = as.numeric(), Notes = as.character())
for(i in seq(1, length(add_nwis))){
  site_number <- add_nwis[i]
  cnrwa_site_data <- cnrwa[cnrwa$site_no==site_number,]
  start.date <- "1984-01-01"
  end.date <- "2021-12-31"
  nwis_lvls <- readNWISgwl(add_nwis[i], startDate = start.date, endDate = end.date)
  nwis_lvls$lev_va <- -1 * nwis_lvls$lev_va
  append_lvls <- full_join(nwis_lvls, cnrwa_site_data, by=c("site_no" = "site_no", "lev_va" = "lvl_below_surf", 
                                                            "lev_dt" = "DATE"))
  annavg <- aggregate(lev_va ~ cut(lev_dt, "1 year"), append_lvls, mean)
  annavg$YEAR <- as.numeric(substr(annavg$`cut(lev_dt, "1 year")`, 1, 4))
  annavg <- annavg %>% dplyr::filter(YEAR >= 1984)
  if (var(annavg$lev_va) > 0 & nrow(annavg) >= 5){
    modMK <- data.frame(t(mmkh(annavg$lev_va)))
    df <- data.frame(SITENO = site_number, New_Pval = modMK$new.P.value, Old_Pval = modMK$old.P.value,
                     n_effective = modMK$N.N., n_raw = nrow(annavg[which(!is.na(annavg$lev_va)),]),
                     Sens_Slope = modMK$Sen.s.slope, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), 
                     Notes = "USGS groundwater levels supplemented by 2021 CNRWA report data")
  } else if (var(annavg$lev_va) == 0 & nrow(annavg) >= 5){
    print(cat(site_number,"Only has one water level value - cannot run modified mann-kendall"))
    df <- data.frame(SITENO = site_number, New_Pval = NA, Old_Pval = NA,
                     n_effective = NA, n_raw = nrow(annavg[which(!is.na(annavg$lev_va)),]),
                     Sens_Slope = 0, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), 
                     Notes = "Only one water level value; trend == 0")
  } else {
    print(cat(site_number,"does not have enough annual readings to calculate trend"))
    df <- data.frame(SITENO = site_number, New_Pval = NA, Old_Pval = NA,
                     n_effective = NA, n_raw = NA,
                     Sens_Slope = NA, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR),
                     Notes = "Not enough annual measurements")
  }
  nwis_new <- rbind(nwis_new, df)
}
nwis_new
nwis_new$SITENO %in% annWL_df_NWIS$SITENO

# Append to other trend data
annWL_df_NWIS <- rbind(annWL_df_NWIS, nwis_new)

#----------------------------------
# Optional, only get points for sites with calculated statistics
#annWL_df_NWIS_stats <- annWL_df_NWIS[!is.na(annWL_df_NWIS$Sens_Slope),]


# Get sites as spatialpointsdf
# Use for-loop - too many to call all at once
site_loc1 <- readNWISsite(annWL_df_NWIS$SITENO[1:5000])
site_loc2 <- readNWISsite(annWL_df_NWIS$SITENO[5001:10000])
site_loc3 <- readNWISsite(annWL_df_NWIS$SITENO[10001:nrow(annWL_df_NWIS)])

test <- rbind(site_loc1, site_loc2, site_loc3)

# # Some sites missing from NWIS database - new ones established by CNRWA?
# which(annWL_df_NWIS$SITENO %notin% test$site_no)
# annWL_df_NWIS[c(12382, 12383, 12384,12385),]
# 
# # Create spdf for these points
# cn_data <- annWL_df_NWIS[c(12382, 12383, 12384,12385),]
# readNWISsite(cn_data$SITENO)
# cn_data$NDWR_Sttn_Nm <- c("078 N27 E27 13AADD1", "178A N30 E63 31BAAA1    ITCAINA WELL",
#                           "174 N14 E60 04DBD 1", "143 S02 E40 10CC 1")
# cn_pts <- data.frame(SITENO = as.character(c(40125111832701, 40262511453001, 390612115134801, 374633117320201)),
#                      lat_va = c(401250.5, 402632.7, 390611.8, 374653.7),
#                      long_va = c(1185326.7, 1145441.7, 1151348.4, 1173202.0))
# 
# 
# cn_pts <- merge(cn_pts, cn_data, by.x='SITENO', by.y='SITENO')
# cn_pts
# 
# 
# cn_sp <- SpatialPointsDataFrame(coords = data.frame(cn_pts$long_va, cn_pts$lat_va),
#                                     data = cn_pts)
# crs(cn_sp) = crs('+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs')
# # Somethings up with the ymin/amax (latitude)
# 
# #x <- spTransform(cn_sp, "+proj=longlat +datum=WGS84")
# x <- spTransform(cn_sp, CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))
# lonlat <- geom(x)[, c("x", "y")]
# head(lonlat, 3)
# plot(x, add=T, col="red")


out_sites <- SpatialPointsDataFrame(coords = data.frame(test$dec_long_va, test$dec_lat_va),
                                    data = test)


crs(out_sites) <- crs(ha)
plot(ha)
plot(out_sites, add=T, col="forestgreen")

# Attribute sites with stats
sites_out <- merge(out_sites, y = annWL_df_NWIS, by.x = "site_no", by.y = "SITENO", all.x= TRUE)
#sites_out <- merge(nwis_sites, y = annWL_df_NWIS_stats, by.x = "site_no", by.y = "SITENO", all.x= TRUE)
nwis_pts <- sites_out


#---------------------------------------------------------------------
#---------------------------------------------------------------------

# Get NDWR data from well level database
filename <- "D:/data_temp/gde_assessment/gde_threats/wellnet_gwlevels_092321.csv"
lines <- readLines(filename)
lines <- gsub('([^,])"([^,])', '\\1""\\2', lines)
data <- read.csv(textConnection(lines))
wellobs <- read.csv(filename)

dim(wellobs)
head(wellobs)

# Unique site names
site_names = as.vector(unique(wellobs$Site_Name))

# Fix site names to remove double-space
site_names_fix <- sub(pattern = "  ", " ", site_names)
head(site_names_fix)
wellobs$Site_Name_Fix <- sub("  ", " ", wellobs$Site_Name)
head(wellobs)
"117 S01 E35 09CC  1" %in% site_names_fix
site_names = as.vector(unique(wellobs$Site_Name_Fix))

# Read in NDWR site locations
gdb <- "D:/data_temp/gde_assessment/gde_threats/gde_threats.gdb"
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(gdb)
ndwr <- readOGR(dsn=gdb,layer="wellnet_gwsites_092321")
# Fix site name again - remove double space
ndwr$Site_Name_Fix <- sub(pattern = "  ", " ", ndwr$Site_Name)
head(ndwr)

annWL <- function(site_name){
  x <- wellobs[wellobs$Site_Name_Fix==site_name,]
  if (nrow(x) < 1){
    print(cat(site_name,"does not have enough annual readings to calculate trend"))
    df <- data.frame(SITENO = site_name, New_Pval = NA, Old_Pval = NA,
                     n_effective = NA, n_raw = NA,
                     Sens_Slope = NA, MinYear = NA, MaxYear = NA,
                     Notes = "Not enough annual measurements")
  } #else if (is.na(sum(x$Water_Level) == TRUE)){
  else if ((sum(is.na(x$Water_Level)) == length(x$Water_Level))) {
    print(cat(site_name,"only has NULL readings"))
    df <- data.frame(SITENO = site_name, New_Pval = NA, Old_Pval = NA,
                     n_effective = NA, n_raw = NA,
                     Sens_Slope = NA, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR),
                     Notes = "Only NULL measurement values")
  }
  else {
    print(cat(site_name,"gets a trend calculation..."))
    #just_dates <- unlist(strsplit(as.character(x$Measure_date), "[ ]"))
    #x$Date_char <- just_dates[-grep(":", just_dates)] # Removes time stamps
    ##x$Date_char <- just_dates[!just_dates == "0:00:00"]
    just_dates <- as.Date(as.character(x$Measure_date), "%m/%d/%Y")
    x$Date_char <- just_dates
    x <- x %>% mutate(DATE = as.Date(Date_char, format = "%m/%d/%Y"))
    annavg <- aggregate(Water_Level ~ cut(DATE, "1 year"), x, mean)
    annavg$YEAR <- as.numeric(substr(annavg$`cut(DATE, "1 year")`, 1, 4))
    if (var(annavg$Water_Level) > 0 & nrow(annavg) >= 5){
      modMK <- data.frame(t(mmkh(annavg$Water_Level)))
      df <- data.frame(SITENO = site_name, New_Pval = modMK$new.P.value, Old_Pval = modMK$old.P.value,
                       n_effective = modMK$N.N., n_raw = nrow(annavg[which(!is.na(annavg$lev_va)),]),
                       Sens_Slope = modMK$Sen.s.slope, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), Notes = NA)
    } else if (var(annavg$Water_Level) == 0 & nrow(annavg) >= 5){
      print(cat(site_name,"Only has one water level value - cannot run modified mann-kendall"))
      df <- data.frame(SITENO = site_name, New_Pval = NA, Old_Pval = NA,
                       n_effective = NA, n_raw = nrow(annavg[which(!is.na(annavg$Water_Level)),]),
                       Sens_Slope = 0, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), 
                       Notes = "Only one water level value; trend == 0")
    } 
    else {
      print(cat(site_name,"does not have enough annual readings to calculate trend"))
      df <- data.frame(SITENO = site_name, New_Pval = NA, Old_Pval = NA,
                       n_effective = NA, n_raw = NA,
                       Sens_Slope = NA, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR),
                       Notes = "Not enough annual measurements")
    }
    return(df)}
}
test <- annWL(site_name = "137B N09 E43 03AAAD1")
test <- annWL(site_name = "143 S02 E40 10CC  1")


#----------------------------------
# Calculate trends and format data

annWL_df_NDWR <- data.frame(SITENO = as.character(), New_Pval = as.numeric(), Old_Pval = as.numeric(),
                       n_effective = as.numeric(), n_raw = as.numeric(), Sens_Slope = as.numeric(),
                       MinYear = as.numeric(), MaxYear = as.numeric(), Notes = as.character())
site_names <- sort(site_names)

for (i in seq(1, length(site_names))){
  sitename <- site_names[i]
  tss_df <- annWL(site_name = sitename)
  annWL_df_NDWR <- rbind(annWL_df_NDWR, tss_df)
}

#----------------------------------
# Append CNRWA data to NDWR sites
# Central NV Regional Water Authority has additional data not-yet integrated in NWIS
cnrwa <- read.csv("E:\\RCF Data Recovery\\Recovered Files\\P00\\GDE_Threats\\Hydrology\\cnrwa_add_wells.csv")
head(cnrwa)
unique(cnrwa$site_no)
cnrwa <- cnrwa %>% mutate(DATE = as.Date(Date, format = "%m/%d/%Y"))
cnrwa$YEAR <- as.numeric(substr(cnrwa$DATE, 1, 4)) 
head(cnrwa)

# See which sites have additional data from the 2021 CNRWA report
# Sites that already exist in the trend calculations
# Calculate new trend with additional data and REPLACE
check_ndwr <- as.vector(unique(cnrwa[which(cnrwa$ndwr_site %in% annWL_df_NDWR$SITENO),]$ndwr_site))
check_ndwr %in% annWL_df_NDWR$SITENO
check_ndwr <- check_ndwr[-c(which(check_ndwr == ""))]
check_ndwr


new_trend2 <- function(site_name){
  og_site_trend <- annWL_df_NDWR[annWL_df_NDWR$SITENO==site_name,]
  cnrwa_site_data <- cnrwa[cnrwa$ndwr_site==site_name,]
  cnrwa_site_data$lvl_below_surf <- -1 * cnrwa_site_data$lvl_below_surf
  #print(og_site_trend$MaxYear)
  #print(max(cnrwa_site_data$YEAR))
  if (og_site_trend$MaxYear < max(cnrwa_site_data$YEAR)){
    print(paste("Re-calculate and replace site trend", site_name, sep=": "))
    ndwr_lvls <- wellobs[wellobs$Site_Name_Fix==site_name,]
    just_dates <- as.Date(as.character(ndwr_lvls$Measure_date), "%m/%d/%Y")
    ndwr_lvls$Date_char <- just_dates
    ndwr_lvls <- ndwr_lvls %>% mutate(DATE = as.Date(Date_char, format = "%m/%d/%Y"))
    append_lvls <- full_join(ndwr_lvls, cnrwa_site_data, by=c("Site_Name" = "site_no",
                                                              "Water_Level" = "lvl_below_surf", 
                                                              "DATE" = "DATE"))
    annavg <- aggregate(Water_Level ~ cut(DATE, "1 year"), append_lvls, mean)
    annavg$YEAR <- as.numeric(substr(annavg$`cut(DATE, "1 year")`, 1, 4))
    if (var(annavg$Water_Level) > 0 & nrow(annavg) >= 5){
      modMK <- data.frame(t(mmkh(annavg$Water_Level)))
      df <- data.frame(SITENO = site_name, New_Pval = modMK$new.P.value, Old_Pval = modMK$old.P.value,
                       n_effective = modMK$N.N., n_raw = nrow(annavg[which(!is.na(annavg$lev_va)),]),
                       Sens_Slope = modMK$Sen.s.slope, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR), 
                       Notes = "NDWR groundwater levels supplemented by 2021 CNRWA report data")
    }
    else {
      print(cat(site_name,"does not have enough annual readings to calculate trend"))
      df <- data.frame(SITENO = site_name, New_Pval = NA, Old_Pval = NA,
                       n_effective = NA, n_raw = NA,
                       Sens_Slope = NA, MinYear = min(annavg$YEAR), MaxYear = max(annavg$YEAR),
                       Notes = "Not enough annual measurements")
    } 
    return(df)
  }     
}

ndwr_replace <- data.frame(SITENO = as.character(), New_Pval = as.numeric(), Old_Pval = as.numeric(),
                           n_effective = as.numeric(), n_raw = as.numeric(), Sens_Slope = as.numeric(),
                           MinYear = as.numeric(), MaxYear = as.numeric(), Notes = as.character())
for(i in seq(1, length(check_ndwr))){
  s <- check_ndwr[i]
  x <- new_trend2(s)
  ndwr_replace <- rbind(ndwr_replace, x)
}
ndwr_replace

for(i in seq(1, nrow(ndwr_replace))){
  x <- as.vector(ndwr_replace$SITENO[i])
  print(x)
  ndwr_index <- which(annWL_df_NDWR$SITENO==x)
  print(ndwr_index)
  annWL_df_NDWR[ndwr_index,] <- ndwr_replace[i,]
}
annWL_df_NDWR[4436,]


#----------------------------------
# Get sites as spatialpoints df
head(annWL_df_NDWR)
#sites_with_stats <- annWL_df_NDWR[!is.na(annWL_df_NDWR$Sens_Slope),] # Only grab sites that had slope stats calculated above
#ndwr_stats <- ndwr[ndwr$Site_Name_Fix %in% sites_with_stats$SITENO,] # Get spatial points df of those sites with stats

sites_out <- merge(ndwr, y = annWL_df_NDWR, by.x = "Site_Name_Fix", by.y = "SITENO", all.x= TRUE)
sites_out
sites_out[sites_out$Site_Name_Fix=="209 S04 E61 28CD  1",]@data
crs(sites_out)
sites_out
ctz_data <- sites_out[5324,] # New site at Cortex gold mine - only 1 year of data
sites_out <- sites_out[-c(5324),]

#ndwr_pts <- spTransform(sites_out, CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))
ndwr_pts <- spTransform(sites_out, CRS("+proj=longlat +datum=NAD83 +no_defs"))
plot(ndwr_pts, add=TRUE, col="purple")


ctz <- ndwr[ndwr$Site_Name_Fix=="054 N27 E47 08BBDA2",]@data
ctz_sp <- SpatialPointsDataFrame(coords = data.frame(-1*ctz$Lon_DD_NAD83, ctz$Lat_DD_NAD83), 
                                 proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs"),
                                 data = ctz)
ctz_full <- merge(ctz_sp, ctz_data, by.x='Site_Name_Fix', by.y='Site_Name_Fix', all.y=TRUE)
plot(ctz_full, add=T, col="blue")
ctz_full@data
ctz_full <- ctz_full[,-(25:26)]

ndwr_pts_all <- rbind(ndwr_pts, ctz_full)
ndwr_pts <- ndwr_pts_all

#---------------------------------------------------------------------
#---------------------------------------------------------------------

# Combine NDWR and NWIS site data
ndwr <- ndwr_pts
nwis <- nwis_pts

# Remove duplicate sites -
# Some NWIS sites already in NDWR Wellnet database - going to remove these from NWIS before combining
# Create new columns with no whitespace in NDWR site names
ndwr$sttn_nm2 <- gsub(" ", "", ndwr$Site_Name_Fix, fixed=TRUE)
nwis$sttn_nm2 <- gsub(" ", "", nwis$station_nm, fixed=TRUE)

dupes <- ndwr[which(ndwr$sttn_nm2 %in% nwis$sttn_nm2),]
plot(dupes, add=T, col="orange")

test <- nwis[grep(dupes$sttn_nm2[5], nwis$sttn_nm2, ignore.case=TRUE),]
x <- ndwr[ndwr$sttn_nm2 == test$sttn_nm2,]@data # NDWR has more data for this well
y <- nwis[nwis$sttn_nm2 == test$sttn_nm2,]@data # NWIS record should be scrapped

nwis_rem <- c()
for (s in seq(1, length(dupes))){
  x <- dupes[s,]
  print(x$sttn_nm2)
  nwis_id <- grep(x$sttn_nm2, nwis$sttn_nm2, ignore.case=TRUE)
  print(nwis[nwis_id,]$sttn_nm)
  nwis_rem <- c(nwis_rem, nwis_id)
}

nwis_nodupes <- nwis[-c(nwis_rem),]
nwis_nodupes
nrow(nwis_pts)

#----------------------------------
# Combine spatial object with stats
# ONLY SITES WITH CALCULATED TREND
# Does not rectify ALL columns, just the important ones for now...
combo <- full_join(nwis_nodupes@data, ndwr@data, by=c("site_no" = "Site_Name_Fix", "station_nm" = "Site_Name", 
                                                      "dec_lat_va" = "Lat_DD_NAD83", "dec_long_va" = "Lon_DD_NAD83", 
                                                      "New_Pval" = "New_Pval", "Old_Pval" = "Old_Pval",
                                                      "n_effective" = "n_effective", "n_raw" = "n_raw", "Sens_Slope" = "Sens_Slope",
                                                      "MinYear" = "MinYear", "MaxYear" = "MaxYear", "Notes" = "Notes"))
combosp <- SpatialPointsDataFrame(coords = data.frame(combo$dec_long_va, combo$dec_lat_va), combo)
crs(combosp) <- crs(nwis)
plot(combosp, add=T, col="blue")
# CNRWA data added ~20 data points and updated others

# Filter to p-value <= 0.05 and negative trend value (significant falling)
# Or at least add as a marker/attribute
combosp$SigFall <- 0
sig_falling <- which(combosp$New_Pval <= 0.05 & combosp$Sens_Slope < 0)
combosp$SigFall[c(sig_falling)] <- 1

# Value of 100 where trend is significantly rising
sig_rising <- which(combosp$New_Pval <= 0.05 & combosp$Sens_Slope > 0)
combosp$SigFall[c(sig_rising)] <- 100
head(combosp)

# Assign HA name and ID to wells
crs(combosp) <- crs(ha)
test <- over(combosp, ha[,c("HYD_AREA", "HYD_AREA_N")])
head(test)
dim(test)
dim(combosp@data)

combosp$HYD_AREA <- test$HYD_AREA
combosp$HYD_AREA_NAME <- test$HYD_AREA_N
head(combosp)


# Export shapefile
writeOGR(combosp, 
         dsn="E:/RCF Data Recovery/Recovered Files/P00/GDE_Threats/Hydrology", 
         layer="gwlevels_stats_110321", 
         driver="ESRI Shapefile", overwrite_layer = TRUE)

# Export CSV
write.csv(combosp@data, "E:/RCF Data Recovery/Recovered Files/P00/GDE_Threats/Hydrology/gwlevels_stats_110321.csv")

#-------------------------------------------------------------
#-------------------------------------------------------------
# Basin summary statistics

# Basin names and ids with:
# number of wells that meet trend criteria ("wells that we looked at")
# number of significantly falling gw wells
# most recent year of gw level observation


x <- read.csv("E:/RCF Data Recovery/Recovered Files/P00/GDE_Threats/Hydrology/gwlevels_stats_110321.csv")
colnames(x)
x$ID <- as.character(x$X)

test <- x %>% group_by(SigFall) %>% summarise(NumSites = n_distinct(ID))
test$PropTrend <- test$NumSites/sum(test$NumSites)

# Number of sites in HA
test <- data.frame(x %>% group_by(HYD_AREA) %>% summarise(count = n_distinct(ID)))

# Number of sites in HA with info on whether trend is significantly rising or falling
test2 <- data.frame(x %>% group_by(HYD_AREA, SigFall) %>% summarise(count = n_distinct(ID)))


# Collapse to one set of attributes for HA
# Number of sites
# Number of sig falling sites
# Most recent year of observation

x1 <- data.frame(x %>% group_by(HYD_AREA) %>% summarise(SiteDataCount=n_distinct(ID)))
x2 <- data.frame(x %>% filter(SigFall==1) %>% group_by(HYD_AREA) %>% summarise(FallCount=n_distinct(ID)))
x3 <- data.frame(x %>% filter(SigFall==0) %>% group_by(HYD_AREA) %>% summarise(NoTrendCount=n_distinct(ID)))
x4 <- data.frame(x %>% filter(SigFall==100) %>% group_by(HYD_AREA) %>% summarise(RiseCount=n_distinct(ID)))

xx <- left_join(x1, x2, by="HYD_AREA")
xx <- left_join(xx, x3, by="HYD_AREA")
xx <- left_join(xx, x4, by="HYD_AREA")
xx <- left_join(xx, ha_year, by="HYD_AREA")
head(xx)

# Read in HA data
gdb <- "E:\\RCF Data Recovery\\Recovered Files\\P00\\GDE_Threats\\Maps\\GDE_Threats.gdb"
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(gdb)
ha <- readOGR(dsn=gdb,layer="hydrographic_basin_boundaries")
#ha <- spTransform(ha, CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs "))
plot(ha)
head(ha)

ha@data <- left_join(ha@data, xx, by="HYD_AREA")
head(ha@data)

# Replace NAs with zeroes where applicable
ha@data[is.na(ha$SiteDataCount),]$SiteDataCount <- 0
ha@data[is.na(ha$FallCount),]$FallCount <- 0
ha@data[is.na(ha$NoTrendCount),]$NoTrendCount <- 0
ha@data[is.na(ha$RiseCount),]$RiseCount <- 0
head(ha)

ha@data[is.na(ha$RecentYear),]$RecentYear <- "No sites"
ha@data[ha$RecentYear==-Inf,]$RecentYear <- "No data recorded for sites"
head(ha)

# Calculate proportion of wells that have falling levels for each HA
ha@data$PropFall <- ha$FallCount/ha$SiteDataCount
ha@data[is.na(ha$PropFall),]$PropFall <- -9999

# write to shapefile
writeOGR(ha, 
         dsn="E:/RCF Data Recovery/Recovered Files/P00/GDE_Threats/Hydrology", 
         layer="hydrographic_area_gwstats", 
         driver="ESRI Shapefile", overwrite_layer = TRUE)

# Write to csv
ha@data[ha$PropFall==-9999,]$PropFall <- NA
write.csv(ha@data, "E:/RCF Data Recovery/Recovered Files/P00/GDE_Threats/Hydrology/hydrographic_area_gwstats.csv")


#-------------------------------------------------------------
#-------------------------------------------------------------
# Example gw level time series showing trends

# NWIS - 364329116402902
# Amargosa basin 1997 - 2021
s <- readNWISgwl(364329116402902)
s$lev_va <- -1 * s$lev_va # Reverse groundwater values so sens-slope trend is negative where depth to gw is increasing
s <- s %>% mutate(DATE = as.Date(lev_dt, format = "%Y-%m-%d"))
s$YEAR <- as.numeric(substr(s$`DATE`, 1, 4))
annavg <- aggregate(lev_va ~ cut(DATE, "1 year"), s, mean)
annavg$YEAR <- as.numeric(substr(annavg$`cut(DATE, "1 year")`, 1, 4))


# GGplot - raw data/ annual measurements / trend line
p <- ggplot(data=s, aes(x=DATE, y=lev_va)) + geom_line(aes(group=site_no), color="#616161") +
  geom_point(data=annavg, aes(x=as.Date(`cut(DATE, "1 year")`), y= lev_va), color="red") +
  geom_line(data=annavg, aes(x=as.Date(`cut(DATE, "1 year")`), y= lev_va), color="red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(color = "black", size = 6, margin=margin(0,0,0,0)),
        plot.subtitle = element_text(size = 5, vjust = -1)) +
  theme_bw() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "USGS Groundwater Well - 364329116402902",
       subtitle = "Amargosa Valley, NV")
p




brteplot <- function(in_df, brte_no){
  x_sub <- subset(in_df, BRTE > brte_no)
  xagg <- aggregate(NDVI ~ Date, data = x_sub, mean, na.rm=TRUE)
  xagg$PLOTNAME <- 0
  p <- ggplot(data = x_sub, aes(x=Date, y=NDVI)) + geom_line(aes(group=PLOTNAME)) + 
    geom_line(data=xagg, aes(x=Date, y=NDVI, group=PLOTNAME), color="red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(color = "black", size = 6, margin=margin(0,0,0,0)),
          plot.subtitle = element_text(size = 5, vjust = -1)) + 
    labs(title = as.character(x_sub$BRTE[1]),
         subtitle = "WY2015 - WY2019 NDVI")
  return(list(plot = p, sites = unique(unlist(x_sub$PLOTNAME)), all_ts = x_sub))
}
es <- data.frame(xmin = as.Date(c("2015-03-25", "2016-03-25", "2017-03-25", "2018-03-25", "2019-03-25")),
                 xmax = as.Date(c("2015-04-25", "2016-04-25", "2017-04-25", "2018-04-25", "2019-04-25")),
                 ymin = -Inf, ymax = Inf)
ms <- data.frame(xmin = as.Date(c("2015-06-15", "2016-06-15", "2017-06-15", "2018-06-15", "2019-06-15")),
                 xmax = as.Date(c("2015-07-15", "2016-07-15", "2017-07-15", "2018-07-15", "2019-07-15")),
                 ymin = -Inf, ymax = Inf)
p <- ggplot() + 
  geom_rect(data = es, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#59a4ff", alpha = 0.3) +
  geom_rect(data = ms, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5) +
  geom_line(data = combo, aes(x=date, y=mean_ndvi, color=SYSXCLA)) + 
  scale_colour_manual(values=c("#0ea800", "#c76000")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(color = "black", size = 10, margin=margin(0,0,0,0)),
        plot.subtitle = element_text(size = 9, vjust = -1),
        legend.direction = "horizontal", 
        legend.position = "bottom") +
  labs(title = "Big Sagebrush - upland with trees",
       subtitle = "2013 - 2019 Water Year NDVI")
p



# NDWR - 153 N21 E53 03BBDD2
# http://water.nv.gov/WaterLevelDataChart.aspx?autoid=1112
# by Winnemucca 1997 - 2021
# Another option: http://water.nv.gov/WaterLevelDataChart.aspx?autoid=931



