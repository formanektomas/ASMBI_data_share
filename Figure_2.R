###############################################################################
#
library(maptools)
library(ggplot2)
library(rgeos)
library(rgdal) 
library(RColorBrewer)
library(scales)
library(ggmap)
library(reshape2)
library(stringr)
library(SmarterPoland)
library(reshape) 
library(spdep)
library(car) 
library(lmtest)
library(tripack) 

###############################################################################
#
###############################################################################
#
# 
#
rm(list=ls())
## Download shapefiles if necessary
##
# download.file("http://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/download/ref-nuts-2013-60m.shp.zip", 
#               destfile = "NUTS_2013_60M_SH.zip")
## NUTS2016
## http://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/download/ref-nuts-2016-60m.shp.zip
## unzip to SpatialPolygonsDataFrame
# unzip("NUTS_2013_60M_SH.zip")
library(rgdal)
ogrListLayers(dsn = "./NUTS_2013_60M_SH")
# Projection systems:
#  EPSG:4258 (ETRS 1989, coordinates in decimal degrees),
#  EPSG:4326 (WGS84, coordinates in decimal degrees),
#  EPSG:3035 (ETRS 1989 in Lambert Azimutal projection with centre in E52N10, coordinates in meters),
#  EPSG:3857 (WGS84 Web Mercator Auxiliary Sphere, coordinates in meters)

map <- readOGR(dsn = "./NUTS_2013_60M_SH", layer = "NUTS_RG_60M_2013_4326")
#
# unique(map$NUTS_ID)
#
# polygon data at all NUTS levels.
AT.rows <- grep("AT", map$NUTS_ID, fixed = T)
CZ.rows <- grep("CZ", map$NUTS_ID, fixed = T)
DE.rows <- grep("DE", map$NUTS_ID, fixed = T)
HU.rows <- grep("HU", map$NUTS_ID, fixed = T)
PL.rows <- grep("PL", map$NUTS_ID, fixed = T)
SK.rows <- grep("SK", map$NUTS_ID, fixed = T)
BE.rows <- grep("BE", map$NUTS_ID, fixed = T) # Belgium
NL.rows <- grep("NL", map$NUTS_ID, fixed = T) # Neth
LU.rows <- grep("LU", map$NUTS_ID, fixed = T) # Lux
DK.rows <- grep("DK", map$NUTS_ID, fixed = T) # Denmark
SI.rows <- grep("SI", map$NUTS_ID, fixed = T) # Slovenia
#
# Sample plots for all selected countries
CE <- c(AT.rows, CZ.rows, SK.rows, PL.rows, HU.rows, DE.rows,BE.rows,NL.rows,LU.rows,DK.rows,SI.rows)
rm(AT.rows)
rm(CZ.rows)
rm(HU.rows)
rm(DE.rows)
rm(PL.rows)
rm(SK.rows)
rm(BE.rows)
rm(NL.rows)
rm(LU.rows)
rm(DK.rows)
rm(SI.rows)
#
CE_map <- map[CE,] # selected countries only
CE_map <- CE_map[CE_map$LEVL_CODE == 2,] # NUTS2 only
#
summary(CE_map)
CE_map$NUTS_ID
#
#
###########################################################
# Import data
CE_data <- read.csv("ASMBIdf.csv")
CE_data$ERR <- ERR*100 # for perc. points (0-100 scale)
CE_data <- CE_data[,c("geo","time","ERR")]
# CE_data$lY <- exp(CE_data$lY)
?cast
dat2 <- cast(CE_data, geo ~ time, value="ERR")
dat2 <- as.data.frame(dat2)
##########################################################
# 
row.names(dat2) <- dat2$geo
# rownames
dat2 <- dat2[order(row.names(dat2)), ]
row.names(CE_map) <- as.character(CE_map$NUTS_ID)
CE_map <- CE_map[order(row.names(CE_map)), ]
# join
shape <- spCbind(CE_map, dat2)
# fortify spatialpolygondataframe into data.frame
shape$id <- rownames(shape@data)
map.points <- fortify(shape, region = "id")
map.df <- merge(map.points, shape, by = "id")
# Make a top layer for "hole" areas (Bremen, Berlin, Prague, Wien)
map.df$TopL <- 0
map.df$TopL[map.df$id == "DE30"] <- 1 #Berlin
map.df$TopL[map.df$id == "DE50"] <- 1 #Bremen
map.df$TopL[map.df$id == "CZ01"] <- 1 #Prague
map.df$TopL[map.df$id == "AT13"] <- 1 #Wien
map.df.o <- map.df[order(map.df$TopL), ]
# Polygons for State boundaries (NUTS0)
CE_map.2 <- map[CE,]
CE_map.2 <- CE_map.2[CE_map.2$LEVL_CODE == 0,]
# 2 in 1 plot - GDP - preparation
map.df.l <- melt(data = map.df, id.vars = c("id", "long", "lat", "group", "TopL"), 
                 measure.vars = c("X2016"))
# year variable (variable) is class string and type X20xx.  Lets remove
# the X and convert it to numerical
map.df.l$variable <- str_replace_all(map.df.l$variable, "X", "")
map.df.l$variable <- factor(map.df.l$variable)
map.df.l$variable <- as.numeric(levels(map.df.l$variable))[map.df.l$variable]
#
ggplot(map.df.l, aes(long,lat,group=group)) +
  geom_polygon(aes(fill = value)) +
  geom_polygon(data = map.df.l[map.df.l$TopL == 1, ], aes(fill = value)) +
  # Top layer for "hole" regions
  scale_fill_gradientn('Relative \nshare \nRenewable \nEnergy \n2016', colours=brewer.pal(8, "RdYlGn"))+
  geom_polygon(data = map.df.l, aes(long,lat), 
               fill=NA, 
               color = "grey50",
               size=0.1) + # Grey NUTS2 borders
  geom_path(data=CE_map.2, color="black", size = 0.3)+ # NUTS0 borders
  coord_map() +
  theme_minimal()
#
#
#