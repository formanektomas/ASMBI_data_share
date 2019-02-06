#
#
#
# 
#
#
rm(list=ls())
DF <- read.csv("ASMBIdf.csv")
DF <- DF[order(DF$time, DF$geo),]
#
# Connectivity matrix W (W = C1 below)
require(splm)
coords <- DF[DF$time==2010,c("geo","long", "lat")]
IDs <- coords$geo
coords <- coordinates(coords[,2:3])
# Conectivity
nb.km <- dnearneigh(coords, d1=0, d2=500, longlat=T, row.names = IDs) # alt d2=380
listw.km <- nb2listw(nb.km)
C1 <- nb2mat(nb.km,style = "B") #
C1[1:10,1:10]
C1 <-unname(C1)
require(spmoran)
# ?meigen
E1 <- meigen(cmat=C1) 
Eig.df <- as.data.frame(E1$sf)
#
###########################################################################
# Plot the eigenvectors... ################################################
###########################################################################

# Step 1 - Preparation of "mappning data"
#
#
library(rgeos) # install.packages("rgeos")
library(rgdal) # install.packages("rgdal")
# download.file("http://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/download/ref-nuts-2013-60m.shp.zip", 
#               destfile = "NUTS_2013_60M_SH.zip")
# unzip("NUTS_2013_60M_SH.zip")
ogrListLayers(dsn = "./NUTS_2013_60M_SH") # 
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
#plot(CE_map)
summary(CE_map)
CE_map$NUTS_ID

#
######################################################
#
#
# Step 3 - Merging spatial and macroeconomic data 
#
#
# 
Eig.df <- cbind(IDs, Eig.df)
# Now, we "join" both objects
CE_map <- CE_map[order(CE_map$NUTS_ID), ]
row.names(Eig.df) <- Eig.df$IDs
row.names(CE_map) <- as.character(CE_map$NUTS_ID)
library(maptools) # install.packages("maptools")
shape <- spCbind(CE_map, Eig.df)
# Finally, we need to use fortify() in order to transform 
# SpatialPolygonDataFrame into data.frame for use with ggplot2
library(ggplot2)
library(rgeos) #install.packages("rgeos")
shape$id <- rownames(shape@data)
map.points <- fortify(shape, region = "id")
map.df <- merge(map.points, shape, by = "id")
#
# To ensure correct plotting of "hole" NUTS2 regions, we have to
# prepare a top layer to the plot. Relevant for: Bremen, Berlin, Praha, Wien.
map.df$TopL <- 0
map.df$TopL[map.df$id == "DE30"] <- 1 #Berlin
map.df$TopL[map.df$id == "DE50"] <- 1 #Bremen
map.df$TopL[map.df$id == "CZ01"] <- 1 #Praha
map.df$TopL[map.df$id == "AT13"] <- 1 #Wien
map.df.o <- map.df[order(map.df$TopL), ]
#
# Polygons for State boundaries (NUTS0)
CE_map.2 <- map[CE,]
CE_map.2 <- CE_map.2[CE_map.2$LEVL_CODE == 0,]
#
######################################################
#
#
# Step 5
#
#
library(RColorBrewer) #install.packages("RColorBrewer")
library(scales) #install.packages("scales")
library(ggmap) #install.packages("ggmap")
library(reshape2) #install.packages("reshape2")
#
colnames(map.df)
map.df.l <- melt(data = map.df, id.vars = c("id", "long", "lat", "group", "TopL"), 
                 measure.vars = c("V1", "V2", "V3",
                                  "V4", "V5", "V6",
                                  "V7", "V8", "V9"))
library(stringr)
map.df.l$variable <- factor(map.df.l$variable)
map.df.l$group <- map.df.l$id

#
# Eig. plot - 
ggplot(map.df.l, aes(long,lat,group=group)) +
  geom_polygon(aes(fill = value)) +
  geom_polygon(data = map.df.l[map.df.l$TopL == 1, ], aes(fill = value)) +
  # Top layer for "hole" regions
  scale_fill_gradientn('Map of \neigenvectors \nV1 - V9', colours=brewer.pal(8, "Reds"))+
  geom_path(data = CE_map, color = "grey50", size=0.1, na.rm = T) + # Grey NUTS2 borders
  geom_path(data=CE_map.2, color="black", size = 0.4, na.rm = T)+ # NUTS0 borders
  coord_map() +
  facet_wrap(~variable) + 
  theme_minimal()
#








