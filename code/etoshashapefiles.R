
setwd('C:\\Users\\lindesay\\Documents\\EtoshaShapefiles')

require(rgdal)
roads<-readOGR( dsn='C:\\Users\\lindesay\\Documents\\EtoshaShapefiles\\Roads', layer="enp roads")
roads_fortify<-fortify(roads)
proj4string(roads) <- CRS("+proj=longlat +datum=WGS84")  ## for example
roadsUTM <- spTransform(roads, CRS("+proj=utm +zone=33 ellps=WGS84"))

water<-readOGR( dsn='C:\\Users\\lindesay\\Documents\\EtoshaShapefiles\\Water', layer="waters2")
proj4string(water) <- CRS("+proj=longlat +datum=WGS84")  ## for example
waterUTM <- spTransform(water, CRS("+proj=utm +zone=33 ellps=WGS84"))

waterh<-readOGR( dsn='C:\\Users\\lindesay\\Documents\\EtoshaShapefiles\\Water', layer="functional water")
proj4string(waterh) <- CRS("+proj=longlat +datum=WGS84")  ## for example
waterholes <- spTransform(waterh, CRS("+proj=utm +zone=33 ellps=WGS84"))


tracks<-readOGR( dsn='C:\\Users\\lindesay\\Documents\\EtoshaShapefiles\\Roads', layer="TRACK")
tracks_fortify<-fortify(tracks)
proj4string(tracks) <- CRS("+proj=longlat +datum=WGS84")  ## for example
tracksUTM <- spTransform(tracks, CRS("+proj=utm +zone=33 ellps=WGS84"))

veg<-readOGR( dsn='Vegetation', layer="VEG96")
#veg_fortify<-fortify(tracks)
proj4string(veg) <- CRS("+proj=longlat +datum=WGS84")  ## for example
vegUTM <- spTransform(veg, CRS("+proj=utm +zone=33 ellps=WGS84"))

mydat<-read.csv('../../Dropbox/Papers/papers/CarcassSALSApaper/data/datcomb.csv')


plot(roadsUTM)
plot(tracksUTM, add=T, col='lightblue')
points(mydat$x.pos[mydat$death==1], mydat$y.pos[mydat$death==1], pch=20, col='red')
points(mydat$x.pos[mydat$death==0], mydat$y.pos[mydat$death==0], pch=20, col='lightgrey')
plot(waterUTM, pch='x', col='blue', add=T)

plot(waterholes, pch='x', col='green', add=T)

rainfall<-read.csv('Rainfall_DailMonth_ENP 2017Complete.csv')[,1:11]
require(dplyr)

seasonalrainfall<-filter(rainfall, is.na(Year)) %>% mutate(yearid=1:n()) %>% filter(yearid>55)
