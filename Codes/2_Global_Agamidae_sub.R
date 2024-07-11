#########The resting status contribution of Global Agamidae+ diurnal+terrestrial######
#########Zhong-wen Jiang 21 Feb 2023#############

library(rgdal)
library(sf)
library(terra)
library(sp)
library(raster)
library(dplyr)
library(rgeos)
library(bfsMaps)


Shai_database<-read.csv("/Users/jiangzhongwen/Dropbox/Egg-laying opportunity2.0/Data/Database from literatrue/Traits of lizard (mass, Tb)/Appendix S1 - Lizard data version 1.0.csv")

head(Shai_database)
table(Shai_database$Family)

Agamidae_database<-Shai_database[Shai_database$Family=="Agamidae",]
range_agamidae<-subset(Re_range, Binomial%in%Agamidae_database$Binomial)
range_merge_all<-CombinePolygons(range_agamidae,range_agamidae$Group)

Agamidae_diurnal<-Agamidae_database[Agamidae_database$Activity.time=="Diurnal",]
Agamidae_ter<-Agamidae_diurnal[Agamidae_diurnal$substrate=="Terrestrial",]
Agamidae_ter<-subset(Agamidae_ter,Binomial!="NA")


Re_range<-readOGR("/Volumes/My Passport/QT/GLOBAL reptile/reptile_dis/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
range_test<-subset(Re_range, Binomial%in%Agamidae_ter$Binomial)
ele_10m<-raster("/Volumes/My Passport/QT/Elevation/wc2.1_10m_elev.tif")
range_merge<-CombinePolygons(range_test,range_test$Group)

Agamidae_ele<-crop(ele_10m,range_test)
Agamidae_ele<-mask(ele_10m,range_test)
Agamidae_point<-data.frame(rasterToPoints(Agamidae_ele))
colnames(Agamidae_point)<-c("lon","lat","ele")

#write.csv(Agamidae_point,file = "/Users/jiangzhongwen/Desktop/Resting ecology/Data/global_agamidae.csv",row.names = FALSE)

Agamidae_point<-read.csv("/Users/jiangzhongwen/Desktop/Resting ecology/Data/global_agamidae.csv")

ID=1
high_elevation_point<-data.frame()

########2°x 2°#####
for (lon in seq(floor(min(Agamidae_point$lon)),(floor(max(Agamidae_point$lon))+1),by=2)) {
  lon_point<-Agamidae_point[Agamidae_point$lon>=lon&Agamidae_point$lon<lon+2,]
  
  if(nrow(lon_point)!=0){
    for (lat in seq(floor(min(lon_point$lat)),(floor(max(lon_point$lat))+1),by=2)) {
      lat_point<-lon_point[lon_point$lat>=lat&lon_point$lat<lat+2,]
      
      if(nrow(lat_point)>1){
        max_ele_point<-lat_point[lat_point$ele==max(lat_point$ele),]
        min_ele_point<-lat_point[lat_point$ele==min(lat_point$ele),]
        ele_point<-rbind(max_ele_point,min_ele_point)
        ele_point["ID"]<-ID
        ID=ID+1
        high_elevation_point<-rbind(high_elevation_point,ele_point)
      }
    }
  }
}

write.csv(high_elevation_point,file = "/Users/jiangzhongwen/Desktop/Resting ecology/Data/highele_region.csv",row.names = FALSE)
high_elevation_point<-read.csv("/Users/jiangzhongwen/Desktop/Resting ecology/Data/highele_region.csv")

