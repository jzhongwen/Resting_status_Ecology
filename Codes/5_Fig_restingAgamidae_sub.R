#######analysis of global agamidae resting status#########
#######Zhong-wen Jiang 2-Mar 2023########################
library(raster)
library(dplyr)
library(rgeos)
library(sf)
library(terra)
library(maps) 
library(maptools)
library(ggplot2)

site_all<-read.csv("/Users/jiangzhongwen/Desktop/Resting ecology/Data/highele_region.csv") 


allyear_data<-data.frame()
for (year in 1986:2015) {
  oneyear<-read.csv(paste0("/Users/jiangzhongwen/Desktop/Resting ecology/Data/global_agamidae/",year,"_physiology.csv"))
  allyear_data<-rbind(oneyear,allyear_data)
}

ecto_result1<-data.frame()

for (i in 1:nrow(site_all)) {
  i_result_year<-site_all[i,]
  oneyear_data<-allyear_data[allyear_data$LON==i_result_year$lon&allyear_data$LAT==i_result_year$lat,]
  oneyear_data<-na.omit(oneyear_data)

  oneyear_act<-data.frame(LON=i_result_year$lon,LAT=i_result_year$lat,
                          FP=mean(oneyear_data[oneyear_data$STATUS=="ACT",]$FP),
                          MR=mean(oneyear_data[oneyear_data$STATUS=="ACT",]$MR),
                          STATUS="ACT",N=nrow(oneyear_data[oneyear_data$STATUS=="ACT",]))
  
  oneyear_rest<-data.frame(LON=i_result_year$lon,LAT=i_result_year$lat,
                          FP=mean(oneyear_data[oneyear_data$STATUS=="REST",]$FP),
                          MR=mean(oneyear_data[oneyear_data$STATUS=="REST",]$MR),
                          STATUS="REST",N=nrow(oneyear_data[oneyear_data$STATUS=="REST",]))
  
  oneyear_result<-rbind(oneyear_act,oneyear_rest)
  oneyear_result["ID"]<-site_all[i,]$ID
  oneyear_result["ELE"]<-site_all[i,]$ele
  ecto_result1<-rbind(ecto_result1,oneyear_result)
  
}
ecto_result2<-ecto_result1[ecto_result1$N!=0,]



dif_result<-data.frame()

for(i in unique(ecto_result2$ID)){
  ID_data<-ecto_result2[ecto_result2$ID==i,]
  FP_differ_rest<-mean(ID_data[ID_data$ELE==min(ID_data$ELE)&ID_data$STATUS=="REST",]$FP)-mean(ID_data[ID_data$ELE==max(ID_data$ELE)&ID_data$STATUS=="REST",]$FP)
  FP_differ_act<-mean(ID_data[ID_data$ELE==min(ID_data$ELE)&ID_data$STATUS=="ACT",]$FP)-mean(ID_data[ID_data$ELE==max(ID_data$ELE)&ID_data$STATUS=="ACT",]$FP)
  FP_differ_all<-FP_differ_rest+FP_differ_act
#abs(FP_contri_rest)>1,rest have more contributionn than act, and if FP_contri_rest<-0,here rest and act have different pattern
  FP_contri_rest<-FP_differ_rest/FP_differ_act
  
  MR_differ_rest<-mean(ID_data[ID_data$ELE==min(ID_data$ELE)&ID_data$STATUS=="REST",]$MR)-mean(ID_data[ID_data$ELE==max(ID_data$ELE)&ID_data$STATUS=="REST",]$MR)
  MR_differ_act<-mean(ID_data[ID_data$ELE==min(ID_data$ELE)&ID_data$STATUS=="ACT",]$MR)-mean(ID_data[ID_data$ELE==max(ID_data$ELE)&ID_data$STATUS=="ACT",]$MR)
  MR_differ_all<-MR_differ_rest+MR_differ_act
  MR_contri_rest<-MR_differ_rest/MR_differ_act
  
  LON<-mean(ID_data$LON)
  LAT<-mean(ID_data$LAT)
  
  ID_result<-data.frame(ID=i,LON=LON,LAT=LAT,
                        dif=c(FP_differ_rest,FP_differ_act,FP_differ_all,FP_contri_rest,MR_differ_rest,MR_differ_act,MR_differ_all,MR_contri_rest),
                        STATUS=c("REST","ACT","ALL","REST_contri","REST","ACT","ALL","REST_contri"),
                        Traits=c("FP","FP","FP","FP","MR","MR","MR","MR"),
                        ELE_dif=mean(max(ID_data$ELE))-mean(min(ID_data$ELE)))
  dif_result<-rbind(dif_result,ID_result)
}



Rest_contri_FP<-dif_result[dif_result$Traits=="FP"&dif_result$STATUS=="REST_contri",]
Rest_contri_MR<-dif_result[dif_result$Traits=="MR"&dif_result$STATUS=="REST_contri",]
Rest_contri_FP<-na.omit(Rest_contri_FP)
Rest_contri_MR<-na.omit(Rest_contri_MR)

Rest_FP_Category<-data.frame()
for (i in 1:nrow(Rest_contri_FP)) {
  i_Category<-Rest_contri_FP[i,]
  if(i_Category$dif>1){
    i_Category["Category"]="A"
  }else{
    if(i_Category$dif<=1&i_Category$dif>0){
      i_Category["Category"]="C"
    }else {
      if(i_Category$dif<=0&abs(i_Category$dif)<1){
        i_Category["Category"]="D"
      }else {
        i_Category["Category"]="B"
      }
    }
  } 
  Rest_FP_Category<-rbind(Rest_FP_Category,i_Category)
}


Rest_MR_Category<-data.frame()
for (i in 1:nrow(Rest_contri_MR)) {
  i_Category<-Rest_contri_MR[i,]
  if(i_Category$dif>1){
    i_Category["Category"]="A"
  }else{
    if(i_Category$dif<=1&i_Category$dif>0){
      i_Category["Category"]="C"
    }else {
      if(i_Category$dif<=0&abs(i_Category$dif)<1){
        i_Category["Category"]="D"
      }else {
        i_Category["Category"]="B"
      }
    }
  } 
  Rest_MR_Category<-rbind(Rest_MR_Category,i_Category)
}


write.csv(Rest_FP_Category,file ="/Users/jiangzhongwen/Desktop/Resting ecology/Data/Rest_FP_Category.csv",row.names = FALSE)
write.csv(Rest_MR_Category,file ="/Users/jiangzhongwen/Desktop/Resting ecology/Data/Rest_MR_Category.csv",row.names = FALSE)

Rest_FP_Category2<-Rest_FP_Category[Rest_FP_Category$ELE_dif>500,]
Rest_MR_Category2<-Rest_MR_Category[Rest_MR_Category$ELE_dif>500,]


ext<-extent(min(Rest_FP_Category2$LON),max(Rest_FP_Category2$LON),min(Rest_FP_Category2$LAT),max(Rest_FP_Category2$LAT))
grid <- st_bbox(ext) %>% 
  st_make_grid(cellsize = 2, what = "polygons") %>%
  st_set_crs("+proj=longlat +datum=WGS84")
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

data_fp<-Rest_FP_Category2[,c(2,3,8)]
data_fp <- data_fp %>% sf::st_as_sf(coords = c("LON","LAT"), crs = "+proj=longlat +datum=WGS84")
grid_fp<-grid %>% 
  sf::st_join(data_fp, left = FALSE)

data_mr<-Rest_MR_Category2[,c(2,3,8)]
data_mr <- data_mr %>% sf::st_as_sf(coords = c("LON","LAT"), crs = "+proj=longlat +datum=WGS84")
grid_mr<-grid %>% 
  sf::st_join(data_mr, left = FALSE)

world <- st_as_sf(map('world', plot = FALSE, fill = TRUE))

Figure_fp<-ggplot()+
  geom_sf(data = world,color="grey90",fill="grey90")+
  geom_sf(data = grid_fp[2],aes(fill=Category),color=NA)+
  scale_fill_manual(values = c("#437199","#5ED1BD","#CB4B58","#FDAE61"),name="Absolute difference in cumulative digestion",
                    labels=c(paste0("Rest>Activity (same sign;",round(nrow(grid_fp[grid_fp$Category=="A",])/nrow(grid_fp),2)*100,"%)"),
                    paste0("Rest>Activity (opposite signs;",round(nrow(grid_fp[grid_fp$Category=="B",])/nrow(grid_fp),2)*100,"%)"),
                    paste0("Rest<Activity (same sign;",round(nrow(grid_fp[grid_fp$Category=="C",])/nrow(grid_fp),2)*100,"%)"),
                    paste0("Rest<Activity (opposite signs;",round(nrow(grid_fp[grid_fp$Category=="D",])/nrow(grid_fp),2)*100,"%)")))+
xlim(-20,150)+
  ylim(-40,60)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


Figure_mr<-ggplot()+
  geom_sf(data = world,color="grey90",fill="grey90")+
  geom_sf(data = grid_mr[2],aes(fill=Category),color=NA)+
  scale_fill_manual(values = c("#437199","#5ED1BD","#CB4B58","#FDAE61"),name="Absolute difference in cumulative metabolism",
                    labels=c(paste0("Rest>Activity (same sign;",round(nrow(grid_mr[grid_mr$Category=="A",])/nrow(grid_mr),2)*100,"%)"),
                             paste0("Rest>Activity (opposite signs;",round(nrow(grid_mr[grid_mr$Category=="B",])/nrow(grid_mr),2)*100,"%)"),
                             paste0("Rest<Activity (same sign;",round(nrow(grid_mr[grid_mr$Category=="C",])/nrow(grid_mr),2)*100,"%)"),
                             paste0("Rest<Activity (opposite signs;",round(nrow(grid_mr[grid_mr$Category=="D",])/nrow(grid_mr),2)*100,"%)")))+
  
  xlim(-20,150)+
  ylim(-40,60)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  

library(ggpubr)

pdf("/Users/jiangzhongwen/Desktop/Resting ecology/Figure/Figure5.pdf",height=15,width=18)

ggarrange(Figure_fp,Figure_mr,nrow=2,labels=c('a','b'),
          font.label = list(size =30, color = "black",family="serif"))

dev.off()





########extract the environmental data to analyze the reason

Category_point_fp<-data.frame()

for (ID in unique(Rest_FP_Category$ID)) {
  ID_point<-ecto_result2[ecto_result2$ID==ID,]
  ID_point["Category"]<-Rest_FP_Category[Rest_FP_Category$ID==ID,]$Category
  Category_point_fp<-rbind(Category_point_fp,ID_point)
}

Category_point_fp<-Category_point_fp[Category_point_fp$STATUS=="ACT",][,c(1,2,7,8,9)]

Category_point_mr<-data.frame()

for (ID in unique(Rest_MR_Category$ID)) {
  ID_point<-ecto_result2[ecto_result2$ID==ID,]
  ID_point["Category"]<-Rest_MR_Category[Rest_MR_Category$ID==ID,]$Category
  Category_point_mr<-rbind(Category_point_mr,ID_point)
}

Category_point_mr<-Category_point_mr[Category_point_mr$STATUS=="ACT",][,c(1,2,7,8,9)]

write.csv(Category_point_fp,file ="/Users/jiangzhongwen/Desktop/Resting ecology/Data/Category_allpoint_fp.csv",row.names = FALSE)
write.csv(Category_point_mr,file ="/Users/jiangzhongwen/Desktop/Resting ecology/Data/Category_allpoint_mr.csv",row.names = FALSE)
