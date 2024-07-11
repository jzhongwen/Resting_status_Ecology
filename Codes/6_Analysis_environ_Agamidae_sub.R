########6. Analyze the environmental data###########

library(NicheMapR)
library(plyr)
library(parallel)
library(R.utils)

ncore=20

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

Category_point<-read.csv("/scratch/gpfs/liangm/Resting_status_Jiang/Data/Category_allpoint.csv")

####load a function for running microclimate model using downloaded TerraClim data
source("/scratch/gpfs/liangm/desert_bird/WD_planB/Customized functions/micro_terra_customized.R")

####Create a function for running the microclimate model for a site (lon, lat), a year and a scenario
cal.micro<-function(lon,lat,year,scenario){
  micro<-try(micro_terra_Liang(
    #Type in parameters that are not default
    scenario=scenario,
    terra_source="/scratch/gpfs/liangm/TerraClim/data", #default='c:/era5_data/era5': specify folder and file prefix with local ERA5 data extracted via the mcera5 package (no trailing forward slash); Specify path of downloaded ERA5 climate data
    loc=c(lon,lat), #default=c(-91.415669, -0.287145): Longitude and latitude (decimal degrees);
    ystart=year, 
    yfinish=year,
    timeinterval = 12
  ))
  metout<-as.data.frame(micro[["metout"]])
  if(lat>0){
    monthdata<-metout[metout$DOY==unique(metout$DOY)[7],]
  }else{
    monthdata<-metout[metout$DOY==unique(metout$DOY)[1],]
  }
  #### TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 2m default)
  T_data<-data.frame(Tmean_2m=mean(monthdata$TAREF),Trange_2m=max(monthdata$TAREF)-min(monthdata$TARE),
                Tmean_1cm=mean(monthdata$TALOC),Trange_1cm=max(monthdata$TALOC)-min(monthdata$TALOC))
  T_data
}

####Run microclimate model
sce=0

micro_point<-data.frame()
for(year_to_run in c(1986:2015)){
  micro_list<-list()
    timer<-vector()
  for(p in 1:floor(nrow(Category_point)/ncore)){
    ptm<-proc.time() 
    result<-mclapply(((p-1)*ncore+1):(p*ncore),function(i){
      i_point<-Category_point[i,]
      result_i<-withTimeout({
        cal.micro(lon=Category_point$LON[i],lat=Category_point$LAT[i],year=year_to_run,scenario=sce)
      },timeout = 60,onTimeout = "silent")
      
      if(class(result_i)=="NULL") {
        result_i<-data.frame(Tmean_2m=NA,Trange_2m=NA,Tmean_1cm=NA,Trange_1cm=NA) 
      }
      i_bind<-cbind(i_point,result_i)
      i_bind
      
      },mc.cores=ncore)
    
    micro_list<-append(micro_list,result)
    timer[p]<-as.numeric((proc.time() - ptm)[3])
    rest<-round(mean(timer,na.rm=T)*(floor(nrow(Category_point)/ncore)-p)/3600,digits=2)
    print(paste0("sce = ",sce," year_to_run = ",year_to_run,"; ","p = ",p,"; Time left for this year: ",rest," hours"))
  }
  
  micro_list1<-append(micro_list,mclapply((floor(nrow(Category_point)/ncore)*ncore+1):nrow(Category_point),function(i){
    i_point<-Category_point[i,]
    result_i<-withTimeout({cal.micro(lon=Category_point$LON[i],lat=Category_point$LAT[i],year=year_to_run,scenario=sce)},
                          timeout = 60,onTimeout = "silent")
    if(class(result_i)=="NULL") {
      result_i<-data.frame(Tmean_2m=NA,Trange_2m=NA,Tmean_1cm=NA,Trange_1cm=NA) 
    }
    i_bind<-cbind(i_point,result_i)
    i_bind
  },mc.cores=ncore))
  
  micro_point_year<-as.data.frame(do.call(rbind,micro_list1))
  micro_point_year["YEAR"]<-year_to_run
  micro_point<-rbind(micro_point,micro_point_year)
}

write.csv(micro_point, file=paste0("/scratch/gpfs/liangm/Resting_status_Jiang/Data/point_micro.csv"),row.names = FALSE)




#####Figure########
library(plyr)
library(ggplot2)

micro_point<-read.csv("/Users/jiangzhongwen/Desktop/Resting ecology/Data/point_micro.csv")
Rest_FP_Category<-read.csv("/Users/jiangzhongwen/Desktop/Resting ecology/Data/Rest_FP_Category.csv")
Rest_FP_Category2<-Rest_FP_Category[Rest_FP_Category$ELE_dif>500,]
Category_point_mr<-read.csv("/Users/jiangzhongwen/Desktop/Resting ecology/Data/Category_allpoint_mr.csv")


micro_point2<-data.frame()
for (ID in unique(micro_point$ID)) {
  ID_point<-micro_point[micro_point$ID==ID,]
  ID_mr<-Category_point_mr[Category_point_mr$ID==ID,]
  ID_high<-ID_point[ID_point$ELE==max(ID_point$ELE),]
  ID_low<-ID_point[ID_point$ELE==min(ID_point$ELE),]
  ID_high["Elevation"]<-"High"
  ID_low["Elevation"]<-"Low"
  ID_point<-rbind(ID_high,ID_low)
  ID_point["Category_mr"]<-ID_mr$Category
  micro_point2<-rbind(micro_point2,ID_point)
}
#####elevational range >500m
micro_point2<-subset(micro_point2,ID%in%Rest_FP_Category2$ID)
#####food passage time
Tmean_2m_sum_fp<-ddply(micro_point2,c("LON","LAT","ID","ELE","Category_fp","Elevation"),summarise,N=sum(!is.na(Tmean_2m)),
                 mean=mean(Tmean_2m,na.rm=T),
                 sd=sd(Tmean_2m,na.rm=T),
                 se=sd/sqrt(N))
Trange_2m_sum_fp<-ddply(micro_point2,c("LON","LAT","ID","ELE","Category_fp","Elevation"),summarise,N=sum(!is.na(Trange_2m)),
                    mean=mean(Trange_2m,na.rm=T),
                    sd=sd(Trange_2m,na.rm=T),
                    se=sd/sqrt(N))
Tmean_1cm_sum_fp<-ddply(micro_point2,c("LON","LAT","ID","ELE","Category_fp","Elevation"),summarise,N=sum(!is.na(Tmean_1cm)),
                    mean=mean(Tmean_1cm,na.rm=T),
                    sd=sd(Tmean_1cm,na.rm=T),
                    se=sd/sqrt(N))
Trange_1cm_sum_fp<-ddply(micro_point2,c("LON","LAT","ID","ELE","Category_fp","Elevation"),summarise,N=sum(!is.na(Trange_1cm)),
                    mean=mean(Trange_1cm,na.rm=T),
                    sd=sd(Trange_1cm,na.rm=T),
                    se=sd/sqrt(N))

#####metabolic rate
Tmean_2m_sum_mr<-ddply(micro_point2,c("LON","LAT","ID","ELE","Category_mr","Elevation"),summarise,N=sum(!is.na(Tmean_2m)),
                       mean=mean(Tmean_2m,na.rm=T),
                       sd=sd(Tmean_2m,na.rm=T),
                       se=sd/sqrt(N))
Trange_2m_sum_mr<-ddply(micro_point2,c("LON","LAT","ID","ELE","Category_mr","Elevation"),summarise,N=sum(!is.na(Trange_2m)),
                        mean=mean(Trange_2m,na.rm=T),
                        sd=sd(Trange_2m,na.rm=T),
                        se=sd/sqrt(N))
Tmean_1cm_sum_mr<-ddply(micro_point2,c("LON","LAT","ID","ELE","Category_mr","Elevation"),summarise,N=sum(!is.na(Tmean_1cm)),
                        mean=mean(Tmean_1cm,na.rm=T),
                        sd=sd(Tmean_1cm,na.rm=T),
                        se=sd/sqrt(N))
Trange_1cm_sum_mr<-ddply(micro_point2,c("LON","LAT","ID","ELE","Category_mr","Elevation"),summarise,N=sum(!is.na(Trange_1cm)),
                         mean=mean(Trange_1cm,na.rm=T),
                         sd=sd(Trange_1cm,na.rm=T),
                         se=sd/sqrt(N))


Figure_Tmean_1cm_fp<-ggplot(Tmean_1cm_sum_fp, aes(x=Category_fp, y=mean, fill=Elevation)) +
  geom_boxplot(size=1)+
  scale_fill_manual(values=c("#437199","#CB4B58"))+
  labs(x="Category", y = "Mean daily temperature (°C)",title = "Cumulative digestion")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(linewidth = 0.8),
        strip.text = element_text(colour = "black",
                                  family = "serif",face = "bold", size = 14), 
        strip.background = element_rect(linetype = 0 ),
        plot.subtitle = element_text(family = "serif", size = 19, face = "bold", colour = "black", hjust = 0.5), 
        plot.caption = element_text(family = "serif", face = "bold", hjust = 0.5), 
        axis.title = element_text(family = "serif",face=c("bold"),size = 24), 
        axis.text = element_text(size = 24,family = "serif",color="black",face=c("bold")), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        legend.background = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5,vjust = -1,face = "bold",size = 24,family = "serif"))

Figure_Trange_1cm_fp<-ggplot(Trange_1cm_sum_fp, aes(x=Category_fp, y=mean, fill=Elevation)) +
  geom_boxplot(size=1)+
  scale_fill_manual(values=c("#437199","#CB4B58"))+
  labs(x="Category", y = "Temperature daily range (°C)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(linewidth = 0.8),
        strip.text = element_text(colour = "black",
                                  family = "serif",face = "bold", size = 14), 
        strip.background = element_rect(linetype = 0 ),
        plot.subtitle = element_text(family = "serif", size = 19, face = "bold", colour = "black", hjust = 0.5), 
        plot.caption = element_text(family = "serif", face = "bold", hjust = 0.5), 
        axis.title = element_text(family = "serif",face=c("bold"),size = 24), 
        axis.text = element_text(size = 24,family = "serif",color="black",face=c("bold")), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        legend.background = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.025,vjust = -5,face = "bold",size = 24,family = "serif"))

Figure_Ele_fp<-ggplot(Trange_1cm_sum_fp, aes(x=Category_fp, y=ELE, fill=Elevation)) +
  geom_boxplot(size=1)+
  scale_fill_manual(values=c("#437199","#CB4B58"))+
  labs(x="Category", y = "Elevation (m)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(linewidth = 0.8),
        strip.text = element_text(colour = "black",
                                  family = "serif",face = "bold", size = 14), 
        strip.background = element_rect(linetype = 0 ),
        plot.subtitle = element_text(family = "serif", size = 19, face = "bold", colour = "black", hjust = 0.5), 
        plot.caption = element_text(family = "serif", face = "bold", hjust = 0.5), 
        axis.title = element_text(family = "serif",face=c("bold"),size = 24), 
        axis.text = element_text(size = 24,family = "serif",color="black",face=c("bold")), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        legend.background = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.025,vjust = -5,face = "bold",size = 24,family = "serif"))


Figure_Tmean_1cm_mr<-ggplot(Tmean_1cm_sum_mr, aes(x=Category_mr, y=mean, fill=Elevation)) +
  geom_boxplot(size=1)+
  scale_fill_manual(values=c("#437199","#CB4B58"))+
  labs(x="Category", y = "Mean daily temperature (°C)",title = "Cumulative metabolism")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(linewidth = 0.8),
        strip.text = element_text(colour = "black",
                                  family = "serif",face = "bold", size = 14), 
        strip.background = element_rect(linetype = 0 ),
        plot.subtitle = element_text(family = "serif", size = 19, face = "bold", colour = "black", hjust = 0.5), 
        plot.caption = element_text(family = "serif", face = "bold", hjust = 0.5), 
        axis.title = element_text(family = "serif",face=c("bold"),size = 24), 
        axis.text = element_text(size = 24,family = "serif",color="black",face=c("bold")), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        legend.background = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5,vjust = -1,face = "bold",size = 24,family = "serif"))

Figure_Trange_1cm_mr<-ggplot(Trange_1cm_sum_mr, aes(x=Category_mr, y=mean, fill=Elevation)) +
  geom_boxplot(size=1)+
  scale_fill_manual(values=c("#437199","#CB4B58"))+
  labs(x="Category", y = "Temperature daily range (°C)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(linewidth = 0.8),
        strip.text = element_text(colour = "black",
                                  family = "serif",face = "bold", size = 14), 
        strip.background = element_rect(linetype = 0 ),
        plot.subtitle = element_text(family = "serif", size = 19, face = "bold", colour = "black", hjust = 0.5), 
        plot.caption = element_text(family = "serif", face = "bold", hjust = 0.5), 
        axis.title = element_text(family = "serif",face=c("bold"),size = 24), 
        axis.text = element_text(size = 24,family = "serif",color="black",face=c("bold")), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        legend.background = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.025,vjust = -5,face = "bold",size = 24,family = "serif"))


Figure_Ele_mr<-ggplot(Trange_1cm_sum_mr, aes(x=Category_mr, y=ELE, fill=Elevation)) +
  geom_boxplot(size=1)+
  scale_fill_manual(values=c("#437199","#CB4B58"))+
  labs(x="Category", y = "Elevation (m)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(linewidth = 0.8),
        strip.text = element_text(colour = "black",
                                  family = "serif",face = "bold", size = 14), 
        strip.background = element_rect(linetype = 0 ),
        plot.subtitle = element_text(family = "serif", size = 19, face = "bold", colour = "black", hjust = 0.5), 
        plot.caption = element_text(family = "serif", face = "bold", hjust = 0.5), 
        axis.title = element_text(family = "serif",face=c("bold"),size = 24), 
        axis.text = element_text(size = 24,family = "serif",color="black",face=c("bold")), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        legend.background = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.025,vjust = -5,face = "bold",size = 24,family = "serif"))


library(ggpubr)

pdf("/Users/jiangzhongwen/Desktop/Resting ecology/Figure/Figure S4.pdf",height=18,width=18)

ggarrange(Figure_Tmean_1cm_fp,Figure_Tmean_1cm_mr,
          Figure_Trange_1cm_fp,Figure_Trange_1cm_mr,
          Figure_Ele_fp,Figure_Ele_mr,
          ncol=2,nrow=3,labels=c('a','b','c','d',"e","f"),
          font.label = list(size =30, color = "black",family="serif"))

dev.off()
