######Eco-physiology of resting lizards######
#####Zhong-wen Jiang 24 Dec 2022#############

library(lmerTest)
library(emmeans)
library(car)
library(plyr)
library(nlme)
library(ggplot2)
library(ggpmisc)
library(ggtext)
library(microclima)
library(NicheMapR)
library(mcera5)####era5 0.25degree=15min
library(raster)
library(ModelMetrics)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(viridis)

Lizard_traits<-read.csv("Lizard_traits.csv")

#########Figure 1, The resting depth#######
RestTb<-read.csv("Burrow.csv")

#Depth of RestTb
bartlett.test(Depth~Elevation,data=RestTb)

mod_Depth<-aov(Depth~SVL+Elevation,data = RestTb)
summary(mod_Depth)
shapiro.test(resid(mod_Depth))#0.68
emmeans(mod_Depth, list(pairwise~Elevation), adjust = "tukey")
# 1             estimate   SE  df t.ratio p.value
#2600m - 3400m    -1.43 1.53 346 -0.932  0.6204 
#2600m - 4200m    11.30 1.59 346  7.124  <.0001 
#3400m - 4200m    12.73 1.49 346  8.540  <.0001 

Restdepth_sum<-ddply(RestTb,c("Elevation"),summarise,N=sum(!is.na(Depth)),
                  mean=mean(Depth,na.rm=T),
                  sd=sd(Depth,na.rm=T),
                  se=sd/sqrt(N))


FigDepth<-ggplot(RestTb, aes(x=Elevation, y=Depth, fill=Elevation)) +
  geom_boxplot(size=1)+
  scale_fill_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  labs(x="Population", y = "Resting depth (cm)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.8,0.8),
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

ggsave(FigDepth,file="Figuredepth.pdf",device = "pdf",height=6,width=6,dpi = 300)



#########Figure 2, body temperature and thermal safety margin#########
#####Active body temperature analysis and Figure 1 b####
ActTb<-read.csv("Tbody.csv")
#ActTb$Year<-as.factor(ActTb$Year)

#Remove the 8 & 19 o'clock#
#ActTb<-subset(ActTb,Tbody>30&Time!=8&Time!=19 )
#statistical analysis#
leveneTest(Tbody~Elevation,data = ActTb)#abnormal data, LeveneTest be used to test the Homogeneity of Variance;normal data, Bartlett.test be used
mod_Tbody<-lme(Tbody~Elevation,random = ~1|Time/Date,data=ActTb)
shapiro.test(resid(mod_Tbody))
qqnorm(resid(mod_Tbody))
hist(resid(mod_Tbody))
summary(mod_Tbody)
anova(mod_Tbody)
emmeans(mod_Tbody, list(pairwise~Elevation), adjust = "tukey")
#1             estimate    SE   df t.ratio p.value
#2600m - 3400m     1.23 0.402 97   3.058  0.0080
#2600m - 4200m     2.62 0.450 97   5.819  <.0001
#3400m - 4200m     1.39 0.403 97   3.438  0.0025 

#Heat tolerance
#Read in data
CTmax<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/CTm.csv")[,c(1,3,6)]
colnames(CTmax)<-c("Elevation","SVL","CTmax")
dat_CTmax<-CTmax[CTmax$SVL>5 ,c(1:3)]



###Figure of active body temperatures & CTmax

ActTb_sum<-ddply(ActTb,c("Elevation","Time"),summarise,N=sum(!is.na(Tbody)),
                 mean=mean(Tbody,na.rm=T),
                 sd=sd(Tbody,na.rm=T),
                 se=sd/sqrt(N))

FigActTb<-ggplot() +
  geom_line(data=ActTb_sum,aes(x = Time, y =mean,group = Elevation,linetype=Elevation,color=Elevation),linewidth=1,position=position_dodge(0.3))+
  geom_point(data=ActTb_sum,aes(x = Time, y =mean,group = Elevation,shape=Elevation,color=Elevation),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=ActTb_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=Elevation), width=1,
                position=position_dodge(0.3))+
  geom_hline(yintercept =mean(dat_CTmax$CTmax),col="#CB4B58",linewidth=1)+
  annotate("text", x = 16, y = 46.5, label = paste0("CTmax=",round(mean(dat_CTmax$CTmax),1)),size=8,family="serif")+
  scale_color_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  scale_x_continuous(breaks = seq(8,20,1))+
 # ylim(10,50)+
  labs(x="Hour", y = "Active body temperature (°C)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.2,0.8),
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
ggsave(FigActTb,file="Figure1_b.pdf",device = "pdf",height=6,width=6,dpi = 300)



#####Resting body temperature analysis and Figure 1 d#####
RestTb<-read.csv("Burrow.csv")
#hypothesis test
shapiro.test(RestTb$Tb)
leveneTest(Tb~Elevation,data=RestTb)####When the data is abnormal, the levenTest be used to test homegeneity of Variance 

mod_RestTb<-gls(Tb~Time+Elevation,weights=varIdent(form=~1|Elevation),data = RestTb)
summary(mod_RestTb)
anova(mod_RestTb)
shapiro.test(resid(mod_RestTb))
emmeans(mod_RestTb, list(pairwise~Elevation), adjust = "tukey")
# 1             estimate    SE  df t.ratio p.value
#2600m - 3400m     5.67 0.290 151 19.559  <.0001 
#2600m - 4200m     8.74 0.352 237 24.819  <.0001 
#3400m - 4200m     3.07 0.251 193 12.237  <.0001 

#Cold tolerance
#Read in data
CTmin<-read.csv("CTm.csv")
colnames(CTmin)<-c("Elevation","SVL","CTmin")



###Figure of active body temperatures & CTmin
RestTb$Time<-factor(RestTb$Time,levels = c("20","21","22","23","0","1","2","3","4","5","6","7"))
RestTb$TimeNew<-factor(RestTb$TimeNew,levels = c("20","22","24","2","4","6","7"))
RestTb_sum<-ddply(RestTb,c("Elevation","Time"),summarise,N=sum(!is.na(Tb)),
                  mean=mean(Tb,na.rm=T),
                  sd=sd(Tb,na.rm=T),
                  se=sd/sqrt(N))

FigRestTb<-ggplot() +
  geom_line(data=RestTb_sum,aes(x = Time, y =mean,group = Elevation,linetype=Elevation,color=Elevation),linewidth=1,position=position_dodge(0.3))+
  geom_point(data=RestTb_sum,aes(x = Time, y =mean,group = Elevation,shape=Elevation,color=Elevation),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=RestTb_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=Elevation), width=1,
                position=position_dodge(0.3))+
  scale_color_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  scale_y_continuous(limits=c(0,29),breaks = seq(0.0,29.0,2.0))+
  labs(x="Hour", y = "Resting body temperature (°C)")+
  geom_hline(yintercept =mean(CTmin[CTmin$Elevation=="4200m",]$CTmin),col="#437199",linewidth=1)+
  annotate("text", x = 8.5, y = 3, label = paste0("CTmin(4200m)=",round(mean(CTmin[CTmin$Elevation=="4200m",]$CTmin),1)),size=8,family="serif")+
  geom_hline(yintercept =mean(CTmin[CTmin$Elevation!="4200m",]$CTmin),col="#FDAE61",linewidth=1)+
  annotate("text", x = 8, y = 5.5, label = paste0("CTmin(2600m,3400m)=",round(mean(CTmin[CTmin$Elevation!="4200m",]$CTmin),1)),size=8,family="serif")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="none",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth  = 0.8),
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
ggsave(FigRestTb,file="Figure1_d.pdf",device = "pdf",height=6,width=6,dpi = 300)



#####histgram of whole day body temperatures#####

head(ActTb)
head(RestTb)
ActTb_rb<-ActTb[,c(1,2,3,4,5,6,7)]
RestTb_rb<-RestTb[,c(1,3,6,4,7,8,9)]
colnames(RestTb_rb)<-c("Year","Date","Elevation","Time","Sex","SVL","Tbody")
Tb_all<-rbind(ActTb_rb,RestTb_rb)

head(Tb_all)

Tball_sum<-ddply(Tb_all,c("Elevation","Time"),summarise,N=sum(!is.na(Tbody)),
                  mean=mean(Tbody,na.rm=T),
                  sd=sd(Tbody,na.rm=T),
                  se=sd/sqrt(N))

Figure2_hist<-ggplot(data=Tball_sum,aes(x=mean,fill=Elevation))+
#geom_rect(data=NULL,aes(xmin=0.5,xmax=26,ymin=-Inf,ymax=Inf),fill="grey")+
#geom_histogram(alpha=0.8,position = "identity") +
geom_histogram(alpha=1) +
scale_fill_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
labs(x="Body temperature (°C)")+
theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.1,0.8),
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
ggsave(Figure2_hist,file="Figure2_hist.pdf",device = "pdf",height=6,width=12,dpi = 300)



#####The validation of operative temperature########

#####Minimum operative temperature#
#Resting environmental temperature
Tsoil_2600m<-read.csv("Tsoil_2600m.csv")
Tsoil_2600m<-na.omit(Tsoil_2600m)
#Tsoil_2600m$Time<-factor(Tsoil_2600m$Time,levels = c("12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5","0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5"))
Tsoil_2600m["POP"]<-"2600m"

Tsoil_3400m<-read.csv("Tsoil_3400m.csv")
Tsoil_3400m<-na.omit(Tsoil_3400m)
#Tsoil_3400m$Time<-factor(Tsoil_3400m$Time,levels = c("12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5","0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5"))
Tsoil_3400m["POP"]<-"3400m"

Tsoil_4200m<-read.csv("Tsoil_4200m.csv")
Tsoil_4200m<-na.omit(Tsoil_4200m)
#Tsoil_4200m$Time<-factor(Tsoil_4200m$Time,levels = c("12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5","0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5"))
Tsoil_4200m["POP"]<-"4200m"

Tsoil<-rbind(Tsoil_2600m,Tsoil_3400m,Tsoil_4200m)
Tsoil<-Tsoil[Tsoil$Time<7.5|Tsoil$Time>19.5,]
Tsoil["DOY"]<-unclass(as.POSIXlt(Tsoil$Date))[["yday"]]
#####10/20/30/50cm in microclimate model and  10cm 20cm 30cm 40cm 50cm 60cm 80cm in copper

year=2020
#pop="2600m"
Tsoil_mod<-data.frame()
for (pop in c("2600m","3400m","4200m")) {
  lon=Lizard_traits[Lizard_traits$Ele==pop,]$lon
  lat=Lizard_traits[Lizard_traits$Ele==pop,]$lat

  micro<-micro_era5( 
    run.gads = 2,
    loc=c(lon,lat), #default=c(-91.415669, -0.287145): Longitude and latitude (decimal degrees);
    dstart=paste0("01/01/",year), #default="01/01/2019": First day to run, date in format "d/m/Y" e.g. "01/01/2015";
    dfinish=paste0("31/12/",year), #default="31/07/2019": Last day to run, date in format "d/m/Y" e.g. "31/12/2015";
    #slope = slope,
    #aspect=aspect,
    REFL=0.35, #default=0.15: Soil solar reflectance, decimal %; Sand, dry white, in Table 11.2 in Campbell and Norman 2012
    Usrhyt=0.01, #default=0.01: Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest; The height of lizards were sitting roughly 0.01 above ground
    runshade=1, #default=1: Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)? For speeding up
    warm=0, #default=0: uniform warming, °C; For simulating climate change
    spatial="QT_Restingstatus", #default='c:/era5_data/era5': specify folder and file prefix with local ERA5 data extracted via the mcera5 package (no trailing forward slash); Specify path of downloaded ERA5 climate data
    RUF=0.0004, #default=0.004: Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m.; Smooth desert
    CMH2O=0.5, #default=1: Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground); Very dry
    minshade = 0,
    maxshade = 1
  )
  soil<-as.data.frame(micro$soil)
  mod_soil<-data.frame(POP=pop,DOY=soil$DOY,TIME=soil$TIME,"10cm"=soil$D10cm,"20cm"=soil$D20cm,"30cm"=soil$D30cm,"50cm"=soil$D50cm)
  Tsoil_mod<-rbind(Tsoil_mod,mod_soil)
}

colnames(Tsoil_mod)<-c("POP","DOY","TIME","10cm","20cm","30cm","50cm")
head(Tsoil_mod)
Tsoil_obs<-Tsoil[,c(6,7,3,4,5)]
head(Tsoil_obs)

Tsoil_mod<-melt(Tsoil_mod,id.vars = c("POP","DOY","TIME"))
colnames(Tsoil_mod)<-c("POP","DOY","TIME","Depth","Temprature")
colnames(Tsoil_obs)<-c("POP","DOY","TIME","Depth","Temprature")
Tsoil_obs<-Tsoil_obs[Tsoil_obs$Depth=="10cm"|Tsoil_obs$Depth=="20cm"|Tsoil_obs$Depth=="30cm"|Tsoil_obs$Depth=="50cm",]
#POP="2600m"
#TIME=6
#DOY=225
#Depth="10cm"

Tsoil_validation<-data.frame()

for (POP in c("2600m","3400m","4200m")) {
  
  Tsoilobs_pop<-Tsoil_obs[Tsoil_obs$POP==POP,]
  
  for (TIME in unique(round(Tsoilobs_pop$TIME,0))) {
    
    Tsoilobs_time<-Tsoilobs_pop[Tsoilobs_pop$TIME==TIME,]
    
    for (DOY in unique(Tsoilobs_time$DOY)) {
      
      Tsoilobs_doy<-Tsoilobs_time[Tsoilobs_time$DOY==DOY,]
      
      if(TIME>19){
        Tsoilmod_doy<-Tsoil_mod[Tsoil_mod$POP==POP&Tsoil_mod$TIME/60==(TIME-8)&Tsoil_mod$DOY==DOY,]
      }else{
        Tsoilmod_doy<-Tsoil_mod[Tsoil_mod$POP==POP&Tsoil_mod$TIME/60==(24+TIME-8)&Tsoil_mod$DOY==DOY-1,]
        
      }
      
      for (Depth in c("10cm","20cm","30cm","50cm")) {
        Tsoilobs_depth<-Tsoilobs_doy[Tsoilobs_doy$Depth==Depth,]
        Tsoilmod_depth<-Tsoilmod_doy[Tsoilmod_doy$Depth==Depth,]
        Tsoil_result<-data.frame(POP=POP,DOY=DOY,TIME=TIME,Depth=Depth,Tsoil_obs=mean(Tsoilobs_depth$Temprature),Tsoil_mod=Tsoilmod_depth$Temprature)
        Tsoil_validation<-rbind(Tsoil_validation,Tsoil_result)
      }
    }
  }
}

Tsoil_validation<-na.omit(Tsoil_validation)

cor<-cor.test(Tsoil_validation$Tsoil_mod,Tsoil_validation$Tsoil_obs,method = c("pearson"))
Rmse_soil<-rmse(Tsoil_validation$Tsoil_mod,Tsoil_validation$Tsoil_obs)#4.39
rho_soil<-cor$estimate#0.87


Figure_S1b<-ggplot()+
  geom_line(aes(x=0:30,y=0:30),linewidth=2,linetype=1,color="black")+
  geom_point(data = Tsoil_validation,aes(x=Tsoil_mod,y=Tsoil_obs),size=4,color="#4169E1",alpha=0.5)+
  annotate("richtext", x = 6, y = 30, label =paste0("**",paste0("R=",round(rho_soil,2),";RMSE=",round(Rmse_soil,2)),"**"),size=6,family="serif",label.colour = NA,fill = NA)+
  theme_classic()+
  xlab("Predicted soil temperatures")+
  ylab("Observed soil temperatures")+
  theme(legend.title=element_blank(),
        legend.position="none",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth  = 0.8),
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
 # facet_grid(.~Depth)

ggsave(Figure_S1b,file=paste("FigureS1b.pdf",sep=""),device = "pdf",height=6,width=6,dpi = 300)

pdf("Figure_S1.pdf",height=6,width=12)
ggarrange(Figure_S1a,Figure_S1b,ncol=2,labels=c('a','b'),
          font.label = list(size =30, color = "black",family="serif"))

dev.off()

######################################################################
#########Figure 3, Digestive accumulation and carbon dioxide production##########

ActTb_2600m<-ActTb[ActTb$Elevation=="2600m",]
ActTb_3400m<-ActTb[ActTb$Elevation=="3400m",]
ActTb_4200m<-ActTb[ActTb$Elevation=="4200m",]

RestTb_2600m<-RestTb[RestTb$Elevation=="2600m",]
RestTb_3400m<-RestTb[RestTb$Elevation=="3400m",]
RestTb_4200m<-RestTb[RestTb$Elevation=="4200m",]


######Food passage time 
######Model built according to Angilletta Jr, Michael J. "Thermal and physiological constraints on energy assimilation in a widespread lizard (Sceloporus undulatus)." Ecology 82.11 (2001): 3044-3056.
FPT<-read.csv("Passagetime.csv")
FPT$Temperature<-as.factor(FPT$Temperature)

leveneTest(passagetime~Ele,data=FPT)
mod_FPT<-aov(passagetime~Temperature*Ele,data = FPT)
summary.lm(mod_FPT)
shapiro.test(resid(mod_FPT))

####The figure of the multiple regression of three populations.

formula<-y~poly(x,2,raw = TRUE)
Figure_S2<-ggplot(FPT,aes(x=Temperature,y=passagetime,colour=Ele)) +
  geom_point()+
  geom_smooth(aes(fill=Ele),method = 'lm',formula =formula)+
  stat_poly_eq(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),rr.digits=2,
    formula = formula, parse = TRUE,hjust=-1
  )+
  scale_fill_manual(values = c("#5ED1BD","#5FA3E1","#562abc"))+
  scale_color_manual(values = c("#5ED1BD","#5FA3E1","#562abc"))+
  labs(x="Body temperature (°C)", y = "Food passage time (h)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.85,0.9),
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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
        plot.title = element_text(hjust = 0.5,face = "bold",size = 24,family = "serif"))
ggsave(Figure_S2,file="/Users/jiangzhongwen/Desktop/Resting ecology/Figure/Figure_S2.pdf",device = "pdf",height=6,width=9,dpi = 300)


#FPT$Temperature<-as.factor(FPT$Temperature)
FPT2600m<-FPT[which(FPT$Ele=="2600m"),]
FPT3400m<-FPT[which(FPT$Ele=="3400m"),]
FPT4200m<-FPT[which(FPT$Ele=="4200m"),]


####high-elevation population
FPT_mod2<-lm(passagetime~Temperature+I(Temperature^2),data=FPT[FPT$Ele!="4200m",])
summary(FPT_mod2)#R²=0.83 y=-38.95Tb+0.57Tb²+702.92

####2600m
FPT2600m_mod2<-lm(passagetime~Temperature+I(Temperature^2),data=FPT2600m)
summary(FPT2600m_mod2)#R²=0.86 y=-38.95Tb+0.57Tb²+701.08

####3400m
FPT3400m_mod2<-lm(passagetime~Temperature+I(Temperature^2),data=FPT3400m)
summary(FPT3400m_mod2)#R²= 0.81，y=-38.89Tb+0.57Tb²+703.83

####4200m
FPT4200m_mod2<-lm(passagetime~Temperature+I(Temperature^2),data=FPT4200m)
summary(FPT4200m_mod2)#R²=0.75,y=-31.45Tb+0.46Tb²+579.35

Rate_digest_2600m<-function(Tb){
  rate<-(1/(-38.95*Tb+0.57*(Tb)^2+701.08))
  return(rate)
}

Rate_digest_3400m<-function(Tb){
  rate<-(1/(-38.89*Tb+0.57*(Tb)^2+703.83))
  return(rate)
}

Rate_digest_4200m<-function(Tb){
  rate<-(1/(-31.45*Tb+0.46*(Tb)^2+579.35))
  return(rate)
}

ActTb_2600m["PassageRate"]<-Rate_digest_2600m(ActTb_2600m$Tbody)
ActTb_3400m["PassageRate"]<-Rate_digest_3400m(ActTb_3400m$Tbody)
ActTb_4200m["PassageRate"]<-Rate_digest_4200m(ActTb_4200m$Tbody)

ActPassageRate<-rbind(ActTb_2600m,ActTb_3400m,ActTb_4200m)

RestTb_2600m["PassageRate"]<-Rate_digest_2600m(RestTb_2600m$Tb)
RestTb_3400m["PassageRate"]<-Rate_digest_3400m(RestTb_3400m$Tb)
RestTb_4200m["PassageRate"]<-Rate_digest_4200m(RestTb_4200m$Tb)

RestPassageRate<-rbind(RestTb_2600m,RestTb_3400m,RestTb_4200m)

ActPassageRate["status"]<-"Activity"
RestPassageRate["status"]<-"Rest"

ActPassageRate<-ActPassageRate[,c(1,2,3,4,5,6,7,8,11,12)]
RestPassageRate<-RestPassageRate[,c(1,3,6,4,7,8,9,10,16,17)]
colnames(RestPassageRate)<-c("Year","Date","Elevation","Time","Sex","SVL","Tbody","Gravid","PassageRate","status")
PassageRate<-rbind(ActPassageRate,RestPassageRate)


Tbody_statistical<-ddply(PassageRate,c("Elevation","status"),summarise,N=sum(!is.na(Tbody)),
                       mean=mean(Tbody,na.rm=T),
                       sd=sd(Tbody,na.rm=T),
                       se=sd/sqrt(N))



PassageRate_statistical<-ddply(PassageRate,c("Elevation","status"),summarise,N=sum(!is.na(PassageRate)),
                       mean=mean(PassageRate,na.rm=T),
                       sd=sd(PassageRate,na.rm=T),
                       se=sd/sqrt(N))

PassageRate_sum<-ddply(PassageRate,c("Time","Elevation","status"),summarise,N=sum(!is.na(PassageRate)),
                           mean=mean(PassageRate,na.rm=T),
                           sd=sd(PassageRate,na.rm=T),
                           se=sd/sqrt(N))

act_pass<-PassageRate_sum[PassageRate_sum$status=="Activity",]
res_pas<-PassageRate_sum[PassageRate_sum$status=="Rest",]


wilcox.test(mean~Elevation,data = act_pass[act_pass$Elevation!="2600m",],paired=T)
wilcox.test(mean~Elevation,data = act_pass[act_pass$Elevation!="3400m",],paired=T)
wilcox.test(mean~Elevation,data = act_pass[act_pass$Elevation!="4200m",],paired=T)

wilcox.test(mean~Elevation,data = res_pas[res_pas$Elevation!="2600m",],paired=T)
wilcox.test(mean~Elevation,data = res_pas[res_pas$Elevation!="3400m",],paired=T)
wilcox.test(mean~Elevation,data = res_pas[res_pas$Elevation!="4200m",],paired=T)


PassageRate_sum1<-data.frame()
for (Ele in c("2600m","3400m","4200m")) {
  
  Pass<-PassageRate_sum[which(PassageRate_sum$Elevation==Ele),]

  for (status in c("Rest","Activity")) {
    
    pass_status<-Pass[Pass$status==status,]
      
    result<-data.frame(Elevation=Ele,status=status,Acum=sum(pass_status$mean))
    
    PassageRate_sum1<-rbind(PassageRate_sum1,result)
  }
}


Pass_compare1<-data.frame(Pair="2600m-3400m",status=c("Rest","Activity"),
                     differences=c(PassageRate_sum1[PassageRate_sum1$Elevation=="2600m"&PassageRate_sum1$status=="Rest",]$Acum-PassageRate_sum1[PassageRate_sum1$Elevation=="3400m"&PassageRate_sum1$status=="Rest",]$Acum,
                     PassageRate_sum1[PassageRate_sum1$Elevation=="2600m"&PassageRate_sum1$status=="Activity",]$Acum-PassageRate_sum1[PassageRate_sum1$Elevation=="3400m"&PassageRate_sum1$status=="Activity",]$Acum))

Pass_compare2<-data.frame(Pair="2600m-4200m",status=c("Rest","Activity"),
                     differences=c(PassageRate_sum1[PassageRate_sum1$Elevation=="2600m"&PassageRate_sum1$status=="Rest",]$Acum-PassageRate_sum1[PassageRate_sum1$Elevation=="4200m"&PassageRate_sum1$status=="Rest",]$Acum,
                                   PassageRate_sum1[PassageRate_sum1$Elevation=="2600m"&PassageRate_sum1$status=="Activity",]$Acum-PassageRate_sum1[PassageRate_sum1$Elevation=="4200m"&PassageRate_sum1$status=="Activity",]$Acum))

Pass_compare3<-data.frame(Pair="3400m-4200m",status=c("Rest","Activity"),
                     differences=c(PassageRate_sum1[PassageRate_sum1$Elevation=="3400m"&PassageRate_sum1$status=="Rest",]$Acum-PassageRate_sum1[PassageRate_sum1$Elevation=="4200m"&PassageRate_sum1$status=="Rest",]$Acum,
                                   PassageRate_sum1[PassageRate_sum1$Elevation=="3400m"&PassageRate_sum1$status=="Activity",]$Acum-PassageRate_sum1[PassageRate_sum1$Elevation=="4200m"&PassageRate_sum1$status=="Activity",]$Acum))

foodpassage_differ<-rbind(Pass_compare1,Pass_compare2,Pass_compare3)

all_fpt<-data.frame()
for (pair in unique(foodpassage_differ$Pair)) {
  pair_fpt<-foodpassage_differ[foodpassage_differ$Pair==pair,]
  result<-data.frame(Pair=pair,status="All",differences=sum(pair_fpt$differences))
  all_fpt<-rbind(all_fpt,result)
}
head(all_fpt)

foodpassage_differ<-rbind(foodpassage_differ,all_fpt)
foodpassage_differ$status<-factor(foodpassage_differ$status,level=c("Activity","Rest","All"))

write.csv(foodpassage_differ,file = "foodpassage_differ.csv",row.names = FALSE)

#####Figrue 3a_hist and 3b_hist food passage (a) and carbon dioxide production (b) distribution histogram
#####Figure 3a,b food passage (a) vary with time and (b) accumulation at activity and rest

Figure3a_hist<-ggplot(data=PassageRate_sum,aes(x=mean,fill=Elevation))+
  #geom_rect(data=NULL,aes(xmin=0.5,xmax=26,ymin=-Inf,ymax=Inf),fill="grey")+
  #geom_histogram(alpha=0.5,position = "identity") +
  geom_histogram(alpha=1) +
  scale_fill_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  labs(x="Digestion rate (%/h)",y="Count (n)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.5,0.8),
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
ggsave(Figure3a_hist,file="Figure3a_hist.pdf",device = "pdf",height=6,width=12,dpi = 300)


Figure3a<-ggplot() +
  geom_rect(data=NULL,aes(xmin=-0.5,xmax=7.5,ymin=-Inf,ymax=Inf),
            fill="grey")+
  geom_rect(data=NULL,aes(xmin=19.5,xmax=23.5,ymin=-Inf,ymax=Inf),
            fill="grey")+
  geom_line(data=PassageRate_sum,aes(x = Time, y =mean,group = Elevation,linetype=Elevation,color=Elevation),linewidth=1,position=position_dodge(0.3))+
  geom_point(data=PassageRate_sum,aes(x = Time, y =mean,group = Elevation,shape=Elevation,color=Elevation),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=PassageRate_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=Elevation), width=1,
                position=position_dodge(0.3))+
   scale_color_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  scale_x_continuous(breaks = seq(0,23,1))+
  labs(x="Hour", y = "Food passage rate (%)")+
  annotate("richtext", x = 3.5, y = 0.03, label = "**Rest**",size=8,family="serif",label.colour = NA,fill = NA)+
  annotate("richtext", x = 21.5, y = 0.03, label = "**Rest**",size=8,family="serif",label.colour = NA,fill = NA)+
  annotate("richtext", x = 13.5, y = 0.03, label = "**Activity**",size=8,family="serif",label.colour = NA,fill = NA)+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.15,0.8),
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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
  
Figure3b_2<-ggplot(data=foodpassage_differ,aes(x=status,y=differences,fill=Pair))+
  geom_bar(stat = 'identity',position = position_dodge(preserve = 'single'),width = 0.5)+
  scale_fill_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  labs(x="Status", y = "Difference in daily cumulative digestion (%)")+
  #geom_text(aes(label=round(differences,2)), 
   #         position = position_dodge2(width = 0.9, preserve = 'single'), 
    #        vjust = -0.2, hjust = 0.5,size=8,family="serif",fontface="bold")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.2,0.8),
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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


Figure3b<-ggplot(data=PassageRate_sum1,aes(x=Elevation,y=Acum,fill=status))+
  geom_bar(stat = 'identity',position = position_dodge(preserve = 'single'))+
  scale_fill_manual(values=c("grey90","grey"))+
  labs(x="Population", y = "Food passage accumulation (%)")+
  theme_classic()+

ggsave(Figure3b_2,file="Figure3b.pdf",device = "pdf",height=6,width=12,dpi = 300)

####MR from field data collection

ActTb_2600m<-ActTb[ActTb$Elevation=="2600m",]
ActTb_3400m<-ActTb[ActTb$Elevation=="3400m",]
ActTb_4200m<-ActTb[ActTb$Elevation=="4200m",]

RestTb_2600m<-RestTb[RestTb$Elevation=="2600m",]
RestTb_3400m<-RestTb[RestTb$Elevation=="3400m",]
RestTb_4200m<-RestTb[RestTb$Elevation=="4200m",]

MR<-read.csv("MRtest.csv")#### the best and final data is MRtest 

MR2600m<-subset(MR,Elevation=="2600m")
MR3400m<-subset(MR,Elevation=="3400m")
MR4200m<-subset(MR,Elevation=="4200m")

#all population
MR_mod<-lm(log10(V)~log10(BM)+Tb,data = MR[MR$Elevation!="4200m",])
summary(MR_mod)#R²=0.77, lgV=0.32lgBM+0.040Tb-1.574

#2600m
MR2600m_mod<-lm(log10(V)~log10(BM)+Tb,data = MR2600m)
summary(MR2600m_mod)#R²=0.7521, lgV=0.059479lgBM+0.040919Tb-1.372556 
#3400m
MR3400m_mod<-lm(log10(V)~log10(BM)+Tb,data = MR3400m)
summary(MR3400m_mod)#R²=0.8017, lgV=0.463171lgBM+0.038788Tb-1.670339 
#4200m
MR4200m_mod<-lm(log10(V)~log10(BM)+Tb,data = MR4200m)
summary(MR4200m_mod)#R²=0.808, lgV=1.020382lgBM+0.031665Tb-1.599509

Rate_MR_2600m<-function(Tb){
  rate<-(((0.042*6.31^0.059))*(10^(0.041*Tb)))
  return(rate)
}

Rate_MR_3400m<-function(Tb){
  rate<-(((0.021*6.31^0.463))*(10^(0.039*Tb)))
  return(rate)
}

Rate_MR_4200m<-function(Tb){
  rate<-((0.025*6.31^1.02)*(10^(0.032*Tb)))
  return(rate)
}

ActTb_2600m["MR"]<-Rate_MR_2600m(ActTb_2600m$Tbody)
ActTb_3400m["MR"]<-Rate_MR_3400m(ActTb_3400m$Tbody)
ActTb_4200m["MR"]<-Rate_MR_4200m(ActTb_4200m$Tbody)

ActMRRate<-rbind(ActTb_2600m,ActTb_3400m,ActTb_4200m)

RestTb_2600m["MR"]<-Rate_MR_2600m(RestTb_2600m$Tb)
RestTb_3400m["MR"]<-Rate_MR_3400m(RestTb_3400m$Tb)
RestTb_4200m["MR"]<-Rate_MR_4200m(RestTb_4200m$Tb)

RestMRRate<-rbind(RestTb_2600m,RestTb_3400m,RestTb_4200m)

ActMRRate["status"]<-"Activity"
RestMRRate["status"]<-"Rest"

ActMRRate<-ActMRRate[,c(1,2,3,4,5,6,7,8,11,12)]
RestMRRate<-RestMRRate[,c(1,3,6,4,7,8,9,10,16,17)]
colnames(RestMRRate)<-c("Year","Date","Elevation","Time","Sex","SVL","Tbody","Gravid","MR","status")
MRRate<-rbind(ActMRRate,RestMRRate)


MRRate_statistical<-ddply(MRRate,c("Elevation","status"),summarise,N=sum(!is.na(MR)),
                               mean=mean(MR,na.rm=T),
                               sd=sd(MR,na.rm=T),
                               se=sd/sqrt(N))


MRRate_sum<-ddply(MRRate,c("Time","Elevation","status"),summarise,N=sum(!is.na(MR)),
                       mean=mean(MR,na.rm=T),
                       sd=sd(MR,na.rm=T),
                       se=sd/sqrt(N))



act_mr<-MRRate_sum[MRRate_sum$status=="Activity",]
res_mr<-MRRate_sum[MRRate_sum$status=="Rest",]


wilcox.test(act_mr[act_mr$Elevation=="2600m",]$mean,act_mr[act_mr$Elevation=="3400m",]$mean,paired=T)
wilcox.test(act_mr[act_mr$Elevation=="2600m",]$mean,act_mr[act_mr$Elevation=="4200m",]$mean,paired=T)
wilcox.test(act_mr[act_mr$Elevation=="3400m",]$mean,act_mr[act_mr$Elevation=="4200m",]$mean,paired=T)


wilcox.test(res_mr[res_mr$Elevation=="2600m",]$mean,res_mr[res_mr$Elevation=="3400m",]$mean,paired=T)
wilcox.test(res_mr[res_mr$Elevation=="2600m",]$mean,res_mr[res_mr$Elevation=="4200m",]$mean,paired=T)
wilcox.test(res_mr[res_mr$Elevation=="3400m",]$mean,res_mr[res_mr$Elevation=="4200m",]$mean,paired=T)



MRRate_sum1<-data.frame()
for (Ele in c("2600m","3400m","4200m")) {
  
  MR<-MRRate_sum[which(MRRate_sum$Elevation==Ele),]
  
  for (status in c("Rest","Activity")) {
    
    MR_status<-MR[MR$status==status,]
    
    result<-data.frame(Elevation=Ele,status=status,Acum=sum(MR_status$mean))
    
    MRRate_sum1<-rbind(MRRate_sum1,result)
  }
  
}


MR_compare1<-data.frame(Pair="2600m-3400m",status=c("Rest","Activity"),
                          differences=c(MRRate_sum1[MRRate_sum1$Elevation=="2600m"&MRRate_sum1$status=="Rest",]$Acum-MRRate_sum1[MRRate_sum1$Elevation=="3400m"&MRRate_sum1$status=="Rest",]$Acum,
                                        MRRate_sum1[MRRate_sum1$Elevation=="2600m"&MRRate_sum1$status=="Activity",]$Acum-MRRate_sum1[MRRate_sum1$Elevation=="3400m"&MRRate_sum1$status=="Activity",]$Acum))

MR_compare2<-data.frame(Pair="2600m-4200m",status=c("Rest","Activity"),
                          differences=c(MRRate_sum1[MRRate_sum1$Elevation=="2600m"&MRRate_sum1$status=="Rest",]$Acum-MRRate_sum1[MRRate_sum1$Elevation=="4200m"&MRRate_sum1$status=="Rest",]$Acum,
                                        MRRate_sum1[MRRate_sum1$Elevation=="2600m"&MRRate_sum1$status=="Activity",]$Acum-MRRate_sum1[MRRate_sum1$Elevation=="4200m"&MRRate_sum1$status=="Activity",]$Acum))

MR_compare3<-data.frame(Pair="3400m-4200m",status=c("Rest","Activity"),
                          differences=c(MRRate_sum1[MRRate_sum1$Elevation=="3400m"&MRRate_sum1$status=="Rest",]$Acum-MRRate_sum1[MRRate_sum1$Elevation=="4200m"&MRRate_sum1$status=="Rest",]$Acum,
                                        MRRate_sum1[MRRate_sum1$Elevation=="3400m"&MRRate_sum1$status=="Activity",]$Acum-MRRate_sum1[MRRate_sum1$Elevation=="4200m"&MRRate_sum1$status=="Activity",]$Acum))

mr_differ<-rbind(MR_compare1,MR_compare2,MR_compare3)

all_mr<-data.frame()
for (pair in unique(mr_differ$Pair)) {
  pair_mr<-mr_differ[mr_differ$Pair==pair,]
  result<-data.frame(Pair=pair,status="All",differences=sum(pair_mr$differences))
  all_mr<-rbind(all_mr,result)
}

head(all_mr)

mr_differ<-rbind(mr_differ,all_mr)

mr_differ$status<-factor(mr_differ$status,level=c("Activity","Rest","All"))

MRRate_sum$Time<-as.numeric(MRRate_sum$Time)

#####Figure 3c,d Carbon dioxide production (c) with time and (d) accumulation in activity and rest


Figure3c_hist<-ggplot(data=MRRate_sum,aes(x=mean,fill=Elevation))+
  #geom_rect(data=NULL,aes(xmin=0.5,xmax=26,ymin=-Inf,ymax=Inf),fill="grey")+
  geom_histogram(alpha=1) +
  scale_fill_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  labs(x="Metabolic rate (VCO2;ml/h)",y="Count (n)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="none",
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
ggsave(Figure3c_hist,file="Figure3c_hist.pdf",device = "pdf",height=6,width=12,dpi = 300)



Figure3c<-ggplot() +
  geom_rect(data=NULL,aes(xmin=-0.5,xmax=7.5,ymin=-Inf,ymax=Inf),
            fill="grey")+
  geom_rect(data=NULL,aes(xmin=19.5,xmax=23.5,ymin=-Inf,ymax=Inf),
            fill="grey")+
  geom_line(data=MRRate_sum,aes(x = Time, y =mean,group = Elevation,linetype=Elevation,color=Elevation),linewidth=1,position=position_dodge(0.3))+
  geom_point(data=MRRate_sum,aes(x = Time, y =mean,group = Elevation,shape=Elevation,color=Elevation),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=MRRate_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=Elevation), width=1,
                position=position_dodge(0.3))+
  scale_color_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  scale_x_continuous(breaks = seq(0,23,1))+
  annotate("richtext", x = 3.5, y = 2.5, label = "**Rest**",size=8,family="serif",label.colour = NA,fill = NA)+
  annotate("richtext", x = 21.5, y = 2.5, label = "**Rest**",size=8,family="serif",label.colour = NA,fill = NA)+
  annotate("richtext", x = 13.5, y = 2.5, label = "**Activity**",size=8,family="serif",label.colour = NA,fill = NA)+
  labs(x="Hour", y = "Metabolic rate (VCO2;ml/h)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="none",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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



Figure3d_2<-ggplot(data=mr_differ,aes(x=status,y=differences,fill=Pair))+
  geom_bar(stat = 'identity',position = position_dodge(preserve = 'single'),width = 0.5)+
  scale_fill_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  labs(x="Status", y = "Difference in daily cumulative metabolism (VCO2;ml)")+
# geom_text(aes(label=round(Contribution,2)), 
#          position = position_dodge2(width = 0.9, preserve = 'single'), 
#          vjust = -0.2, hjust = 0.5,size=8,family="serif",fontface="bold")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="none",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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



Figure3d<-ggplot(data=MRRate_sum1,aes(x=Elevation,y=Acum,fill=status))+
  geom_bar(stat = 'identity',position = position_dodge(preserve = 'single'))+
  scale_fill_manual(values=c("grey90","grey"))+
  labs(x="Population", y = "Carbon dioxide production (ml)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="none",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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

pdf("Figure_3.pdf",height=14,width=20)
ggarrange(Figure3a_hist,Figure3b_2,Figure3c_hist,Figure3d_2,ncol=2,nrow=2,labels=c('a','b','c','d'),hjust = -5,vjust = 1.5,
          font.label = list(size =30, color = "black",family="serif"))

dev.off()


MR<-read.csv("MRtest.csv")[,c(1,2,4,7,14)]#### the best and final data is MRtest 
colnames(MR)<-c("ID","Ele","Tb","BM","MR")


MR_sum<-ddply(MR,c("Ele","Tb"),summarise,N=sum(!is.na(MR)),
              mean=mean(MR,na.rm=T),
              sd=sd(MR,na.rm=T),
              se=sd/sqrt(N))

Figure_S3<-ggplot() +
  geom_line(data=MR_sum,aes(x = Tb, y =mean,group = Ele,color=Ele),linewidth=1,position=position_dodge(0.3))+
  geom_point(data=MR_sum,aes(x = Tb, y =mean,group = Ele,shape=Ele,color=Ele),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=MR_sum,size=1,aes(x = Tb,ymin=mean-se, ymax=mean+se,color=Ele), width=1,
                position=position_dodge(0.3))+
  scale_color_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  scale_x_continuous(breaks = seq(10,40,5))+
  labs(x="Body temperature (°C)", y = "Metabolic rate (VCO2;ml/h/g)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.15,0.9),
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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
        plot.title = element_text(hjust = 0.5,face = "bold",size = 24,family = "serif"))

ggsave(Figure_S3,file="Figure_S3.pdf",device = "pdf",height=6,width=9,dpi = 300)







#######Figure 4, warming =6, after a hundred with ssp5-8.5####
Lizard_traits<-read.csv("Lizard_traits.csv")

#Cox DTC, Maclean IMD, Gardner AS, Gaston KJ. Global variation in diurnal asymmetry in temperature, cloud cover, specific humidity and precipitation and its association with leaf area index. Glob Change Biol. 2020;26:7099–7111. https://doi.org/10.1111/gcb.15336
#The Tibetan Plateau has experienced median night-time minimum temperature increases of 1.2°C (range: 0.21°C; 3.44°C) more than daytime maximum temperatures
#Tnight_w-Tday_w=1.2,(Tnight_w+Tday_w)/2=Tmean_w
#Tnight_w=Tmean_w+0.6,Tday_w=Tmean-0.6


######NicheMapR ACT=0 or 1&2

All_Tb<-data.frame()

for (warm in c(0,6,5.4,6.6)) {
  for (pop in c("2600m","3400m","4200m")) {
    lon=Lizard_traits[Lizard_traits$Ele==pop,]$lon
    lat=Lizard_traits[Lizard_traits$Ele==pop,]$lat
    trait_list<-Lizard_traits[Lizard_traits$Ele==pop,]
    for (year in 2020:2021) {
      
      micro<-micro_era5( 
        run.gads = 2,
        loc=c(lon,lat), #default=c(-91.415669, -0.287145): Longitude and latitude (decimal degrees);
        dstart=paste0("01/01/",year), #default="01/01/2019": First day to run, date in format "d/m/Y" e.g. "01/01/2015";
        dfinish=paste0("31/12/",year), #default="31/07/2019": Last day to run, date in format "d/m/Y" e.g. "31/12/2015";
        #slope = slope,
        #aspect=aspect,
        REFL=0.35, #default=0.15: Soil solar reflectance, decimal %; Sand, dry white, in Table 11.2 in Campbell and Norman 2012
        Usrhyt=0.01, #default=0.01: Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest; The height of lizards were sitting roughly 0.01 above ground
        runshade=1, #default=1: Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)? For speeding up
        warm=warm, #default=0: uniform warming, °C; For simulating climate change
        spatial="QT_Restingstatus", #default='c:/era5_data/era5': specify folder and file prefix with local ERA5 data extracted via the mcera5 package (no trailing forward slash); Specify path of downloaded ERA5 climate data
        RUF=0.0004, #default=0.004: Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m.; Smooth desert
        CMH2O=0.5, #default=1: Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground); Very dry
        minshade = 0,
        maxshade = 50
      )
      
      lizard_loc<-ectotherm(
        live=1,
        Ww_g	= trait_list$Ww_g, #default=40; Wet weight of animal (g), note this model is 'steady state' so no lags in heating/cooling due to mass
        alpha_max	= 0.95, #default=0.85; Maximum solar absorptivity, 0-1; keep this variable for sensitivity analysis
        alpha_min	= 0.95, #default=0.85; Minimum solar absorptivity, 0-1; keep this variable for sensitivity analysis
        T_F_min	= trait_list$T_F_min, #default=24; Minimum foraging temperature, °C (also affects burrow depth selection)
        T_F_max	= trait_list$T_F_max, #default=34; Maximum foraging temperature, °C
        T_B_min = trait_list$T_B_min, #default=17.5; Minimum basking temperature, °C
        T_RB_min	= trait_list$T_B_min, #default=17.5; Minimum temperature at which animal will move from retreat to basking site, °C
        T_pref	= trait_list$T_pref, #default=30; Preferred body temperature, °C
        CT_max	= trait_list$CT_max, #default=40; Critical thermal maximum, °C (affects burrow depth selection and may be used to impose death from heat stress)
        CT_min	= trait_list$CT_min, #default=6; Critical thermal minimum, °C (affects burrow depth selection and may be used to impose death from cold stress)
        M_1	= trait_list$M_1, #Metabolic rate parameter 1 V_O2=M_1*M^M_2*10^(M_3*Tb), in ml O2 / h, default parameters for lizards based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
        M_2	= trait_list$M_2, #Metabolic rate parameter 2
        M_3	= trait_list$M_3, #Metabolic rate parameter 3
        pct_wet=1,#%of surface area acting as a free-water exchanger, for computing cutaneous water loss,calculate by get_p_wet/ data from (Marco 2018 Zoology)
        mindepth=1,
        maxdepth=9,
        aestdepth	= 9, #Depth (soil node #) to which animal retreats if burrowing and aestivating due to desiccation
        nyears = 1, #micro$nyears, #Number of years the simulation runs for - must be consistent with dimensions of environmental input data
        shade_seek =1,
        burrow = 1
      )
      
      lizard_environ_out<- as.data.frame(lizard_loc[["environ"]])
      Tb_result<-data.frame(lizard_environ_out$DOY,lizard_environ_out$TC,lizard_environ_out$ACT)
      Tb_result["YEAR"]<-year
      Tb_result["POP"]<-pop
      Tb_result["TIME"]<-lizard_environ_out$TIME
      Tb_result["WARM"]<-warm
      colnames(Tb_result)<-c("DOY","Tb","ACT","YEAR","POP","TIME","WARM")
      All_Tb<-rbind(All_Tb,Tb_result)
      
    }
  }
}

All_Tb_2600m<-All_Tb[All_Tb$POP=="2600m",]
All_Tb_3400m<-All_Tb[All_Tb$POP=="3400m",]
All_Tb_4200m<-All_Tb[All_Tb$POP=="4200m",]

All_Tb_2600m["Passagerate"]<-Rate_digest_2600m(All_Tb_2600m$Tb)
All_Tb_3400m["Passagerate"]<-Rate_digest_3400m(All_Tb_3400m$Tb)
All_Tb_4200m["Passagerate"]<-Rate_digest_4200m(All_Tb_4200m$Tb)

All_Tb_2600m["MR"]<-Rate_MR_2600m(All_Tb_2600m$Tb)
All_Tb_3400m["MR"]<-Rate_MR_3400m(All_Tb_3400m$Tb)
All_Tb_4200m["MR"]<-Rate_MR_4200m(All_Tb_4200m$Tb)

All_phys<-rbind(All_Tb_2600m,All_Tb_3400m,All_Tb_4200m)
#July
All_phys<-All_phys[All_phys$DOY>180&All_phys$DOY<212,]

phys_change<-data.frame()

for (POP in unique(All_phys$POP)) {
  symchange_act_fpr<-sum(All_phys[All_phys$POP==POP&All_phys$ACT!="0"&All_phys$WARM=="6",]$Passagerate)-sum(All_phys[All_phys$POP==POP&All_phys$ACT!="0"&All_phys$WARM=="0",]$Passagerate)
  symchange_act_mr<-sum(All_phys[All_phys$POP==POP&All_phys$ACT!="0"&All_phys$WARM=="6",]$MR)-sum(All_phys[All_phys$POP==POP&All_phys$ACT!="0"&All_phys$WARM=="0",]$MR)
  symchange_rest_fpr<-sum(All_phys[All_phys$POP==POP&All_phys$ACT=="0"&All_phys$WARM=="6",]$Passagerate)-sum(All_phys[All_phys$POP==POP&All_phys$ACT=="0"&All_phys$WARM=="0",]$Passagerate)
  symchange_rest_mr<-sum(All_phys[All_phys$POP==POP&All_phys$ACT=="0"&All_phys$WARM=="6",]$MR)-sum(All_phys[All_phys$POP==POP&All_phys$ACT=="0"&All_phys$WARM=="0",]$MR)
  
  symchange_fpr<-sum(All_phys[All_phys$POP==POP&All_phys$WARM=="6",]$Passagerate)-sum(All_phys[All_phys$POP==POP&All_phys$WARM=="0",]$Passagerate)
  symchange_mr<-sum(All_phys[All_phys$POP==POP&All_phys$WARM=="6",]$MR)-sum(All_phys[All_phys$POP==POP&All_phys$WARM=="0",]$MR)
  
  
  asychange_act_fpr<-sum(All_phys[All_phys$POP==POP&All_phys$ACT!="0"&All_phys$WARM=="5.4",]$Passagerate)-sum(All_phys[All_phys$POP==POP&All_phys$ACT!="0"&All_phys$WARM=="0",]$Passagerate)
  asychange_act_mr<-sum(All_phys[All_phys$POP==POP&All_phys$ACT!="0"&All_phys$WARM=="5.4",]$MR)-sum(All_phys[All_phys$POP==POP&All_phys$ACT!="0"&All_phys$WARM=="0",]$MR)
  asychange_rest_fpr<-sum(All_phys[All_phys$POP==POP&All_phys$ACT=="0"&All_phys$WARM=="6.6",]$Passagerate)-sum(All_phys[All_phys$POP==POP&All_phys$ACT=="0"&All_phys$WARM=="0",]$Passagerate)
  asychange_rest_mr<-sum(All_phys[All_phys$POP==POP&All_phys$ACT=="0"&All_phys$WARM=="6.6",]$MR)-sum(All_phys[All_phys$POP==POP&All_phys$ACT=="0"&All_phys$WARM=="0",]$MR)
  
  asychange_fpr<-(sum(All_phys[All_phys$POP==POP&All_phys$WARM=="6.6"&All_phys$ACT=="0",]$Passagerate)+sum(All_phys[All_phys$POP==POP&All_phys$WARM=="5.4"&All_phys$ACT!="0",]$Passagerate))-sum(All_phys[All_phys$POP==POP&All_phys$WARM=="0",]$Passagerate)
  asychange_mr<-(sum(All_phys[All_phys$POP==POP&All_phys$WARM=="6.6"&All_phys$ACT=="0",]$MR)+sum(All_phys[All_phys$POP==POP&All_phys$WARM=="5.4"&All_phys$ACT!="0",]$MR))-sum(All_phys[All_phys$POP==POP&All_phys$WARM=="0",]$MR)
  
  change_result<-data.frame(POP=POP,
                            fpr_change=c(symchange_act_fpr,symchange_rest_fpr,symchange_fpr,asychange_act_fpr,asychange_rest_fpr,asychange_fpr),
                            mr_change=c(symchange_act_mr,symchange_rest_mr,symchange_mr,asychange_act_mr,asychange_rest_mr,asychange_mr),
                            Status=c("ACT","REST","ALL","ACT","REST","ALL"),
                            WARM=c("Symmetry","Symmetry","Symmetry","Asymmetry","Asymmetry","Asymmetry"))
  phys_change<-rbind(change_result,phys_change)
}


phys_change$WARM<-factor(phys_change$WARM,level=c("Symmetry","Asymmetry"))
phys_change$Status<-factor(phys_change$Status,level=c("ACT","REST","ALL"))

write.csv(phys_change,file = "/Users/jiangzhongwen/Desktop/Resting ecology/Data/phys_change.csv",row.names = FALSE)

phys_change<-read.csv("/Users/jiangzhongwen/Desktop/Resting ecology/Data/phys_change.csv")
phys_change$WARM<-factor(phys_change$WARM,level=c("Symmetry","Asymmetry"))
phys_change$Status<-factor(phys_change$Status,level=c("ACT","REST","ALL"))

Figure4_ab<-ggplot(data=phys_change,aes(x=Status,y=fpr_change,fill=POP))+
  geom_bar(stat = 'identity',position = position_dodge(preserve = 'single'),width = 0.5)+
  scale_fill_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  #geom_text(aes(label=round(fpr_change,1)), 
   #         position = position_dodge2(width = 0.9, preserve = 'single'), 
    #        vjust = -0.2, hjust = 0.5,size=8,family="serif",fontface="bold")+
 # annotate("richtext", x = 2, y = 0.08, label = "**Asymmetry warming**",size=8,family="serif",label.colour = NA,fill = NA)+
  labs(x="Status", y = "Change in cumulative digestion (%)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.1,0.9),
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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
        plot.title = element_text(hjust = 0.025,vjust = -5,face = "bold",size = 24,family = "serif"))+
  facet_grid(.~WARM)

Figure4_cd<-ggplot(data=phys_change,aes(x=Status,y=mr_change,fill=POP))+
  geom_bar(stat = 'identity',position = position_dodge(preserve = 'single'),width = 0.5)+
  scale_fill_manual(values=c("#5ED1BD","#5FA3E1","#562abc"))+
  #geom_text(aes(label=round(mr_change,0)), 
   #         position = position_dodge2(width = 0.9, preserve = 'single'), 
    #        vjust = -0.2, hjust = 0.5,size=8,family="serif",fontface="bold")+
  #annotate("richtext", x = 2, y = 0.08, label = "**Asymmetry warming**",size=8,family="serif",label.colour = NA,fill = NA)+
  labs(x="Status", y = "Change in cumulative metabolism (CO2;ml)")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="none",
        legend.text=element_text(size = 24,face=c("bold"),family = "serif"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8),
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
        plot.title = element_text(hjust = 0.025,vjust = -5,face = "bold",size = 24,family = "serif"))+
  facet_grid(.~WARM)



pdf("Figure_4.pdf",height=10,width=14)

ggarrange(Figure4_ab,Figure4_cd,nrow=2,labels=c('a','c'),hjust = -5.5,vjust = 1.5,
          font.label = list(size =30, color = "black",family="serif"))

dev.off()
