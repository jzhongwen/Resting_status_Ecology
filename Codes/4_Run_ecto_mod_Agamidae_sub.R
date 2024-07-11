#####-----2.Run development model and ectotherm model-----#####
#Note: Run ectotherm model to get body temperature and activity state of juveniles for each species at all sites in their distribution range
#for years 1986-2015 and sudo 1986-2015 considering the global mean air temperature becomes 2 degree warmer than pre-industrial period

####Load or install&load packages
library(NicheMapR)
library(plyr)
library(parallel)
ncore=20



loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


######We used the physiological traits of high-elevation population to apply on Agamidae####
#Food passage time:R²=0.83 y=-38.95Tb+0.57Tb²+702.92
#MR:R²=0.77, lgV=0.32lgBM+0.040Tb-1.574

Rate_digest<-function(Tb){
  rate<-(1/(-36.89*Tb+0.54*(Tb)^2+668.57))
  return(rate)
}

Rate_MR<-function(Tb){
  rate<-(((0.027*6.31^0.32))*(10^(0.04*Tb)))
  return(rate)
}

Lizard_traits<-read.csv("/scratch/gpfs/liangm/Resting_status_Jiang/Data/Lizard_traits.csv")
traits_spe<-Lizard_traits[Lizard_traits$Ele=="2600m",]

site_all<-read.csv("/scratch/gpfs/liangm/Resting_status_Jiang/Data/highele_region.csv") 

cal.physiology<-function(micro,traits_spe,year_to_run){
  
    ###Run ectotherm model
    ecto<-try(ectotherm(
        Ww_g	= 6.31, #Wet weight of animal (g), note this model is 'steady state' so no lags in heating/cooling due to mass
        alpha_max	= 0.95, #Maximum solar absorptivity, 0-1
        alpha_min	= 0.95, #Maximum solar absorptivity, 0-1
        T_F_min	= traits_spe$T_F_min, #Minimum foraging temperature, °C (also affects burrow depth selection)
        T_F_max	= traits_spe$T_F_max, #Maximum foraging temperature, °C
        T_B_min	= traits_spe$T_RB_min, #Minimum basking temperature, °C
        T_RB_min	= traits_spe$T_RB_min, #Minimum temperature at which animal will move from retreat to basking site, °C
        T_pref	= traits_spe$T_pref, #Preferred body temperature, °C
        CT_max	= traits_spe$CT_max, #Critical thermal maximum, °C (affects burrow depth selection and may be used to impose death from heat stress)
        CT_min	= traits_spe$CT_min, #Critical thermal minimum, °C (affects burrow depth selection and may be used to impose death from cold stress)
        shade_seek	= 1, #Shade seeking allowed? 1=yes, 0=no
        burrow	= 1, #Shelter in burrow allowed? 1=yes, 0=no
        pct_wet=1,
        mindepth = 1,
        maxdepth = 9,
        aestdepth = 9,
        
        #Environmental inputs:
        minshades =micro[["minshade"]], #Vector of daily minimum shade values - can be different to value used in microclimate model (e.g. to simulate sunspot tracking on a forest floor) (%)
        maxshades =micro[["maxshade"]], #Vector of daily maximum shade values - can be different to value used in microclimate model (e.g. to simulate use of fine-scale shade in a generally unshaded habitat) (%)
        
        alpha_sub = 1-micro[["REFL"]], #Vector of daily substrate reflectances (0-1)?ù̶?ɳ?ӵ?
        DEP = micro$DEP, #Depths available from the microclimate model simulation
        metout = micro$metout, #Microclimate model output for above ground, minimum shade conditions
        shadmet = micro$shadmet, #Microclimate model output for above ground, maximum shade conditions
        soil = micro$soil, #Microclimate model output for soil temperature, minimum shade conditions
        shadsoil = micro$shadsoil, #Microclimate model output for soil temperature, maximum shade conditions
        soilmoist = micro$soilmoist, #Microclimate model output for soil moisture, minimum shade conditions
        shadmoist = micro$shadmoist, #Microclimate model output for soil moisture, maximum shade conditions
        humid = micro$humid, #Microclimate model output for soil humidity, minimum shade conditions
        shadhumid = micro$shadhumid, #Microclimate model output for soil humidity, maximum shade conditions
        soilpot = micro$soilpot, #Microclimate model output for soil water potential, minimum shade conditions
        shadpot = micro$shadpot, #Microclimate model output for soil water potential, maximum shade conditions
        rainfall = as.integer(micro$RAINFALL), #Vector of daily rainfall (mm)
        rainhr = rep(-1,length(micro$metout)), #Vector of hourly rainfall (mm), overwrites rainfall if not negative
        elev = as.numeric(micro$elev), #Elevation of simulation (m), obtained from microclimate model output by default
        longitude = micro$longlat[1], #Longitude (decimal degrees), obtained from microclimate model output by default
        latitude = micro$longlat[2] #Latitude (decimal degrees), obtained from microclimate model output by default
        
        
    ))

    environ<-as.data.frame(ecto$environ)
    
    phy_result<-data.frame(Lon=micro[["longlat"]][1],Lat=micro[["longlat"]][2],year_to_run,environ$DOY,environ$TIME,environ$TC,environ$ACT)
    
    colnames(phy_result)<-c("LON","LAT","YEAR","DOY","TIME","TB","ACT")
    
    
    if(micro[["longlat"]][2]>0){
      ####July for Northern Hemisphere
      phy_result<-phy_result[phy_result$DOY>180&phy_result$DOY<212,]
    }else{
      ####January for Southern Hemisphere
      phy_result<-phy_result[phy_result$DOY>=0&phy_result$DOY<32,]
    }
    
    phy_result["FPR"]<-Rate_digest(phy_result$TB)
    phy_result["MR"]<-Rate_MR(phy_result$TB)
    
    phy_result_year<-data.frame(YEAR=year_to_run,LON=unique(phy_result$LON),LAT=unique(phy_result$LAT),
                                FP=c(sum(phy_result[phy_result$ACT!=0,]$FPR),sum(phy_result[phy_result$ACT==0,]$FPR)),
                                MR=c(sum(phy_result[phy_result$ACT!=0,]$MR),sum(phy_result[phy_result$ACT==0,]$MR)),
                                STATUS=c("ACT","REST"))
    phy_result_year
    
}


  for (year_to_run in 1989:2015) {
    ecto_result<-list()
    micro_spe<-loadRData(paste0("/scratch/gpfs/liangm/Resting_status_Jiang/WD/Intermediate results/Agamidae/0/",year_to_run,".RData"))
    micro<-micro_spe[[1]]
    timer<-vector()
    for (p in 1:floor(length(micro_spe)/ncore)) {
      ptm<-proc.time() 
        result<-mclapply(((p-1)*ncore+1):(p*ncore),function(x){
        micro<-micro_spe[[x]]
        result_i<-try(cal.physiology(micro,traits_spe,year_to_run))
        if(class(result_i)=="try-error") {
          result_i<-NA 
        }
        result_i
      },mc.cores=ncore)
      
      ecto_result<-append(ecto_result,result) 
      timer[p]<-as.numeric((proc.time() - ptm)[3])
      rest<-round(mean(timer,na.rm=T)*(floor(length(micro_spe)/ncore)-p)/3600,digits=2)
      print(paste0(" year_to_run = ",year_to_run,"; ","p = ",p,"; Time left for this year: ",rest," hours"))
    }
    
    ecto_result1<-append(ecto_result,mclapply((floor(length(micro_spe)/ncore)*ncore+1):length(micro_spe),function(i){
      result_i<-try(cal.physiology(x,micro,traits_spe,year_to_run))
      if(class(result_i)=="try-error") {
        result_i<-NA 
      }
      result_i
    },mc.cores=ncore))
    
    ecto_result1<-do.call(rbind,ecto_result1)
    write.csv(ecto_result1, file=paste0("/scratch/gpfs/liangm/Resting_status_Jiang/WD/Intermediate results/Agamidae/",year_to_run,"_physiology.csv"),row.names = FALSE)
    rm(ecto_result)
    rm(ecto_result1) 
    rm(micro_spe)
    rm(micro)
  }

  
  ecto_result1<-read.csv("/scratch/gpfs/liangm/Resting_status_Jiang/WD/Intermediate results/Agamidae/1986_physiology.csv")

  ecto_result2<-data.frame()
  
  for (i in 1:nrow(na.omit(ecto_result1))) {
    
    i_result<-na.omit(ecto_result1)[i,]
    i_result["ID"]<-site_all[site_all$lon==i_result$LON&site_all$lat==i_result$LAT,]$ID
    i_result["ELE"]<-site_all[site_all$lon==i_result$LON&site_all$lat==i_result$LAT,]$ele
    ecto_result2<-rbind(ecto_result2,i_result)
  }  
  
dif_result<-data.frame()

  for(i in unique(ecto_result2$ID)){
    ID_data<-ecto_result2[ecto_result2$ID==i,]
    FP_differ_rest<-mean(abs(ID_data[ID_data$ELE==min(ID_data$ELE)&ID_data$STATUS=="REST",]$FP))-mean(abs(ID_data[ID_data$ELE==max(ID_data$ELE)&ID_data$STATUS=="REST",]$FP))
    FP_differ_act<-mean(abs(ID_data[ID_data$ELE==min(ID_data$ELE)&ID_data$STATUS=="ACT",]$FP))-mean(abs(ID_data[ID_data$ELE==max(ID_data$ELE)&ID_data$STATUS=="ACT",]$FP))
    FP_differ_all<-FP_differ_rest+FP_differ_act
    FP_contri_rest<-FP_differ_rest/FP_differ_act
    
    MR_differ_rest<-mean(abs(ID_data[ID_data$ELE==min(ID_data$ELE)&ID_data$STATUS=="REST",]$MR))-mean(abs(ID_data[ID_data$ELE==max(ID_data$ELE)&ID_data$STATUS=="REST",]$MR))
    MR_differ_act<-mean(abs(ID_data[ID_data$ELE==min(ID_data$ELE)&ID_data$STATUS=="ACT",]$MR))-mean(abs(ID_data[ID_data$ELE==max(ID_data$ELE)&ID_data$STATUS=="ACT",]$MR))
    MR_differ_all<-MR_differ_rest+MR_differ_act
    MR_contri_rest<-MR_differ_rest/MR_differ_act
    
    LON<-mean(ID_data$LON)
    LAT<-mean(ID_data$LAT)
      
    ID_result<-data.frame(YEAR=unique(ID_data$YEAR),ID=i,LON=LON,LAT=LAT,
                          dif=c(FP_differ_rest,FP_differ_act,FP_differ_all,FP_contri_rest,MR_differ_rest,MR_differ_act,MR_differ_all,MR_contri_rest),
                          STATUS=c("REST","ACT","ALL","REST_contri","REST","ACT","ALL","REST_contri"),
                          Traits=c("FP","FP","FP","FP","MR","MR","MR","MR"),
                          ELE_dif=mean(max(ID_data$ELE))-mean(min(ID_data$ELE)))
    dif_result<-rbind(dif_result,ID_result)
  }

library(ggplot2)
Rest_contri_FP<-dif_result[dif_result$Traits=="FP"&dif_result$STATUS=="REST_contri",]
Rest_contri_MR<-dif_result[dif_result$Traits=="MR"&dif_result$STATUS=="REST_contri",]

Rest_contri_FP<-na.omit(Rest_contri_FP)

Rest_FP_Category<-data.frame()


for (i in 1:nrow(Rest_contri_FP)) {
  i_Category<-Rest_contri_FP[i,]
  if(i_Category$dif>1){
    i_Category["Category"]="A"
  }else{
    if(i_Category$dif<=1&i_Category$dif>0){
      i_Category["Category"]="B"
    }else {
      if(i_Category$dif<=0&abs(i_Category$dif)<1){
        i_Category["Category"]="C"
      }else {
        i_Category["Category"]="D"
      }
    }
  } 
  Rest_FP_Category<-rbind(Rest_FP_Category,i_Category)
  
}


Rest_contri_MR<-na.omit(Rest_contri_MR)
Rest_contri_MR<-Rest_contri_MR[Rest_contri_MR$dif<100,]



ggplot()+
  geom_point(data=Rest_FP_Category,aes(LON,LAT,color=Category))+
  scale_colour_manual(values = c("blue","orange","red","green"))


ggplot()+
  geom_point(data=Rest_contri_FP,aes(LON,LAT,color=dif),size)+
  scale_colour_gradient(low = "green", high = "red")

ggplot()+
  geom_point(data=Rest_contri_MR,aes(LON,LAT,color=dif))+
  scale_colour_gradient(low = "green", high = "red")