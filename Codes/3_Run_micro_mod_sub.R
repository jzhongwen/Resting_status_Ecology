#####-----3.Run microclimate model-----#####
#Note: Run microcliamte model to get microclimates for all sites that contains lizards (lizards considered in this study)
#for years 1986-2015 and sudo 1986-2015 considering the global mean air temperature becomes 2 degree warmer than pre-industrial period
#For resting ecology

####load or install&load packages
library(NicheMapR)
library(parallel)
library(R.utils)
ncore=20


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
    timeinterval = 365
  ))
  micro[c(1:10,15,16,21,23:25,26,28:30,34,37,38)]
}

####Import a dataframe containing all sites that contains lizards considered in this study
site_all<-read.csv("/scratch/gpfs/liangm/Resting_status_Jiang/Data/highele_region.csv") #ID means the unique ID for each site; lon is the longitude; lat is the latitude

####Run microclimate model

        sce=0
        for(year_to_run in c(1986:2015)){

            micro_list<-list()
            timer<-vector()
            for(p in 1:floor(nrow(site_all)/ncore)){
                ptm<-proc.time() 

                      result<-mclapply(((p-1)*ncore+1):(p*ncore),function(i){
                        
                        result_i<-withTimeout({
                          cal.micro(lon=site_all$lon[i],lat=site_all$lat[i],year=year_to_run,scenario=sce)
                          
                        },timeout = 60,onTimeout = "silent")
                        if(class(result_i)=="NULL") {
                          result_i<-NA 
                        }
                        result_i
                      },mc.cores=ncore)
                 
                micro_list<-append(micro_list,result)
                timer[p]<-as.numeric((proc.time() - ptm)[3])
                rest<-round(mean(timer,na.rm=T)*(floor(nrow(site_all)/ncore)-p)/3600,digits=2)
                print(paste0("sce = ",sce," year_to_run = ",year_to_run,"; ","p = ",p,"; Time left for this year: ",rest," hours"))
            }
      
            micro_list1<-append(micro_list,mclapply((floor(nrow(site_all)/ncore)*ncore+1):nrow(site_all),function(i){
                result_i<-withTimeout({cal.micro(lon=site_all$lon[i],lat=site_all$lat[i],year=year_to_run,scenario=sce)},
                                      timeout = 60,onTimeout = "silent")
                if(class(result_i)=="NULL") {
                    result_i<-NA 
                }
                result_i
            },mc.cores=ncore))
      
            save(micro_list1, file=paste0("/scratch/gpfs/liangm/Resting_status_Jiang/WD/Intermediate results/Agamidae/",sce,"/",year_to_run,".RData"))
            rm(micro_list)
            rm(micro_list1)
        }
  