#' @title Reading field measured data from BKG station Wettzell
#'
#' @description Reading measurement data from field various field equipment from the BKG Observatory Wettzell (Bavaria).
#' Available dataset include:
#' 
#' \itemize{
#' \item soil data (5 TDR clusters each with approx. 40 sensors)
#' \item snow data
#' \item lysimeter data (potential evaporation, weight data from the lysimeter)
#' \item precipition data from 4 locations ARROUND the observatory
#' \item groundwater data from 10 groundwater wells
#' \item clima data Wettzell (official data from the BKG)
#' \item clima data GFZ (data collected locally operated by GFZ)
#' \item runoff data from 4 locations arround the observatory
#' }
#'
#' @param param insert one of the following options to receive the newest of this data:
#' "Soil", "Snow", "LysiEto", "LysiWeight", "Precip", "GW", "ClimaBKG", ClimaGFZ", "Runoff"

#' @details Data represents NON filtered raw-data directly from each individual measurement equipment loggers.
#' If not other stated in SI-units.
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples readWettzell("Soil")
#' 
read_wettzell <- function(
  param = "Soil" #"Snow", "LysiEto", "LysiWeight", "Precip", "GW", "ClimaBKG", ClimaGFZ", "Runoff"
  )
  {
  #home directory
  home="/home/mreich/server/hygra/DataWettzell/"
  ## insert path to SERVER !!
  #load packages
# library(zoo)
  Sys.setenv(TZ = "GMT") #set to get no errors when merging zoo-objects !!
  
  #distinguish parameter
  switch(param,
         Soil={ #read soil data TS
            setwd(paste(home,"SoilMoisture/Cluster_Data/",sep="")) #working directory soil moisture
            name1=load(file="./SG_complete/SGold_raw.rdata",.GlobalEnv)
            name2=load(file="./Hang_complete/Hang_raw.rdata",.GlobalEnv)
            name3=load(file="./Senke_complete/Senke_raw.rdata",.GlobalEnv)
            name4=load(file="./Mast_complete/Mast_raw.rdata",.GlobalEnv)
            name5=load(file="./SGnew_complete/SGnew_raw.rdata",.GlobalEnv)
            dataread= c("Data read:",name1,name2,name3,name4,name5)
            return(dataread)
         },# end Soil
         Snow={
           setwd(paste(home,"Snow/data_complete/",sep=""))
           name1=load(file=".rdata",.GlobalEnv)
           dataread= c("Data read:",name1)
           return(dataread)
         },
         LysiEto={
           setwd(paste(home,"Lysimeter/data_complete/",sep=""))
           name1=load(file="eto_lysimeter_agg.rdata",.GlobalEnv)
           dataread= c("Data read:",name1)
           return(dataread)
         },
         LysiWeight={
           setwd(paste(home,"Lysimeter/data_complete/",sep=""))
           name1=load(file="weight_lysimeter_raw.rdata",.GlobalEnv)
           dataread= c("Data read:",name1)
           return(dataread)
         },
         Precip={
           setwd(paste(home,"Precipitation/",sep=""))
           name1=load(file=".rdata",.GlobalEnv)
           dataread= c("Data read:",name1)
           return(dataread)
         },
         GW={
           setwd(paste(home,"Groundwater/GW_Pegel_complete/",sep=""))
           #NON pivot-format data!?! name1=load(file="groundwater_corrected_pivot.rdata",.GlobalEnv)
           name1=load(file="groundwater_corrected_pivot.rdata",.GlobalEnv)
           dataread= c("Data read:",name1)
           return(dataread)
         },
         ClimaBKG={
           setwd(paste(home,"Climate/",sep=""))
           name1=load(file="./30min_wettzell/clima30min_raw.rdata",.GlobalEnv)
           dataread= c("Data read:",name1)
           return(dataread)
         },
         ClimaGFZ={
           setwd(paste(home,"Lysimeter/data_complete/",sep=""))
           name1=load(file="clima_lysimeter_raw.rdata",.GlobalEnv)
           dataread= c("Data read:",name1)
           return(dataread)
         },
         Runoff={
           setwd(paste(home,"Runoff/data_complete/",sep=""))
           name1=load(file="augraben_raw.rdata",.GlobalEnv)
           name2=load(file="pegelHoellensteinSEBA_raw.rdata",.GlobalEnv)
           name3=load(file="wehrKELLER_raw.rdata",.GlobalEnv)
           name4=load(file="wehrSEBA_raw.rdata",.GlobalEnv)
           dataread= c("Data read:",name1,name2,name3,name4)
           return(dataread)
         }
         )  
}

### end readWettzell ###

