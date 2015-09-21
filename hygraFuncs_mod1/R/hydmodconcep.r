#' @title Calculate Hydrological Model
#'
#' @description Various single storage buckets conected in 1D with a conceptional approach
#'
#' @param data.in input data as a data.frame with columns structured: $time, $snowH, $precip, $ETo
#' @param parameters set of parameters classifying the different soil layers
#' @param output set output to "absolut": absolute values of each storage at each timestep, "dif": values with differences to value at t0, "meanval": values with diffeence to mean value of each storage compartment.
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 

hydmodconcep <- function(data.in, parameters, output="absolut"){

resultsMOD = data.frame(time=data.in$time, IF=0) #prepare result-file

#### SNOW ####
#SnowHeightT$SS = apply(SnowHeightT, 1, function(x,VOR){if(x[2]==0)0 else{ if(x[4]<0){x[2]*VOR[5]/VOR[2]} else {VOR[5]*x[3]*Pcorr}} }, VOR=SnowHeightT[-1,])
SnowHeight.df = data.frame(h=data.in$snowH) #create temporary dataset for SnowHeight to storage sub-calculation-data
SnowHeight.df$dif = c(0,SnowHeight.df$h[2:length(SnowHeight.df$h)]-SnowHeight.df$h[1:(length(SnowHeight.df$h)-1)]) #calculate SnowHeight decline or incline per timestep
SnowHeight.df$SS = 0; SnowHeight.df$SSdif = 0 

for(i in 2:length(data.in$time)){
  if(SnowHeight.df$h[i] == 0) SnowHeight.df$SS[i] = 0
  else {
    if(SnowHeight.df$dif[i] < 0) SnowHeight.df$SS[i] = SnowHeight.df$h[i] * SnowHeight.df$SS[i-1] / SnowHeight.df$h[i-1]
    else SnowHeight.df$SS[i] = SnowHeight.df$SS[i-1] + data.in$precip[i] *parameters$Pcorr
  }
  SnowHeight.df$SSdif[i] = SnowHeight.df$SS[i] - SnowHeight.df$SS[i-1]
}
#SnowHeight.df$SSdif = c(0,SnowHeight.df$SS[2:length(SnowHeight.df$SS)]-SnowHeight.df$SS[1:(length(SnowHeight.df$SS)-1)]) #calculate SnowHeight decline or incline per timestep
#write to result-file
resultsMOD$snow = SnowHeight.df$SS
resultsMOD$snow[which(is.na(resultsMOD$snow)==T)] = 0 #correct NA-values due to division by "0", in order to keep TS on track
resultsMOD$snowmelt = SnowHeight.df$SSdif
resultsMOD$snowmelt[which(is.na(resultsMOD$snowmelt)==T)] = 0 #correct NA-values due to division by "0", in order to keep TS on track

#### RAIN and SNOWMELT ####
#input
resultsMOD$precipIN = data.in$precip * parameters$Pcorr - resultsMOD$snowmelt
resultsMOD$precipIN[which(resultsMOD$precipIN <0)] = 0 #only use positive values of input
#### SoilMoisture ####
#initial conditions
resultsMOD$SM = parameters$FC
resultsMOD$ETa = 0

for(i in 1:(length(data.in$time)-1)){
  #ETa
  if((resultsMOD$SM[i]/parameters$FC) > parameters$LP){
	  resultsMOD$ETa[i+1] = data.in$ETo[i+1]
	  qexcess = resultsMOD$SM[i]-parameters$FC
	  resultsMOD$SM[i] = parameters$FC
  }
  else{ resultsMOD$ETa[i+1] = data.in$ETo[i+1] * resultsMOD$SM[i]/(parameters$FC*parameters$LP)
  	resultsMOD$SM[i+1] = resultsMOD$SM[i] + resultsMOD$precipIN[i+1] * (1-parameters$alpha) - resultsMOD$ETa[i+1]
 	qexcess=0 
  }
  #infiltration into Soil
  resultsMOD$IF[i+1] = resultsMOD$precipIN[i+1]*parameters$alpha + qexcess #route firstly always part of alpha of input into VsoilS #SoilMoisture if(resultsMOD$SM[i]+resultsMOD$precipIN[i+1]*(1-parameters$alpha) > parameters$FC){ resultsMOD$SM[i+1] = parameters$FC #soil moisture = field capacity resultsMOD$IF[i+1] = resultsMOD$IF[i+1]+resultsMOD$precipIN[i+1]*(1-parameters$alpha) - parameters$FC + resultsMOD$SM[i]
}

#### SOIL ####
#initial conditions
resultsMOD$Q0 = 0
resultsMOD$Qsoil = 0

#First storage (Soil)
for(i in 1:(length(data.in$time)-1)){
  #deltaT = as.numeric(difftime(data.in$time[i+1],data.in$time[i],units="min")) #get time different in between timesteps, format=POSIXct
  deltaT = data.in$time[i+1] - data.in$time[i] #get time different in between timesteps, format=n_days
  if(resultsMOD$Qsoil[i] * parameters$ksoil > parameters$UZL){ #possible check for "full storage" (limit = UZL)
    resultsMOD$Q0[i+1] = resultsMOD$Qsoil[i] - parameters&UZL/parameters$ksoil + resultsMOD$IF[i+1]*(1-exp(-deltaT/parameters$ksoil)) #possibly overland flow !? -> NOT routed anywhere in the hyd.System!
    resultsMOD$Qsoil[i+1] = parameters$UZL/parameters$ksoil * exp(-deltaT/parameters$ksoil)
  }else {resultsMOD$Qsoil[i+1] = resultsMOD$Qsoil[i]*exp(-deltaT/parameters$ksoil) + resultsMOD$IF[i+1]*(1- exp(-deltaT/parameters$ksoil)) 
 	 }
}
#soil storage
resultsMOD$Vsoil = resultsMOD$Qsoil * parameters$ksoil
#initial conditions
resultsMOD$Qsap = 0
resultsMOD$Qgw = 0
#### SAPROLITE $ GROUNDWATER ####
#Other storages (saprolite,groundwater)
for(i in 1:(length(data.in$time)-1)){
  #deltaT = as.numeric(difftime(data.in$time[i+1],data.in$time[i],units="min")) #get time different in between timesteps
  deltaT = data.in$time[i+1] - data.in$time[i] #get time different in between timesteps, format=n_days
  resultsMOD$Qsap[i+1] = resultsMOD$Qsap[i]*exp(-deltaT/parameters$ksap) + resultsMOD$Qsoil[i+1]*(1-exp(-deltaT/parameters$ksap))
  resultsMOD$Qgw[i+1] = resultsMOD$Qgw[i]*exp(-deltaT/parameters$kgw) + resultsMOD$Qsap[i+1]*(1-exp(-deltaT/parameters$kgw))
}

#saprolite storage & Groundwater
resultsMOD$Vsap = resultsMOD$Qsap * parameters$ksap
resultsMOD$Vgw = resultsMOD$Qgw * parameters$kgw

#for (i in 1:(length(data.in$time)-1)){
#  resultsMOD$GWdH[i] = data.in$GWlevel[i] - data.in$GWlevel[1] #delta GWlevel compared to FIRST VALUE OF DATASET
#}

#### storage compartments: WSC and outputs####
# library(reshape2)
# library(ggplot2)
# library(plyr)
# source("/home/mreich/R-gravity-package/hygraFuncs/R/convertzootodataframe.r")
switch(output,
       absolut = {
WaterStorage = data.frame(time=resultsMOD$time, snow=resultsMOD$snow, SM=resultsMOD$SM, soil=resultsMOD$Vsoil, saprolite=resultsMOD$Vsap, GW=resultsMOD$Vgw)
data.plot = mutate(WaterStorage[,-1], time=index(data.timeseries), type="model")
data.melt = melt(data.plot, id.vars=c("time","type") ,measure.vars=c("snow","SM","soil","saprolite","GW"))
WSC.plot = ggplot(data.melt, aes(x=time, y=value)) + geom_line() + facet_grid(variable ~ .) + xlab("") + ylab("Waterstorage in mm")
plot(WSC.plot)
return(WaterStorage)},
       dif = {
WSC = data.frame(time=resultsMOD$time-resultsMOD$time[1], snow=resultsMOD$snow-resultsMOD$snow[1], SM=resultsMOD$SM-resultsMOD$SM[1], soil=resultsMOD$Vsoil-resultsMOD$Vsoil[1], saprolite=resultsMOD$Vsap-resultsMOD$Vsap[1], GW=resultsMOD$Vgw-resultsMOD$Vgw[1])
data.mean.plot = mutate(WSC[,-1], time=index(data.timeseries), type="model")
data.melt = melt(data.mean.plot, id.vars=c("time","type") ,measure.vars=c("snow","SM","soil","saprolite","GW"))
ggplot(data.melt, aes(x=time, y=value)) + geom_line() + facet_grid(variable ~ .) + xlab("") + ylab("Waterstorage in mm")
return(WSC)},
       meanval =  {
WSC.mean = data.frame(time=resultsMOD$time-resultsMOD$time[1], snow=resultsMOD$snow-mean(resultsMOD$snow, na.rm=T), SM=resultsMOD$SM-mean(resultsMOD$SM, na.rm=T), soil=resultsMOD$Vsoil-mean(resultsMOD$Vsoil, na.rm=T), saprolite=resultsMOD$Vsap-mean(resultsMOD$Vsap, na.rm=T), GW=resultsMOD$Vgw-mean(resultsMOD$Vgw, na.rm=T))
data.mean.plot = mutate(WSC.mean[,-1], time=index(data.timeseries), type="model")
data.melt = melt(data.mean.plot, id.vars=c("time","type") ,measure.vars=c("snow","SM","soil","saprolite","GW"))
ggplot(data.melt, aes(x=time, y=value)) + geom_line() + facet_grid(variable ~ .) + xlab("") + ylab("Waterstorage in mm")
return(WSC.mean)}
       ) #end switch
}
