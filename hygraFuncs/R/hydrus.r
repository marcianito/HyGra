######################################################################
### functions concerning hydrological software "hydrus" (1D,2D,3D) ###
######################################################################

#' @title Data preparation for HYDRUS
#'
#' @description dem TS 
#'
#' @param inputpath DEM-data (realtive or absolute paths)
#' @param inputfile data.frame containing information about the DEM (row wise): columns, rows, starting x, starting y, lengths of one cellsize
#' @param dz vector containing the thickness of each subsurface layer
#' @param limit.area vector of coordinates declaring the area of interest within the inputfile (DEM). Structure: X(left), Y(top), X(right), Y(bottom).
#' 
#' @details Both DEM-data and its information file can be created using readDEM(). 
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing 
#' @export

hydrus_dem <- function(inputpath,inputfile,dz,limit.area){
  ##loading until package status..
      library(reshape2)
  #adjust DEM for hydrus input (local coordinates, vertical layer adjustments, etc...)
      readDEM(inputpath, inputfile)
      discre(dem.info[1,],dem.info[2,],dem.info[3,],dem.info[4,],dem.info[5,])
      #change DEM cols & rows to coordinates
      colnames(dem) = xp
      rownames(dem) = yp
      YXZ=melt(as.matrix(dem))
      XYZ = cbind(X=YXZ[,2],Y=YXZ[,1],Z=YXZ[,3])
      #XYZ.hydrus = cbind(X=(YXZ[,2]-min(YXZ[,2])),Y=(YXZ[,1]-min(YXZ[,1])),Z=YXZ[,3]) #real XYZ coordinates
      #XYZlayers.hydrus = cbind(X=(YXZ[,2]-min(YXZ[,2])),Y=(YXZ[,1]-min(YXZ[,1])),Z0=(min(YXZ[,3])-dzmax),Z1=(min(YXZ[,3])-dzmax), Zmax=YXZ[,3]) #modified coordinates for hydrus input
      outside.area=which((XYZ[,1]<limit.area[1] & XYZ[,2]>limit.area[2]) | (XYZ[,1]>limit.area[3] & XYZ[,2]<limit.area[4]))#, arr.ind=T)
      outside.area=which((XYZ[,1]<limit.area[1] | XYZ[,2]>limit.area[2]) | (XYZ[,1]>limit.area[3] | XYZ[,2]<limit.area[4]))#, arr.ind=T)
      XYZ.filter = XYZ[-outside.area,]
      #outside.X=which(XYZ[,1]<limit.area[1] | XYZ[,1]>limit.area[3])
      #outside.Y=which(XYZ[,2]>limit.area[2] | XYZ[,2]<limit.area[4])
      dzmax=sum(dz)
      #XYZlayers.hydrus = cbind(X=(XYZ.filter[,1]-min(XYZ[,1])),Y=(XYZ.filter[,2]-min(XYZ[,2])),Zbase=(min(XYZ.filter[,3])-dzmax),ZGW=(min(XYZ.filter[,3])-dzmax)+dz[4],Zsap=(min(XYZ.filter[,3])-dzmax)+sum(dz[3:4]),Zsoil=(min(XYZ.filter[,3])-dzmax)+dz[2:4],Zmax=XYZ.filter[,3]) #modified coordinates for hydrus input
      XYZlayers.hydrus = cbind(X=(XYZ.filter[,1]-min(XYZ[,1])),Y=(XYZ.filter[,2]-min(XYZ[,2])),Zbase=(min(XYZ.filter[,3])-dzmax))
      for(i in length(dz):1){
       XYZlayers.hydrus = cbind(XYZlayers.hydrus, (XYZ.filter[,3]-sum(dz[1:i])))
        #,ZGWlow=XYZ.filter[,3]-sum(dz[1:4]),ZGW=XYZ.filter[,3]-sum(dz[1:3]),Zsap=XYZ.filter[,3]-sum(dz[1:2]),Zsoil=XYZ.filter[,3]-dz[1],Zmax=XYZ.filter[,3]) #modified coordinates for hydrus input 
      }
      XYZlayers.hydrus = cbind(XYZlayers.hydrus, Zmax=XYZ.filter[,3])
      return(XYZlayers.hydrus)      
}

#' @title Data preparation for HYDRUS
#'
#' @description Create evapotranspiration time series 
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_eto <- function(date.in, date.out, freq="hour", aszoo=F){
      ## später mit readwettzell daten laden/kombinieren !!!
      #ETO from lysimeter data
      library(zoo); Sys.setenv(TZ = "GMT")
      load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/eto_lysimeter_agg.rdata")
      #load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/eto_lysimeter_raw.rdata")
      eto = eto.agg[which(index(eto.agg)==date.in):which(index(eto.agg)==date.out)]
      eto.agg = aggregate(eto, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),sum, na.rm=T) #autocorrection of non-rounded timestampts
      if(aszoo==F){
      return(zootodf(eto.agg))
      }
      else {return(eto.agg)}
}

#' @title Data preparation for HYDRUS
#'
#' @description Create lysimeter drainage time series
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param area bottom surface area of the lysimeter. In m².
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#'
#' @details The outout will be given in meters drained water column per timestep (differences).
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_lysidrain_ben <- function(date.in, date.out, freq="hour", area=1,winlength=10, aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT")
# load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/weight_lysimeter_raw.rdata")
      load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/weight_lysimeter20082011_raw.rdata")
      #adjust script to drainage from lysimeter
      lysimeter = weight_lysimeter_20082011[which(index(weight_lysimeter_20082011)==date.in):which(index(weight_lysimeter_20082011)==date.out)]
      # drain_kg = lysimeter$SiwaWaage #weight in[kg] of drained water
      drain_kg = lysimeter$SiwaWaageCUM #weight in[kg] of drained water
      #change to data.frame
      drain_kg_df = zootodf(drain_kg)
      #####
      #set parameters
      # lag = laglength/length(drain_kg_df$value)
      # robustness = 5
      win=rep(1/(winlength+1),(winlength+1))
      ## remove NA-values
      drain_kg_df$value[which(is.na(drain_kg_df$value)==T)] = 0
      ## smooth using rlowess
      # drain_kg_df$siwa_smooth = lowess(drain_kg_df$value, f=lag, iter=robustness)
      # drain_kg_df$siwa_smooth = lowess(drain_kg_df$time, drain_kg_df$value, f=lag, iter=robustness)$y
      # drain_kg_df$siwa_mv = rollmean(drain_kg_df$value, f=lag, iter=robustness)
      drain_kg_df$siwa_mv = stats::filter(coredata(drain_kg_df$value),win,sides=2)
      ## interpolate missing time steps
      # drain_kg_df$siwa = na.approx(drain_kg_df$siwa_mv, na.rm=T)
      ######
      #shift data one timestep
      drain_kg_df$siwaSHIFT = c(drain_kg_df$siwa_mv[-1],NA) #create new column with values shifted one timestep
      #calculate difference (t+1) - t
      drain_kg_df$siwaDIF = -1*(drain_kg_df$siwaSHIFT - drain_kg_df$siwa_mv)
      #eliminate NA-values
      # drain_kg_df$SiwaWaageDIF[which(is.na(drain_kg_df$SiwaWaageDIF)==T)] = 0
      #correct the times of pumping out the water tank
			# drain_kg_df[which(drain_kg_df$SiwaWaageDIF < -0.5)] = 0
      #change timeseries back to zoo
      drain_dif_kg = zoo(drain_kg_df$siwaDIF, order.by=drain_kg_df$time)
      #convert to [m] water column
      density_water = 999.9720 #[kg/m³]
      lysi_drain = (drain_dif_kg/lysi_area)*(1/density_water)
      # drain.agg = lysi_drain
      drain.agg = aggregate(lysi_drain, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),sum, na.rm=T) #autocorrection of non-rounded timestampts
      if(aszoo==F){
      return(zootodf(drain.agg))
      }
      else {return(drain.agg)}
}
#' @title Data preparation for HYDRUS
#'
#' @description Create lysimeter drainage time series
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param area bottom surface area of the lysimeter. In m².
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#'
#' @details The outout will be given in meters drained water column per timestep (differences).
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_lysidrain <- function(date.in, date.out, freq="hour", area=1, aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT")
      load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/weight_lysimeter_raw.rdata")
      #adjust script to drainage from lysimeter
      lysimeter = weight.lysimeter[which(index(weight.lysimeter)==date.in):which(index(weight.lysimeter)==date.out)]
      drain_kg = lysimeter$SiwaWaage #weight in[kg] of drained water
      #aggregate to hourly data before processing
      # drain_kg = aggregate(drain_kg, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),mean, na.rm=T) #autocorrection of non-rounded timestampts
      #change to data.frame
      drain_kg_df = zootodf(drain_kg)
      #shift data one timestep
      drain_kg_df$SiwaWaageSHIFT = c(drain_kg_df$value[-1],NA) #create new column with values shifted one timestep
      #calculate difference (t+1) - t
      drain_kg_df$SiwaWaageDIF = -1*(drain_kg_df$SiwaWaageSHIFT - drain_kg_df$value)
      #eliminate tank pumping times 
      drain_kg_df$SiwaWaageDIF[which(drain_kg_df$SiwaWaageDIF > 0.22)] = 0 #0.5; opt=0.22
      #eliminate unrealistic tank increases (too much drainage?) 
      # check value for reasonability !!!
      drain_kg_df$SiwaWaageDIF[which(drain_kg_df$SiwaWaageDIF < -0.22)] = 0 #-0.1; opt=-0.22
      #eliminate NA-values
      drain_kg_df$SiwaWaageDIF[which(is.na(drain_kg_df$SiwaWaageDIF)==T)] = 0
      #correct the times of pumping out the water tank
			# drain_kg_df[which(drain_kg_df$SiwaWaageDIF < -0.5)] = 0
      #change timeseries back to zoo
      drain_dif_kg = zoo(drain_kg_df$SiwaWaageDIF, order.by=drain_kg_df$time)
      #convert to [m] water column
      density_water = 999.9720 #[kg/m³]
      lysi_drain = (drain_dif_kg/lysi_area)*(1/density_water)
      # drain.agg = lysi_drain
      drain.agg = aggregate(lysi_drain, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),sum, na.rm=T) #autocorrection of non-rounded timestampts
      if(aszoo==F){
      return(zootodf(drain.agg))
      }
      else {return(drain.agg)}
}
#' @title Data preparation for HYDRUS
#'
#' @description Create precipitation time series
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#'
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_precip <- function(date.in, date.out, freq="hour", aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT")
      load(file="/home/mreich/server/hygra/DataWettzell/Climate/30min_wettzell/clima30min_raw.rdata")
      clima.all = clima.raw[which(index(clima.raw)==date.in):which(index(clima.raw)==date.out)]
      precip = (clima.all$Prec_Sen1)
      precip.agg = aggregate(precip, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),sum, na.rm=T) #autocorrection of non-rounded timestampts
      if(aszoo==F){
      return(zootodf(precip.agg))
      }
      else {return(precip.agg)}
}

#' @title Data preparation for HYDRUS
#'
#' @description groundwater TS 
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param name name of groundwater well to use. Possible so far "BK1", "BK2", "BK3", "BK14".
#' @param depthLB depths of lower model boundary. All groundwater level values will be subtracted from this value.
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#' @param aprox if the time series has some NA-values, they can be approximized. Default value is TRUE.
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_gw <- function(date.in, date.out, name, depthLB, freq="hour",aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT") #;library(xts)
      gw.name = load(file=paste("/home/mreich/server/hygra/DataWettzell/Groundwater/GW_Pegel_complete/",name,"_corrected.rdata",sep=""))
      gw.agg = aggregate(get(gw.name), function(x) as.POSIXct(trunc(x, freq),ts="GMT"),mean, na.rm=T) #autocorrection of non-rounded timestampts
      #gw.agg = period.apply(get(gw.name), endpoints(get(gw.name), freq), mean, na.rm=T)
      date.in.freq = as.POSIXct(trunc(date.in, freq),ts="GMT")
      date.out.freq = as.POSIXct(trunc(date.out, freq),ts="GMT")
      #gw = gw.agg[which(index(gw.agg)==date.in):which(index(gw.agg)==date.out)]
      gw = gw.agg[which(index(gw.agg)==date.in.freq):which(index(gw.agg)==date.out.freq)]
      gw.corrected = depthLB - gw #recalculate so output TS is in GW-table in reference to lower model boundary (LB)
      if(aszoo==F){
      return(zootodf(gw.corrected))
      }
      else {return(gw.corrected)}
}

#' @title Data preparation for HYDRUS
#'
#' @description soil moisture TS 
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param name name of groundwater well to use. Possible so far "BK1", "BK2", "BK3", "BK14".
#' @param depthLB depths of lower model boundary. All groundwater level values will be subtracted from this value.
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#' @param aprox if the time series has some NA-values, they can be approximized. Default value is TRUE.
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_sm <- function(date.in, date.out, name, datacol, freq="hour",aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT") #;library(xts)
      sm.name = load(file=paste("/home/mreich/server/hygra/DataWettzell/SoilMoisture/Cluster_Data/Data_filtered/",name,"_filtered_6hourmean.rdata",sep=""))
      sm.agg = aggregate(get(sm.name)[,datacol], function(x) as.POSIXct(trunc(x, freq),ts="GMT"),mean, na.rm=T) #autocorrection of non-rounded timestampts
      date.in.freq = as.POSIXct(trunc(date.in, freq),ts="GMT")
      date.out.freq = as.POSIXct(trunc(date.out, freq),ts="GMT")
      #gw = gw.agg[which(index(gw.agg)==date.in):which(index(gw.agg)==date.out)]
      sm = sm.agg[which(index(sm.agg)==date.in.freq):which(index(sm.agg)==date.out.freq)]
      if(aszoo==F){
      return(zootodf(sm))
      }
      else {return(sm)}
}
#' @title Read Hydrus 1D output data
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param vertical_nodes mumber of nodes used in the vertical model discretisation
#' @param timesteps number of modeled timesteps
#' @param t_printout number of model printout times
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param datetime vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 

read_hydrus1d <- function(folder, timesteps, vertical_nodes, t_printout, plotting=F, timestamps){
	# library(dplyr)
	# library(dplyrExtras)

 path_hydrus = "/home/mreich/.wine/drive_c/users/Public/Documents/PC-Progress/Hydrus-1D 4.xx/Examples/Direct/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read output files
 num_lines = length(readLines("T_Level.out"))
 tlevel = read.table(file="T_Level.out", header=F, skip=9, sep="",nrows=(num_lines-10), dec=".")
 colnames(tlevel) = c("Time","rTop","rRoot","vTop","vRoot","vBot","sum(rTop)","sum(rRoot)","sum(vTop)","sum(vRoot)","sum(vBot)","hTop","hRoot","hBot","RunOff","sum(RunOff)","Volume","sum(Infil)","sum(Evap)","TLevel","Cum(WTrans)","SnowLayer")
 #remove duplicated data lines (interpolation artefacts probably)
 tlevel$Timecut=round(tlevel$Time,0)
 tlevel = tlevel[!duplicated(tlevel$Timecut),]

 nodal_infoTIME = read.table(file="Nod_Inf.out", header=F, skip=7, sep="", nrows=1, dec=".")
 nodal_infoIN = read.table(file="Nod_Inf.out", header=F, skip=13, sep="", nrows=vertical_nodes, dec=".")
 nodal_info = cbind(rep(nodal_infoTIME[,2],times=vertical_nodes), nodal_infoIN)
 balanceTIME = read.table(file="Balance.out", header=F, skip=20, sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=30, sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=31, sep="", nrows=1, dec=".")
 balance = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])

 for(i in 2:(t_printout+1)){
 nodal_infoTIME = read.table(file="Nod_Inf.out", header=F, skip=((i-1)*(vertical_nodes+9) + 7), sep="", nrows=1, dec=".")
 nodal_infoIN = read.table(file="Nod_Inf.out", header=F, skip=((i-1)*(vertical_nodes+9) + 13), sep="", nrows=vertical_nodes, dec=".")
 nodal_infoDATA = cbind(rep(nodal_infoTIME[,2],times=vertical_nodes), nodal_infoIN)
 nodal_info = rbind(nodal_info, nodal_infoDATA)
 if(i == (t_printout+1)) break
 balanceTIME = read.table(file="Balance.out", header=F, skip=((i-1)*15+20), sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=((i-1)*15+30), sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=((i-1)*15+31), sep="", nrows=1, dec=".")
 balanceDATA = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
 balance = rbind(balance, balanceDATA)
 }
 colnames(nodal_info) = c("Time","Node","Depth","Head","Moisture","K","C","Flux","Sink","Kappa","v/KsTop","Temp")
 colnames(balance) = c("Time","WaterBalanceVolume","WaterBalancePercent")
 tlevel = mutate(tlevel, datetime=timestamps)
 tlevel <<- tlevel
 nodal_info <<- nodal_info
 balance <<- balance

 #plotting
 if(plotting){
# library(ggplot2)
# library(reshape2)
# library(dplyrExtras)
# library(grid)
# library(gridExtra)
# library(scales)
 tlevel.plot = melt(tlevel, id.vars=c("datetime","Time"), measure.vars=colnames(tlevel)[2:22], variable.name = "Parameters")
 nodal.plot = melt(nodal_info, id.vars=c("Time", "Depth"), measure.vars=colnames(nodal_info)[-1:-3], variable.name = "Parameters")

 printtimes = unique(nodal.plot$Time) 
 printdates = c(datetime[1],datetime[printtimes])
 dates_times = data.frame(Time=printtimes, Dates=format(printdates, "%Y-%m-%d"))

 xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(t_printout+1)
 #  xcolors = colorRampPalette(c("blue","green"))(t_printout+1)

 #  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 #subsetting to choose standard parameters
 choose.params.tlevel = c("vTop", "vRoot", "vBot", "hTop", "hRoot", "hBot")
 units = data.frame(type=c("Flux","Head"), Unit=c("[m/h]","[m]"))
 cols = data.frame(Parameters = choose.params.tlevel, Boundary = c("Top","Root","Bottom"))
 type = data.frame(Parameters = c("vTop","vBot","hTop","hBot"), type=c("Flux","Flux","Head","Head"))
 # tlevel.plot.filter = 	filter(tlevel.plot, grepl("vTop|vBot|hTop|hBot",Parameters)) %>%
 #                         filter(Parameters != "sum(vTop)" & Parameters != "sum(vRoot)" & Parameters != "sum(vBot)") %>%
 #                         mutate(type="Flux") %>%
 #                         mutate_if(Parameters == "hTop" | Parameters == "hRoot" | Parameters == "hBot", type=="Head") %>%
 #                         inner_join(units,by="type") %>%
 #                        mutate(typeUnits = paste(type,Unit)) %>%
 #                        inner_join(cols, by="Parameters")
 # browser()			
 tlevel.plot.filter = filter(tlevel.plot, grepl("vTop|vBot|hTop|hBot",Parameters)) %>%
 			inner_join(type,by="Parameters") %>%
 			inner_join(units,by="type") %>%
			mutate(typeUnits = paste(type,Unit)) %>%
			inner_join(cols, by="Parameters")

			# tlevel.plot.filter$ParamUnits = factor(tlevel.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 tlevel.plot.filter$Boundary = factor(tlevel.plot.filter$Boundary,levels=c("Top","Root","Bottom"))
 #  tlevel.gg = ggplot(tlevel.plot.filter, aes(x=datetime, y=value, colour=Boundary)) + geom_line() + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates), colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #  tlevel.gg = ggplot(tlevel.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 tlevel.gg = ggplot(tlevel.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour=xcolors, size=2, alpha=0.5) + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 choose.params.nodal = c("Head", "Moisture", "K", "Flux")
 units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]", "[m/h]", "[m/h]"))
 #  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture|K|Flux",Parameters)) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
 			inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") 
			
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))

 #  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 #merge both plots together
 grid.arrange(tlevel.gg,nodal.gg,main=folder)
}
 return(print("Output stored in variables tlevel, nodal_info and balance"))
}

#' @title Read Hydrus 2D output data (rectangular mesh)
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param vertical_nodes number of nodes used in the vertical model discretisation
#' @param horizontal_nodes number of nodes used in the horizontal model discretisation
#' @param timesteps number of modeled timesteps
#' @param t_printout number of model printout times
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' @param px number of horizontal node, where the vertical profile should be analyzed
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 

read_hydrus2d_rect <- function(folder, timesteps, vertical_nodes, horizontal_nodes, t_printout, plotting=F, timestamps, px){
	# library(dplyr)
	# library(dplyrExtras)

 path_hydrus = "/home/mreich/server/marvin_reich/hydrus2d/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read hydrus 2d output files
 cord_depthIN = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=vertical_nodes*horizontal_nodes, dec=".") #read z coordinates of grid
 cord_depth = unique(cord_depthIN[,3]);
 cord_depth = cord_depth - max(cord_depth)
 h_num_lines = length(readLines("h_Mean.out"))
 v_num_lines = length(readLines("v_Mean.out"))
 h_Mean = read.table(file="h_Mean.out", header=F, skip=6, sep="",nrows=(h_num_lines-7), dec=".")
 v_Mean = read.table(file="v_Mean.out", header=F, skip=13, sep="",nrows=(v_num_lines-14), dec=".")
 colnames(h_Mean) = c("Time","hAtm","hRoot","hKode3","hKode1","hSeep","hKode5","hKode6","hKode7","hKode8","hKode9")
 colnames(v_Mean) = c("Time","rAtm","rRoot","vAtm","vRoot","vKode3","vKode1","vSeep","vDrain","vBottom","vKode7","vKode8","vKode9","RunOff","Evapor","Infiltr","SnowLayer")
 #remove duplicated data lines (interpolation artefacts probably)
 h_Mean$Timecut=round(h_Mean$Time,0)
 v_Mean$Timecut=round(v_Mean$Time,0)
 h_Mean = h_Mean[!duplicated(h_Mean$Timecut),]
 v_Mean = v_Mean[!duplicated(v_Mean$Timecut),]
 #joint datasets h_Mean and v_Mean
 h_v_Mean = inner_join(h_Mean, v_Mean, by="Time")
 h_v_Mean = mutate(h_v_Mean, datetime=timestamps) #add real dates along with timesteps

 #generate sequence of vertical profiles to read out; each vertical profile needs their own sequence!!
 vert1 = seq(px,by=horizontal_nodes,length.out=vertical_nodes)
 
 th_h_profiles = data.frame()
 for(i in 0:t_printout){
 TH_infoTIME = read.table(file="TH.TXT", header=F, skip=1 + i*(ceiling(horizontal_nodes*vertical_nodes/10) + 3), sep="", nrows=1, dec=".")[,3]
 TH_infoIN = scan(file="TH.TXT", skip=3+ i*(ceiling(horizontal_nodes*vertical_nodes/10) + 3), sep="", nmax=vertical_nodes*horizontal_nodes, dec=".")
 H_infoTIME = read.table(file="H.TXT", header=F, skip=1+ i*(ceiling(horizontal_nodes*vertical_nodes/10) + 3), sep="", nrows=1, dec=".")[,3]
 H_infoIN = scan(file="H.TXT", skip=3+ i*(ceiling(horizontal_nodes*vertical_nodes/10) + 3), sep="", nmax=vertical_nodes*horizontal_nodes, dec=".")
 #fetch times and values of each corresponding vertical profile! 
 TH_profile1 = TH_infoIN[vert1]
 H_profile1 = H_infoIN[vert1]
 profiles = data.frame(Time=H_infoTIME, x= px, z= (1:vertical_nodes), Depth= cord_depth, Head= H_profile1, Moisture=TH_profile1)
 th_h_profiles = rbind(th_h_profiles,profiles)
 }

 balanceTIME = read.table(file="Balance.out", header=F, skip=21, sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=29, sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=30, sep="", nrows=1, dec=".")
 balance = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])

 for(i in 2:(t_printout+1)){
 if(i == (t_printout+1)) break
 balanceTIME = read.table(file="Balance.out", header=F, skip=((i-1)*13+21), sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=((i-1)*13+29), sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=((i-1)*13+30), sep="", nrows=1, dec=".")
 balanceDATA = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
 balance = rbind(balance, balanceDATA)
 }
 colnames(balance) = c("Time","WaterBalanceVolume","WaterBalancePercent")
 h_Mean <<- h_Mean
 v_Mean <<- v_Mean
 th_h_profiles <<- th_h_profiles 
 balance <<- balance

 #plotting
 if(plotting){
# library(ggplot2)
# library(reshape2)
# library(dplyrExtras)
# library(grid)
# library(gridExtra)
# library(scales)
 h_v_Mean.plot = melt(h_v_Mean, id.vars=c("datetime","Time"), measure.vars=colnames(h_v_Mean)[2:22], variable.name = "Parameters")
 nodal.plot = melt(th_h_profiles, id.vars=c("Time", "Depth"), measure.vars=colnames(th_h_profiles)[-1:-4], variable.name = "Parameters")

 printtimes = unique(th_h_profiles$Time)
 #  printtimes = unique(balance[,1]) 
 printdates = c(datetime[1],datetime[printtimes])
 #  printdates = c(datetime[printtimes+1], datetime[timesteps])
 dates_times = data.frame(Time=printtimes, Dates=format(printdates, "%Y-%m-%d:%H"))

 xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(t_printout+1)
 #  xcolors = colorRampPalette(c("blue","green"))(t_printout+1)

 #  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 #subsetting to choose standard parameters
 choose.params.h_v_Mean = c("vAtm", "vRoot", "vKode3", "hAtm", "hRoot", "hKode3")
 units = data.frame(type=c("Flux","Head"), Unit=c("[m²/h]","[m]"))
 type = data.frame(Parameters = c("vAtm","vKode3","hAtm","hKode3"), type=c("Flux","Flux","Head","Head"))
 cols = data.frame(Parameters = choose.params.h_v_Mean, Boundary = c("Top","Root","Bottom"))
 #  h_v_Mean.plot.filter = 	filter(h_v_Mean.plot, Parameters == choose.params.h_v_Mean) %>%
 # browser()
 # h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 #                         filter(Parameters != "sum(vAtm)" & Parameters != "sum(vRoot)" & Parameters != "sum(vKode3)") %>%
 #                         mutate(type="Flux") %>%
 #                         mutate_if(Parameters == "hAtm" | Parameters == "hRoot" | Parameters == "hKode3", type=="Head") %>%
 #                         inner_join(units,by="type") %>%
 #                        mutate(typeUnits = paste(type,Unit)) %>%
 #                        inner_join(cols, by="Parameters")
			
 h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 			inner_join(type,by="Parameters") %>%
 			inner_join(units,by="type") %>%
			mutate(typeUnits = paste(type,Unit)) %>%
			inner_join(cols, by="Parameters")

			# h_v_Mean.plot.filter$ParamUnits = factor(h_v_Mean.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 h_v_Mean.plot.filter$Boundary = factor(h_v_Mean.plot.filter$Boundary,levels=c("Top","Root","Bottom"))

 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value, colour=Boundary)) + geom_line() + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates), colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour=xcolors, size=2, alpha=0.5) + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 
 choose.params.nodal = c("Head", "Moisture")
 units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]"))
 #  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture",Parameters)) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
 			inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") 
			
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]"))

 #  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 #merge both plots together
 grid.arrange(h_v_Mean.gg,nodal.gg,main=paste(folder,": vertical profile at x =", px, sep=" "))
}
 return(print("Output stored in variables h_Mean, v_Mean, th_h_profiles and balance"))
}


#' @title Read Hydrus 2D output data
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param timesteps number of modeled timesteps
#' @param t_printout string of dates of time series where informations are printed
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' @param px number of horizontal node, where the vertical profile should be analyzed
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 

# read_hydrus2d_irregular <- function(folder, timesteps, vertical_nodes, horizontal_nodes, t_printout, plotting=F, timestamps, px){
read_hydrus2d <- function(folder, timesteps, factorforprintout, plotting=F, timestamps, profile_loc, vertical=T){
	# library(dplyr)
	# library(dplyrExtras)

 #path_hydrus = "/home/mreich/server/marvin_reich/hydrus2d/experiments/" 
 #path_hydrus = "/home/mreich/server/hydro72/hydrus2d/experiments/" 
 path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus2d/experiments/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read hydrus 2d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="", nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,2]
 cord_depthIN = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid
 #order_x = cord_depthIN[,2]
 #order_z = cord_depthIN[,3]
 grid_cords = data.frame(x=cord_depthIN[,2] ,z=cord_depthIN[,3])
 # cord_depth = unique(cord_depthIN[,3]);
 # cord_depth = cord_depth - max(cord_depth)
 h_num_lines = length(readLines("h_Mean.out"))
 v_num_lines = length(readLines("v_Mean.out"))
 h_Mean = read.table(file="h_Mean.out", header=F, skip=6, sep="",nrows=(h_num_lines-7), dec=".")
 v_Mean = read.table(file="v_Mean.out", header=F, skip=13, sep="",nrows=(v_num_lines-14), dec=".")
 colnames(h_Mean) = c("Time","hAtm","hRoot","hKode3","hKode1","hSeep","hKode5","hKode6","hKode7","hKode8","hKode9")
 colnames(v_Mean) = c("Time","rAtm","rRoot","vAtm","vRoot","vKode3","vKode1","vSeep","vDrain","vBottom","vKode7","vKode8","vKode9","RunOff","Evapor","Infiltr","SnowLayer")
 #remove duplicated data lines (interpolation artefacts probably)
 h_Mean$Timecut=round(h_Mean$Time,0)
 v_Mean$Timecut=round(v_Mean$Time,0)
 h_Mean = h_Mean[!duplicated(h_Mean$Timecut),]
 v_Mean = v_Mean[!duplicated(v_Mean$Timecut),]
 #joint datasets h_Mean and v_Mean
 h_v_Mean = inner_join(h_Mean[-1,], v_Mean[-1,], by="Time")
 h_v_Mean = mutate(h_v_Mean, datetime=timestamps) #add real dates along with timesteps

 #generate sequence of vertical profiles to read out; each vertical profile needs their own sequence!!
 #vertiklales profil bei x=2 & x=10m
 #vert1 = seq(profile_loc,by=horizontal_nodes,length.out=vertical_nodes)
 
system("sed -n '/^[[:space:]]*[T]/p' TH.TXT > bash_time.out")
system("sed -n '/^[[:space:]]*[0-9]/p' TH.TXT > bash_dataTH.out")
system("sed -n '/^[[:space:]]*-*[0-9]/p' H.TXT > bash_dataH.out")

printout = read.table(file="bash_time.out", dec=".")
data.th = scan(file="bash_dataTH.out", dec=".")
data.h = scan(file="bash_dataH.out", dec=".")
time_print = printout[,3]
th_h_profiles = data.frame(Time=rep(time_print, each=nodes_max),
		      x=rep(grid_cords$x, each=length(time_print)),
		      #y=rep(grid_cords$y, each=length(time_print)),
		      Depth=rep(grid_cords$z, each=length(time_print)),
		      Moisture=data.th,Head=data.h)
##interpolate hydrus mesh output to regular grid
#library(fields)
library(gstat)
#generate regular-spaced grid
grid.x <- seq(min(grid_cords$x), max(grid_cords$x), by=1)
#grid.y <- seq(min(grid_cords$y), max(grid_cords$y), by=1)
grid.z <- seq(min(grid_cords$z), max(grid_cords$z), by=0.05)
#grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)
grid.xz <- expand.grid(x=grid.x, Depth=grid.z)

#interpolate and "stack" for each timestep
nodal.plot=data.frame()
for(i in time_print){
th_h_data = filter(th_h_profiles, Time==i) #filter for one timestep
#theta
idw.gstat <- gstat(formula = Moisture ~ 1, locations = ~ x + Depth, data = th_h_data, nmax = 15, set = list(idp = 2))
theta_interpolated <- predict(idw.gstat, grid.xz)
data_interpolated = cbind(Time = i, theta_interpolated[,-4])
colnames(data_interpolated)[4] = "Moisture"
theta_t = melt(data_interpolated, id.vars=c("Time","Depth","x"), measure.vars=c("Moisture"), variable.name = "Parameters")
#head
idw.gstat <- gstat(formula = Head ~ 1, locations = ~ x + Depth, data = th_h_data, nmax = 15, set = list(idp = 2))
head_interpolated <- predict(idw.gstat, grid.xz)
data_interpolated = cbind(Time = i, head_interpolated[,-4])
colnames(data_interpolated)[4] = "Head"
head_t = melt(data_interpolated, id.vars=c("Time","Depth","x"), measure.vars=c("Head"), variable.name = "Parameters")
#join data together
nodal.plot = rbind(nodal.plot, theta_t, head_t)
}

 balanceTIME = read.table(file="Balance.out", header=F, skip=21, sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=29, sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=30, sep="", nrows=1, dec=".")
 balance = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
for(i in 2:length(time_print)){
 if(i == length(time_print)) break
 balanceTIME = read.table(file="Balance.out", header=F, skip=((i-1)*13+21), sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=((i-1)*13+29), sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=((i-1)*13+30), sep="", nrows=1, dec=".")
 balanceDATA = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
 balance = rbind(balance, balanceDATA)
 }
 colnames(balance) = c("Time","WaterBalanceVolume","WaterBalancePercent")
 h_Mean <<- h_Mean
 v_Mean <<- v_Mean
 th_h_profiles <<- th_h_profiles 
 balance <<- balance

 #plotting
 if(plotting){
# library(ggplot2)
# library(reshape2)
# library(dplyrExtras)
# library(grid)
# library(gridExtra)
# library(scales)
 h_v_Mean.plot = melt(h_v_Mean, id.vars=c("datetime","Time"), measure.vars=colnames(h_v_Mean)[2:22], variable.name = "Parameters")
 #nodal.plot = melt(th_h_profiles, id.vars=c("Time", "Depth","x"), measure.vars=c("Moisture","Head"), variable.name = "Parameters")

 #printtimes = unique(th_h_profiles$Time)
 #printdates = c(datetime[1],datetime[printtimes])
 printdates = c(datetime[1],datetime[time_print])
 #  printdates = c(datetime[printtimes+1], datetime[timesteps])
 #dates cut to "days", better for displaying in x-axis
 #but problematic if 2 dates are selected within one day,then use second option
 #dates_times = data.frame(Time=time_print, Dates=as.POSIXct(format(printdates, "%Y-%m-%d:%H:%M")))
 dates_times = data.frame(Time=time_print, Dates=printdates)
 t_printout = time_print[seq(1, length(time_print), factorforprintout)]
 dates_times_print = filter(dates_times, Time%in%t_printout)

 xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(length(t_printout))

 #  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 #subsetting to choose standard parameters
 choose.params.h_v_Mean = c("vAtm", "vRoot", "vKode3", "hAtm", "hRoot", "hKode3")
 units = data.frame(type=c("Flux","Head"), Unit=c("[m²/h]","[m]"))
 type = data.frame(Parameters = c("vAtm","vKode3","hAtm","hKode3"), type=c("Flux","Flux","Head","Head"))
 cols = data.frame(Parameters = choose.params.h_v_Mean, Boundary = c("Top","Root","Bottom"))
 #  h_v_Mean.plot.filter = 	filter(h_v_Mean.plot, Parameters == choose.params.h_v_Mean) %>%
 # browser()
 # h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 #                         filter(Parameters != "sum(vAtm)" & Parameters != "sum(vRoot)" & Parameters != "sum(vKode3)") %>%
 #                         mutate(type="Flux") %>%
 #                         mutate_if(Parameters == "hAtm" | Parameters == "hRoot" | Parameters == "hKode3", type=="Head") %>%
 #                         inner_join(units,by="type") %>%
 #                        mutate(typeUnits = paste(type,Unit)) %>%
 #                        inner_join(cols, by="Parameters")
			
 h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 			inner_join(type,by="Parameters") %>%
 			inner_join(units,by="type") %>%
			mutate(typeUnits = paste(type,Unit)) %>%
			inner_join(cols, by="Parameters")

			# h_v_Mean.plot.filter$ParamUnits = factor(h_v_Mean.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 h_v_Mean.plot.filter$Boundary = factor(h_v_Mean.plot.filter$Boundary,levels=c("Top","Root","Bottom"))

 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value, colour=Boundary)) + geom_line() + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates), colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour=xcolors, size=2, alpha=0.5) + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + 
	 geom_vline(xintercept = as.numeric(dates_times_print$Dates),color=rep(xcolors,2), size=2, alpha=0.5) +
	 scale_x_datetime(breaks=dates_times_print$Dates, labels=date_format("%d-%m-%Y")) 
 
 choose.params.nodal = c("Head", "Moisture")
 units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]"))
 #  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 if(vertical==T){
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture",Parameters)) %>%
			filter(x==profile_loc) %>%
			filter(Time%in%t_printout) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
 			inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") %>%
			arrange(Time, Parameters, Depth, x)
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]"))
 nodal.plot.filter$Dates = as.factor(nodal.plot.filter$Dates)
 profile_type="vertical"
 nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 }else{
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture",Parameters)) %>%
			filter(Depth==profile_loc) %>%
			filter(Time%in%t_printout) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
 			inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") %>%
			arrange(Time, Parameters, Depth, x)
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]"))
 nodal.plot.filter$Dates = as.factor(nodal.plot.filter$Dates)
 profile_type="horizontal"
 nodal.gg = ggplot(nodal.plot.filter, aes(x=x, y=value, colour=Dates, group=Dates)) + geom_path() + xlab("Cross-section [m]") + ylab ("") + facet_grid(ParamUnits ~ ., scale="free_y") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 }
 #  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 #nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 #merge both plots together
 grid.arrange(h_v_Mean.gg,nodal.gg,top=textGrob(paste(folder,":", profile_type," profile at: ", profile_loc,"[m]", sep="")))
}
 return(print("Output stored in variables h_Mean, v_Mean, th_h_profiles and balance"))
}


#' @title Read Hydrus 3D output data
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param timesteps number of modeled timesteps
#' @param t_printout number of model printout times
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' @param px number of horizontal node, where the vertical profile should be analyzed
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 

read_hydrus3d <- function(folder, timesteps, t_printout, plotting=F, timestamps, px){
	# library(dplyr)
	# library(dplyrExtras)

 #path_hydrus = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus3d_out/" 
 path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus3d/hydrus3D_experiments/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read hydrus 2d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="",skip=5, nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,1]
 grid_cords = read.table("MESHTRIA.TXT", header=T, skip=6 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

 h_num_lines = length(readLines("h_Mean.out"))
 v_num_lines = length(readLines("v_Mean.out"))
 h_Mean = read.table(file="h_Mean.out", header=F, skip=6, sep="",nrows=(h_num_lines-7), dec=".")
 v_Mean = read.table(file="v_Mean.out", header=F, skip=13, sep="",nrows=(v_num_lines-14), dec=".")
 colnames(h_Mean) = c("Time","hAtm","hRoot","hKode3","hKode1","hSeep","hKode5","hKode6","hKode7","hKode8","hKode9")
 colnames(v_Mean) = c("Time","rAtm","rRoot","vAtm","vRoot","vKode3","vKode1","vSeep","vDrain","vBottom","vKode7","vKode8","vKode9","RunOff","Evapor","Infiltr","SnowLayer")
 #remove duplicated data lines (interpolation artefacts probably)
 h_Mean$Timecut=round(h_Mean$Time,0)
 v_Mean$Timecut=round(v_Mean$Time,0)
 h_Mean = h_Mean[!duplicated(h_Mean$Timecut),]
 v_Mean = v_Mean[!duplicated(v_Mean$Timecut),]
 #joint datasets h_Mean and v_Mean
 h_v_Mean = inner_join(h_Mean, v_Mean, by="Time")
 h_v_Mean = mutate(h_v_Mean, datetime=timestamps) #add real dates along with timesteps

 #generate sequence of vertical profiles to read out; each vertical profile needs their own sequence!!
 #vertiklales profil bei x=2 & x=10m
 #vert1 = seq(px,by=horizontal_nodes,length.out=vertical_nodes)
 
system("sed -n '/^[[:space:]]*[T]/p' TH.TXT > bash_time.out")
system("sed -n '/^[[:space:]]*[0-9]/p' TH.TXT > bash_dataTH.out")
system("sed -n '/^[[:space:]]*-*[0-9]/p' H.TXT > bash_dataH.out")

printout = read.table(file="bash_time.out", dec=".")
data.th = scan(file="bash_dataTH.out", dec=".")
data.h = scan(file="bash_dataH.out", dec=".")
time_print = printout[,3]
th_h_profiles = data.frame(Time=rep(time_print, each=nodes_max),
		      x=rep(grid_cords$x, each=length(time_print)),
		      y=rep(grid_cords$y, each=length(time_print)),
		      Depth=rep(grid_cords$z, each=length(time_print)),
		      Moisture=data.th,Head=data.h)

##interpolate hydrus mesh output to regular grid
#library(fields)
library(gstat)

#generate regular-spaced grid
grid.x <- seq(min(grid_cords$x), max(grid_cords$x), by=1)
grid.y <- seq(min(grid_cords$y), max(grid_cords$y), by=1)
grid.z <- seq(min(grid_cords$z), max(grid_cords$z), by=0.05)
grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)

#interpolate and "stack" for each timestep
nodal.plot=data.frame()
for(i in time_print){
th_h_data = filter(th_h_profiles, Time==i) #filter for one timestep
#theta
idw.gstat <- gstat(formula = Moisture ~ 1, locations = ~ x + y + Depth, data = th_h_data, nmax = 15, set = list(idp = 2))
theta_interpolated <- predict(idw.gstat, grid.xyz)
data_interpolated = cbind(Time = i, theta_interpolated[,-5])
colnames(data_interpolated)[5] = "Moisture"
theta_t = melt(data_interpolated, id.vars=c("Time","Depth","x","y"), measure.vars=c("Moisture"), variable.name = "Parameters")
#head
idw.gstat <- gstat(formula = Head ~ 1, locations = ~ x + y + Depth, data = th_h_data, nmax = 15, set = list(idp = 2))
head_interpolated <- predict(idw.gstat, grid.xyz)
data_interpolated = cbind(Time = i, head_interpolated[,-5])
colnames(data_interpolated)[5] = "Head"
head_t = melt(data_interpolated, id.vars=c("Time","Depth","x","y"), measure.vars=c("Head"), variable.name = "Parameters")
#join data together
nodal.plot = rbind(nodal.plot, theta_t, head_t)
}

#th_h_profiles = data.frame()
 #for(i in 0:t_printout){
 #TH_infoTIME = read.table(file="TH.TXT", header=F, skip=1 + i*(ceiling(nodes_max/10) + 3), sep="", nrows=1, dec=".")[,3]
 #TH_infoIN = scan(file="TH.TXT", skip=3+ i*(ceiling(nodes_max/10) + 3), sep="", nmax=nodes_max, dec=".")
 #H_infoTIME = read.table(file="H.TXT", header=F, skip=1+ i*(ceiling(nodes_max/10) + 3), sep="", nrows=1, dec=".")[,3]
 #H_infoIN = scan(file="H.TXT", skip=3+ i*(ceiling(nodes_max/10) + 3), sep="", nmax=nodes_max, dec=".")
 #profiles = data.frame(x=order_x ,Depth=(order_z-max(order_z)), Moisture=TH_infoIN, Head=H_infoIN,Time=TH_infoTIME)
 #th_h_profiles = rbind(th_h_profiles,profiles)
 #}

 balanceTIME = read.table(file="Balance.out", header=F, skip=21, sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=29, sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=30, sep="", nrows=1, dec=".")
 balance = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
for(i in 2:(t_printout+1)){
 if(i == (t_printout+1)) break
 balanceTIME = read.table(file="Balance.out", header=F, skip=((i-1)*13+21), sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=((i-1)*13+29), sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=((i-1)*13+30), sep="", nrows=1, dec=".")
 balanceDATA = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
 balance = rbind(balance, balanceDATA)
 }
 colnames(balance) = c("Time","WaterBalanceVolume","WaterBalancePercent")
 h_Mean <<- h_Mean
 v_Mean <<- v_Mean
 th_h_profiles <<- th_h_profiles 
 balance <<- balance

 #plotting
 if(plotting){
# library(ggplot2)
# library(reshape2)
# library(dplyrExtras)
# library(grid)
# library(gridExtra)
# library(scales)
 h_v_Mean.plot = melt(h_v_Mean, id.vars=c("datetime","Time"), measure.vars=colnames(h_v_Mean)[2:22], variable.name = "Parameters")
 #nodal.plot = melt(th_h_profiles, id.vars=c("Time", "Depth","x","y"), measure.vars=c("Moisture","Head"), variable.name = "Parameters")

 printtimes = unique(th_h_profiles$Time)
 #printdates = c(datetime[1],datetime[printtimes])
 printdates = c(timestamps[1],timestamps[printtimes])
 dates_times = data.frame(Time=printtimes, Dates=format(printdates, "%Y-%m-%d:%H"))

 xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(t_printout+1)
 #  xcolors = colorRampPalette(c("blue","green"))(t_printout+1)

 #  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 #subsetting to choose standard parameters
 choose.params.h_v_Mean = c("vAtm", "vRoot", "vKode3", "hAtm", "hRoot", "hKode3")
 units = data.frame(type=c("Flux","Head"), Unit=c("[m²/h]","[m]"))
 type = data.frame(Parameters = c("vAtm","vKode3","hAtm","hKode3"), type=c("Flux","Flux","Head","Head"))
 cols = data.frame(Parameters = choose.params.h_v_Mean, Boundary = c("Top","Root","Bottom"))
 #  h_v_Mean.plot.filter = 	filter(h_v_Mean.plot, Parameters == choose.params.h_v_Mean) %>%
 # browser()
 # h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 #                         filter(Parameters != "sum(vAtm)" & Parameters != "sum(vRoot)" & Parameters != "sum(vKode3)") %>%
 #                         mutate(type="Flux") %>%
 #                         mutate_if(Parameters == "hAtm" | Parameters == "hRoot" | Parameters == "hKode3", type=="Head") %>%
 #                         inner_join(units,by="type") %>%
 #                        mutate(typeUnits = paste(type,Unit)) %>%
 #                        inner_join(cols, by="Parameters")
			
 h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 			inner_join(type,by="Parameters") %>%
 			inner_join(units,by="type") %>%
			mutate(typeUnits = paste(type,Unit)) %>%
			inner_join(cols, by="Parameters")

			# h_v_Mean.plot.filter$ParamUnits = factor(h_v_Mean.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 h_v_Mean.plot.filter$Boundary = factor(h_v_Mean.plot.filter$Boundary,levels=c("Top","Root","Bottom"))

 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value, colour=Boundary)) + geom_line() + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates), colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour=xcolors, size=2, alpha=0.5) + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 
 choose.params.nodal = c("Head", "Moisture")
 units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]"))
 units$Parameters = factor(units$Parameters,levels=c("Moisture","Head"))

 #  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture",Parameters)) %>%
			filter(x==5) %>%
			filter(y==5) %>%
			#filter(x==px[1]) %>%
			#filter(y==px[1]) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
		        inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") %>%
			arrange(Time, Parameters, Depth, x)
			
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]"))

 #  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 #merge both plots together
 grid.arrange(h_v_Mean.gg,nodal.gg,main=paste(folder,": vertical profile at x =", px, sep=" "))
}
 return(print("Output stored in variables h_Mean, v_Mean, th_h_profiles and balance"))
}

#' @title Select time series of one modeled node
#'
#' @description test
#'
#' @param nodal_info_in input datasets. should be a data.frame
#' @param loc_hor horizontal coordinate (in model coordintaes)
#' @param loc_vert vertical coordinate (in model coordinates)
#' @param sensorname name of sensor (for plotting and structuring)
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 
obsNode_hydrus = function(nodal_info_in, loc_hor, loc_vert,sensorname){
nodal.filter = 	filter(nodal_info_in, grepl("Moisture",Parameters)) %>%
			filter(x==loc_hor) %>%
			filter(Depth==loc_vert) %>%
			mutate(sensor = sensorname) %>%
			mutate(type = "modeled") %>%
 			inner_join(dates_times,by="Time") #%>%
			#arrange(Dates, sensor, type, value)
		node_out = data.frame(datetime = nodal.filter$Dates, sensor = nodal.filter$sensor, value = nodal.filter$value, type = nodal.filter$type)
		return(node_out)
}



#' @title Read Hydrus 2D observation node data 
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param timesteps number of modeled timesteps
#' @param t_printout string of dates of time series where informations are printed
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' @param px number of horizontal node, where the vertical profile should be analyzed
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 

# read_hydrus2d_irregular <- function(folder, timesteps, vertical_nodes, horizontal_nodes, t_printout, plotting=F, timestamps, px){
read_obsNode2d <- function(folder,realTime=T,startdate,plotting=F){
 library(stringr)
 #path_hydrus = "/home/mreich/server/hydro72/hydrus2d/experiments/" 
 path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus2d/experiments/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read mesh
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="", nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,2]
 nodes_cords = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

 #read obsveration node output file
 num_lines = length(readLines("obsNod.out"))
 obsNodeData = read.table(file="obsNod.out", header=T, skip=5, sep="",nrows=(num_lines-7), dec=".")
 nodeNames.in = read.table(file="obsNod.out", header=F, skip=3, sep="",nrows=1)
 nodeNames = vector()
 for(i in 1:(length(nodeNames.in)/2)){
 nodeNames[i] = as.numeric(str_extract(nodeNames.in[,i*2], "[0-9]+"))
 }
 #select columns
 theta_data = select(obsNodeData, contains("theta"))
 #generate zoo-TS
 if(realTime){
 obsNodeTheta = zoo(theta_data, order.by=as.POSIXct(obsNodeData$time*3600, format="%H", origin=startdate))
 }
 else{obsNodeTheta = zoo(theta_data, order.by=obsNodeData$time)}
 colnames(obsNodeTheta) = nodeNames
 #offer possibility to get output as pivot-table with node coordinates!
 obsNode.melt = melt(zootodf(obsNodeTheta), id="time", variable.name="node", value.name="theta")
 #extract node cordinates from cords
 ## dont know if this line works, but its something like this..!
 #obsNode_cords = match(nodeNames, nodes_cords)
 #obsNode.melt = inner_join(obsNode_cords, by=)
 #...
 if(plotting){
	ggplot(obsNode.melt, aes(x=time, y=theta, colour=node)) + geom_line() + facet_grid(node~.)
 }
 return(obsNodeTheta)
} # end function read obsNode data






#' @title Plot modeled nodes and observed soil moisture sensors
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param timesteps number of modeled timesteps
#' @param t_printout string of dates of time series where informations are printed
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' @param px number of horizontal node, where the vertical profile should be analyzed
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

obsNodePlot = function(filename, normdata=T, dim2d=T){
#load SM; filtered at 6 hourly intervalls
load(file="/home/mreich/server/hygra/DataWettzell/SoilMoisture/Cluster_Data/Data_filtered/SGnew_filtered_6hourmean.rdata")
beneathBuilding = SGnew.filter[,11:18] #get sensor beneath SG building
besidesBuilding = SGnew.filter[,19:25] #get closest sensors outside SG building
colnames(beneathBuilding) = c("c shallow","b shallow","d shallow","c middle high","c middle low","b deep","c deep","d deep")
colnames(besidesBuilding) = c("a02","a03","a04","a06","a10","a14","a18")

SMsensors = merge.zoo(beneathBuilding,besidesBuilding, all=T, fill=NA)
if(normdata){#normalize theta data
SMsensors.norm = zoo(apply(coredata(SMsensors),2,normalize), order.by=index(SMsensors))
SMsensors.melt = cbind(melt(zootodf(SMsensors.norm), id="time",variable.name="sensor"),type="observed")
}else{
SMsensors.melt = cbind(melt(zootodf(SMsensors), id="time",variable.name="sensor"),type="observed")
}

#load precipitation TS
load("/home/mreich/server/hygra/DataWettzell/Climate/30min_wettzell/clima30min_raw.rdata")
precip = clima.raw$Prec_Sen1
precip.df = zootodf(precip); colnames(precip.df)[2] = "Precip"
precip.melt = cbind(melt(precip.df, id="time", variable.name="sensor"), type="observed")
#adjust time series length to SM observations
precip.melt = filter(precip.melt, time > as.POSIXct("2010-03-10 16:00:00"))

#read hydrus observation nodes
if(dim2d) folpath = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus2_out/"
else folpath = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus3d_out/"
nodes=load(file=paste(folpath, filename, sep=""))
if(normdata){#normalize theta data
nodes.norm = zoo(apply(coredata(get(nodes)),2,normalize), order.by=index(get(nodes)))
nodes.melt = cbind(melt(zootodf(nodes.norm), id="time", variable.name="node"), type="modeled")
}else{
nodes.melt = cbind(melt(zootodf(get(nodes)), id="time", variable.name="node"), type="modeled")
}
#IMPORTANT!!
#sensor names have to be in the ORDER of node number names from hydrus
#this is somewhat caotic due to internal hydrus numbering
#this is for 2D only!!
if(dim2d) sensors_nodes = data.frame(sensor = c("b shallow","b deep","c shallow","c middle high","c middle low","d shallow","d deep","a02","a03","a04","a06","a10","a18","a14","c deep"), node = colnames(get(nodes)))
#3d
else sensors_nodes = data.frame(sensor = c("a04","a02","a03","a06","a10","a14","b shallow","b deep","c shallow","c middle high","c middle low","d shallow","a deep","c deep","d deep"), node = colnames(get(nodes)))

nodes.melt = inner_join(nodes.melt, sensors_nodes, by="node") %>%
		select(-node) %>%
		select(time, sensor, value, type)
#without precipitation
#SMdata = rbind(SMsensors.melt, nodes.melt)
#with precipitation
SMdata = rbind(SMsensors.melt, nodes.melt, precip.melt)
SMdata.gg = ggplot(SMdata, aes(x=time, y=value, colour=type)) + geom_line() + ylab("Theta [%VWC]") + xlab("") +
		facet_grid(sensor~., scale="free_y")
#save plot
#png(file="/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/compare/NORMALIZED_SMobs-mod_ts49025_2d.png", width=1800, height=1200, res=100)
return(SMdata.gg)
#dev.off()
}
