#' @title Time lag estimation via cross-correlation of two timeseries
#'
#' @description This function is based on standard time series decomposition and afterwards cross-correlation estimation between the two input datasets.
#' @description The main aim is to find dominant time shifts between two time series (e.g. soil moisture, climate data, precipitation stations, etc).
#'
#' @param data1,data2 Timeseries of type .zoo.
#' @param nmax Number of correlation values output.
#' @param norm logical. Use normalized data (default) or not.
#' @param stl logical. Use decomposed data for cross-correlation (default) or not.
#' @param plotting logical. Output cross-correlation plots. Default is F, own output plots are provided. Should be kept FALSE.
#' @param swin Time window used for seasonality estimation. Only used if stl=T.
#' @param twin Time window used for trend estimation. Only used if stl=T.
#' @param ... Further parameters passed to internal functions.

#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples outside.cor = ccf.zoo(besidesBuilding$mux43_04,besidesBuilding$mux43_08,5,T,T,T,24,17521)
#' @export
#' 

timeshift = function(data1,data2,nmax,norm=T,stl=T,plotting=T,swin=24,twin,sensor1="data1", sensor2="data2",...){
# library(zoo);
  Sys.setenv(TZ = "GMT")
# library(xts)
# library(gridExtra)
# if(plotting==TRUE){
# library(reshape2)
# library(ggplot2)}
  #get both datasets to use the same timestamp-frequency
  data1_sameF = aggregate(data1, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  data2_sameF = aggregate(data2, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  #merging all data together can result in problems when na.approx doesn't find 2 non-NA values for interpolation
  #therefor merge has to be used with option all=F
  #data.merge = merge(data1_sameF, data2_sameF, all=T, fill= NA)
  data.merge = merge(data1_sameF, data2_sameF, all=F)
  #merge the datasets, depending on the analysis time
  #this is necessary to avoid na.approx-problems with too many NA's
  #and on the other hand don't loose any data if stl=F
  #if(stl==T) data.merge=merge.zoo(na.approx(data1),na.approx(data2), all=T, fill=NA) 
  #else data.merge=merge.zoo(data1,data2, all=F) 
  #aggregate to hourly time
  #   data1.agg = period.apply(data1, endpoints(data1,"hours"), mean, na.rm=T)
  #   data2.agg = period.apply(data2, endpoints(data2,"hours"), mean, na.rm=T)
  #split in yearly intervalls
  ts.years = endpoints(data.merge, "years") 
  data1.year=list()
  data2.year=list()
  for(i in 1:(length(ts.years)-1)){
	data1.year[[i]]=data.merge[ts.years[i]:ts.years[i+1],1]
	data2.year[[i]]=data.merge[ts.years[i]:ts.years[i+1],2]
  }
  if(norm==T){ #use normalize data as input
  data1.year.norm = lapply(data1.year, function(x) (x - mean(x, na.rm=T))/sd(x,na.rm=T))
  data2.year.norm = lapply(data2.year, function(x) (x - mean(x, na.rm=T))/sd(x,na.rm=T))
  }
  else{ #use original data as input
  data1.year.norm = data1.year
  data2.year.norm = data2.year
  }
  #cross-correlation for each year
  correlation=data.frame(matrix(ncol=(length(ts.years)-1)*2,nrow=nmax)) 
  res.plots = list(); data.plots= list(); data1.stl = list(); data2.stl = list()
  j=1;k=2; cornames=NA
  for(i in 1:(length(correlation)/2)){
  	p1.flag=1;p2.flag=1 #setting standard plotting parameter
    	timedif.data = as.numeric(difftime(index(data1.year.norm[[i]][2]),index(data1.year.norm[[i]][1]), units="hours"))
	freq = 24/timedif.data #frequency in hours
	year=unique(format(index(data1.year.norm[[i]]), "%Y")) #get year from current TS
	if(stl==TRUE){ #decompose TS
	#print(i) #only for debugging
	check1=try(stl(ts(data1.year.norm[[i]],frequency=freq),s.window=swin,t.window=twin),TRUE)
	if(class(check1)=="try-error"){
		data1.na= na.approx(data1.year.norm[[i]])
		#if(i==2){browser()}
		check1.2=try(stl(ts(data1.na,frequency=freq),s.window=swin,t.window=twin),TRUE)
		if(class(check1.2)=="try-error"){
		correlation[,j] = NA; cornames[j]= paste("cor_",max(year),sep="") ;j=j+2
		correlation[,k] = NA; cornames[k]=paste("lagHOUR_",max(year),sep="") ;k=k+2
		print(paste(deparse(substitute(data1)),"Year",max(year), ": calculation not possible. Too much NA data.", sep=" "))
	       	next #skip year if stl is not possible because of missing data
		}
		data1.stl[[i]]= stl(ts(data1.na,frequency=freq),s.window=swin,t.window=twin)
		data1.in= data1.stl[[i]]$time.series[,3] #residuals of stl()
		print(paste(deparse(substitute(data1)),"Year",max(year), ": some NA data was approximized.", sep=" "))
		p1.flag=NA #set flag for choosing what data to plot
	}
	else{
	data1.stl[[i]]= stl(ts(data1.year.norm[[i]],frequency=freq),s.window=swin,t.window=twin)
	data1.in= data1.stl[[i]]$time.series[,3] #residuals of stl()
	} #end check-run data1
	check2=try(stl(ts(data2.year.norm[[i]],frequency=freq),s.window=swin,t.window=twin),TRUE)
	if(class(check2)=="try-error"){ #skip year if stl is not possible because of missing data
		data2.na= na.approx(data2.year.norm[[i]])
		check2.2=try(stl(ts(data2.na,frequency=freq),s.window=swin,t.window=twin),TRUE)
		if(class(check2.2)=="try-error"){
		correlation[,j] = NA; cornames[j]= paste("cor_",max(year),sep="") ;j=j+2
		correlation[,k] = NA; cornames[k]=paste("lagHOUR_",max(year),sep="") ;k=k+2
		print(paste(deparse(substitute(data2)),"Year",max(year), ": calculation not possible. Too much NA data.", sep=" "))
	       	next #skip year if stl is not possible because of missing data
		}
		data2.stl[[i]]= stl(ts(data2.na,frequency=freq),s.window=swin,t.window=twin)
		data2.in= data2.stl[[i]]$time.series[,3] #residuals of stl()
		print(paste(deparse(substitute(data2)),"Year",max(year), ": some NA data was approximized.", sep=" "))
		p2.flag=NA #set flag for choosing what data to plot
	}
	else{
	data2.stl[[i]]= stl(ts(data2.year.norm[[i]],frequency=freq),s.window=swin,t.window=twin)
	data2.in= data2.stl[[i]]$time.series[,3] #residuals of stl()
	} #end check-run data2
	}
	else{ #use NON-decomposed TS
	data1.in = coredata(data1.year.norm[[i]])
	data2.in = coredata(data2.year.norm[[i]])
	}
	# cross-correlation
	#browser()
	res = ccf(data1.in,data2.in,na.action=na.pass,plot=F,main=max(year),...)
	# browser()
    	acf.maxs=NA; cornmax=NA
    	for(ii in 1:nmax){ #finding nth max values
		 if(is.na(sort(res$acf, TRUE)[ii])==T){ cornmax[ii]=NA; acf.maxs[ii]=NA; next}
		 cornmax[ii] = sort(res$acf, TRUE)[ii]
     		 acf.maxs[ii] = which(res$acf == cornmax[ii])
   	}
    	lags.dom = res$lag[acf.maxs]
	#convert output from frequency to time units[hours]
	if(stl==TRUE){lags.real = lags.dom*freq*timedif.data}
	#if(stl==TRUE){lags.real = lags.dom/(timedif.data)} #old
	#else{lags.real = lags.dom*frequency(data1.year.norm[[i]])/(timedif.data)} #old
	else{lags.real = lags.dom*timedif.data}
	correlation[,j] = cornmax; cornames[j]= paste("cor_",max(year),sep="") ;j=j+2
	correlation[,k] = lags.real; cornames[k]=paste("lagHOUR_",max(year),sep="") ;k=k+2
	#plotting
	if(plotting==TRUE){
		if(is.na(p1.flag)==T){ #plot approximized data1
			dataplot1.ts = data.frame(Time=index(data1.na),Timeseries = coredata(data1.na), Residuals = data1.in) #use already processed (na.aproxx) data to show in as timeseries
			# dataplot1.ts = data.frame(Time=index(data1.year[[i]]),Timeseries = coredata(data1.year[[i]]), Residuals = data1.in)
		}
		else{ #plot original normalized data1
			dataplot1.ts = data.frame(Time=index(data1.year[[i]]),Timeseries = coredata(data1.year[[i]]), Residuals = data1.in)
		}
		if(is.na(p2.flag)==T){ #plot approximized data2
			dataplot2.ts = data.frame(Time=index(data2.na),Timeseries = coredata(data2.na), Residuals = data2.in) #use already processed (na.aproxx) data to show in as timeseries
			# dataplot2.ts = data.frame(Time=index(data2.year[[i]]),Timeseries = coredata(data2.year[[i]]), Residuals = data2.in)
		}
		else{ #plot original normalized data2
			dataplot2.ts = data.frame(Time=index(data2.year[[i]]),Timeseries = coredata(data2.year[[i]]), Residuals = data2.in)
		}
		# browser()
		# dataplot1.ts = merge(data1.year[[i]],data1.in, all=T, fill=NA); colnames(dataplot1.ts)=c("Time","Timeseries","Residuals")
		# dataplot2.ts = merge(data2.year[[i]],data2.in, all=T, fill=NA); colnames(dataplot2.ts)=c("Time","Timeseries","Residuals")
		# dataplot1.ts = data.frame(Time=index(data1.year[[i]]),Timeseries = coredata(data1.year[[i]]), Residuals = data1.in)
		# dataplot2.ts = data.frame(Time=index(data2.year[[i]]),Timeseries = coredata(data2.year[[i]]), Residuals = data2.in)
	data1.ts.resh = melt(dataplot1.ts, id="Time")
	data2.ts.resh = melt(dataplot2.ts, id="Time")
	# datasets = rbind(cbind(data1.ts.resh, Sensor=factor(deparse(substitute(data1)))),cbind(data2.ts.resh, Sensor=factor(deparse(substitute(data2)))))
	datasets = rbind(cbind(data1.ts.resh, Sensor=factor(sensor1)),cbind(data2.ts.resh, Sensor=factor(sensor2)))
	TS.plot= ggplot(datasets, aes(x=Time, y=value, colour=Sensor)) + geom_line() + facet_grid(variable ~ ., scale="free_y") + xlab("") + ylab("Signal") # + labs(title=paste(deparse(substitute(data1)),"&",deparse(substitute(data2)), ":", max(year), sep=" "))
	dataccf = data.frame(Lag=res$lag[,,1], Correlation=res$acf[,,1])
	ccf.plot= ggplot(dataccf, aes(x=Lag, y=Correlation)) + geom_bar(stat = "identity")
	#combining plots
	# grid.arrange(TS.plot, ccf.plot,main=textGrob(paste(deparse(substitute(data1)),"&",deparse(substitute(data2)), ":", max(year), sep=" "), gp=gpar(cex=1.5)))
	#grid.arrange(TS.plot, ccf.plot,main=textGrob(paste(sensor1,"&", sensor2, ":", max(year), sep=" "), gp=gpar(cex=1.5)))
	grid.arrange(TS.plot, ccf.plot,top=textGrob(paste(sensor1,"&", sensor2, ":", max(year), sep=" "), gp=gpar(cex=1.5)))
	data.plots[[i]] = recordPlot() #store plot in variable
	}
   }
    colnames(correlation) = cornames #set col names as specific years 
    data.plots <<- data.plots
    data1.stl <<- data1.stl
    data2.stl <<- data2.stl
    ccf.results <<-res 
    return(correlation)
}
### end function timeshift


#' @title Event-based calulation of correlation coefficients, lagtimes and signal to noise ratios
#'
#' @description 
#' @description So far hard-coded for SM data of wettzell SGnew TS
#'
#' @param event vector of start and end date of event.
#' @param event_name character, name of event (mainly used for the plot).
#' @param mv_win length of window for moving average to smooth data in order to calculate signal to noise ratios.
#' @param plots logical, should the results be plotted? (default is TRUE).
#'
#' @details Output is a data.frame with correlation coeficients, lagtime and signal to soise ratios for each sensor.
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples
#' @export
#' 

ccf_events = function(event, event_name, mv_win, plots=T, ...){
#load data
#soil moisture
load(file="/home/mreich/server/hygra/DataWettzell/SoilMoisture/Cluster_Data/Data_filtered/SGnew_filtered_6hourmean.rdata")
beneathBuilding = SGnew.filter[,11:18] #get sensor beneath SG building
besidesBuilding = SGnew.filter[,19:25] #get closest sensors outside SG building
#GW data
load(file="/home/mreich/server/hygra/DataWettzell/Groundwater/GW_Pegel_complete/BK14_corrected.rdata") #GW 14
BK14.agg = aggregate(BK14, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
load(file="/home/mreich/server/hygra/DataWettzell/Groundwater/GW_Pegel_complete/BK3_corrected.rdata") #GW 14
BK3.agg = aggregate(BK3, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
gw = merge(BK14.agg,BK3.agg=(BK3.agg-2.01),all=T,fill=NA)
BK14_na = which(is.na(gw$BK14)==T)
gw$BK14[BK14_na] = gw$BK3[BK14_na]
BK14 = gw$BK14 * -1
#precipitation
load("/home/mreich/server/hygra/DataWettzell/Climate/30min_wettzell/clima30min_raw.rdata")
precip = clima.raw$Prec_Sen1
#trim ts to event period
a06_event = besidesBuilding$mux43_05[which(index(besidesBuilding$mux43_05)==event[1]):which(index(besidesBuilding$mux43_05)==event[2])]
a10_event = besidesBuilding$mux43_06[which(index(besidesBuilding$mux43_06)==event[1]):which(index(besidesBuilding$mux43_06)==event[2])]
a14_event = besidesBuilding$mux43_07[which(index(besidesBuilding$mux43_07)==event[1]):which(index(besidesBuilding$mux43_07)==event[2])]
a18_event = besidesBuilding$mux43_08[which(index(besidesBuilding$mux43_08)==event[1]):which(index(besidesBuilding$mux43_08)==event[2])]
b06_event = beneathBuilding$mux44_02[which(index(beneathBuilding$mux44_02)==event[1]):which(index(beneathBuilding$mux44_02)==event[2])]
b16_event = beneathBuilding$mux44_06[which(index(beneathBuilding$mux44_06)==event[1]):which(index(beneathBuilding$mux44_06)==event[2])]
c06_event = beneathBuilding$mux44_01[which(index(beneathBuilding$mux44_01)==event[1]):which(index(beneathBuilding$mux44_01)==event[2])]
c08_event = beneathBuilding$mux44_04[which(index(beneathBuilding$mux44_04)==event[1]):which(index(beneathBuilding$mux44_04)==event[2])]
c14_event = beneathBuilding$mux44_05[which(index(beneathBuilding$mux44_05)==event[1]):which(index(beneathBuilding$mux44_05)==event[2])]
c17_event = beneathBuilding$mux44_07[which(index(beneathBuilding$mux44_07)==event[1]):which(index(beneathBuilding$mux44_07)==event[2])]
d07_event = beneathBuilding$mux44_03[which(index(beneathBuilding$mux44_03)==event[1]):which(index(beneathBuilding$mux44_03)==event[2])]
d18_event = beneathBuilding$mux44_08[which(index(beneathBuilding$mux44_08)==event[1]):which(index(beneathBuilding$mux44_08)==event[2])]
gw_event = BK14[which(index(BK14)==event[1]):which(index(BK14)==event[2])]
precip_event = precip[which(index(precip)==event[1]):which(index(precip)==event[2])]
#merging all dataset for a visual check
ts_event = merge(a06_event,a10_event,a14_event,a18_event,b06_event,b16_event,c06_event,c08_event,c14_event,c17_event,d07_event,d18_event,gw_event,precip_event,all=T,fill=NA)
ts_event$gw_event = na.approx(ts_event$gw_event)
ts_event$precip_event = na.approx(ts_event$precip_event)
colnames(ts_event) = c("a shallow","a middle high","a middle low","a deep","b shallow","b deep","c shallow","c middle high","c middle low","c deep","d shallow","d deep","gw","precip")
if(plots==T) plot(ts_event, xlab="", main=event_name)
#correlation analysis and snr-calculation
stl_flag=F; plotting=F
swin=17521 ;twin=17521
#shallow
a06_b06_event.cor = timeshift(a06_event,b06_event,1,T,stl_flag,plotting,swin,twin,"",""); b06_event.snr = round(snr(b06_event,mv_win,T,plotting),2)
a06_c06_event.cor = timeshift(a06_event,c06_event,1,T,stl_flag,plotting,swin,twin,"",""); c06_event.snr = round(snr(c06_event,mv_win,T,plotting),2)
a06_d07_event.cor = timeshift(a06_event,d07_event,1,T,stl_flag,plotting,swin,twin,"",""); d07_event.snr = round(snr(d07_event,mv_win,T,plotting),2)
#in between sensors (only 2 from profile c)
a10_c08_event.cor = timeshift(a10_event,c08_event,1,T,stl_flag,plotting,swin,twin,"",""); c08_event.snr = round(snr(c08_event,mv_win,T,plotting),2)
a14_c14_event.cor = timeshift(a14_event,c14_event,1,T,stl_flag,plotting,swin,twin,"",""); c14_event.snr = round(snr(c14_event,mv_win,T,plotting),2)
#deep
a18_b16_event.cor = timeshift(a18_event,b16_event,1,T,stl_flag,plotting,swin,twin,"",""); b16_event.snr = round(snr(b16_event,mv_win,T,plotting),2)
a18_c17_event.cor = timeshift(a18_event,c17_event,1,T,stl_flag,plotting,swin,twin,"",""); c17_event.snr = round(snr(c17_event,mv_win,T,plotting),2)
a18_d18_event.cor = timeshift(a18_event,d18_event,1,T,stl_flag,plotting,swin,twin,"",""); d18_event.snr = round(snr(d18_event,mv_win,T,plotting),2)
#correlation results
results_event = rbind(cbind(Profile="shallow",a06_b06_event.cor,sensor="b", depth=0.6, event=event_name, snr=b06_event.snr),
		cbind(Profile="shallow",a06_c06_event.cor,sensor="c", depth=0.6, event=event_name, snr=c06_event.snr),
		cbind(Profile="shallow",a06_d07_event.cor,sensor="d", depth=0.7, event=event_name, snr=d07_event.snr),
		cbind(Profile="middle high",a10_c08_event.cor,sensor="c", depth=0.9, event=event_name, snr=c06_event.snr),
		cbind(Profile="middle low",a14_c14_event.cor,sensor="c", depth=1.4, event=event_name, snr=c14_event.snr),
		cbind(Profile="deep",a18_b16_event.cor,sensor="b", depth=1.6, event=event_name, snr=b16_event.snr),
		cbind(Profile="deep",a18_c17_event.cor,sensor="c", depth=1.7, event=event_name, snr=c17_event.snr),
		cbind(Profile="deep",a18_d18_event.cor,sensor="d", depth=1.8, event=event_name, snr=d18_event.snr))
colnames(results_event)[2:3] = c("correlation","lagtime")

return(results_event)
} #end of function

#' @title Multi-event based ccf approach, including lagtimes and SNR
#'
#' @description 
#' @description So far hard-coded for SM data of wettzell SGnew TS
#'
#' @param events data.frame of events; each row limits one event, event[,1] = start dates, event[,2] = end dates, event[,3] = event name.
#' @param mv_win length of window for moving average to smooth data in order to calculate signal to noise ratios.
#' @param limit_snr threshold for signal to noise ratio values; default value is 1 and should not be changed. see details. 
#' @param limit_cor threshold for correlation coefficient values; all values below limit_cor will be discarted.
#' @param limit_lag threshold for lagtime values; all values above limit_lag will be discarted.
#' @param plotting logical, should the results be plotted? (default is TRUE).
#'
#' @details Output is a data.frame with mean values for correlation coeficients, lagtime and quality of signal to soise ratio.
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples
#' @export
#' 

ccf_multievent = function(events, mv_win, limit_snr=1, limit_cor, limit_lag, plotting=T){

	numberofevents = length(events[,1])
	results_all = data.frame()
	#creating results for each event
	for(i in 1:numberofevents){
		event = c(events[i,1:2])
		name=as.character(events[i,3])
		res = ccf_events(event, name, mv_win, F)
		#filter results for bad snr AND count occurances of bad snr
		res$snrtest = 1; res$cortest = 1; res$lagtest = 1
		#filter bad SNR
		snr_out = which(res$snr < limit_snr)
		#filter bad correlation
		cor_out = which(res$correlation < limit_cor)
		#filter bad lagtime
		lag_out = which(res$lagtime > limit_lag)
		#apply cor and lag filter
		res$correlation[snr_out] = NA
		res$lagtime[snr_out] = NA
		res$snrtest[snr_out] = 0
		res$correlation[cor_out] = NA
		res$lagtime[cor_out] = NA
		res$cortest[cor_out] = 0 
		res$correlation[lag_out] = NA
		res$lagtime[lag_out] = NA
		res$lagtest[lag_out] = 0
		results_all = rbind(results_all, res)
		#result_all[,,i] = as.matrix(cbind(res$correlation,res$lagtime,res$snr))
		sensor_metadata = data.frame(Profile=res$Profile,sensor=res$sensor,depth=res$depth)
		sensor_metadata = mutate(sensor_metadata, sensorprofile = paste(sensor,Profile,sep=" "))
	}
	#statisical merging
	#adding unique site-identification (profile + sensor)
	results_ordered = mutate(results_all, sensorprofile = paste(sensor,Profile,sep=" ")) %>%
	 	  group_by(sensorprofile)
	#calculate means and check stats of bad data quality (=snr > limit_snr)
	cor_mean = summarise(results_ordered, cor_mean=mean(correlation, na.rm=T),cor_min=min(correlation, na.rm=T),cor_max=max(correlation, na.rm=T))
	lag_mean = summarise(results_ordered, lag_mean=mean(lagtime, na.rm=T),lag_min=min(lagtime, na.rm=T),lag_max=max(lagtime, na.rm=T))
	snrpass = summarise(results_ordered, SNRpass = sum(snrtest, na.rm=T))
	corpass = summarise(results_ordered, CORpass = sum(cortest, na.rm=T))
	lagpass = summarise(results_ordered, LAGpass = sum(lagtest, na.rm=T))
	#putting together with sensor meta data
	results_mean = cbind(cor_mean[,1:2], lag_mean[,2], snrpass[,2], corpass[,2], lagpass[,2])
	valuerange = cbind(cor_mean[,-2], lag_mean[,3:4])
	results_mean = inner_join(results_mean, sensor_metadata, by="sensorprofile") %>%
		       mutate(SNR = SNRpass / numberofevents) %>%
		       mutate(COR = CORpass / numberofevents) %>%
		       mutate(LAG = LAGpass / numberofevents)
	results_table = data.frame(Sensor = results_mean$sensorprofile, Depth = results_mean$depth, Correlation = results_mean$cor_mean, Lagtime = results_mean$lag_mean, SNRpass = results_mean$SNR, CORpass = results_mean$COR,LAGpass = results_mean$LAG)
	if(plotting ==T){
	#prepare for plotting
	results_mean_plot = melt(results_mean, id=c("Profile","sensor","depth","sensorprofile"))
	#actual plotting
	#correlation
	results_cor = mutate(results_mean_plot, sensor = factor(sensor, levels=c("d","c","b"))) %>%
		      filter(variable == "cor_mean") %>%
		      inner_join(valuerange, by="sensorprofile") %>%
		      group_by(sensor, depth)
	correlation.gg = ggplot(data=results_cor, aes(x=value, y=(depth*-1),colour=value)) + ylab("Depth [m]") + xlab("Correlation coefficient [%]") + geom_point(size=4.5) + theme_bw() + ggtitle("") +
	geom_errorbarh(aes(xmax=cor_max,xmin=cor_min)) + geom_path() +
	scale_y_continuous(breaks = c(-0.6,-1.0,-1.5,-1.8)) +
	facet_grid(.~sensor) + 
	scale_colour_gradientn("",colours=c("black","blue","red"))#, breaks=breaks_val)
	#lagtime
	results_lag = mutate(results_mean_plot, sensor = factor(sensor, levels=c("d","c","b"))) %>%
		      filter(variable == "lag_mean") %>%
		      inner_join(valuerange, by="sensorprofile")
	lagtime.gg = ggplot(data=results_lag, aes(x=value, y=(depth*-1),colour=value)) + ylab("Depth [m]") + xlab("Lagtime [h]") + geom_point(size=4.5) + theme_bw() + ggtitle("") +
	geom_errorbarh(aes(xmax=lag_max,xmin=lag_min)) + geom_path() +
	scale_y_continuous(breaks = c(-0.6,-1.0,-1.5,-1.8)) +
	facet_grid(.~sensor) + 
	scale_colour_gradientn("",colours=c("black","blue","red"))#, breaks=breaks_val)
	#snrtest failure
	results_snr = mutate(results_mean_plot, sensor = factor(sensor, levels=c("d","c","b"))) %>%
		      filter(variable == "SNR" | variable == "COR" | variable == "LAG")
	snrpass.gg = ggplot(data=results_snr, aes(x=value, y=(depth*-1),colour=variable)) + ylab("Depth [m]") + xlab("events passed quality measurement[%]") + geom_point(size=4.5, position = position_jitter(w=0,h=0.05)) + theme_bw() + ggtitle("") +
	scale_y_continuous(breaks = c(-0.6,-1.0,-1.5,-1.8)) +
	facet_grid(.~sensor)# + 
	#scale_colour_gradientn("",colours=c("black","blue","red"))#, breaks=breaks_val)
	#all plots together
	grid.arrange(correlation.gg,lagtime.gg,snrpass.gg)
	}
	#
	#return table of mean result
	return(results_table)
}



