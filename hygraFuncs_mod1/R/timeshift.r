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
  #merge the datasets
  data.merge=merge.zoo(data1,data2, all=F) 
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
	grid.arrange(TS.plot, ccf.plot,main=textGrob(paste(sensor1,"&", sensor2, ":", max(year), sep=" "), gp=gpar(cex=1.5)))
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
