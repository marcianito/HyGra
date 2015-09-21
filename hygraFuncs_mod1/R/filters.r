#
#' @title Filtering raw soil moisture data
#'
#' @description Filtering raw soil moisture data with adjustable filter criterias.
#'
#' @param file.input input dataset (zoo)
#' @param val.max/min cuts all data greater/smaller than this value
#' @param val.incr/decr maximum increase/decrease in absolute m³/m³ within timesteps lag.max
#' @param lag.max number of timesteps to consider for maximal increase / decrease
#' @param per.max/min value in percent for maximum / minimum deviation between actual value and a mean value of the time series
#' @param lag.mean period for moving mean calculation
#' @param k splits datasets at every k-th element of lag.mean
#' @details The output will be as zoo-object.
#' @details Default values for val.max=1 and val.min=0.
#' @details Possible values for lag.mean are "hours", "days", "weeks", "months", "quarters" or "years".
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export
#' 
filter_soildata <- function(file.input,val.max=1,val.min=0,val.incr,val.decr,lag.max, per.max, per.min, lag.mean,k){
# source("/home/mreich/R-gravity-package/package09/hygra_test/R/convertzootodataframe.r")
# require(xts)
file.filter = zootodf(file.input)
#filtering
high=which(file.filter[,-1] > val.max, arr.ind=T); high[,2] = high[,2]+1 #unrealistic values upper boundary
low=which(file.filter[,-1] < val.min, arr.ind=T);low[,2] = low[,2]+1 #unrealistic values lower boundary
file.filter[high] = NA
file.filter[low] = NA
num.high=length(high[,1])
num.low=length(low[,1])

# lag.max=10
# num.increase=0; num.decrease=0
#for(comlag in 1:lag.max){
#decrease=which(sapply(file.filter[,-1],diff,comlag) < -val.decr, arr.ind=T); decrease[,1] = decrease[,1]+comlag; decrease[,2] = decrease[,2]+1
#increase=which(sapply(file.filter[,-1],diff,comlag) > val.incr, arr.ind=T); increase[,1] = increase[,1]+comlag; increase[,2] = increase[,2]+1
#file.filter[decrease] = NA
#file.filter[increase] = NA
#num.increase=num.increase + length(increase[,1])
#num.decrease=num.decrease + length(decrease[,1])
#}  

breaks=endpoints(file.filter$time,lag.mean,k) #create breaks in the chosen intervall for splitting up data-analysis (mean calculation)
num.rollingmean = 0 #initialize counting number of corrections
for(i in 2:dim(file.filter)[2]){ #mux (column)-wise
  for(j in 1:(length(breaks)-1)){ #week-wise
    mean.actual = mean(file.filter[breaks[j]:breaks[j+1],i], na.rm=T)
    filter.out = which(file.filter[breaks[j]:breaks[j+1],i]/mean.actual > per.max | file.filter[breaks[j]:breaks[j+1],i]/mean.actual < per.min, arr.ind=T)
    file.filter[(filter.out+breaks[j]-1),i] = NA
    num.rollingmean = num.rollingmean + length(filter.out)
  }
}

# stats = data.frame(number=c(num.high,num.low,num.increase,num.decrease,num.rollingmean)); rownames(stats) = c("greater max", "lower min", "increasing too fast", "decreasing too fast", "deviation from moving mean")
stats = data.frame(number=c(num.high,num.low,num.rollingmean)); rownames(stats) = c("greater max", "lower min", "deviation from moving mean")
print(stats)

file.out = read.zoo(file.filter,format="%Y-%m-%d %H:%M:%S",tz="GMT",index.column=1) #convert back to .zoo object
return(file.out)
}

### end filtering soil data ###
