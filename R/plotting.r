##########################
### plotting functions ###
##########################

#' @title Plot a DEM 
#'
#' @description Plots a Digital Elevation Model / Matrix using ggplot2-style. 
#'
#' @param dem x * y-grid containing z-coordinates 
#' @param dem.info data.frame containing information about the DEM (row wise): columns, rows, starting x, starting y, lengths of one cellsize
#' @param locs data.frame with column strucure $x (coordinate), $y (coordinate), $name (name to be printed in plot)
#' ...
#' @details missing 
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing 
#' 
plot_dem <- function(dem,dem.info,locs){
# require(reshape2)
# require(ggplot2)
  
dem.col = topo.colors(10)
locs$value=c(rep(1,times=length(locs[,1])))
#dem parameters
R = dem.info[2,1]
C = dem.info[1,1]
x0 = dem.info[3,1]
y0 = dem.info[4,1]
dxdy = dem.info[5,1]
#calculate dem-parameters
discre(C,R,x0,y0,dxdy)
#change dem cols & rows
colnames(dem) = xp
rownames(dem) = yp
dem.plot = ggplot(melt(as.matrix(dem)), aes(x=Var2,y=Var1, fill=value)) + geom_raster() + xlab("x") + ylab("y") + scale_fill_gradientn(limits=range(dem),colours=dem.col, name="height mNN", na.value = "grey50")  + labs(title="DGM") + geom_point(data=locs, aes(x,y), colour="black", size=3) + geom_text(data=locs, aes(x,y,label=name),hjust=0.5, vjust=-0.5)
  
return(dem.plot)
  
}

### end plot DEM ###

#' @title Plot soil data
#'
#' @description Generally designed for plotting soil moisture data of different sensor profiles
#'
#' @param file dataset of measured values as a zoo object
#' @param show profile (data column) to plot
#' @param top limiting data values at top
#' @param bottom limiting data values at bottom
#' @param plottogether plot all data in one single plot (TRUE) or subplots for each timeseries (FALSE; default)
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 
plot_soildata <- function(file,show,top,bottom,plottogether=F){
# require(ggplot2)
# require(reshape2)
# require(plyr)
data.prep = cbind(index(file),as.data.frame(file,row.names=NULL))
colnames(data.prep)[1]="time"
data.melt = melt(data.prep, id.vars="time",measure.vars=colnames(file))#data.melt = melt(t1, id.vars="index(SGold.agg)",measure.vars=colnames(file))
data.melt$mux=NA; data.melt$sensor=NA
muxes = unique(sapply(strsplit(colnames(file), "_"), "[[", 1)) #c("mux00","mux01","mux02","mux03","mux04","mux05","mux06")
sensors = sort(unique(sapply(strsplit(colnames(file), "_"), "[[", 2)))#c("1","2","3","4","5","6","7","8")
for(i in 1:length(muxes)) data.melt$mux[grep(muxes[i], as.character(data.melt$variable))] = muxes[i] #matching=pmatch("mux", as.character(data.melt$mux))
for(i in 1:length(sensors)) data.melt$sensor[grep(paste("_",sensors[i],sep=""), as.character(data.melt$variable))] = sensors[i]
data.melt$sensor = factor(data.melt$sensor)

data.plot = subset(data.melt,mux==show & value < top & value > bottom)# & sensor=="7") #ddply(subset(data.melt,data.melt$mux==plotvar),.(layer),summarise,values=mean(value))
if(plottogether==T){
	plot.final=ggplot(data.plot, aes(x=time, y=value, colour=sensor)) + geom_line() + ylab("soil moisture [%]") + labs(title=show)
}
else{
	plot.final=ggplot(data.plot, aes(x=time, y=value, colour=sensor)) + geom_line() + facet_grid(sensor~.) + ylab("soil moisture [%]") + labs(title=show)
}

print("mux shown:")
print(show)
#print(muxes)
#print(sensors)
return(plot.final)
  
}

### end plot soil data ###

#' @title Displaying gravity results
#'
#' @description Visualizes gravity data array (2D)
#'
#' @param grav.array array of values [rows, columns, layers] containing data (usually gravity components or gravity signals), both must have the same dimensions
#' @param grav.name name of device used in graphic
#' @param vertsum logical, whether values will show vertical sums (TRUE) or only one layer (FALSE), default FALSE
#' @param layer layer to show, only affective when vertsum=FALSE
#' @param cuts determinating breaks categories of input values
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' 
gplot_array <- function(grav.array, #gravity array
                  grav.name, #name of device
                  vertsum=F, #layerwise or vertical sums
                  layer, #layer to display
                  cuts = c(-10,-2,-1,-0.5,-0.2,-0.1,-0.01,-0.001,0,0.001,0.01,0.1,0.2,0.5,1,2,10)){
  #dem.show=F;dem
  #load needed packages
# require(reshape2)
# require(ggplot2)

#calculations
g1.sum = apply(grav.array, c(1,2), sum, na.rm=T) #grid

#preparing data
#layerwise or vertical sums
if(vertsum==F){ #plot only indicated layer
#data perparation
  type="Layer"
  dev=grav.name; g1 = cbind(melt(as.matrix(grav.array[,,layer])),dev,layer,type)
  plot.data = g1
  #   plot.data = rbind(g1,g2,gdif)
#automatic boundaries
  #data.range = 
  #cuts = cutnum
  cols = colorRampPalette(c("darkgreen", "orange", "darkblue"))(length(cuts))
  plot.data$dimfill = cut(plot.data$value,cuts)
  #facetting
  plot.final = ggplot(plot.data, aes(Var2,rev(Var1), fill=dimfill)) + geom_raster() + xlab("x") + ylab("y") + labs(title=paste(type,layer,sep=" ")) + facet_wrap(~ dev, ncol=3)  + scale_fill_manual(values=cols, name="nm/s²")
}
else{ #plot vertical sums per grid-cell
  type="Vertical layer sums"
  dev=grav.name; g1 = cbind(melt(g1.sum),dev,type)
  plot.data = g1
  #   plot.data = rbind(g1,g2,gdif)
#automatic boundaries
  #data.range = 
  #cuts = cutnum
  cols = colorRampPalette(c("darkgreen", "orange", "darkblue"))(length(cuts))
  plot.data$dimfill = cut(plot.data$value,cuts)
  #facetting
  plot.final = ggplot(plot.data, aes(Var2,rev(Var1), fill=dimfill)) + geom_raster() + xlab("x") + ylab("y") + labs(title=type) + facet_wrap(~ dev, ncol=3)  + scale_fill_manual(values=cols, name="nm/s²")
}
  
#plotting
return(plot.final)
}

### end displaying gravity output ###
