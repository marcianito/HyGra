#' @title Convert DEM to data.frame
#'
#' @description interpolate DEM to necessary grid of corresponding layer
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

convert_demtodf = function(dem, dem.info){
# convert (melt) dem
dem_x = seq(dem.info[3,1],length.out=length(dem[1,]),by=dem.info[5,1])
dem_y = seq(dem.info[4,1],length.out=length(dem[,1]),by=dem.info[5,1])
colnames(dem) = dem_x
rownames(dem) = rev(dem_y)
dem.melt =  melt(as.matrix(dem))#, id.vars=
colnames(dem.melt) = c("y","x","z")
return(dem.melt)
}

#' @title Convert zoo-object to data frame
#'
#' @description Converting zoo-object to data frame with proper indexing of time column and column names
#'
#' @param value must be a zoo-object
#' @references Marvin Reich (2017), mreich@@posteo.de
# ' @examples examples MISSING
#' 
zootodf <- function(value) {
    df <- as.data.frame(value)
    df$time <- index(value) #create a Date column
    rownames(df) <- NULL #so row names not filled with dates
    df <- df[,c(ncol(df), 1:(ncol(df)-1))] #reorder columns so Date first
    return(df)
}

#' @title Calculate number of decimal places in a number
#'
#' @description 
#'
#' @param x Numeric, the number to check decimal places.
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}

#' @title Construct vertices for SG building foundation
#'
#' @description Following ptinpoly-packages, this function creates vertices
#' which can be used to correct gravity component grids for building structures
# (e.g. SG pillar, building baseplate, building walls, etc.).
#'
#' @param x,y,z_cords Vector, holding min and max values of the x,y and z points which fix the polyheder.
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

construct_vertices = function(
                    x_cords,
                    y_cords,
                    z_cords
){
    vertice = as.matrix(data.frame(
            x = c(min(x_cords), min(x_cords), min(x_cords), min(x_cords), max(x_cords), max(x_cords), max(x_cords), max(x_cords)),
            y = c(min(y_cords), min(y_cords), max(y_cords), max(y_cords), min(y_cords), min(y_cords), max(y_cords), max(y_cords)),
            z = c(min(z_cords), max(z_cords), min(z_cords), max(z_cords), min(z_cords), max(z_cords), min(z_cords), max(z_cords))
    ))
    return(vertice)
}


#' @title Normalize data 
#'
#' @description Normalize data in subtracting the mean and dividing by the standard deviation.
#'
#' @param x input dataset
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export
#' 

normalize = function(x){
	data_norm = (x - mean(x, na.rm=T))/sd(x,na.rm=T)
	return(data_norm)
}

#' @title Normalize data only sbutracting mean
#'
#' @description Normalize data in subtracting the mean and dividing by the standard deviation.
#'
#' @param x input dataset
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export
#' 

normalize_mean = function(x){
	data_norm = (x - mean(x, na.rm=T))
	return(data_norm)
}

#' @title Calculate differences between last and first value of a vector
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2018), mreich@@posteo.de
#' @import dplyr
#' @examples missing

dif_firstLast = function(
   data_in
){
   dif_val = dplyr::last(data_in) - dplyr::first(data_in)
   return(dif_val)
}

#' @title Read Digital Elevation Models 
#'
#' @description Read DEM's and output information files for further usage. 
#'
#' @param dempath Path of input DEM-file. 
#' @param filename DEM-file in .acs-format (see details)
#' 
#' @details So far only it is only possible to read DEM-file in .asc arquitecture format.
#' @details Outputs files: data.DEM, info.DEM (as data.frames).
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing
 
read_dem = function(dempath, filename){
    info.DEM = read.table(file=paste(dempath,filename,sep=""),nrows=6, row.names=1, header=F, sep="", stringsAsFactors=F, colClasses=c("character","numeric"))
    nodata = info.DEM[6,1]
    data.DEM = read.table(file=paste(dempath,filename,sep=""), sep="", dec=".", skip=6, stringsAsFactors=F, na.strings=nodata, header=F)
    dem.info <<- info.DEM; dem<<- data.DEM
    sucess = "output: dem.info and dem"
    return(sucess)
}

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
plot_dem <- function(dem,dem.info,locs = NA){
## DEBUGGING
# locs = SGlocs  
##
dem.col = topo.colors(10)
if(!is.na(locs)){
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
dem.plot = ggplot(melt(as.matrix(dem)), aes(x=Var2,y=Var1, fill=value)) +
    geom_raster() + xlab("x") + ylab("y") + scale_fill_gradientn(limits=range(dem),colours=dem.col, name="height mNN", na.value = "grey50")  +
    labs(title="DGM") + 
    geom_point(data=locs, aes(x,y), colour="black", size=3) + 
    geom_text(data=locs, aes(x,y,label=name),hjust=0.5, vjust=-0.5)
}else{
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
dem.plot = ggplot(melt(as.matrix(dem)), aes(x=Var2,y=Var1, fill=value)) +
    geom_raster() + xlab("x") + ylab("y") + scale_fill_gradientn(limits=range(dem),colours=dem.col, name="height mNN", na.value = "grey50")  +
    labs(title="DGM")
}
# return data  
return(dem.plot)
}

### discretization ###
discre <- function(C,R,x0,y0,dxdy){
  
  xl=matrix(0, nrow=1, ncol=C) #=zeros(1,C);
  xl[1]=x0;
  for (c1 in 1:(C-1)){
    xl[c1+1]=xl[c1]+dxdy
  }
  xr=xl+dxdy
  xp=(xl+xr)/2
  
  yl=matrix(0, nrow=1, ncol=R) #=zeros(1,R);
  yl[R]=y0
  for (c2 in seq(R,2, by=-1)){ #R:-1:2
    yl[c2-1]=yl[c2]+dxdy
  }
  yr=yl+dxdy
  yp=(yl+yr)/2
  
  #making variables globally accessable
  xl <<- xl; xr <<- xr; xp <<- xp
  yl <<- yl; yr <<- yr; yp <<- yp
  #res <- list()
  #print("xl, xr, xp"); print("yl, yr, yp")
  #return(print("discretization parameters calculated"))
}

#' @title RMSE of simulated and observed vectors
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2018), mreich@@posteo.de
#' @import 
#' @examples missing

rmse = function(
    ts_sim,
    ts_obs
){
  ts_df = data.frame(ts_sim,ts_obs)
  ts_df$sqr_dif = (ts_df$ts_sim - ts_df$ts_obs)^2
  rmse_val = sqrt(sum(ts_df$sqr_dif) / length(ts_sim))
  return(rmse_val)
}

#' @title tif_to_ascii: converts geoTiff fromat raster data into arc ascii format
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2018), mreich@@posteo.de
#' @examples missing

tif_to_ascii = function(
            file_name_in,
            file_name_out,
            input_dir,
            output_dir
){
    #get name of file
    # file_name = 
    # 
    DEM_input = raster::raster(paste0(input_dir, file_name_in, ".tif"))
    # raster::writeRaster(DEM_input, filename = paste0(output_dir, file_name_out, "_DEM.asc"), format = "ascii") #, datatype="INT4S")
    raster::writeRaster(DEM_input, filename = paste0(output_dir, file_name_in, "/", file_name_out, ".asc"), format = "ascii") #, datatype="INT4S")
    return(NULL)
}

#' @title Aggregate data.frame time series to other time periods
#'
#' @description test
#'
#' @param timeseries_data data.frame, formatted with columns [datetime, values].
#' @param newPeriod charater string, options are "hourly", "daily", "weekly".
#' @param fun character string, indicating aggregation function.
#' @param time_offset integer, indicating the offset to include in the aggregation.
#' Units are with respect to declared new period.
#' 
#' @return test
#' 
#' @details Convention for setting dates corresponding to new period: POSIXct is
#' always used from the LAST hour, day, etc. of the aggregation period.
#' @references Marvin Reich (2018), mreich@@posteo.de
#' @import dplyr
#' @examples missing

aggTS = function(
    timeseries_data,
    newPeriod,
    fun = "sum",
    time_offset = 0,
    conserve_columns = NA
){
    ## DEBUGGING
    # timeseries_data = as.data.frame(ET_hourly)
    # timeseries_data = as.data.frame(data_radar_mod)
    # timeseries_data = tt_data_radar
    # newPeriod = "hourly"
    # fun = "sum"
    # time_offset = 0
    # conserve_columns = "cell_index"
    # conserve_columns = c("x", "y")
    ##
    ## force input to be a data.frame
    timeseries_data = as.data.frame(timeseries_data)
    # construct list of columns to conserve during summarizing
    if(!is.na(conserve_columns)){
    columns = as.list(c("datetime", conserve_columns))
    }else{
    columns = as.list(c("datetime"))
    }
    # format after new period
    # and include offset
    switch(newPeriod,
           daily = {
        ts_newPeriod = timeseries_data %>%
        dplyr::mutate(datetime = as.POSIXct(trunc(datetime, "days"))) %>%
        dplyr::mutate(datetime = datetime + 3600 * 24 * time_offset)
           },
           hourly = {
        ts_newPeriod = timeseries_data %>%
        dplyr::mutate(datetime = as.POSIXct(trunc(datetime, "hours"))) %>%
        dplyr::mutate(datetime = datetime + 3600 * time_offset)
           },
           weekly = {
               ## still not sure how to do this !! ?
               ## now, grouping should be over "datetime_weeks"
        # week needed for correct time setting:
        # convention: end of week
        week = 3600 * 24 * 7
        ts_newPeriod = timeseries_data %>%
        dplyr::mutate(datetime_weeks = paste0(format(datetime, "%Y-%W"), "-0")) %>%
        dplyr::mutate(datetime = as.POSIXct(strptime(datetime_weeks, "%Y-%W-%w")) + week) %>%
        dplyr::mutate(datetime = datetime + 3600 * 24 * 7 * time_offset) %>%
        dplyr::select(-datetime_weeks)
           }
           )
    #
    # format after aggregation function
    switch(fun,
           sum = {
        ts_newPeriod = ts_newPeriod %>%
        # dplyr::group_by(datetime) %>%
        dplyr::group_by_(.dots = columns) %>%
        # dplyr::summarize(observed = sum(value, na.rm = T))
        dplyr::summarize(value = sum(value, na.rm = T))
           },
           mean = {
        ts_newPeriod = ts_newPeriod %>%
        # dplyr::group_by(datetime) %>%
        dplyr::group_by_(.dots = columns) %>%
        # dplyr::summarize(observed = mean(value, na.rm = T))
        dplyr::summarize(value = mean(value, na.rm = T))
           }
           )
    #
    # convert time pack to POSIXct
    ts_newPeriod$datetime = as.POSIXct(ts_newPeriod$datetime)
    # return new time series data
    return(ts_newPeriod)
}




