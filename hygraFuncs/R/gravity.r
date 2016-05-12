#################################################
### functions concerning gravity calculations ###
#################################################

#' @title Calculate gravity component
#'
#' @description Calculates the gravity components (100percent values) for a given DEM using a nested approach
#'
#' @param dem DEM in UTM
#' @param dem.info information to the used DEM, see details
#' @param dz vector of thickness of each layer horizon, number of values equals number of layers used for calculation
#' @param sgx,sgy coordinates of gravimeter location
#' @param sgdz high of gravimeter sensor above DEM
#' @param umbrella logical, wether to compute excluding the umbrella effect, default TRUE
#' @param sghouse coordinates of the gravimeter house, used as umbrella surface
#' @param umbrella.max maximal depth until which umbrella effect will be included
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export

gcomp <- function(dem,dem.info,dz,sgx,sgy,sgdz,umbrella=F,sghouse,umbrella.max){
                  #w=T #change units
                  #output=c(array,grid,layers,value)
#constants
gama <- 6.673e-11 #m³/(kg*s²)
rho <- 1000 #kg/m³
#w <- 1e8 #gravity units µGal
w <- 1e9 #gravity units nm/s²
#read DEM information  
Row = dem.info[2,1]
Col = dem.info[1,1]
x0 = dem.info[3,1]
y0 = dem.info[4,1]
dxdy = dem.info[5,1]
dx = dxdy
dy = dxdy
#calculate DEM parameters
discre(Col,Row,x0,y0,dxdy)
#change DEM cols & rows to coordinates
colnames(dem) = xp
rownames(dem) = yp
#z coordinate SG
sgz=dem[which(match(yp,sgy)==T),which(match(xp,sgx)==T)]+sgdz
#number of z layer
L=length(dz)
#calculate absolute z-coord from each layer, based on dz and L
absZint = vector()
absZend = vector()
for(l in 1:L){
  if(l == 1){absZint[l] = 0; absZend[l] = absZint[l] + dz[l]; next} #start of vertical depth
  absZint[l] = absZint[l-1] + dz[l-1]
  absZend[l] = absZint[l] + dz[l]
}
#SG house / umbrella effect area
if(umbrella==T){
lo = sghouse[1,]
lu = sghouse[2,]
ro = sghouse[3,]
ru = sghouse[4,]
}
#radius criteria
r2exac=2^2 #bei 2.5 x 2.5 cellsizes ~=   50 m
r2macm=9^2 #bei 2.5 x 2.5 cellsizes ~= 1000 m
#prepare SG-data-file
sg_gridRCL=array(0, dim=c(Row,Col,L))
# missing nested approach with forsberg and pointmass!!
##gravity after macmillan
print("calculating gravity components....patience..")
for (l in 1:L){
      for (cc in 1:Col){
        for (rr in 1:Row){
##debugging
# browser()
##debugging
        #z.values  
        Zint=dem[rr,cc] - absZint[l]
        Zend= dem[rr,cc] - absZend[l]
        zp=(Zint+Zend)/2 #z midpoint of cell
        if(is.na(zp) == TRUE){next} #next cell(x,y) for NA on DEM-information
        #umbrella effect
        if(umbrella==T){
        if(xp[cc] >= lo[1] & xp[cc] >= lu[1] & xp[cc] <= ro[1] & xp[cc] <= ru[1] & yp[rr] >= lu[2] & yp[rr] >= ru[2] & yp[rr] <= lo[2] & yp[rr] <= ro[2] & Zend >= (dem[rr,cc]-umbrella.max)){
        print(c(cc,rr,l)); sg_gridRCL[rr,cc,l]=0 ;next}
        }
        #distances
        rad = sqrt((xp[cc]-sgx)^2+(yp[rr]-sgy)^2+(zp-sgz)^2) #radial distance to SG
        r2=rad^2
        dr2=dx^2+dy^2+(dz[l])^2 #radial "size" of DEM / coordinate-data system
        f2=r2/dr2 #abstand zelle-SG / diagonale aktueller berechnungs-quader
        
        # different methods after the distance from mass to SG
                if (f2<=r2exac){ #very close to SG
                  sg_gridRCL[rr,cc,l]=forsberg(gama,w,rr,cc,xl,xr,yl,yr,Zint,Zend,sgx,sgy,sgz,rho) #unit depends on w
                }
                if(f2>r2macm){ #very far from SG
                  sg_gridRCL[rr,cc,l]=pointmass(gama,zp,sgz,dx,dy,dz[l],rad,w,rho) #unit depends on w
                }
                if(f2>r2exac & f2<r2macm){ #in the "middlle"
                  sg_gridRCL[rr,cc,l]=macmillan(gama,rr,cc,xp,yp,zp,sgx,sgy,sgz,dx,dy,dz[l],rad,w,rho) #unit depends on w
                }
        
        #zp=dem[r,c]+z[l]
        #rad = sqrt((xp[c]-sgx)^2+(yp[r]-sgy)^2+(zp-sgz)^2) #radial distance to SG
        #sg_gridRCL[r,c,l]=macmillan(gama,r,c,xp,yp,zp,sgx,sgy,sgz,dx,dy,dz,rad,w,rho) #in µGal
        }
      }
    }  
print("..finished!")
##generating output
#comp.array = sg_gridRCL #output: array
comp.layer = vector() #one value per layer
comp.grid = matrix(0, ncol=Col, nrow=Row) #one layer in spatial extend R x C
for(l in 1:L){
comp.layer[l] = sum(sg_gridRCL[,,l], na.rm=T) #output: layer
comp.grid= comp.grid + sg_gridRCL[,,l] #output: grid
}
comp.value = sum(sg_gridRCL,na.rm=T) #output: value

return(sg_gridRCL) #standard output = array

}

### end gcomp ###

#' @title Displaying gravity results
#'
#' @description Visualizes gravity data from two devices in varies posibilities
#'
#' @param grav1/2.array array of values [rows, columns, layers] containing data (usually gravity components or gravity signals), both must have the same dimensions
#' @param grav1/2.name name of devices used in graphic
#' @param vertsum logical, whether values will show vertical sums (TRUE) or only one layer (FALSE), default FALSE
#' @param layer layer to show, only affective when vertsum=FALSE
#' @param onlydif logical, should absolute values be shown of layer x grav1 and grav2 (FALSE) or should each layer be shown as difference grav1 - grav2, default FALSE
#' @param cuts determinating breaks categories of input values
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export
#' 
gdisp <- function(grav1.array, #gravity array 1
                  grav1.name, #name of device
                  grav2.array, #gravity array 2
                  grav2.name, #name of device
                  vertsum=F, #layerwise or vertical sums
                  layer, #layer to display
                  onlydif=F, #show all layers but with differences
                  cuts = c(-10,-2,-1,-0.5,-0.2,-0.1,-0.01,-0.001,0,0.001,0.01,0.1,0.2,0.5,1,2,10)){
  #dem.show=F;dem
  #load needed packages
# require(reshape2)
# require(ggplot2)

#calculations
gdif.array = grav1.array - grav1.array #array
g1.sum = apply(grav1.array, c(1,2), sum, na.rm=T) #grid
g2.sum = apply(grav2.array, c(1,2), sum, na.rm=T) #grid
gdif.sum = apply(gdif.array, c(1,2), sum, na.rm=T) #grid

#preparing data
#individual or only sg-differences
if(onlydif==F){ #individual
#layerwise or vertical sums
if(vertsum==F){ #plot only indicated layer
#data perparation
  type="Layer"
  dev=grav1.name; g1 = cbind(melt(grav1.array[,,layer]),dev,layer,type)
  dev=grav2.name; g2 = cbind(melt(grav2.array[,,layer]),dev,layer,type)
  dev="difference"; gdif = cbind(melt(gdif.array[,,layer]),dev,layer,type)
  plot.data = rbind(g1,g2,gdif)
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
  dev=grav1.name; g1 = cbind(melt(g1.sum),dev,type)
  dev=grav2.name; g2 = cbind(melt(g2.sum),dev,type)
  dev="difference"; gdif = cbind(melt(gdif.sum),dev,type)
  plot.data = rbind(g1,g2,gdif)
#automatic boundaries
  #data.range = 
  #cuts = cutnum
  cols = colorRampPalette(c("darkgreen", "orange", "darkblue"))(length(cuts))
  plot.data$dimfill = cut(plot.data$value,cuts)
  #facetting
  plot.final = ggplot(plot.data, aes(Var2,rev(Var1), fill=dimfill)) + geom_raster() + xlab("x") + ylab("y") + labs(title=type) + facet_wrap(~ dev, ncol=3)  + scale_fill_manual(values=cols, name="nm/s²")
}
}
else{ #sg-differences per layer
  type="Signal difference per layer"
  #i = layer;... -> automatically adjusting to number of layers
  dev="SGdif layer 1";layer="1"; gdif_l1 = cbind(melt(gdif.array[,,1]),dev,type,layer)
  dev="SGdif layer 2";layer="2"; gdif_l2 = cbind(melt(gdif.array[,,2]),dev,type,layer)
  dev="SGdif layer 3";layer="3"; gdif_l3 = cbind(melt(gdif.array[,,3]),dev,type,layer)
  dev="SGdif layer 4";layer="4"; gdif_l4 = cbind(melt(gdif.array[,,4]),dev,type,layer)
  dev="SGdif layer 5";layer="5"; gdif_l5 = cbind(melt(gdif.array[,,5]),dev,type,layer)
  dev="SGdif all layers";layer="all"; gdif = cbind(melt(gdif.sum),dev,type,layer)
  plot.data = rbind(gdif_l1,gdif_l2,gdif_l3,gdif_l4,gdif_l5,gdif)
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

#' @title Calulate gravity signals 
#'
#' @description After calculating the gravity components (gcomp()), this function can be used to get real values using soil moisture data (theta). The output of the gravity signal is possible in different dimensions. 
#'
#' @param gcompfiles vector containing the names (complete paths or relative) of the gravity component files to be used for calculations. These files are generated using gcomp().
#' @param theta data.frame with column structure $time, $timestep, $theta (value), $layer
#' @param output Defines the output to be a singel value, layer or grid (default is value)
#' ...
#' @details Usually one needs relative gravity signals, therefore the inputed theta timeseries should be acutally delta theta vaules per timestep. The
#' other option is calculating the differences in values per timestep of the output of this function.
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing 
#' @export
#' 
gsignal <- function(gcompfiles,theta,output="value"){
                    #theta strucutre: time-timestep-layer-x-y-value ;so later the same structure can be used for grids !!
library(plyr)
#get information
ts = max(theta$timestep)#number of timesteps
L = max(theta$layer) #number of layers
gcompnum = length(gcompfiles)
#create output-files
totalg=data.frame(time=unique(theta$time),timestep=unique(theta$timestep), g=0) #one value per TIMESTEP
totalg_layers=matrix(0, nrow=L, ncol=ts) #one value per TIMESTEP per LAYER
#totalg_grid= array....
theta$layer = factor(theta$layer)

for(i in 1:gcompnum){
  filename = load(file=gcompfiles[i]) #load gravity component files
  input.array = get(filename) #need to remove old object for space and memory!! ;rm()
  comp.layer = apply(input.array, 3, sum, na.rm=T) #layer
  #comp.grid = apply(input.array, c(1,2), sum, na.rm=T) #grid
  totalg_grid=array(data=0,dim=c(dim(input.array)[1:2],ts)) #set to automatic dim-detection !!!
    for (TS in 1:ts){
	    # browser()
	    # this line seemed to calculate wrongfully ALL theta values to 0.1 !?! reason: summarise?
	    # totalg$g[TS] = totalg$g[TS] + sum(comp.layer*ddply(subset(theta,theta$timestep==TS),.(layer),summarise,values=mean(value))$values, na.rm=T)
	    theta.calc = theta
    totalg$g[TS] = totalg$g[TS] + sum(comp.layer*subset(theta.calc, theta.calc$timestep==TS)$theta, na.rm=T)
    # totalg$g[TS] = totalg$g[TS] + sum(comp.layer*subset(theta,theta$timestep==TS)$theta, na.rm=T)
    # totalg_layers[,TS] = totalg_layers[,TS] + (comp.layer*ddply(subset(theta,theta$timestep==TS),.(layer),summarise,values=mean(value))$values)
    totalg_layers[,TS] = totalg_layers[,TS] + (comp.layer*subset(theta.calc,theta.calc$timestep==TS)$theta)
      if(output=="grid"){
      for(l in 1:L){
      totalg_grid[,,TS] = totalg_grid[,,TS] + (input.array[,,l]*ddply(subset(theta,theta$timestep==TS),.(layer),summarise,values=mean(value))$values[l])          
      }      
      }
    }
}
switch(output,
       value=return(totalg),
       layer=return(totalg_layers),
       grid=return(totalg_grid))
}

### end gravity signal ###

#' @title Calculate gravity signal from snow storage
#'
#'
#' @param snow.input data.frame with columns $time, $SWE (snow water equivalent)
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 
gsnow <- function(snow.input){
  
snow.g = -3.9e-6 * (snow.input[,2])^2 - 9e-4*snow.input[,2]
#snow.out = as.data.frame(cbind(snow.input[,1], snow.g))
snow.out = cbind(snow.input[,1], snow.g)
colnames(snow.out) = c("time", "gsnow")
return(snow.out)
}

### end gsignal snow

#' @title Calculate gravity signal from water storage changes
#'
#' @description Outputs the gravity signal of each storage deparment where water mass changes occured
#'
#' @param WSC.input data.frame with information about $time and further columns with the timeseries of water storage changes in [mm]
#' @param comp.layer vector of gravity component of each layer(100\% values), corresponding to a water storage
#' 
#' @details The number of components has to match the number of columns in WSC.input (+ column $time)
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 
gWSC <- function(WSC.input, comp.layer){ #input: comp.names
  
if((dim(WSC.input)[2]-1) != length(comp.layer)) print("number WSC compartments and layers are not equal")
n = length(comp.layer)
WSC.out = WSC.input[,1]

for(i in 1:n){
WSC.out = cbind(WSC.out,WSC.input[,i+1]*comp.layer[i])
colnames(WSC.out)[i+1]=paste("gP",i, sep="")
}
return(WSC.out)
}

### end gsignal WSC ###


#' @title Interpolate DEM to layer grid (3D)
#'
#' @description interpolate DEM to necessary grid of corresponding layer
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

demtogrid = function(dem, dem.info, grid_discr, depth_split){
library(gstat)
# convert (melt) dem
dem_x = seq(dem.info[3,1],length.out=length(dem[1,]),by=dem.info[5,1])
dem_y = seq(dem.info[4,1],length.out=length(dem[,1]),by=dem.info[5,1])
colnames(dem) = dem_x
rownames(dem) = dem_y
dem.melt =  melt(as.matrix(dem))#, id.vars=
colnames(dem.melt) = c("y","x","z")
#generate regular-spaced grid
# temporarily hardgecoded hydrus 3d Area
grid.x <- seq(4564029.26662,4564049.35606 , by=grid_discr[1])
grid.y <- seq(5445650.34147,5445670.49548 , by=grid_discr[2])
#grid.x <- seq(min(dem.melt$x), max(dem.melt$x), by=grid_discr[1])
#grid.y <- seq(min(dem.melt$y), max(dem.melt$y), by=grid_discr[2])
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr[3])
grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)
grid.surface <- expand.grid(x=grid.x, y=grid.y)

#interpolate dem data to new grid
idw.gstat = gstat(formula = z ~ 1, locations = ~ x + y, data = dem.melt, nmax = 4, set = list(idp = 2))
surface = predict(idw.gstat, grid.surface)
gcomp_grid = cbind(grid.xyz[,1:2],z=surface[,3]-grid.xyz$Depth,zgrid=grid.xyz$Depth, layer=rep(seq(1,length(grid.z), by=1),each=(length(grid.x)*length(grid.y))))
return(gcomp_grid)
} #end function

#' @title Interpolate DEM to layer grid (2D)
#'
#' @description interpolate DEM to necessary grid of corresponding layer
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

demtogrid_2d = function(dem, dem.info, grid_discr, depth_split, grid.y){
library(gstat)
# convert (melt) dem
dem_x = seq(dem.info[3,1],length.out=length(dem[1,]),by=dem.info[5,1])
dem_y = seq(dem.info[4,1],length.out=length(dem[,1]),by=dem.info[5,1])
colnames(dem) = dem_x
rownames(dem) = dem_y
dem.melt =  melt(as.matrix(dem))#, id.vars=
colnames(dem.melt) = c("y","x","z")
#generate regular-spaced grid
# temporarily hardgecoded hydrus 3d Area
grid.x <- seq(4564029.26662,4564049.35606 , by=grid_discr[1])
#grid.y <- seq(5445650.34147,5445670.49548 , by=grid_discr[2])
#grid.x <- seq(min(dem.melt$x), max(dem.melt$x), by=grid_discr[1])
#grid.y <- seq(min(dem.melt$y), max(dem.melt$y), by=grid_discr[2])
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr[3])
grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)
grid.surface <- expand.grid(x=grid.x, y=grid.y)

#interpolate dem data to new grid
#idw.gstat = gstat(formula = z ~ 1, locations = ~ x + y, data = dem.melt, nmax = 4, set = list(idp = 2))
idw.gstat = gstat(formula = z ~ 1, locations = ~ x, data = select(dem.melt,x,z), nmax = 4, set = list(idp = 2))
surface = predict(idw.gstat, grid.surface)
#gcomp_grid = cbind(grid.xyz[,1:2],z=surface[,3]-grid.xyz$Depth, layer=rep(seq(1,length(grid.z), by=1),each=(length(grid.x)*length(grid.y))))
gcomp_grid = cbind(grid.xyz[,1],z=surface[,3]-grid.xyz$Depth, layer=rep(seq(1,length(grid.z), by=1),each=(length(grid.x)*length(grid.y))))
return(gcomp_grid)
} #end function

#' @title Interpolate DEM to layer SURFACE grid (3D)
#'
#' @description interpolate DEM to necessary grid of corresponding layer
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

demtoSurfacegrid = function(dem, dem.info, grid_discr, depth_split){
library(gstat)
# convert (melt) dem
dem_x = seq(dem.info[3,1],length.out=length(dem[1,]),by=dem.info[5,1])
dem_y = seq(dem.info[4,1],length.out=length(dem[,1]),by=dem.info[5,1])
colnames(dem) = dem_x
rownames(dem) = dem_y
dem.melt =  melt(as.matrix(dem))#, id.vars=
colnames(dem.melt) = c("y","x","z")
#generate regular-spaced grid
grid.x <- seq(min(dem.melt$x), max(dem.melt$x), by=grid_discr[1])
grid.y <- seq(min(dem.melt$y), max(dem.melt$y), by=grid_discr[2])
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr[3])
grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)
grid.surface <- expand.grid(x=grid.x, y=grid.y)

#interpolate dem data to new grid
idw.gstat = gstat(formula = z ~ 1, locations = ~ x + y, data = dem.melt, nmax = 4, set = list(idp = 2))
surface = predict(idw.gstat, grid.surface)
return(surface)
} #end function

#' @title Calculate gravity component, WITHOUT umbrella effect
#'
#' @description Calculates the gravity components (100percent values) for a given DEM using a nested approach
#'
#' @param g_grid grid for which to calculate gravity components. (generally derived from DEM of region in combination with hydrological modeling mesh/grid).
#' @param senloc coordinates of gravity sensor location. data.frame with x,y,z columns.
#' @param g_discr discretization of g_grid in differences in x,y,and z-direction. they have to be uniform steps!
#' @param edge if grid to calulate is located on a domain edge (top or bottom). for top set edge='first', for bottom set edge='last'.
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export

gcomp_raw <- function(g_grid, senloc, g_discr, edge){
#constants
gama <- 6.673e-11 #m³/(kg*s²)
rho <- 1000 #kg/m³
#w <- 1e8 #gravity units µGal
w <- 1e9 #gravity units nm/s²
#radius criteria
#r2exac=2^2 #bei 2.5 x 2.5 cellsizes ~=   50 m
#r2macm=9^2 #bei 2.5 x 2.5 cellsizes ~= 1000 m
r_inner = 50 #[m]
r_outer = 1000 #[m]
## function: heart of gravity component calculation
# estiamte distances of actual point and gravity sensor
# calculate needed gravity cell (volumina) information / locations..
# cell midpoint = actual grid (mesh) points (x,y,z)
# choose gravity estimation method after distance of cell to gravity sensor
select_gMethod = function(xloc, yloc, zloc, gloc, gdiscr, edge, layer, layermax){
        #distances
        rad = sqrt((xloc-gloc$x)^2+(yloc-gloc$y)^2+(zloc-gloc$z)^2) #radial distance to SG
        r2=rad^2
        dr2=gdiscr$x^2+gdiscr$y^2+gdiscr$z^2 #radial "size" of DEM / coordinate-data system
        f2=r2/dr2 #abstand zelle-SG / diagonale aktueller berechnungs-quader
        # different methods after the distance from mass to SG
        #if (f2<=r2exac){ #very close to SG
        if (rad <= r_inner){ #very close to SG
		xl=xloc-(0.5*gdiscr$x);xr=xloc+(0.5*gdiscr$x)
		yl=yloc-(0.5*gdiscr$y);yr=yloc+(0.5*gdiscr$y)
		if(layer==1 | layer == layermax){
			if(edge=="first" & layer==1){ # first z-layer
			 Zint =zloc 
			 Zend = zloc-(0.5*gdiscr$z)
			}
			if(edge =="last" & layer==layermax){ # last z-layer
			 Zint = zloc+(0.5*gdiscr$z)
			 Zend = zloc
			}
			else{
			 Zint = zloc+(0.5*gdiscr$z)
			 Zend = zloc-(0.5*gdiscr$z)
			}
		}
		else{ # all other z-layers
		Zint = zloc+(0.5*gdiscr$z)
		Zend = zloc-(0.5*gdiscr$z)
		}
           gcomp_cell=forsberg_raw(gama,w,xl,xr,yl,yr,Zint,Zend,gloc$x,gloc$y,gloc$z,rho) #unit depends on w
        }
         #if(f2>r2macm){ #very far from SG
         if(rad >= r_outer){ #very far from SG
		if(layer==1 | layer == layermax){
			if(edge=="first" & layer==1){ # first z-layer
		  	 zdiscr = 0.5*gdiscr$z
		 	 zloc_mid = zloc+0.5*zdiscr
			}
			if(edge =="last" & layer==layermax){ # last z-layer
			 zdiscr = 0.5*gdiscr$z
		  	 zloc_mid = zloc-0.5*zdiscr
			}
			else{
		 	zloc_mid = zloc
		 	zdiscr = gdiscr$z
			}
		}
		else{ # all other z-layers
		 zloc_mid = zloc
		 zdiscr = gdiscr$z
		}
           gcomp_cell=pointmass(gama,zloc_mid,gloc$z,gdiscr$x,gdiscr$y,zdiscr,rad,w,rho) #unit depends on w
        }
        #if(f2>r2exac & f2<r2macm){ #in the "middlle"
        if(rad > r_inner & rad < r_outer){ #in the "middlle"
		if(layer==1 | layer == layermax){
			if(edge=="first" & layer==1){ # first z-layer
		  	 zdiscr = 0.5*gdiscr$z
		 	 zloc_mid = zloc+0.5*zdiscr
			}
			if(edge =="last" & layer==layermax){ # last z-layer
			 zdiscr = 0.5*gdiscr$z
		  	 zloc_mid = zloc-0.5*zdiscr
			}
			else{
		 	zloc_mid = zloc
		 	zdiscr = gdiscr$z
			}
		}
		else{ # all other z-layers
		 zloc_mid = zloc
		 zdiscr = gdiscr$z
		}
           gcomp_cell=macmillan_raw(gama,xloc,yloc,zloc_mid,gloc$x,gloc$y,gloc$z,gdiscr$x,gdiscr$y,zdiscr,rad,w,rho) #unit depends on w
        }
# output one value: gravity component for corresponding input cell
return (gcomp_cell)
}
#prepare SG-data-file
# gravity component is just a new column on g_grid
# containing the g_comp, conditionally to some causes, calculated using one of the 3 methods
# !! important
# due to maintaining original grid, gravity component cells (volumini) of the first row (surface) and last row (lower boundary),
# have to be adjusted to only HALF SIZE (value), due to cut off
# advantage of this approach:
# maintaining org. grid facilitates trasnformation of (hydrological) model output
# because no grid_transformation is needed
# if NOT, use layer "would be lost"
# due to midpoint of g component cell (volumina)
# here the calculations for each line (grid point)
## ! problem with dplyr-way: 
#g_grid_result = mutate(g_grid, gcomp = select_gMethod(x,y,z,senloc,g_discr))

for(i in 1:length(g_grid[,1])){
g_grid$gcomp[i] = select_gMethod(g_grid$x[i],g_grid$y[i],g_grid$z[i],senloc,g_discr,edge, g_grid$layer[i], max(g_grid$layer))
}

## for first and last row: cut in halt g_component
## depending on y..
## first row
#switch(edge,
       #first={
	#g_grid$gcomp[which(g_grid$z == min(unique(g_grid$z)))] = 0.5 * g_grid$gcomp[which(g_grid$z == min(unique(g_grid$z)))]
       #},
       #last={
	#g_grid$gcomp[which(g_grid$z == max(unique(g_grid$z)))] = 0.5 * g_grid$gcomp[which(g_grid$z == max(unique(g_grid$z)))]
       #},
       #both={
	#g_grid$gcomp[which(g_grid$z == min(unique(g_grid$z)))] = 0.5 * g_grid$gcomp[which(g_grid$z == min(unique(g_grid$z)))]
	#g_grid$gcomp[which(g_grid$z == max(unique(g_grid$z)))] = 0.5 * g_grid$gcomp[which(g_grid$z == max(unique(g_grid$z)))]
       #}
       #)
# return results
#return(g_grid_result)
return(g_grid)
}

### end gcomp_raw ###

#' @title subtract umbrella EFFECT from gcomp (raw) grid
#'
#' @description Calculates the gravity components (100percent values) for a given DEM using a nested approach
#'
#' @param g_grid grid for which to calculate gravity components. (generally derived from DEM of region in combination with hydrological modeling mesh/grid).
#' @param senloc coordinates of gravity sensor location. data.frame with x,y,z columns.
#' @param g_discr discretization of g_grid in differences in x,y,and z-direction. they have to be uniform steps!
#' @param edge if grid to calulate is located on a domain edge (top or bottom). for top set edge='first', for bottom set edge='last'.
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export

subtract_umbrellaEffect = function(grid_input, u_polygon, u_depth, sgz){
	library(SDMTools)
	grid_xy = data.frame(x=grid_input$x, y=grid_input$y)
	in_poly = pnt.in.poly(grid_xy, u_polygon)
	z_depth = sgz - u_depth
	grid_xy_pip = cbind(grid_input,pip = in_poly$pip)
	grid_xy_pip$gcomp[which(grid_xy_pip$pip == 1 & grid_xy_pip$z > z_depth)] = 0
	#grid_output = filter(cbind(grid_input,pip = in_poly$pip), pip == 0 & z > z_depth) %>%
	grid_output = select(grid_xy_pip, x,y,z,gcomp)
# return result
return(grid_output)
}

#' @title subtract umbrella SPACe from gcomp (raw) grid
#'
#' @description Calculates the gravity components (100percent values) for a given DEM using a nested approach
#'
#' @param g_grid grid for which to calculate gravity components. (generally derived from DEM of region in combination with hydrological modeling mesh/grid).
#' @param senloc coordinates of gravity sensor location. data.frame with x,y,z columns.
#' @param g_discr discretization of g_grid in differences in x,y,and z-direction. they have to be uniform steps!
#' @param edge if grid to calulate is located on a domain edge (top or bottom). for top set edge='first', for bottom set edge='last'.
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export

subtract_umbrellaSpace = function(grid_input, u_space){
	# 


return(grid_output)
}

#' @title subtract everything out of certain radius from gcomp (raw) grid
#'
#' @description 
#'
#' @param grid_input grid from where to limit gcomp to a certain radius. Has to be data.frame with columns (x,y,z).
#' @param gloc location of gravity sensor (x,y,z).
#' @param max_rad radius of exclusion in [meters].
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export

subtract_maxrad = function(grid_input, gloc, max_rad){
	# limit maximum radius of gravity grid
	# if max_rad = NA, no limitation will be realized
	grid_output = mutate(grid_input, distance_rad = sqrt((x-gloc$x)^2+(y-gloc$y)^2+(z-gloc$z)^2)) %>%
			mutate(excluse = ifelse(distance_rad > max_rad, TRUE, FALSE)) %>%
			# delete all rows with exclude == TRUE
			filter(excluse == FALSE) %>%
			select(x,y,z,gcomp)
# return result
return(grid_output)
}

#' @title Calulate gravity signals 
#'
#' @description After calculating the gravity components (gcomp()), this function can be used to get real values using soil moisture data (theta). The output of the gravity signal is possible in different dimensions. 
#'
#' @param gcompfiles vector containing the names (complete paths or relative) of the gravity component files to be used for calculations. These files are generated using gcomp().
#' @param theta data.frame with column structure $time, $timestep, $theta (value), $layer
#' @param output Defines the output to be a singel value, layer or grid (default is value)
#' ...
#' @details Usually one needs relative gravity signals, therefore the inputed theta timeseries should be acutally delta theta vaules per timestep. The
#' other option is calculating the differences in values per timestep of the output of this function.
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing 
#' @export
#' 
gsignal_grids <- function(gcompfile,SMdif){

SMdif_mod = SMdif
SMdif_mod$x = SMdif_mod$x + min(gcompfile$x)
SMdif_mod$y = SMdif_mod$y + min(gcompfile$y)

gcomp_sort = group_by(gcompfile, x,y)
SMdif_sort = group_by(SMdif_mod, x,y)
SMdif_comp = data.frame(Timestep = SMdif_sort$Timestep, SMobs=SMdif_sort$value, gcomp=gcomp_sort$gcomp)
SMdif_comp$gsignal = SMdif_comp$SMobs * SMdif_comp$gcomp / 100
#SMdif_comp = data.frame(Timestep = SMdif_sort$Timestep, difValue=SMdif_sort$DifValue, gcomp=gcomp_sort$gcomp)
#SMdif_comp$gsignal = SMdif_comp$difValue * SMdif_comp$gcomp / 100

#SMtest =  group_by(SMdif_mod, x,y) %>%
#t= data.frame(ti= SMtest$Timestep,tt = SMtest$value*SMtest$DifValue/100)
#tgr = group_by(t,ti) %>%
	#summarize(sum(tt))
#gsignals_allTS = data.frame(Timestep= SMdif_sort$Timestep, gsig_allTS = SMdif_sort$DifValue * SMdif_sort$gcomp/100)
gsignals = group_by(SMdif_comp, Timestep) %>%
		summarize(gmod = sum(gsignal))

return(gsignals)
}

### end gravity signal ###

#' @title Calulate gravity signals 2D
#'
#' @description After calculating the gravity components (gcomp()), this function can be used to get real values using soil moisture data (theta). The output of the gravity signal is possible in different dimensions. 
#'
#' @param gcompfiles vector containing the names (complete paths or relative) of the gravity component files to be used for calculations. These files are generated using gcomp().
#' @param theta data.frame with column structure $time, $timestep, $theta (value), $layer
#' @param output Defines the output to be a singel value, layer or grid (default is value)
#' ...
#' @details Usually one needs relative gravity signals, therefore the inputed theta timeseries should be acutally delta theta vaules per timestep. The
#' other option is calculating the differences in values per timestep of the output of this function.
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing 
#' @export
#' 
gsignal_grids_2d <- function(gcompfile,SMinput){

SMdif_mod = SMinput
#SMdif_mod$xDEM = SMdif_mod$xDEM-1
#SMdif_mod$x = SMdif_mod$x + min(gcompfile$x)
#SMdif_mod$y = SMdif_mod$y + min(gcompfile$y)

#gcomp_sort = group_by(gcompfile, x,zgrid)
#SMdif_sort = group_by(SMdif_mod, xDEM,Depth)
gcomp_sort = arrange(gcompfile, x,zgrid)
SMdif_sort = arrange(SMdif_mod, xDEM,Depth)
SMdif_comp = data.frame(Timestep = SMdif_sort$Timestep, SMobs=SMdif_sort$value, gcomp=gcomp_sort$gcomp)
SMdif_comp$gsignal = SMdif_comp$SMobs * SMdif_comp$gcomp
#SMdif_comp = data.frame(Timestep = SMdif_sort$Timestep, difValue=SMdif_sort$DifValue, gcomp=gcomp_sort$gcomp)
#SMdif_comp$gsignal = SMdif_comp$difValue * SMdif_comp$gcomp / 100

#SMtest =  group_by(SMdif_mod, x,y) %>%
#t= data.frame(ti= SMtest$Timestep,tt = SMtest$value*SMtest$DifValue/100)
#tgr = group_by(t,ti) %>%
	#summarize(sum(tt))
#gsignals_allTS = data.frame(Timestep= SMdif_sort$Timestep, gsig_allTS = SMdif_sort$DifValue * SMdif_sort$gcomp/100)
gsignals = group_by(SMdif_comp, Timestep) %>%
		summarize(gmod = sum(gsignal))

return(gsignals)
}

