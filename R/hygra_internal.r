##########################
### internal functions ###
##########################

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

### end discretization ###

### layers ###
layers <- function(C,R,Z,L,dz){
  
  lytop=array(0, dim=c(R,C,L))
  lybot=array(0, dim=c(R,C,L))
  lymid=array(0, dim=c(R,C,L))
  lytop[,,1]=as.matrix(Z)
  lybot[,,1]=as.matrix(Z-dz)
  lymid[,,1]=as.matrix((lytop[,,1]+lybot[,,1])/2)
  for (c1 in 1:(L-1)){
    lytop[,,c1+1]=as.matrix(lytop[,,c1]-dz)
    lybot[,,c1+1]=as.matrix(lybot[,,c1]-dz)
    lymid[,,c1+1]=as.matrix((lytop[,,c1+1]+lybot[,,c1+1])/2)
  }
  
  
  #making variables globally accessable
  lytop <<- lytop; lybot <<- lybot; lymid <<- lymid
  #res <- list()
  #print("lytop, lybot, lymid")
  #return(print("layer parameters calculated"))
}

### end layers ###



### forsberg gravity component ###
forsberg <- function(gama,w,r,c,xl,xr,yl,yr,Zint,Zend,xs,ys,zs,rho){

  #distance mass to SG
  x=0
  x[1]=xl[c]-xs
  x[2]=xr[c]-xs
  y=0
  y[1]=yl[r]-ys
  y[2]=yr[r]-ys
  z=0
  z[1]=Zint-zs
  z[2]=Zend-zs

  sum=0
  for (i in 1:2){
    for (ii in 1:2){
      for (iii in 1:2){
        rf=sqrt(x[i]^2+y[ii]^2+z[iii]^2)
        sum=sum+(-1)^(i+ii+iii)*(x[i]*log(y[ii]+rf)+y[ii]*log(x[i]+rf)-z[iii]*atan(x[i]*y[ii]/z[iii]/rf))
      } #end iii
    } #end ii
  } #end i
  d_forsberg=-w*gama*rho*sum #in µGal, if w = 1e8 
  return(d_forsberg)

}

### end forsberg gravity component ###
### forsberg RAW gravity component ###
forsberg_raw <- function(gama,w,xl,xr,yl,yr,Zint,Zend,xs,ys,zs,rho){

  #distance mass to SG
  x=0
  x[1]=xl-xs
  x[2]=xr-xs
  y=0
  y[1]=yl-ys
  y[2]=yr-ys
  z=0
  z[1]=Zint-zs
  z[2]=Zend-zs

  sum=0
  for (i in 1:2){
    for (ii in 1:2){
      for (iii in 1:2){
        rf=sqrt(x[i]^2+y[ii]^2+z[iii]^2)
        sum=sum+(-1)^(i+ii+iii)*(x[i]*log(y[ii]+rf)+y[ii]*log(x[i]+rf)-z[iii]*atan(x[i]*y[ii]/z[iii]/rf))
      } #end iii
    } #end ii
  } #end i
  d_forsberg=-w*gama*rho*sum #in µGal, if w = 1e8 
  return(d_forsberg)

}

### end forsberg RAW gravity component ###

### macmillan gravity component ###
macmillan <- function(gama,r,c,xp,yp,zp,xs,ys,zs,dx,dy,dz,rad,w,rho){
  
  # calculations for distances
  alfa=2*dx^2-dy^2-dz^2
  beta=-dx^2+2*dy^2-dz^2
  ome=-dx^2-dy^2+2*dz^2
  abg=alfa*(xp[c]-xs)^2+beta*(yp[r]-ys)^2+ome*(zp-zs)^2
  # 3 different macmillan terms
  tm1=-((zp-zs)/rad^3)
  tm2=-5/24*(zp-zs)*abg/rad^7
  tm3=ome/12*(zp-zs)/rad^5 ##12!?! benjamin hatte hier mal eine 24 stehen...warum??
  # multiply together for final result for this layer, spacial extent R&C and for one step in time
  d_macmillan=w*gama*rho*dx*dy*dz*(tm1+tm2+tm3) #in µGal, if w = 1e8 #NEGATIVE SIGN REMOVED!!
  return(d_macmillan)
  
}

### end macmillan gravity component ###
### macmillan gravity RAW component ###
macmillan_raw <- function(gama,xp,yp,zp,xs,ys,zs,dx,dy,dz,rad,w,rho){
  
  # calculations for distances
  alfa=2*dx^2-dy^2-dz^2
  beta=-dx^2+2*dy^2-dz^2
  ome=-dx^2-dy^2+2*dz^2
  abg=alfa*(xp-xs)^2+beta*(yp-ys)^2+ome*(zp-zs)^2
  # 3 different macmillan terms
  tm1=-((zp-zs)/rad^3)
  tm2=-5/24*(zp-zs)*abg/rad^7
  tm3=ome/12*(zp-zs)/rad^5 ##12!?! benjamin hatte hier mal eine 24 stehen...warum??
  # multiply together for final result for this layer, spacial extent R&C and for one step in time
  d_macmillan=w*gama*rho*dx*dy*dz*(tm1+tm2+tm3) #in µGal, if w = 1e8 #NEGATIVE SIGN REMOVED!!
  return(d_macmillan)
  
}

### end macmillan RAW gravity component ###

### pointmass gravity component ###
pointmass <- function(gama,zp,zs,dx,dy,dz,rad,w,rho){
  
  d_pointmass=-w*gama*rho*dx*dy*dz*(zp-zs)/rad^3 #in µGal, if w = 1e8 
  return(d_pointmass)
  
}

### end pointmass gravity component ###

#' @title Convert zoo-object to data frame
#'
#' @description Converting zoo-object to data frame with proper indexing of time column and column names
#'
#' @param value must be a zoo-object
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
# ' @examples examples MISSING
#' 
zootodf <- function(value) {
    df <- as.data.frame(value)
    df$time <- index(value) #create a Date column
    rownames(df) <- NULL #so row names not filled with dates
    df <- df[,c(ncol(df), 1:(ncol(df)-1))] #reorder columns so Date first
    return(df)
}
### end convert too to data.frame ###

#' @title Split timeseries
#'
#' @description Splits a timeseries into differently long pieces of the data, according to the split-indices used.
#' By default it is intended to split at NA-values. Afterwards splits can be sorted out which don't meet a specified length.
#'
#' @param data.in Dataset to be split (must be zoo-object)
#' @param n indicates minimum length of splits
#' ...
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 

split_ts <- function(data.in,n){
  ts.cut=list(NA); i.ts=1
  ind=which(is.na(data.in)==T)
  ind=c(0,ind,length(data.in)+1)
  for(i in 1:length(ind)){
    if(is.na(data.in[ind[i]+1])==T) next #goes to next index to avoid long NA-sequences
    if(ind[i]==ind[length(ind)]) break #terminates for-loop at end of indices
    data.cut=data.in[(ind[i]+1):(ind[i+1]-1)]
    if(length(data.cut)<n)next #cut data-splits shorter than n  
    ts.cut[[i.ts]]=data.cut
    i.ts=i.ts+1
  }
  plot(data.in, main="split up timeseries", ylab=as.character(substitute(data.in)), xlab="")
  for(k in 1:length(ts.cut)){
    lines(ts.cut[[k]], col=rainbow(length(ts.cut))[k])
  }
  return(ts.cut)  
}


#' @title Convert pf-values to suction potential (cm water column)
#'
#' @description 
#'
#' @param PF input PF-values. List, vector.
#' ...
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing
#' 
PFtoSuc = function(PF){
	suction = round(10^PF,2)
	return(suction)
}


#' @title find number of decimal places
#'
#' @description find number of decimal places
#' 
#' @param x numerical input.
#' @details missing
#' @references Marvin Reich (2017), mreich@@gfz-potsdam.de
#' @examples missing

decimalplaces = function(x) {
    if ((x %% 1) != 0) {
        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}
