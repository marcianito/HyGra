###########################################################################
### calculations concerning hydrological parameters and soil parameters ###
###########################################################################

#' @title Calculate K unsaturated
#'
#' @description Get unsaturated hydraulic conductivities for any soil moisture value based on the van Genuchten parameters.
#'
#' @param theta soil moisture value which to calculate K unsaturated
#' @param Ksat saturated hydraulic conductivity
#' @param thr residual moisture content
#' @param ths saturated moisture content
#' @param alpha inverse of air entry point
#' @param n pore size distribution
#' 
#' @details missing 
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

kvadose <- function(theta,vgparam,ksat){
	thr = vgparam[1,2]
	ths = vgparam[2,2]
	alpha = vgparam[3,2]
	n = vgparam[4,2]
	if(theta > 1) theta=theta/100 #make sure theta is in the right unit
	if(theta < thr) print("Entered soil moisture number < residual soil moisture."); stop
	if(theta > ths) print("Entered soil moisture number > saturated soil moisture."); stop
	Se = (theta-thr)/(ths-thr)
	Kun = Ksat*Se^0.5*(1-(1-Se^(n/(n-1)))^(1-1/n))^2
	return(Kun)
}

#' @title Calculate soil moisture content from suction potential 
#'
##' @description test
#'
#' @param suction Suction pressure value (in cm) to be converted to soil moisture
#' @param vgparam data.frame of van Genuchten parameters with the column structure: $thr, $ths, $alpha, $n
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' @export
#' 

theta <- function(suction,vgparam){
	thr = vgparam[1,2]
	ths = vgparam[2,2]
	alpha = vgparam[3,2]
	n = vgparam[4,2]
	theta = thr + (ths-thr)*(1/(1+abs(alpha*suction)^n))^(1-1/n)
	return(theta)
}

#' @title Calculate suction potential from soil moisture content
#'
##' @description test
#'
#' @param theta Soil moisture value (decimal number) to be converted to suction pressure
#' @param vgparam data.frame of van Genuchten parameters with the column structure: $thr, $ths, $alpha, $n 
#' 
#' @details missing 
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 

suction <- function(theta,vgparam){
	thr = vgparam[1,2]
	ths = vgparam[2,2]
	alpha = vgparam[3,2]
	n = vgparam[4,2]
	Se = (theta - thr)/(ths-thr)
	suction = (1/-alpha)*((1/(Se)^(1/(1-1/n)))-1)^1/n
	return(suction)
}

#' @title Calculate van Genuchten parameters 
#'
#' @description From a given retention curve dataset van Genuchten parameters thr, ths, alpha and n are calculated.
#'
#' @param data.in Values for water retention curve. Column structure: $theta, $suction
#' @param use.method string of method name. possible aer VG (vanGenuchten), brook, gardner, kosugi.
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' @export
#' 

soilparam <- function(data.in,use.method="VG"){
# library(HydroMe) #load van Genuchten model
if(use.method == "VG"){ #standard method
	param.estimate = nlsLM(theta ~ SSvgm(suction,thr,ths,alp,nscal,mscal),data=data.in,control=nls.lm.control(maxiter=500))
}
if(use.method == "brook"){
	data.brook = data.frame(x=data.in$suction, y=data.in$theta)
	param.estimate = nlsLM(theta ~ Brook(suction,thr,ths,alp,nscal),data=data.in,control=nls.lm.control(maxiter=200),start=c(thr=Dstart(data.brook)[1],ths=Dstart(data.brook)[2], alp=Dstart(data.brook)[3], nscal=Dstart(data.brook)[4]-1))
}
if(use.method == "gardner"){
	param.estimate = nlsLM(theta ~ SSgardner(suction,thr,ths,alp,nscal),data=data.in,control=nls.lm.control(maxiter=200))
}
if(use.method == "kosugi"){
	param.estimate = nlsLM(theta ~ SSkosugi(suction,thr,ths,alp,nscal),data=data.in,control=nls.lm.control(maxiter=200,options(warn=-1)))
}
thr=summary(param.estimate)$parameters[1,1]
ths=summary(param.estimate)$parameters[2,1]
alpha=summary(param.estimate)$parameters[3,1]
n=summary(param.estimate)$parameters[4,1]
output = data.frame(variable = c("thr","ths","alpha","n"),
		    value = c(round(thr,3),round(ths,3),round(alpha,3),round(n,3)))
return(output)  
  
}


