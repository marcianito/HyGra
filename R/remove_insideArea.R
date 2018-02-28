#' @title removes all DEM points inside a certain area
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

remove_insideArea = function(grid_input, area_poly){
	#library(SDMTools)
	##area_poly = data.frame(x=rem_area[1:2],y=rem_area[3:4])
	#points_in=data.frame(x=grid_input$x, y=grid_input$y)
	#grid_pips = pnt.in.poly(points_in, area_poly)
	##grid_pips = pnt.in.poly(grid_input[,1:2], area_poly)
	#grid_output = dplyr::mutate(grid_pips, z=grid_input$z) %>%
		      #dplyr::filter(pip == 0) %>%
		      #dplyr::select(x,y,z)
	
	grid_output = dplyr::mutate(grid_input, pip = ifelse(
					x > area_poly$x[1] &
					x < area_poly$x[2] &
					y > area_poly$y[3] &
					y < area_poly$y[1]  
					,1,0)) %>%
		      dplyr::filter(pip == 0) %>%
		      dplyr::select(x,y,z,Depth,layer)
return(grid_output)
}
