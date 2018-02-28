#' @title removes all DEM points inside of a certain radius
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

remove_insideRadius = function(grid_input, gloc, max_rad){
	# limit maximum radius of gravity grid
	# if max_rad = NA, no limitation will be realized
	grid_output = dplyr::mutate(grid_input, distance_rad = sqrt((x-gloc$x)^2+(y-gloc$y)^2+(z-gloc$z)^2)) %>%
			dplyr::mutate(excluse = ifelse(distance_rad > max_rad, TRUE, FALSE)) %>%
			# delete all rows with exclude == TRUE
			dplyr::filter(excluse == TRUE) %>%
			dplyr::select(x,y,z,Depth,layer)
# return result
return(grid_output)
}

