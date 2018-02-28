#' @title Calculate topographic gravity factor
#'
#' @description This is used to transform units: [mm] to [nm/s²]
#'
#' @param gravity_component_grid3d data.frame, supplying the gravity grid.
#' Needs columns x, y, z, Depth, gcomp (gravity component value).
#' 
#' @return Returns the topographic conversion factor which represents
#' the gravity effect resulting from 1 mm of vertical water storage change.
#' This factor is used to convert mm of water from precipitation or evapotranspiration
#' measurements, etc. to gravity units [nm/s²].
#' 
#' @details missing
#' @references Marvin Reich (2018), mreich@@posteo.de
#' @import 
#' @examples missing

facTopo_from_grid = function(
    gravity_component_grid3d,
    ...
){
    ############
    ## DEBUGGING
    # gravity_component_grid3d = 
    ##
    ## aggregating all components for one topography factor:
    # sum horizontally
    # then average vertically
    # this results in one vertically averaged gravity component per radar cell
    TopoFactor = gravity_component_grid3d %>%
        dplyr::group_by(Depth) %>%
        # sum horizontally
        dplyr::summarize(gcomp_z = sum(gcomp, na.rm = T)) %>%
        # average vertically
        # 100% value (filled gravity cubes)
        # convert [m] to [mm] to get units:
        # 1 mm change of water causes HOW MUCH gravity change in [nm / s²]
        dplyr::summarize(0.001 * mean(gcomp_z, na.rm = T))
    # convert to correct data type (NOT data.frame !)
    TopoFactor = as.numeric(TopoFactor)
    # return value 
    return(TopoFactor)
}

