#' @title Create a gravity component grid out of 3 spatially different discretized grids
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

create_singleGravityGrid = function(
    DEM_input_file,
    DEM_dir,
    SGloc,
    grid_discretizations,
    grid_vertical,
    range_coords_x = NA,
    range_coords_y = NA,
    grid_radius = NA,
    correct_SGpillar = NA 
){
    ## DEBUGGING
    # DEM_input_files = DEM_file
    # DEM_dir = dir_DEM
    # SGloc = iGrav_locs
    # grid_discretizations = grid3d_discrezitation
    # grid_vertical = grid3d_vertDepth
    # threshold_radi = radi_limits
    # correct_SGpillar = SGpillar_data
    ## new
    # DEM_input_file = DEM_input_file
    # DEM_dir = dir_input
    # SGloc = SGloc
    # grid_discretizations = grid3d_discr
    # grid_vertical = grid3d_depth
    # range_coords_x = sprinklingArea_x
    # range_coords_y = sprinklingArea_y
    # correct_SGpillar = c(thres_radius, thres_depth)
    ## 
    if(is.na(grid_radius)){
    gravity_grid = create_gravityGrid(
                DEM_input_file = DEM_input_file,
                dir_input_DEM = DEM_dir,
                SG_coordinates = SGloc,
                grid_discretization = grid_discretizations,
                grid_depth = grid_vertical,
                range_coords_x = range_coords_x,
                range_coords_y = range_coords_y
                )
    }else{
    gravity_grid = create_gravityGrid(
                DEM_input_file = DEM_input_file,
                dir_input_DEM = DEM_dir,
                SG_coordinates = SGloc,
                grid_discretization = grid_discretizations,
                grid_depth = grid_vertical,
                radius_inner = 0,
                radius_outer = grid_radius
                )
    }
    #
    ## remove SG pillar (only for inner grid)
    # cylinder shape
    if(i == 1 & !is.na(correct_SGpillar[1]) & length(correct_SGpillar) == 2){
      # read threshold values for correction
      thres_radius = correct_SGpillar[1]
      thres_depth = correct_SGpillar[2]
      gravity_grid = correct_SGpillar(
                  gravity_comp3d = gravity_grid,
                  correct_radius = thres_radius,
                  correct_depth = thres_depth,
                  SG_X = as.numeric(SGloc$x),
                  SG_Y = as.numeric(SGloc$y)
      )
    # rectangular shape
    }else{
      pillar_x = correct_SGpillar[1:2]
      pillar_y = correct_SGpillar[3:4]
      pillar_z = correct_SGpillar[5:6]
      gravity_grid = correct_SGpillar(
                  gravity_comp3d = gravity_grid,
                  Pillar_x = pillar_x,
                  Pillar_y = pillar_y,
                  Pillar_z = pillar_z,
                  SG_X = as.numeric(SGloc$x),
                  SG_Y = as.numeric(SGloc$y)
      )
    }
    #
    #
    # return nested grid data
    return(gravity_grid)
}

