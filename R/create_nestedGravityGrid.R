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

create_nestedGravityGrid = function(
    DEM_input_files,
    DEM_dir,
    SGloc,
    grid_discretizations,
    grid_vertical,
    threshold_radi,
    correct_SGpillar = NA
){
    ## DEBUGGING
    # DEM_input_files = DEM_input_file
    # DEM_dir = dir_input
    # SGloc = SGloc
    # grid_discretizations = grid3d_discr
    # grid_vertical = grid3d_depth
    # threshold_radi = limits_grids
    # correct_SGpillar = NA
    ## 
    # check how much grids (=iteration) are supplied
    n_grid = length(grid_discretizations[,1])
    # create result variable
    nested_gGrid = data.frame()
    #
    # run iteratively over all nested grids
    for(i in 1:n_grid){
      # i = 1
      # construct parameter for this nested grid run
      rad_inner = threshold_radi[i]
      rad_outer = threshold_radi[i + 1]
      grid_discr = grid_discretizations[i,]
      DEM_in = DEM_input_files[i]
      #
      # construct gravity grid
    # Rprof("prof.out") #start a profiling session and put the output to a file "prof.out"
      gravity_grid = create_gravityGrid(
                  DEM_input_file = DEM_in,
                  dir_input_DEM = DEM_dir,
                  SG_coordinates = SGloc,
                  grid_discretization = grid_discr,
                  grid_depth = grid_vertical,
                  radius_inner = rad_inner,
                  radius_outer = rad_outer
                  )
      print(paste0("gravity grid created", i))
    # ggrids = create_gravityGrid()
    # Rprof(NULL) #stop the profiling session
    # summaryRprof(filename = "prof.out")
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
                    SG_Y = as.numeric(SGloc$y),
                    grid_discretization = grid_discr
        )
      }
      # rectangular shape
      if(i == 1 & !is.na(correct_SGpillar[1]) & length(correct_SGpillar) == 6){
        pillar_x = correct_SGpillar[1:2]
        pillar_y = correct_SGpillar[3:4]
        pillar_z = correct_SGpillar[5:6]
        gravity_grid = correct_SGpillar(
                    gravity_comp3d = gravity_grid,
                    Pillar_x = pillar_x,
                    Pillar_y = pillar_y,
                    Pillar_z = pillar_z,
                    SG_X = as.numeric(SGloc$x),
                    SG_Y = as.numeric(SGloc$y),
                    grid_discretization = grid_discr
        )
      }
      #
      # join grid to overall grid-data.frame
      nested_gGrid = rbind(nested_gGrid, gravity_grid)
      #
    }
    #
    # return nested grid data
    return(nested_gGrid)
}

