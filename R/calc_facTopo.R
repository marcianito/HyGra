#' @title Calculate topographic gravity factor
#'
#' @description This is used to transform units: [mm] to [nm/s²]
#'
#' @param grid_extent vector, giving as first argument x and as second y grid extent in [m]
#' @param grid3d_discr data.frame(x = 100, y = 100, z = 0.5)
#' @param grid3d_depth c(-10, 0) # min, max
#' @param dir_input directory of DEM file
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

calc_facTopo = function(
    site_name,
    grid_extent,
    grid3d_discr,
    grid3d_depth,
    topo_flat = F,
    pillar_correction = T,
    dir_input_DEM,
    DEM_fileType = "ascii"

){
    ############
    ## DEBUGGING
    # site_name = "SU"
    # grid_extent = c(getStationParam(site, "TopoFactor_gridExtent_x"), getStationParam(site, "TopoFactor_gridExtent_y"))
    # grid3d_discr = data.frame(x = getStationParam(site, "TopoFactor_gridDiscretization_x"),
    #                                  y = getStationParam(site, "TopoFactor_gridDiscretization_y"),
    #                                  z = getStationParam(site, "TopoFactor_gridDiscretization_z"))
    # grid3d_depth = c(getStationParam(site, "TopoFactor_gridExtent_z"), 0)
    # topo_flat = F
    # dir_input_DEM = paste0(dir_input, "DEM/")
    # DEM_fileType = "ascii"
    # # test for bougher
    # grid_extent = c(50000,50000)
    # grid3d_discr = data.frame(x = 1000,
    #                           y = 1000,
    #                           z = 1)
    # grid3d_depth = c(-3, 0)
    # topo_flat = T
    # 
    ############
    ## Gravimeter location
    # necessarily in UTM coordinate system
    # or relative (then so far wont work with surface DEM)
    # SG_x = getStationParam(site_name, "Longitude")
    # SG_y = getStationParam(site_name, "Latitude")
    SG_x = getStationParam(site_name, "UTM_x")
    SG_y = getStationParam(site_name, "UTM_y")
    SG_Z = getStationParam(site_name, "Height")
    # unifrom for ALL sites !?
    # units: [m]
    # same height as in "our setups"
    SG_SensorHeight = 1.05
    # SG pillar specifications
    # units [m]
    thres_radius = 0.5
    thres_depth = 1.2
    ## Model domain
    # extent in [m]
    grid_ext_x = grid_extent[1]
    grid_ext_y =  grid_extent[2]
    # local grid or UTM, depending on the coordinates of the SG !
    # UTM
    grid_x = c(SG_x - (grid_ext_x / 2), SG_x + (grid_ext_x / 2)) # min, max
    grid_y = c(SG_y - (grid_ext_y / 2), SG_y + (grid_ext_y / 2)) # min, max
    ## Model discretization
    # in [m]
    # grid3d_discr = data.frame(x = 100, y = 100, z = 0.5)
    # grid3d_depth = c(-10, 0) # min, max
    # ## SG pillar
    # # in [m]
    # thres_radius = 1.0
    # thres_depth = 0.7
    ## DEM input file
    # if left empty, a flat topographie will be assumed
    if(topo_flat){
        DEM_input_file = ""
    }else{
      switch(DEM_fileType,
        ascii = { DEM_input_file = paste0(site_name, "_DEM.asc") },
        GeoTiff = { DEM_input_file = paste0(site_name, "_DEM.tif") }
        )
        # DEM_input_file = paste0(getStationParam(site_name, "SRTM_file"), ".asc")
    }
    # DEM_input_file = paste0(site, "_dem_up_to_10km_res10m.acs")
    
    ####################
    ## calculating: surface grid, gravity_component grid and then topography factor
    ## Gravimeter location
    SG_z = SG_Z + SG_SensorHeight
    SGloc = data.frame(x=SG_x, y=SG_y, z=SG_z)
    #
    ## Generate cropped DEM and surface grid
    surface_grid = surface_grid(
                DEM = DEM_input_file,
                grid_domain_x = grid_x,
                grid_domain_y = grid_y,
                grid_discretization = grid3d_discr,
                input_dir = dir_input_DEM,
                output_dir = dir_input_DEM
                # , sep = "a", etc.
    )
    # # plot to check
    # ggplot(data = melt(surface_grid, id = c("x", "y")), aes(x=x, y=y)) + 
    #   geom_raster(aes(fill = value))
    #
    ## Generate 3d gravity component grid 
    gravity_component_grid3d = gravity_comp_grid(
                surface = surface_grid,
                SG_coordinates = SGloc,
                grid_discretization = grid3d_discr,
                grid_depth = grid3d_depth,
                range_coords_x = grid_x,
                range_coords_y = grid_y,
                grid_edges = "regular"
    )
    # # plot to check
    # layer_num = 1
    # ggplot(data = dplyr::filter(melt(gravity_component_grid3d, id = c("x", "y", "z", "Depth", "layer")), layer == layer_num), aes(x=x, y=y)) + 
    #   geom_raster(aes(fill = value))
    # 
    ## Correct gravity component grid for SG pillar 
    if(pillar_correction){
        gravity_component_grid3d = correct_SGpillar(
                    gravity_comp3d = gravity_component_grid3d,
                    correct_radius = thres_radius,
                    correct_depth = thres_depth,
                    SG_X = SG_x,
                    SG_Y = SG_y
        )
    }
    # 
    # plot to check
    #  plot_gcomp_grid(
    #     grid_input = gravity_component_grid3d,
    #     yloc = SG_y,
    #     output_dir = "~/temp/",
    #     grid_discretization = grid3d_discr
    # )
    # 
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
    
    return(TopoFactor)
}

