#' @title Surface grid
#'
#' @description Generates a 2d grid of the surface of topography at the observatory building.
#' This is done on basis of a DEM or assuming a flat surface.
#'
#' @param DEM Character string, containing the name of the DEM-file to use.
#' If left empty, a flat topography with value 0 (zero) will be assumed.
#' @param grid_domain_x test
#' @param input_dir Character string, specifying the path of the directory where the DEM-file is located.
#' @param output_dir Character string, specifying the path of the directory where output should be stored.
#' 
#' @return Returns a data.frame, which holds topographical information, including x,y coordinates,
#' from the area of the SG building.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

DEM_to_grid3d = function(
            DEM,
            grid_domain_x,
            grid_domain_y,
            grid_discr,
            depth_split,
            loc_z = NA,
            input_dir,
            output_dir,
            ...
){
    # grid_domain_x = Building_x
    # grid_domain_y = Building_y
    # grid_domain_x = sprinklingArea_x
    # grid_domain_y = sprinklingArea_y
    # DEM = DEM_input_file
    # grid_discr = grid3d_discr
    # input_dir = dir_input
    # DEM = DEM_input_file
    # grid_domain_x = grid_x
    # grid_domain_y = grid_y
    # grid_discr = grid3d_discr
    # input_dir = dir_input_DEM
    # output_dir = dir_input_DEM
    # DEM = DEM_input_file
    # grid_domain_x = grid_x
    # grid_domain_y = grid_y
    # grid_discr = grid_discr[1,]
    # depth_split = grid_depth
    # input_dir = dir_input_DEM
    # output_dir = dir_input_DEM
    # DEM = DEM_file[1]
    # DEM_dir = dir_DEM
    # radius_inner = 0
    # radius_outer = 100
    # grid_domain_x = c(iGrav_locs$x - (radius_outer), iGrav_locs$x + (radius_outer)) # min, max
    # grid_domain_y = c(iGrav_locs$y - (radius_outer), iGrav_locs$y + (radius_outer)) # min, max
    # SGloc = iGrav_locs
    # grid_discr = grid3d_discrezitation[1,]
    # depth_split = grid3d_vertDepth
    # input_dir = dir_DEM
    # output_dir = dir_DEM

    ## check if DEM is non-empty, if not supplied: construct grid
    if(DEM == ""){
        # create raster of new grid extent
        dem_grid_domain = raster::raster(xmn = min(grid_domain_x),
                                 xmx = max(grid_domain_x),
                                 ymn = min(grid_domain_y),
                                 ymx = max(grid_domain_y),
                                 nrows = length(seq(min(grid_domain_x),max(grid_domain_x),by=grid_discr$x)),
                                 ncols = length(seq(min(grid_domain_y),max(grid_domain_y),by=grid_discr$y)),
                                 vals = loc_z
                                 # resolution = .5
                                 )
    }else{
        # read DEM as raster
        dem_raster = raster::raster(paste0(input_dir,DEM))
        ## check visually
        # read_dem(input_dir, DEM)
        # plot_dem(dem, dem.info, data.frame(x=max(grid_domain_x),y=max(grid_domain_y),name="dummy"))

        # create raster of new grid extent
        grid_domain_new = raster::raster(xmn = min(grid_domain_x),
                                 xmx = max(grid_domain_x),
                                 ymn = min(grid_domain_y),
                                 ymx = max(grid_domain_y),
                                 nrows = length(seq(min(grid_domain_x),max(grid_domain_x),by=grid_discr$x)),
                                 ncols = length(seq(min(grid_domain_y),max(grid_domain_y),by=grid_discr$y))
                                 # resolution = .5
                                 )
        # if(cropped){
        #     # crop raster
        #     # if not done, getting a high discr,
        #     # from a big DEM would take forever!
        #     # format: extent(row1, row2, col1, col2)
        #     extent_around_coords = raster::extent(
        #         min(grid_domain_x),
        #         max(grid_domain_x),
        #        min(grid_domain_y) ,
        #        max(grid_domain_y))
        #     message("cropping DEM raster to smaller extent")
        #     dem_raster = raster::crop(dem_raster, extent_around_coords, snap = "out")
        # }
        # resampel DEM no new extent and resolution
        dem_grid_domain = raster::resample(dem_raster, grid_domain_new)
        
        # raster::writeRaster(dem_grid_domain, filename = paste0(output_dir, "dem_grid"), format="ascii", NAflag=-9999, overwrite=T)
        # 
        # # reload grid for verification
        # # in order to create the surface grid
        # ## ! here has to be found a more straight forward way,
        # # going directly from dem_grid_domain to surface_grid
        # # WITHOUT writing a physical file !!
        # dem_grid_file = "dem_grid.asc" 
        # read_dem(output_dir, dem_grid_file)
        # # # plot DEM 
        # # plot_dem(dem, dem.info, locations)
        # # convert DEM grid to data.frame
        # surface_grid = convert_demtodf(dem, dem.info)
    }
    ## create vertical domain of constructed or resampled grid
    # convert raster to data.frame
    # print(str(dem_grid_domain))
    # print(class(dem_grid_domain))
    dem_grid = raster::as.data.frame(dem_grid_domain, xy = T, long = T)
    # find out number of vertical layers
    z_layers = seq(max(depth_split), min(depth_split), by = -1 * as.numeric(grid_discr$z))
    num_z_layers = length(z_layers)
    # construct vertical grid
    grid3d = data.frame()
    for(i in 1:num_z_layers){
        # create grid for ONE layer
        grid_z = data.frame(x = dem_grid$x,
                            y = dem_grid$y,
                            z = (dem_grid$value + z_layers[i]),
                            Depth = z_layers[i],
                            layer = i
                            )
        # combine grids
        grid3d = rbind(grid3d, grid_z)
    }
    # grid3d = cbind(grid.xyz[,1:2],
    #                z=surface[,3] + grid.xyz$z,
    #                Depth=grid.xyz$z, 
    #                layer=rep(seq(1,length(grid.z), by=1),each=(length(unique(grid.xyz$x)) * length(unique(grid.xyz$y))))
    #                )

    ## return data: grid in 3d as data.frame
    return(grid3d)
}

