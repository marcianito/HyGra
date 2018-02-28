#' @title test
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

create_gravityGrid = function(
            DEM_input_file,
            dir_input_DEM,
            SG_coordinates,
            grid_discretization,
            grid_depth,
            range_coords_x = NA,
            range_coords_y = NA,
            radius_inner = NA,
            radius_outer = NA,
            ...
){
    ## DEBUGGING
    # DEM_input_file = DEM_file
    # # DEM_input_file = DEM_in
    # dir_input_DEM = dir_DEM
    # # # SG_coordinates = SGloc
    # SG_coordinates = SG_locs
    # grid_discretization = grid3d_discrezitation
    # # grid_discretization = grid3d_discretization[1,]
    # grid_depth = grid3d_vertDepth
    # # # radius_inner = rad_inner
    # # # radius_outer = rad_outer
    # # radius_inner = 0
    # # radius_outer = 100
    # range_coords_x = grid_domain_x
    # range_coords_y = grid_domain_y
    ## testing UmbrellaEffect package:
            # DEM_input_file = DEM_input_file
            # dir_input_DEM = ""
            # SG_coordinates = SGloc
            # grid_discretization = grid3d_discr
            # grid_depth = grid3d_depth
            # range_coords_x = Building_x
            # range_coords_y = Building_y
            # radius_inner = NA
            # radius_outer = NA
    ##
    if(is.na(radius_inner)){
        # rectangle
        # decide between relative or absolute coordinates
        if(length(range_coords_x) == 1){
        # subtracting / adding  half of complete extent
        grid_x = c(SG_coordinates$x - (range_coords_x / 2), SG_coordinates$x + (range_coords_x / 2)) # min, max
        grid_y = c(SG_coordinates$y - (range_coords_y / 2), SG_coordinates$y + (range_coords_y / 2)) # min, max
        }else{
        # setting absolute grid domain borders
        grid_x = c(min(range_coords_x), max(range_coords_x)) # min, max
        grid_y = c(min(range_coords_y), max(range_coords_y)) # min, max
        }
    }else{
        # cirle
        # in fact, here it is still a rectangle
        # but prepared with extentions for making it a cut-out-circle later on
        grid_x = c(SG_coordinates$x - (radius_outer), SG_coordinates$x + (radius_outer)) # min, max
        grid_y = c(SG_coordinates$y - (radius_outer), SG_coordinates$y + (radius_outer)) # min, max
    }
    # # Generate cropped DEM and surface grid
    # message("creating adjusted surface grid")
    # grid_surface = surface_grid(
    #             DEM = DEM_input_file,
    #             grid_domain_x = grid_x,
    #             grid_domain_y = grid_y,
    #             grid_discretization = grid_discretization,
    #             input_dir = dir_input_DEM,
    #             output_dir = dir_input_DEM
    #             # , sep = "a", etc.
    # )
    # # ggplot(data = melt(grid_surface, id = c("x", "y")), aes(x=x, y=y)) + 
    # #   geom_raster(aes(fill = value))
    # # generate 3d grid
    # message("surface to 3d grid")
    # grid3d = surface_to_grid3d(
    #         surface_grid = grid_surface,
    #         grid_discr = grid_discretization,
    #         depth_split = grid_depth,
    #         Bd_x = grid_x,
    #         Bd_y = grid_y 
    # )
    # generate 3d grid
    message("DEM to 3d grid")
    grid3d = DEM_to_grid3d(
                DEM = DEM_input_file,
                grid_domain_x = grid_x,
                grid_domain_y = grid_y,
                grid_discr = grid_discretization,
                depth_split = grid_depth,
                loc_z = SG_coordinates$z,
                input_dir = dir_input_DEM,
                output_dir = dir_input_DEM
    )

    # plot to check
    # ggplot(data = melt(grid3d, id = c("x", "y")), aes(x=x, y=y)) + 
    #   geom_raster(aes(fill = value))

    # define geometry of grid and apply provided limits
    if(!is.na(radius_inner)){
        message("limiting data to OUTER radius..")
        grid_limitOuterRadius = remove_outsideRadius(grid3d, SG_coordinates, radius_outer)
        message("limiting data to INNER radius..")
        grid3d = remove_insideRadius(grid_limitOuterRadius, SG_coordinates, radius_inner)
        #
        # cutting out inner rectangle area
        # .. leave for implementation later on !!
        # grid_radius50m = remove_insideArea(grid_radius50mOuter, hydrus_domain)
    }

    # generate gravity component grid
    message("calculating gravity components per grid cell")
    gcomp_grid = fill_gcompgrid(
            g_grid = grid3d,
            senloc = SG_coordinates,
            g_discr = grid_discretization,
            edge = "regular"
    )
    # plot to check
    # ggplot(data = melt(gcomp_grid, id = c("x", "y", "z")), aes(x=x, y=y)) + 
    #   geom_raster(aes(fill = value))

    # # take out column "layer", which was needed only for internal calculations
    # gcomp_grid = dplyr::select(gcomp_grid, -layer)

    # round all x,y,z to same decimal places
    # if not, joining might not be complete!
    message("rounding data to necessary precision")
    round_x = decimalplaces(grid_discretization$x)
    round_y = decimalplaces(grid_discretization$y)
    round_z = decimalplaces(grid_discretization$z)
    gcomp_grid$x = round(gcomp_grid$x, round_x)
    gcomp_grid$y = round(gcomp_grid$y, round_y)
    gcomp_grid$z = round(gcomp_grid$z, (round_z + 1))
    gcomp_grid$Depth = round(gcomp_grid$Depth, round_z)

    # return data
    return(gcomp_grid)
}

####################
## TESTING
####################
# 
# ####################
# ## radius model domain to 50 m
# filename = "WE_UP_TO_300m_05m.asc" #r=300m, dxdy=0,5m
# read_dem(dempath,filename)
# #convert from dem (asc) to data.frame
# grid_50m = convert_demtodf(dem, dem.info)
# # limit outer and inner boundaries
# grid_radius50mOuter = remove_outsideRadius(grid_50m, SGloc, 50)
# grid_radius50m = remove_insideArea(grid_radius50mOuter, hydrus_domain)
# #ggplot(grid_radius50m, aes(x=x,y=y)) + geom_tile(aes(fill=z)) + geom_point(data=SGloc, aes(x=x,y=y))
# 
# grid3d_radius50m_overhead = demgrid_to_gcompgrid_Edges(grid_radius50m,grid_discr_r50m, depth_r50m)
# # correct AGAIN for inner and outer radius / area due to interpolation method necessary !!!
# grid3d_radius50mOuter = remove_outsideRadius(grid3d_radius50m_overhead, SGloc, 50,T)
# grid3d_radius50m = remove_insideArea(grid3d_radius50mOuter, hydrus_domain,T)
# #ggplot(filter(grid3d_radius50m, layer==50), aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z))
# # overlap correction with inner area
# library(SDMTools)
# hydrus_r50 = pnt.in.poly(cbind(x=grid3d_radius50m$x, y=grid3d_radius50m$y), hydrus_domain)
# olay = which(hydrus_r50$pip == 1)
# grid3d_radius50m[olay,]=NA
# grid3d_radius50m = filter(grid3d_radius50m, is.na(x)==F | is.na(y)==F | is.na(z)==F)
# save(grid3d_radius50m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid3d_radius50m.rdata")
# ##
# ## ???
# ##
# ## r50,300 & 5000 probably have to be calculated on glic !!!
# gcomp_radius50m = gcomp_raw_Edges(grid3d_radius50m, SGloc, grid_discr_r50m, grid_edge_r50m)
# 
# 
# ####################
# ## radius 50 m to 300 m
# filename = "WE_UP_TO_300m_05m.asc" #r=300m, dxdy=0,5m
# read_dem(dempath,filename)
# # convert from dem (asc) to data.frame
# grid_300m = convert_demtodf(dem, dem.info)
# # limit outer and inner boundaries
# grid_radius300mOuter = remove_outsideRadius(grid_300m, SGloc, 300)
# grid_radius300m = remove_insideRadius(grid_radius300mOuter, SGloc, 50)
# #ggplot(grid_radius300m, aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z))
# save(grid_radius300m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius300m.rdata")
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius300m.rdata")
# 
# ####################
# ## radius 300 m - 5560 m
# filename = "WE_UP_TO_12km_10m.asc" #r=12km, dydy=10,0m
# read_dem(dempath,filename)
# # convert from dem (asc) to data.frame
# grid_5000m = convert_demtodf(dem, dem.info)
# # limit outer and inner boundaries
# grid_radius5000mOuter = remove_outsideRadius(grid_5000m, SGloc, 5560)
# grid_radius5000m = remove_insideRadius(grid_radius5000mOuter, SGloc, 300)
# #ggplot(grid_radius5000m, aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z))
# save(grid_radius5000m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius5000m.rdata")
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius5000m.rdata")
# 
# ## combine ALL grids to check for OVERLAPS and make a visual control !!
# #grid_complete_area = rbind(grid_hydrus_domain, grid_radius50mInner, grid_radius300mInner,grid_radius5000mInner)
# #ggplot(grid_complete_area, aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z)) + geom_point(data=SGloc, aes(x=x,y=y))
# 
# #grid_complete_area = rbind(cbind(grid_hydrus_domain, name="hydrus",col="green"),
#                #cbind(grid_radius50m,name="radius 50m",col="orange"),
#                ##cbind(grid_radius300m, name="radius 300m",col="blue"),
#                ##cbind(grid_radius5000m, name="radius 5560m",col="grey"),
#                #cbind(SGloc, name="SGnew",col="red"))
# #ggplot(grid_complete_area, aes(x=x,y=y, colour=name)) + geom_point()# + geom_point(data=SGloc, aes(x=x,y=y))
# 
# # check if no OVERLAPS of adjacent grids exist !!
# # point in polygon test
# library(SDMTools)
# hydrus_border = data.frame(x=c(range(grid_hydrus_domain$x),range(grid_hydrus_domain$x)),y=c(rep(min(grid_hydrus_domain$y),2),rep(max(grid_hydrus_domain$y),2)))
# hydrus_r50 = pnt.in.poly(cbind(x=grid_radius50m$x, y=grid_radius50m$y), hydrus_border)
# olay = which(hydrus_r50$pip == 1)
# overlay_hydrus_r50 = hydrus_r50[olay,]
# ggplot(overlay_hydrus_r50, aes(x=x,y=y)) + geom_point()
# grid_radius50m[olay,]=NA
# t=which(is.na(grid_radius50m)==T)
# grid_radius50m[t,]
# ggplot(grid_radius50m, aes(x=x,y=y)) + geom_point()
# ggplot(rbind(cbind(grid_radius50m,name="r50"),cbind(grid_hydrus_domain,name="hydrus")), aes(x=x,y=y, colour=name)) + geom_point()
# grid_radius50m = filter(grid_radius50m, is.na(x)==F | is.na(y)==F | is.na(z)==F)
# save(grid_radius50m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius50m.rdata")
# 
# # point in polygon method with a rectangle DOES NOT WORK when comparing two discs: min, max extensd pretend going from disc shape to rectangular
# #r50m_border = data.frame(x=c(range(grid_radius50m$x, na.rm=T),range(grid_radius50m$x, na.rm=T)),y=c(rep(min(grid_radius50m$y, na.rm=T),2),rep(max(grid_radius50m$y, na.rm=T),2)))
# #r50_r300 = pnt.in.poly(cbind(x=grid_radius300m$x, y=grid_radius300m$y), r50m_border)
# 
# ##! check with COMPLETE pip-r50 in r300 !!
# 
# r50_r300 = pnt.in.poly(cbind(x=grid_radius300m$x, y=grid_radius300m$y),cbind(x=grid_radius50m$x, y=grid_radius50m$y) )
# olay = which(r50_r300$pip == 1)
# overlay_r50_r300 = r50_r300[olay,]
# overlay_r50 = grid_radius50m[olay,]
# ggplot(overlay_r50_r300, aes(x=x,y=y)) + geom_point()
# grid_radius300m[olay,]=NA
# ggplot(grid_radius300m, aes(x=x,y=y)) + geom_point()
# ggplot(rbind(cbind(grid_radius50m,name="r50"),cbind(grid_radius300m,name="r300")), aes(x=x,y=y, colour=name)) + geom_point()
# save(grid_radius300m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius300m.rdata")
# 
# r300m_border = data.frame(x=c(range(grid_radius300m$x),range(grid_radius300m$x)),y=c(rep(min(grid_radius300m$y),2),rep(max(grid_radius300m$y),2)))
# r300_r5000 = pnt.in.poly(cbind(x=grid_radius5000m$x, y=grid_radius5000m$y), r300m_border)
# # .. 
# save(grid_radius5000m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius5000m.rdata")
# 
# 
# ###############################################
# # generate one 3d grid with individual discretization for each grid
# # load grid data:
# # SGbuilding domain
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_SGbuilding_domain.rdata")
# # hydrus domain
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_hydrus_domain.rdata")
# # hydrus domain - r50m
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius50m.rdata")
# # r50m - r300m
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius300m.rdata")
# # r300m - r5000m
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid_radius5000m.rdata")
# 
# ## discretizations
# # SGbuilding domain
# grid_discr_SGbuilding_domain = data.frame(x=0.1,y=0.1,z=0.1); depth_SGbuilding_domain = c(0,5)
# grid_edge_SGbuilding_domain = "both"
# # hydrus domain
# grid_discr_hydrus_domain = data.frame(x=0.1,y=0.1,z=0.1); depth_hydrus_domain = c(0,5)
# grid_edge_hydrus_domain = "both"
# # r50m
# grid_discr_r50m = data.frame(x=0.1,y=0.1,z=0.1); depth_r50m = c(0,5)
# grid_edge_r50m = "both"
# # r300m
# grid_discr_r300m = data.frame(x=1,y=1,z=0.1); depth_r300m = c(0,5)
# grid_edge_r300m = "both"
# # r5000m
# # use this one if layer (depth) 0 is desired
# #grid_discr_r5000m = data.frame(x=10,y=10,z=1); depth_r5000m = c(0,5)
# # use this one, starting layers at 1st bounary (here: 1m)
# grid_discr_r5000m = data.frame(x=10,y=10,z=1); depth_r5000m = c(1,5)
# #grid_edge_r5000m = "both"
# 
# ####################
# ## SGbuilding domain
# grid3d_SGbuilding_domain = demgrid_to_gcompgrid_Edges(grid_SGbuilding_domain,grid_discr_SGbuilding_domain, depth_SGbuilding_domain)
# # ggplot(filter(grid3d_SGbuilding_domain, layer==10), aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z))
# save(grid3d_SGbuilding_domain, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid3d_SGbuilding_domain.rdata")
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid3d_SGbuilding_domain.rdata")
# gcomp_SGbuilding_domain = gcomp_raw_Edges(grid3d_SGbuilding_domain, SGloc, grid_discr_SGbuilding_domain, grid_edge_SGbuilding_domain)
# save(gcomp_SGbuilding_domain, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/gcomp_SGbuilding_domain_withoutUmbrella.rdata")
# 
# ####################
# ## hydrus domain
# grid3d_hydrus_domain = demgrid_to_gcompgrid_Edges(grid_hydrus_domain,grid_discr_hydrus_domain, depth_hydrus_domain)
# #ggplot(filter(grid3d_hydrus_domain, layer==10), aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z))
# save(grid3d_hydrus_domain, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid3d_hydrus_domain.rdata")
# load(file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid3d_hydrus_domain.rdata")
# gcomp_hydrus_domain = gcomp_raw_Edges(grid3d_hydrus_domain, SGloc, grid_discr_hydrus_domain, grid_edge_hydrus_domain)
# save(gcomp_hydrus_domain, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/gcomp_hydrus_domain_withoutUmbrella.rdata")
# 
# ####################
# ## r50m
# grid3d_radius50m_overhead = demgrid_to_gcompgrid_Edges(grid_radius50m,grid_discr_r50m, depth_r50m)
# # correct AGAIN for inner and outer radius / area due to interpolation method necessary !!!
# grid3d_radius50mOuter = remove_outsideRadius(grid3d_radius50m_overhead, SGloc, 50,T)
# grid3d_radius50m = remove_insideArea(grid3d_radius50mOuter, hydrus_domain,T)
# #ggplot(filter(grid3d_radius50m, layer==50), aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z))
# # overlap correction with inner area
# library(SDMTools)
# hydrus_r50 = pnt.in.poly(cbind(x=grid3d_radius50m$x, y=grid3d_radius50m$y), hydrus_domain)
# olay = which(hydrus_r50$pip == 1)
# grid3d_radius50m[olay,]=NA
# grid3d_radius50m = filter(grid3d_radius50m, is.na(x)==F | is.na(y)==F | is.na(z)==F)
# save(grid3d_radius50m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid3d_radius50m.rdata")
# ##
# ## ???
# ##
# ## r50,300 & 5000 probably have to be calculated on glic !!!
# gcomp_radius50m = gcomp_raw_Edges(grid3d_radius50m, SGloc, grid_discr_r50m, grid_edge_r50m)
# save(gcomp_radius50m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/gcomp_radius50m_withoutUmbrella.rdata")
# 
# ####################
# ## r300m
# grid3d_radius300m = demgrid_to_gcompgrid_Edges(grid_radius300m,grid_discr_r300m, depth_r300m)
# grid3d_radius300mOuter = remove_outsideRadius(grid3d_radius300m, SGloc, 300,T)
# grid3d_radius300m = remove_insideRadius(grid3d_radius300mOuter, SGloc, 50,T)
# #ggplot(filter(grid3d_radius300m, layer==10), aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z))
# save(grid3d_radius300m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid3d_radius300m.rdata")
# gcomp_radius300m = gcomp_raw_Edges(grid3d_radius300m, SGloc, grid_discr_r300m, grid_edge_r300m)
# save(gcomp_radius300m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/gcomp_radius300m_withoutUmbrella.rdata")
# 
# ####################
# ## r5000m
# grid3d_radius5000m = demgrid_to_gcompgrid(grid_radius5000m,grid_discr_r5000m, depth_r5000m)
# grid3d_radius5000mOuter = remove_outsideRadius(grid3d_radius5000m, SGloc, 5560,T)
# grid3d_radius5000m = remove_insideRadius(grid3d_radius5000mOuter, SGloc, 300,T)
# #ggplot(filter(grid3d_radius5000m, layer==10), aes(x=x,y=y, colour=z)) + geom_tile(aes(fill=z))
# save(grid3d_radius5000m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/grid3d_radius5000m_new.rdata")
# gcomp_radius5000m = gcomp_raw(grid3d_radius5000m, SGloc, grid_discr_r5000m)
# save(gcomp_radius5000m, file="/home/mreich/Dokumente/Wettzell/gravitymodelling/data/HydrusModels/gravity_component_grids/gcomp_radius5000m_withoutUmbrella_new.rdata")
# 
# # combine grids into one
# grid_complete = rbind(grid_hydrus_domain,grid_radius50m,grid_radius300m,grid_radius5000m)
#     return()
# }
# 
