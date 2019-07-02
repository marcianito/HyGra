########################################################
### example script: setting up everything for calculation of local hydrological gravity component responses ###
########################################################

####################
## Instructions
####################

##
##
##

####################
## load libraries
# library(HyGra)
library(devtools)
load_all("/home/mreich/Dokumente/written/ResearchRepos/HyGra")
library(zoo); Sys.setenv(TZ = "GMT")
library(xts)
library(dplyr)
library(raster)
library(reshape2)
library(ggplot2)
library(viridis)
library(gstat)
library(ptinpoly)
library(plot3D)
####################

#########################################
## SETUP
####################
## Output and input settings
# Directory
# path should be absolute
# (if not, it will be relative to the packages library path)
dir_input = "~/Dokumente/GravityProcessing/Local_hydrology/testdata/"
dir_output = "~/Dokumente/GravityProcessing/Local_hydrology/testdata/"

## Gravimeter location
# in [m]
# relativ, local coordinate sytem
SG_x = 3
SG_y = 3
SG_Z = 0
SG_SensorHeight = 0
# UTM coordinate system
SG_x = 4564041.87 
SG_y = 5445662.88 
SG_Z = 606.471
SG_SensorHeight = 1.05 

## Model domain
# in [m]
# local grid or UTM, depending on the coordinates of the SG !
modelDomain_x = c(0, 6) # min, max
modelDomain_y = c(0, 6) # min, max
# grid3d_depth = c(-3, 0) # min, max
# UTM
modelDomain_x = c(SG_x - 3, SG_x + 3) # min, max
modelDomain_y = c(SG_y - 3, SG_y + 3) # min, max

## Model discretization
# [m]
# single grid
grid3d_discr = data.frame(x = .1, y = .1, z = .1)
# nested grid
# one discretization per nested grid (=3)
grid3d_discr = data.frame(x = c(0.5, 1, 10),
                          y = c(0.5, 1, 10),
                          z = c(0.1, 1, 1)
                          )
# vertical extension of grids
# this is equal for all grids
# in [m], negativ pointing downwards
grid3d_depth = c(-2, 0) # min, max

# for nested grids
# boundarys between grids
# [m]
radi_limits = c(0, 50, 300, 1000) 

## building parameters
## Parameters for foundation of building
# these include baseplate, walls, SG pillar(s)
# please use same units as in DEM and model domain
Building_walls_x = .6 # extension
Building_walls_y = .6 # extension 
Building_walls_z = 1.5 # extension
# local grid
Building_baseplate_z = c(-.5, 0) # min, max
# UTM
Building_baseplate_z = c(SG_Z - .5, SG_Z) # min, max

# Round pillar
# [m]
thres_radius = 0.5
thres_depth = -1.2
# square pillar
# [m]
# local grid
Building_SGpillar_x = c(2, 4) # min, max
Building_SGpillar_y = c(2, 4) # min, max
Building_SGpillar_z = c(-1, 0) # min, max
# UTM
Building_SGpillar_x = c(SG_x - 1, SG_x + 1) # min, max
Building_SGpillar_y = c(SG_y - 1, SG_y + 1) # min, max
Building_SGpillar_z = c(SG_Z - 1, SG_Z) # min, max

## Input files
## DEM input file
# if left empty, a flat topographie will be assumed
# single grid
DEM_input_file = ""
DEM_input_file = "WE_UP_TO_300m_05m.asc"
# nested grid
DEM_input_files = c("","","")
DEM_input_files = c("ju_dem_up_to_100m_res01m.acs",
             "ju_dem_up_to_300m_res05m.acs",
             "ju_dem_up_to_10km_res10m.acs"
             )

## Soil moisture data time series (observed or modelled)
# fomrat in .rData, .csv or .tsf
soilMoisture_input_file = "SMdata_TS_1d.rData"
# dummy data
soilMoisture_input_file = "SMdata_dummy.csv"
#
# end setup
#########################################
## CALCULATIONS
####################

#########################################
## Gravimeter location
#########################################
SG_z = SG_Z + SG_SensorHeight
SGloc = data.frame(x=SG_x, y=SG_y, z=SG_z)

#########################################
# Generate 3d gravity component grid : single
#########################################
gravity_component_grid3d_single = create_singleGravityGrid(
            DEM_input_file = DEM_input_file,
            DEM_dir = "",
            SGloc = SGloc,
            grid_discretizations = grid3d_discr,
            grid_vertical = grid3d_depth,
            range_coords_x = modelDomain_x,
            range_coords_y = modelDomain_y,
            # correct_SGpillar =  NA
            correct_SGpillar =  c(thres_radius, thres_depth)
)
# plot
# horizontal
plot_gcomp_grid(grid_input = gravity_component_grid3d_single,
                layer_num = 1,
                grid_discretization = grid3d_discr,
                plane = "horizontal"
)
# vertical
plot_gcomp_grid(grid_input = gravity_component_grid3d_single,
                yloc = SG_y,
                grid_discretization = grid3d_discr,
                plane = "vertical"
)

#########################################
# Generate 3d gravity component grid : nested
#########################################
gravity_grid_nested = create_nestedGravityGrid(
            DEM_input_files = DEM_input_files,
            DEM_dir = dir_input,
            SGloc = SGloc,
            grid_discretizations = grid3d_discr,
            grid_vertical = grid3d_depth,
            threshold_radi = radi_limits,
            correct_SGpillar = c(thres_radius, thres_depth)
            )
# plot
# horizontal
plot_gcomp_grid(grid_input = gravity_grid_nested,
                layer_num = 1,
                grid_discretization = grid3d_discr,
                plane = "horizontal"
)
# vertical
plot_gcomp_grid(grid_input = gravity_grid_nested,
                yloc = SG_y,
                grid_discretization = grid3d_discr,
                plane = "vertical"
)

#########################################
## Correct gravity component grid for SG pillar
#########################################
gravity_component_grid3d_single_pillarCorrected = correct_SGpillar(
            gravity_comp3d = gravity_component_grid3d_single,
            # Pillar_x = NA,
            # Pillar_y = NA,
            # Pillar_z = NA,
            correct_radius = thres_radius,
            correct_depth = thres_depth,
            SG_X = SG_x,
            SG_Y = SG_y #,
            # grid_discretization = NA
)
# plotting grid
plot_gcomp_grid(  grid_input = gravity_component_grid3d_single_pillarCorrected,
                  yloc = SG_y,
                  grid_discretization = grid3d_discr
)
#########################################
## Correct gravity component grid for SG building foundation
#########################################
gravity_component_grid3d_single_BuildingCorrected = correct_SGbuilding_foundation(
            gravity_gcomp = gravity_component_grid3d_single,
            Bdwall_ext_x = Building_walls_x,
            Bdwall_ext_y = Building_walls_y,
            Bdwall_ext_z = Building_walls_z,
            Bdbase_x = Building_x,
            Bdbase_y = Building_y,
            Bdbase_z = Building_baseplate_z,
            Pillar_x = Building_SGpillar_x,
            Pillar_y = Building_SGpillar_y,
            Pillar_z = Building_SGpillar_z,
            grid_discretization = grid3d_discr
)
# plotting grid
plot_gcomp_grid(
                  grid_input = gravity_grid_nested,
                  yloc = SG_y,
                  output_dir = dir_output,
                  grid_discretization = grid3d_discrezitation
)

#########################################
## Extrapolate soil moisture time series data (observed or modelled) to gravity grid domain
#########################################
SMgrid3d_dummy = SoilMoisture_grid3d(
            grid_domain = gravity_component_grid3d_single,
            soilMoisture_input = soilMoisture_input_file,
            grid_discretization = grid3d_discr,
            grid_depth = grid3d_depth,
            input_dir = dir_input
            # , sep = "a", etc..
)

#########################################
## Calculate gravity response (from outside of building)
#########################################
gravity_response = calculate_gravity_response(
            gcomp_grid = gravity_component_grid3d_single,
            mass_input = SMgrid3d_dummy
)
# plot ts
ggplot(gravity_response, aes(x = datetime, y = value)) + 
  geom_line() + geom_point() + 
  xlab("time") + ylab("Gravity [nm/s²]")

#########################################
