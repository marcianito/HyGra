#' @title plot irrigated water volume per layer and timestep
#'
#' @description test
#'
#' @param test
#' @param test
#' @param test
#' ...
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

plotInfWater <- function(filein, showlegend = T){
  library(viridis)
  data_plot = dplyr::group_by(filein, Timestep, zgrid) %>%
            dplyr::summarize(value = sum(value, na.rm = T))
  # legend !?
  if(showlegend){ leg = "right"
  }else{ leg = "none"}

  data.gg = ggplot(data_plot, aes(x = Timestep, y = zgrid)) + geom_raster(aes(fill = value)) + ylim(2,0) + 
            theme(legend.position=leg) +
            ylab("Depth [m]") + xlab("Timestep [min]") + 
            scale_fill_gradientn(colours = viridis(7))
  return(data.gg) 
}

#' @title check infiltrated water volumen
#'
#' @description test
#'
#' @param test
#' @param test
#' @param test
#' ...
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

checkWaterVolumes <- function(filein, whichCheck, inf_time, cell_volume, infWater_timestep, ts_check=1, plotting=F){
  switch(whichCheck,
         totalwater = {
                        modelwater = dplyr::summarize(filein, modelsum = sum(value, na.rm=T)) 
                        comparison = data.frame(water_model = as.integer(modelwater) * cell_volume,
                                                water_real = sum(seq(1,inf_time)) * infWater_timestep
                                               )
                      },
         waterpertimestep = {
                        modelwater = dplyr::group_by(filein, Timestep) %>%
                                     dplyr::summarize(modelsum = sum(value/Timestep, na.rm=T)) 
                        comparison = data.frame(timestep = modelwater$Timestep,
                                                water_model = modelwater$modelsum * cell_volume,
                                                water_real = infWater_timestep
                                               )
                        if(plotting){ print(ggplot(modelwater, aes(x=Timestep, y=modelsum)) + geom_line()) }
                      },
         tswater = {
                        modelwater = dplyr::filter(filein, Timestep == ts_check) %>%
                                     dplyr::summarize(modelsum = sum(value, na.rm=T)) 
                        comparison = data.frame(water_model = as.integer(modelwater) * cell_volume,
                                                water_real = ts_check * infWater_timestep
                                               )
                      }
  )
  return(comparison)
}
