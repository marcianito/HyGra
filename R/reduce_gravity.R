#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2019), mreich@@posteo.de
#' @export
#' @examples missing

reduce_gravity = function(
        gravity_obs,
        tides = NA,
        atmo = NA,
        polar = NA,
        drift = NA,
        ghe = NA,
        ntol = NA,
        timePeriod_start,
        timePeriod_end,
        plotting = FALSE
){
    # tides
    if(is.na(tides) == FALSE){
        cor_tides = tides
    }else{
        cor_tides = data.frame(
                        datetime = gravity_obs$datetime,
                        tides = 0
                        )
    }
    # atmo
    if(is.na(atmo) == FALSE){
        cor_atmo = atmo
    }else{
        cor_atmo = data.frame(
                        datetime = gravity_obs$datetime,
                        atmo = 0
                        )
    }
    # polar
    if(is.na(polar) == FALSE){
        cor_polar = polar
    }else{
        cor_polar = data.frame(
                        datetime = gravity_obs$datetime,
                        polar = 0
                        )
    }
    # drift
    if(is.na(drift) == FALSE){
        cor_drift = drift
    }else{
        cor_drift = data.frame(
                        datetime = gravity_obs$datetime,
                        drift = 0
                        )
    }
    # ghe
    if(is.na(ghe) == FALSE){
        cor_ghe = ghe
    }else{
        cor_ghe = data.frame(
                        datetime = gravity_obs$datetime,
                        ghe = 0
                        )
    }
    # ntol
    if(is.na(ntol) == FALSE){
        cor_ntol = ntol
    }else{
        cor_ntol = data.frame(
                        datetime = gravity_obs$datetime,
                        ntol = 0
                        )
    }
    # re-format time period, if necessary
    if(!is(timePeriod_start, "POSIXct")){
        timePeriod_start = as.POSIXct(timePeriod_start)
    }
    if(!is(timePeriod_end, "POSIXct")){
        timePeriod_end = as.POSIXct(timePeriod_end)
    }
    # reduce observed gravity signal
    gravity_residual = gravity_obs %>%
        dplyr::left_join(cor_tides) %>%
        dplyr::left_join(cor_atmo) %>%
        dplyr::left_join(cor_polar) %>%
        dplyr::left_join(cor_drift) %>%
        dplyr::left_join(cor_ghe) %>%
        dplyr::left_join(cor_ntol) %>%
        dplyr::filter(datetime > timePeriod_start) %>%
        dplyr::filter(datetime < timePeriod_end) %>%
        dplyr::mutate(value = gravity_obs - tides - atmo - polar - drift - ghe - ntol)
    if(plotting){
        gravity_plot = gravity_residual
        colnames(gravity_plot)[9] = "gravity_residuals"
        gravity_plot = reshape2::melt(gravity_plot, id.vars = "datetime")
        gravity_plot.pl = ggplot(gravity_plot, aes(x = datetime, y = value, color = variable)) + 
            geom_line() + 
            ylab("Gravity [nm/s²]") + xlab("")
        print(gravity_plot.pl)
    }
    # reduce output columns
    gravity_residual = gravity_residual %>%
        dplyr::select(datetime, value)
    # return time series
    return(gravity_residual)
}

