#' @title Interpolate missing time steps of data time series
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

interpTS = function(
    data_in,
    freq_in,
    freq_out,
    aggFunc = "mean",
    data_col_name
){
    # get start and end time stamps
    ts_start = min(data_in$datetime)
    ts_end = max(data_in$datetime)
    # construct time series dates: hourly
    datetime_new = data.frame(datetime = seq(ts_start, ts_end, by=freq_in))
    # select data column
    gw_selected = data_in %>%
        dplyr::select_("datetime", data_col_name)
    colnames(gw_selected)[2] = "value"
    # join old into new data time stamps
    data_new = datetime_new %>%
        dplyr::left_join(gw_selected)
    # approximate (interpolate, linear) missing data 
    data_new$value = na.approx(data_new$value)
    data_new_out = aggTS(data_new, freq_out, aggFunc)
    # return data
    return(data_new_out)
}

