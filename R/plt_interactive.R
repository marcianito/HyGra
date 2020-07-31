#' @title pass a ggplot plot object to ggplotly for interactive plotting
#'
#' @description 
#'
#' @param ggplot_plt ggplot plot object
#' @param set_browser string; set the browser to use. Default is "google-chrome-stable".
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2020), mreich@@posteo.de
#' @import 
#' @examples missing

plt_interactive = function(ggplot_plt, set_browser = "google-chrome-stable"){
    ## load library
    library(plotly)
    # set browser settings
    options(browser = set_browser)
    ## plot in browser with:
    # dynamic ticks
    # rangeslider
    ggplotly(ggplot_plt, dynamicTicks = TRUE) %>%
    rangeslider()
}

