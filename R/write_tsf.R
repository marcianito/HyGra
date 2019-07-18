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

write_tsf = function(
                file,
                header,
                data
                ){
    # write data
    write.table(header, file=file, row.names=F, col.names=F, quote=F)
    write.table(data, file=file, row.names=F, col.names=F, append=T, quote=F)
    # returns nothing
    return(paste0("Data written to: ", output_file))
}
