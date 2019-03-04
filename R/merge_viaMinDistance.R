#' @title Match data.frames via minimal distance of its (x,y) coordinates
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
#' @examples demo
#' @examples track <- data.frame(xt=runif(10,0,360), yt=rnorm(10,-90, 90))
#' @examples classif <- data.frame(xc=runif(10,0,360), yc=rnorm(10,-90, 90), v1=letters[1:20], v2=1:20)
#' @examples dist.merge(track, classif, 'xt', 'yt', 'xc', 'yc')

dist <- function(x1, y1, x2, y2) {
((x1-x2)^2 + (y1-y2)^2)^0.5
}

merge_viaMinDistance <- function(x, y, xeast, xnorth, yeast, ynorth){
tmp <- t(apply(x[,c(xeast, xnorth)], 1, function(x, y){
dists <- apply(y, 1, function(x, y) dist(x[2],
x[1], y[2], y[1]), x)
cbind(1:nrow(y), dists)[dists == min(dists),,drop=F][1,]
}
, y[,c(yeast, ynorth)]))
tmp <- cbind(x, min.dist=tmp[,2], y[tmp[,1],-match(c(yeast,
ynorth), names(y))])
row.names(tmp) <- NULL
tmp
}

