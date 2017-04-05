## Module for utility functions used internally
## These functions should not be used by the client

#' Function to transform x and y coordinates
#' 
#' @keywords internal
#'
#' @param x x coordinate
#' @param y y coordinate
#' @param pan.row number of panel rows
#' @param pan.col number of panel columns
#' @param prop.fill proportion of how much of the panel should be filled
#'
#' @return A list of transformed x and y variables
utils.trafo_pan <- function(x, y, pan.row, pan.col, prop.fill=0.85){
  
  #shrink to fit into a square panel with range 1...
  xt <- x/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))*(prop.fill)
  yt <- y/(max(y,na.rm=TRUE)-min(y,na.rm=TRUE))*(prop.fill)
  #...and move to right place
  xt <- xt-min(xt, na.rm=TRUE)+pan.col-prop.fill/2
  yt <- yt-min(yt, na.rm=TRUE)-pan.row-prop.fill/2
  
  return(list('xt'=xt, 'yt'=yt))
}


#function to dictate the exact number of bars in a histogram

#' Function to dictate the exact number of bars in a histogram
#' 
#' @keywords internal
#'
#' @param data numeriv vector (Required)
#' @param nbars number of bars (Required)
#'
#' @return A numeric vector of breaks
utils.exBar <- function(data, nbars){
  #create breaks
  return(seq(min(data), max(data), length.out=nbars+1))
}