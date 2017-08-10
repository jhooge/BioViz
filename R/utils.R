## Module for utility functions used internally
## These functions should not be used by the client

#' Function to transform x and y coordinates
#' 
#' @author Sebastian Voss, Adam Skubala
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


#' Function to dictate the exact number of bars in a histogram
#' 
#' @author Sebastian Voss, Adam Skubala
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


#' Function for position dodging in a plot
#' 
#' @author Sebastian Voss, Adam Skubala
#'  
#' @keywords internal
#'
#' @param x a numeric vector with repeated entries (each with frequency n) (Required)
#' @param width width for the distribution of equal entries in 'x' (Default: 0.1)
#'
#' @return a numeric vector of the same length with no equal entries, equal entries
#          of the original vector are distributed evenly over an interval with length
utils.pos_dodge <- function(x, width=0.1){
  x <- as.numeric(x)
  x.freq <- table(x)
  dodge <- seq(-width/2,width/2, length.out=x.freq[1])
  for (i in unique(x)){
    x[x==i] <- x[x==i]+dodge
  }
  return(x)
}


#' Gets the label for a caret model
#' 
#' @param model A trained caret model (train or rfe)
#'
#' @return A human readable model label (character)
utils.get_label <- function(model) {
  label <- model$fit$modelInfo$label
  if (is.null(label)) {
    label <- model$modelInfo$label
  }
  return(label)
}