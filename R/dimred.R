## This module encapsules all functions regarding dimensionality reduction

#' Function to plot components from PCA (Principal Component Analysis)
#' 
#' @TODO Figure out how to print a single legend with class labels
#'       Currently the legend is positioned for every plot in the 
#'       grid.
#' 
#' @description The function returns a pairs plot of a user defined
#' number of components, that can also be colored by a label if available.
#' Be aware that the function also centers and scales the data before 
#' computing the principal component analysis.
#' 
#' \itemize{
#'   \item{"Lower Diagonal"}{The lower diagonal shows the data points 
#'                           projected to the respective two dimensional 
#'                           component space.}
#'   \item{"Upper Diagonal"}{The upper diagonal shows a contour plot of the data
#'                           projected the respective two dimensional component space.}
#'   \item{"Diagonal"}{The diagonal shows a density plot of the projection of variables
#'                     to a principal component axis.}
#'  }
#' The lower diagonal fo the pais plot shows the projected data points to the
#' respective component space.
#'
#' @export
#'
#' @param df The principal component matrix variable by component (data.frame)
#' @param label Vector by which the data points should be labeled (numeric, character or factor)
#' @param ncomp The number of components theat should be plotted against each other (integer)
#'
#' @import GGally
#' @import ggplot2
#'
#' @return A pairs plot (ggplot2) of principal components
#' 
#' @examples 
#' \dontrun{
#' data(fiveClass)
#' Y <- fiveClass$Class
#' X <- fiveClass[, 2:ncol(fiveClass)]
#' 
#' ## Default Pairs Plot of the First 3 Principal Components
#' dimred.plot_pca(X)
#' 
## Default Pairs Plot of the First 3 Principal Components Labeled by Class Vector
#' dimred.plot_pca(X, Y)
#' 
#' ## Varying the Number of Components
#' dimred.plot_pca(X, Y, ncomp=5)
#' }
dimred.plot_pca <- function(df, label=NULL, ncomp=3) {
  
  stopifnot(ncomp >= 2)
  
  color=NULL
  df_projected <- prcomp(df,
                        center = TRUE,
                        scale = TRUE) 
  components <- data.frame(df_projected$x[,1:ncomp])
  stopifnot(ncomp <= ncol(components))
  
  fig <- ggpairs(data=components,
                 lower=list(continuous=wrap("points", alpha= .5)),
                 upper=list(continuous=wrap("density", alpha = .5)),
                 diag=list(continuous=wrap("densityDiag", alpha = .5)))
  
  if(!is.null(label)) {
    if(!is.factor(label)) {
      label <- as.factor(label) 
    }
    
    components <- cbind(label, components)
    colnames(components)[1] <- "Label"
    fig <- ggpairs(data=components, aes_string(color="Label"),
                   columns=2:ncol(components),
                   lower=list(continuous=wrap("points", alpha= .5)),
                   upper=list(continuous=wrap("density", alpha = .5)),
                   diag=list(continuous=wrap("densityDiag", alpha = .5)))
  }
  return(fig)
}

#' Function to plot components from MCA (Multiple Correspondence Analysis)
#'
#' @return Not implemented yet!
#' @export

dimred.plot_mca <- function() {
  stop("Not implemented yet!")
}

#' Function to plot components from PLS (Partial Least Squares)
#'
#' @return Not implemented yet!
#' @export
dimred.plot_pls <- function() {
  stop("Not implemented yet!")
}

#' Function to plot components from t-SNE (t-Stochastic Neighbor Embedding)
#'
#' @return Not implemented yet!
#' @export
dimred.plot_tsne <- function() {
  stop("Not implemented yet!")
}
