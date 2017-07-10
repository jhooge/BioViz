## This module encapsules all functions regarding dimensionality reduction

#' Function to plot components from PCA (Principal Component Analysis)
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
#' @param label Nector by which the data points should be labeled (numeric, character or factor)
#' @param ncomp The number of components theat should be plotted against each other (integer)
#'
#' @import GGally
#' @import ggplot2
#'
#' @return A pairs plot (ggplot2) of principal components
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
    fig <- ggpairs(data=components, aes(color=Label),
                   legends=TRUE,
                   columns=2:ncol(components),
                   lower=list(continuous=wrap("points", alpha= .5)),
                   upper=list(continuous=wrap("density", alpha = .5)),
                   diag=list(continuous=wrap("densityDiag", alpha = .5)))
  }
  
  
  fig <- fig + 
    theme_bw() +
    theme(axis.title = element_text(size = 10), 
          axis.text = element_text(size = 8), 
          legend.title = element_blank())
  
  return(fig)
}

pca2dPlot <- function(X, label) {
  Y <- X[, label]
  X <- X[, -which(names(X) == label)]
  
  X_projected <- pcaTransform(X)
  # preProc   = preProcess(X, method=c("center", "scale", "pca"))
  # X_projected = predict(preProc, X)[, 1:3] # PCA projection
  
  projection <- data.frame(Label=Y,
                           PC1=X_projected[, 1], 
                           PC2=X_projected[, 2],
                           PC3=X_projected[, 3])
  
  tools <- c("pan", "resize", 
             "wheel_zoom", "box_zoom", 
             "box_select", "lasso_select", 
             "reset", "save")
  
  cols <- 2:ncol(projection)
  nms <- expand.grid(names(projection)[cols], 
                     rev(names(projection)[cols]), 
                     stringsAsFactors = FALSE)
  
  splom_list <- vector("list", 9)
  for(ii in seq_len(nrow(nms))) {
    splom_list[[ii]] <- figure(width = 200, height = 200, tools = tools,
                               xlab = nms$Var1[ii], ylab = nms$Var2[ii]) %>%
      ly_points(nms$Var1[ii], nms$Var2[ii], 
                data = projection,
                color = Label, 
                size = 6, 
                hover=list(Label),
                legend = FALSE)
  }
  grid_plot(splom_list, ncol = 3, same_axes = TRUE, link_data = TRUE)
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
