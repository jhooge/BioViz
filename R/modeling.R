## This module encapsules all functions regarding modelling and model performance estimation
source("utils.R")

## TODO: This has to be reimplemented, as Sebastian's ROC curve is bugged

#' Function to plot ROC (Receiver Operating Characteristic) Curve
#'
#' @return Not implemented yet!
#' @export
modeling.roc_curve_plot <- function() {
  stop("Not yet implemented")
}

#' Function to plot the validation performance over a parameter domain
#' 
#' @param model A trained caret model (should be of class "train" or "rfe")
#' 
#' @import caret
#' @export
#'
#' @return A ggplot figure
#'
#' @examples
#' \dontrun{
#' library(caret)
#' 
#' X <- iris[1:ncol(iris)-1]
#' Y <- iris$Species
#' tuneGrid <- expand.grid(k=1:20)
#' model <- train(X, Y, method="knn", tuneGrid=tuneGrid)
#' 
#' modeling.parameter_performance_plot(model)
#' }
modeling.parameter_performance_plot <- function(model) {
    label <- utils.get_label(model)
    fig <- NULL
    
    switch(class(model),
           train={fig <- ggplot(model) + ggtitle(label) + theme_bw()},
           rfe={fig <- ggplot(model$fit) + ggtitle(label) + theme_bw()})
    return(fig)
}

