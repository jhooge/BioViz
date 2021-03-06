% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimred.R
\name{dimred.plot_pca}
\alias{dimred.plot_pca}
\title{Function to plot components from PCA (Principal Component Analysis)}
\usage{
dimred.plot_pca(df, label = NULL, ncomp = 3)
}
\arguments{
\item{df}{The principal component matrix variable by component (data.frame)}

\item{label}{Vector by which the data points should be labeled (numeric, character or factor)}

\item{ncomp}{The number of components theat should be plotted against each other (integer)}
}
\value{
A pairs plot (ggplot2) of principal components
}
\description{
The function returns a pairs plot of a user defined
number of components, that can also be colored by a label if available.
Be aware that the function also centers and scales the data before 
computing the principal component analysis.

TODO: Figure out how to print a single legend with class labels
      Currently the legend is positioned for every plot in the 
      grid.

\itemize{
  \item{"Lower Diagonal"}{The lower diagonal shows the data points 
                          projected to the respective two dimensional 
                          component space.}
  \item{"Upper Diagonal"}{The upper diagonal shows a contour plot of the data
                          projected the respective two dimensional component space.}
  \item{"Diagonal"}{The diagonal shows a density plot of the projection of variables
                    to a principal component axis.}
 }
The lower diagonal fo the pais plot shows the projected data points to the
respective component space.
}
\examples{
\dontrun{
data(fiveClass)
Y <- fiveClass$Class
X <- fiveClass[, 2:ncol(fiveClass)]

## Default Pairs Plot of the First 3 Principal Components
dimred.plot_pca(X)

dimred.plot_pca(X, Y)

## Varying the Number of Components
dimred.plot_pca(X, Y, ncomp=5)
}
}
