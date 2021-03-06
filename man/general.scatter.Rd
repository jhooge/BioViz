% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{general.scatter}
\alias{general.scatter}
\title{Function factory to generate a scatter plots}
\usage{
general.scatter(data, by = NULL, fun = "mean", smooth.fun = "auto",
  smooth.formula = y ~ x, density = FALSE, title = NULL, xlab = NULL,
  ylab = NULL, legend.pos = "none")
}
\arguments{
\item{data}{A dataframe with at least three columns (x,y and a facetting variable) (Required)}

\item{by}{Grouping variable (Required)}

\item{fun}{Function for statistical aggregation function (Default: "mean")}

\item{smooth.fun}{Smoothing method (function) to use, eg. "lm", "glm", "gam", "loess", "rlm". (Default: "auto")
For method = "auto" the smoothing method is chosen based on the size of the largest group 
(across all panels). loess is used for less than 1,000 observations; otherwise gam is used with 
formula = y ~ s(x, bs = "cs"). Somewhat anecdotally, loess gives a better appearance, 
but is O(n^2) in memory, so does not work for larger datasets.}

\item{smooth.formula}{Formula to use in smoothing function, eg. y ~ x, y ~ poly(x, 2), y ~ log(x) (Default: y ~ x)}

\item{density}{Display density contour (Default: FALSE)}

\item{title}{Plot Title (Default: NULL)}

\item{xlab}{x-axis label (Default: NULL)}

\item{ylab}{y-axis label (Default: NULL)}

\item{legend.pos}{(Default: "none")}
}
\value{
ggplot2 figure
}
\description{
Depending on the grouping variable 'by', this function factory generates
either an ungrouped (by=NULL) or a grouped scatterplot (by=<grouping_variable>).
This function provides the option for different regression fits, based on an inbuilt
smooth function or a custom model, defined by smooth.formula.
}
\examples{
\dontrun{
x <- rnorm(500)
y <- rnorm(500)+x^2
df <- data.frame(x=x, y=y, 
                 A=sample(c("a1", "a2"), 500, replace = TRUE),
                 B=sample(c("b1", "b2"), 500, replace = TRUE))

## A simple ungrouped scatter plot with a loess fit
general.scatter(df)

## A grouped scatter plot with a loess fit
general.scatter(df, by="A")

## An ungrouped scatter plot and a marker indicating the data variance
general.scatter(df, fun="var")

## Scatterplots with different regression fits
# Linear Fit
general.scatter(df, by="A", smooth.fun="lm")

## Robust Linear Fit
library(MASS)
general.scatter(df, by="A", smooth.fun="rlm")

## User defined polynomial fit
# Degree = 1 resulting in a linear fit
general.scatter(df, smooth.fun="glm", smooth.formula=y ~ poly(x, 1)) 
# Degree = 2 resulting in a quadratic fit
general.scatter(df, smooth.fun="glm", smooth.formula=y ~ poly(x, 2))

## Equivalently for grouped scatterplots
# Degree = 1 resulting in a linear fit
general.scatter(df, by="A", smooth.fun="glm", smooth.formula=y ~ poly(x, 1)) 
# Degree = 2 resulting in a quadratic fit
general.scatter(df, by="B", smooth.fun="glm", smooth.formula=y ~ poly(x, 2))

## Including a density contour
general.scatter(df, by="B", smooth.fun="glm", smooth.formula=y ~ poly(x, 2),
                density=T, legend.pos="bottom")
}
}
\author{
Jens Hooge
}
