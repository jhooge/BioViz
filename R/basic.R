#' Correlation Coefficients Panel
#'
#' @description
#' Plot absolute correlation coefficients (pearson) into panel
#' and adjust the text size according to the correlation
#'
#' @import graphics
#' @import stats
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param digits Number of digits the correlation will be displayed with (Default: 2)
#' @param prefix Not yet documented
#' @param cex.cor Not yet documented
#' @param ... Additional parameters
#'
#' @return plot text element
panel_cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


#' Histogram Panel
#'
#' @description
#' Plot histogram into diagonal panel of a numeric vector
#'
#' @import graphics
#'
#' @param x numeric vector
#' @param ... Additional parameters
#'
#' @return histogram with colored bars
panel_hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "#67000D", ...)
}

#' Smooth Scatter Panel
#'
#' @description
#' Plot a smoothScatter (color density of a scatter plot)
#' with a loess fit.
#'
#' @import graphics
#' @import grDevices
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param bg Background (Default: NA)
#' @param cex Not yet documented
#' @param col.smooth Not yet documented
#' @param span Not yet documented
#' @param iter Not yet documented
#' @param ... Additional parameters
#'
#' @return smoothed scatter plot
panel_smoothScatter <- function (x, y, bg = NA,
                                 cex = 1, col.smooth = "red",
                                 span = 2/3, iter = 3, ...) {
  # colors for the density
  palette <- colorRampPalette(c("blue", "orange", "red"))
  s <- smoothScatter(x, y, colramp = palette, bg = bg, cex = cex, add=T)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = col.smooth, ...)
}


#' Pairs Plot
#'
#' @description
#' A matrix of plots is produced. The upper diagonal matrix
#' consists of the correlation between pairs of variables. The
#' lower diagonal matrix shows the scatter of each sample including
#' a density overlay and a loess fit. Finally the diagonal shows a histogram
#' of each of the input variables.
#'
#' @import graphics
#' @import grDevices
#'
#' @param mat A matrix or a dataframe with numerical columns
#' @param ... Additional parameters
#'
#' @return Plots a pairs plot on screen
#'
#' @examples
#' data(iris)
#' plot_pairs(iris[, 1:4])
#'
#' @export
plot_pairs <- function(mat, ...) {
  pairs(mat,
       lower.panel = panel_smoothScatter,
       upper.panel = panel_cor,
       diag.panel  = panel_hist,
       ...)
}


#' Plots a volcano plot.
#'
#' @description
#' The plot does not aim to show any details about which genes or
#' pathways are over/under represented, but is meant to show only
#' the overall directionality of the genes between groups. If there
#' is a certain perturbation in a system, the volcano plot may reveal
#' a response in which there is a dramatic increase or decreased in
#' transcription. In either case, we can get a better idea of the overall
#' system by showing a volcano plot. Each gene is represented with its p-value
#' as the x-axis and its (log2) fold change as the y-axis.
#'
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#'
#' @param data A data.frame which includes columns for -log10(p-values) and log2(foldchange)
#' @param title The figure title
#' @param groups A character vector with the treatment groups
#' @param top_labeled A numeric value that controls how many samples should be highlighted based on their pvalue
#' @param xlabel X-Axis Label
#' @param ylabel Y-Axis Label
#' @param xcutoff A numeric value for the -log10(p-value) cutoff
#' @param ycutoff A numeric value for the log2(foldchange) cutoff
#' @param log1p A flag to apply a log(1+x) transformation on the y axis to account for large values
#'
#' @return A ggplot2 figure
#' @export
#'
#' @examples
#' data(iPSC)
#' volcano_plot(iPSC)
#' volcano_plot(iPSC, groups=c("A", "B"))
#' volcano_plot(iPSC, top_labeled=10, xcutoff=c(-log2(2), log2(4)))
#' volcano_plot(iPSC, top_labeled=10, log1p=TRUE)
volcano_plot <- function(data, title="Volcano Plot",
                                 groups=c("Untreated", "Treated"),
                                 top_labeled=0,
                                 xlabel=expression(log[2]("Foldchange")),
                                 ylabel=expression(-log[10]("p-value")),
                                 xcutoff=c(-log2(4), log2(4)),
                                 ycutoff=-log10(.05),
                                 log1p=FALSE) {
  if(is.null(groups)) {
    groups <- ""
  }
  if(length(groups)!=2) {
    warning("Number of group labels should be 2. No group colors applied!")
  }

  data$fc <- log2(data$fc)
  data$pvalue <- -log10(data$pvalue)

  # Color corresponds to fold change directionality
  fig <- ggplot() +
    ggtitle(label = title) + # Add a title
    xlab(xlabel) +  # x-axis label
    ylab(ylabel) +  # y-axis label
    geom_vline(xintercept = 0, colour = "black") + # add line at 0
    geom_hline(yintercept = ycutoff, colour = "black") + # p(0.05) = 1.3
    scale_x_continuous(limits = c(xcutoff[1]-2, xcutoff[2]+2)) +
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position = "none") # remove legend

  if (length(groups)==2) {
    ## Set group labels
    data <- data %>%
      mutate(color = ifelse(data$fc > 0 & data$pvalue > ycutoff,
                            yes = groups[1],
                            no = ifelse(data$fc < 0 & data$pvalue > ycutoff,
                                        yes = groups[2],
                                        no = "none")))

    ## Color points by group
    colors <- c("#E64B35", "#3182bd", "#636363")
    names(colors) <- c(groups[1], groups[2], "none")

    fig <- fig +
      geom_point(data=data, aes(x = fc, y = pvalue, color = factor(color)),
                 size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
      annotate(geom = "text",
               label = groups[1],
               x = xcutoff[1], y = max(data$pvalue)+1,
               size = 5, colour = "black") + # add Untreated text
      annotate(geom = "text",
               label = groups[2],
               x = xcutoff[2], y = max(data$pvalue)+1,
               size = 5, colour = "black") + # add Treated text
      scale_color_manual(values = colors) # change colors
  } else {
    fig <- fig +
      geom_point(data=data, aes(x = fc, y = pvalue),
                 size = 1.75, alpha = 0.8, na.rm = T) + # Make dots bigger
      scale_colour_gradient(low = "black", high = "black", guide = FALSE) # Color black
  }

  if(log1p) {
    fig <- fig + scale_y_continuous(trans = "log1p") ## Apply log(1+x) transformation
  }

  ## Label top pvalue points by their id
  if (top_labeled > 0) {
    top_labeled <- top_n(data, n = top_labeled, wt = pvalue)
    fig <- fig +
      geom_text_repel(data = top_labeled,
                      mapping = aes(x = fc, y = pvalue, label = id),
                      size = 3,
                      fontface = 'bold',
                      color = 'black',
                      box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.5, "lines"))
  }
  return(fig)
}


