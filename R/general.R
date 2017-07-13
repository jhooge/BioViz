source("R/utils.R")

#' Correlation Coefficients Panel
#' 
#' @author Jens Hooge
#'
#' @description
#' Plot absolute correlation coefficients (pearson) into panel
#' and adjust the text size according to the correlation
#'
#' @import graphics
#' @import stats
#' 
#' @keywords internal
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
#' @author Jens Hooge
#'
#' @description
#' Plot histogram into diagonal panel of a numeric vector
#'
#' @import graphics
#' 
#' @keywords internal
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
#' @author Jens Hooge
#'
#' @description
#' Plot a smoothScatter (color density of a scatter plot)
#' with a loess fit.
#'
#' @import graphics
#' @import grDevices
#' 
#' @keywords internal
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
#' @author Jens Hooge
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
#' @seealso pairs
#'
#' @return Plots a pairs plot on screen
#'
#' @examples
#' \dontrun{
#' data(iris)
#' general.plot_pairs(iris[, 1:4])
#' }
#' 
#' @export
general.plot_pairs <- function(mat, ...) {
  pairs(mat,
       lower.panel = panel_smoothScatter,
       upper.panel = panel_cor,
       diag.panel  = panel_hist,
       ...)
}


#' Plots a volcano plot.
#' 
#' @author Jens Hooge
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
#' \dontrun{
#' data(iPSC)
#' ## Plot Simple Volcano Plot
#' general.volcano_plot(iPSC)
#' 
#' ## Plot Simple Volcano Plot with Group Labels
#' general.volcano_plot(iPSC, groups=c("A", "B"))
#' 
#' ## Plot Volcano Plot and add top10 most significant samples
#' general.volcano_plot(iPSC, top_labeled=10, xcutoff=c(-log2(2), log2(4)))
#' 
#' ## Plot Volcano Plot and add Top10 most Significant Samples and Transform y-axis by log(1+x)
#' general.volcano_plot(iPSC, top_labeled=10, log1p=TRUE)
#' }
general.volcano_plot <- function(data, title="Volcano Plot",
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
    data$color <- as.factor(data$color)
    ## This hack is needed for top_n to work properly during package build
    data <- cbind(data$color, data)
    data[,ncol(data)] <- NULL
    colnames(data) <- c("color", "id", "fc", "pvalue")

    ## Color points by group
    colors <- c("#E64B35", "#3182bd", "#636363")
    names(colors) <- c(groups[1], groups[2], "none")

    fig <- fig +
      geom_point(data=data, aes_string(x = "fc", y = "pvalue", color = "color"),
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
      geom_point(data=data, aes_string(x = "fc", y = "pvalue"),
                 size = 1.75, alpha = 0.8, na.rm = T) + # Make dots bigger
      scale_colour_gradient(low = "black", high = "black", guide = FALSE) # Color black
  }

  if(log1p) {
    fig <- fig + scale_y_continuous(trans = "log1p") ## Apply log(1+x) transformation
  }

  ## Label top pvalue points by their id
  if (top_labeled > 0) {
    # top_labeled <- top_n(data, n = top_labeled, wt = pvalue)
    top_labeled <- top_n(data, n = top_labeled)
    fig <- fig +
      geom_text_repel(data = top_labeled,
                      mapping = aes_string(x = "fc", y = "pvalue", label = "id"),
                      size = 3,
                      fontface = 'bold',
                      color = 'black',
                      box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.5, "lines"))
  }
  return(fig)
}


#' Function to create a barplot, with bars split by a second variable
#' 
#' @author Sebastian Voss, Adam Skubala
#'
#' @param freq vector or matrix with the bar heights, e.g. a frequency table ('table()' output), 
#'     names will be used for bar labels if 'labels' is not specified, if 'freq' is a matrix, 
#'     entries of each column are represented by separate bars with no space in between
#' @param labels bar labels (default=NULL)
#' @param labels.ab additional labels, which are placed above the bars (default=NULL), should 
#'     have the same dimension as 'freq'
#' @param file.name name for the PNG output file (without '.png' ending)
#' @param col color for the bars (default='tomato3'). can be a single value or an object of 
#'     the same dimension than freq, specifying the color for each bar
#' @param cex.lab size of the bar labels (default=1)
#' @param cex.lab.ab size of the bar labels above the bars (default=1) 
#' @param cex.x size of the x-labels (default=1)
#' @param xlab label for the x-axis (default='')
#' @param ylab label for the y-axis (default='Frequency')
#' @param add.legend TRUE/FALSE should a legend be added underneath the plot (default=TRUE), 
#'     only works if legend parameters are specified
#' @param leg.lab labels for the legend (default=NULL)
#' @param leg.col colors for the legend (default=NULL)
#' @param leg.title title for the legend (default=NULL)
#' @param mar plot margin (default=c(5,4,3,1))
#' @param png.width width of the PNG-file in inches (default=6)
#' @param png.height height of the PNG-file in inches (default=5)
#' @param png.out TRUE/FALSE if true, a png file is created as output, otherwise the plot 
#'     is created in the R graphic device
#'
#' @return PNG file with the barplot ('file.name.png') numeric vector with size of the PNG-file 
#'     in inches for rescaling in RTF
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- sample(1:15, size=500, replace=TRUE)
#' by <- sample(c('A','B'), size=500, replace=TRUE)
#' freq <- table(by, x)
#' col <- matrix(ncol=ncol(freq), nrow=nrow(freq))
#' 
#' col1 <- colorRampPalette(c('tomato3','grey95'))
#' col[1,] <- col1(9)[1]
#' col[2,] <- col1(9)[7]
#' 
#' col2 <- colorRampPalette(c('skyblue2','grey95'))
#' col[1,c(5,11)] <- col2(9)[1]
#' col[2,c(5,11)] <- col2(9)[7]
#' 
#' cex.lab <- 0.7
#' cex.lab.ab <- 0.5
#' cex.x <- 0.9
#' labels <- NULL
#' labels.ab <- rep(paste0('X', 1:15), each=2)
#' labels.ab[which(1:length(labels.ab) %% 2 == 0)] <- ''
#' mar <- c(7,4,3,1)
#' xlab <- 'numbers'
#' ylab <- 'frequencies'
#' file.name <- 'test'
#' png.width <- 7
#' png.height <- 6
#' add.legend <- TRUE
#' leg.title <- 'color'
#' leg.col <- c(col1(9)[1],col1(9)[7],col2(9)[1],col2(9)[7])
#' leg.lab <- c('red','faded red','blue','faded blue')
#' 
#' general.bar_plot_by(freq=freq, labels=NULL, labels.ab=labels.ab, 
#'     file.name='test1_barplot', col=col, cex.lab=0.7, cex.lab.ab=0.5, cex.x=0.8,
#'       xlab='Number', ylab='Frequency', add.legend=TRUE, 
#'       leg.lab=leg.lab, leg.col=leg.col, leg.title='color',
#'       mar=c(7,4,3,1), png.width=7, png.height=6)
#' general.bar_plot_by(freq=freq, labels=paste0('X', 1:15), labels.ab=NULL, 
#'     file.name='test2_barplot', 
#'     col=col, cex.lab=0.7, cex.lab.ab=0.5, cex.x=0.8,
#'     xlab='Variable', ylab='Frequency', add.legend=TRUE, 
#'     leg.lab=leg.lab, leg.col=leg.col, leg.title='color',
#'     mar=c(7,4,3,1), png.width=7, png.height=6)
#' general.bar_plot_by(freq=table(x), labels.ab=paste0('X', 1:15), 
#'     file.name='test3_barplot', col=col[2,5],
#'       xlab='Number', ylab='Frequency')
#' general.bar_plot_by(freq=table(x), labels.ab=paste0('X', 1:15), 
#'     file.name='test4_barplot', col=col[2,5], mar=c(7,4,3,1),
#'     xlab='Number', ylab='Frequency', add.legend=TRUE, 
#'     leg.title='color', leg.col=leg.col[4], leg.lab=leg.lab[4])
#' }
general.bar_plot_by <- function(freq, labels=NULL, labels.ab=NULL, file.name, col='tomato3', cex.lab=1, cex.lab.ab=1, cex.x=1,
                        xlab='', ylab='Frequency', add.legend=FALSE, leg.lab=NULL, leg.col=NULL, leg.title=NULL,
                        mar=c(5,4,3,1), png.width=6, png.height=5, png.out=TRUE){ 
  
  #extract labels from 'freq' if 'labels' is not specified
  if(is.null(labels)){
    if(is.matrix(freq)){
      labels <- colnames(freq)
    }else{
      labels <- names(freq)
    }
  }
  
  #define length of y-axis
  if(is.null(labels.ab)){
    ymax <- max(freq)
  }else{
    ymax <- 1.05*max(freq)
  }
  
  #create bar plot
  create.plot <- function(){
    
    #middle points of the bars, which will be added to the plot
    bar.mp <- barplot(freq, plot=FALSE, beside=TRUE)
    ##if 'freq' is a vector (no by variable specified), middle points are given as a matrix with one column,
    ##whereas they are given as a matrix with one row for each bylevel, when 'freq' is a matrix
    if (ncol(bar.mp)==1){
      bar.mp <- t(bar.mp)
    }
    bar.mp.mp <- colMeans(bar.mp)
    #set plot margins
    par(mar=mar)
    #empty plot
    plot(0,0, xlim=c(min(bar.mp)-0.5,max(bar.mp)+0.5), ylim=c(0,ymax),
         type='n', xlab='', ylab=ylab, main='', axes=FALSE)
    #color background of the plotting region
    rect(par('usr')[1], par('usr')[3],par('usr')[2],par('usr')[4], col='grey95')
    #add customized axes with 45 degree labels on the x-axis
    box()
    axis(1, at=bar.mp.mp, labels=rep('',length(bar.mp.mp)), tcl=-0.3)
    axis(2)
    #add labels to the x-axis (45-degree angle)
    y.text <- grconvertY(-0.04, from = 'npc', to = 'user') #y-coordinate for the labels
    text(bar.mp.mp, y.text, labels=labels, xpd=TRUE, srt=-45, cex=cex.lab, adj=c(0,0.5))
    #add bars
    barplot(freq, add=TRUE, axes=FALSE, axisnames=FALSE, border=NA, col=col, beside=TRUE)
    
    #add labels above the bars
    if(!is.null(labels.ab)){
      text(bar.mp, freq+0.02*max(freq), labels=labels.ab, adj=c(0.5,0), cex=cex.lab.ab)
    }
    
    ##add legend
    if(add.legend==TRUE){
      x.leg <- (par('usr')[1]+par('usr')[2])/2
      y.leg <- grconvertY(0, from='ndc', to='user')
      l <- legend(x.leg, y.leg, legend=leg.lab, fill=leg.col, border='white',
                  xjust=0.5, yjust=0, title=leg.title, xpd=TRUE, horiz=FALSE, ncol=2,
                  bty='n', cex=0.8)
    }
    
    #add x-label
    if(add.legend==TRUE){
      text((par('usr')[1]+par('usr')[2])/2, l$rect$top, paste(xlab,'\n\n',sep=''), xpd=TRUE, cex=cex.x)
    }else{
      mtext(xlab, cex=cex.x, side=1, line=par('mar')[1]-2)
    }
    
  }
  
  #PNG or standard R graphic device
  if (png.out==TRUE){
    size <- c(png.height,png.width)
    names(size) <- c('height','width')
    #create png-file
    png(gsub(' ','_',paste(file.name,'.png',sep='')), height=size[1], width=size[2], units='in', res=300)
    create.plot()
    dev.off()
    #return PNG size (useful for rescaling)
    return(size)
  }else{
    create.plot()
  }
  
}


#' Function to create a simple barplot
#'
#' @author Sebastian Voss, Adam Skubala
#'
#' @param freq vector with the bar heights, e.g. a frequency table ('table()' output),
#'     names will be used for bar labels if 'labels' is not specified
#' @param labels bar labels (default=NULL)
#' @param labels.ab additional labels, which are placed above the bars (default=NULL)
#' @param file.name name for the PNG output file (without '.png' ending)
#' @param col color for the bars (default='tomato3'). if more than one color is specified,
#                 each bar is filled with the repsective color
#' @param cex.lab size of the bar labels (default=1)
#' @param cex.lab.ab size of the bar labels above the bars (default=1)
#' @param cex.x size of the x-labels (default=1)
#' @param xlab label for the x-axis (default='')
#' @param ylab label for the y-axis (default='Frequency')
#' @param add.legend TRUE/FALSE should a legend be added underneath the plot 
#'                   only works if legend parameters are specified (default=TRUE)
#' @param leg.lab labels for the legend (default=NULL)
#' @param leg.col colors for the legend (default=NULL)
#' @param leg.title title for the legend (default=NULL)
#' @param mar plot margin (default=c(5,4,3,1))
#' @param png.width width of the PNG-file in inches (default=6)
#' @param png.height height of the PNG-file in inches (default=5)
#'
#' @return PNG file with the barplot ('file.name.png') and a 
#'     numeric vector with size of the PNG-file in inches for rescaling in RTF
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- sample(1:30, size=500, replace=TRUE)
#' freq <- table(x)
#' 
#' general.bar_plot_simple(freq, labels.ab=1:30, file.name="test", 
#'                 col=rep(c('tomato3','tomato3','skyblue2'), 10), 
#'                 cex.lab=.7, cex.lab.ab=.8, cex.x=.9,
#'                 ylab="Bar Height", 
#'                 add.legend=TRUE, leg.title="Color",
#'                 leg.lab=c('red','blue'), leg.col=c('tomato3','skyblue2'),
#'                 mar=c(4,4,3,1))
#' }
general.bar_plot_simple <- function(freq, labels=NULL, labels.ab=NULL, file.name, col='tomato3', cex.lab=1, cex.lab.ab=1, cex.x=1,
                            xlab='', ylab='Frequency', add.legend=FALSE, leg.lab=NULL, leg.col=NULL, leg.title=NULL,
                            mar=c(5,4,3,1), png.width=6, png.height=5){ 
  
  #extract labels from 'freq' if 'labels' is not specified
  if(is.null(labels)){
    labels <- names(freq)
  }
  
  #define length of y-axis
  if(is.null(labels.ab)){
    ymax <- max(freq)
  }else{
    ymax <- 1.05*max(freq)
  }
  
  #create bar plot
  png(paste0(file.name,'.PNG'), width=png.width, height=png.height, units='in', res=300)
  
  #middle points of the bars, which will be added to the plot
  bar.mp <- barplot(freq, plot=FALSE)
  #set plot margins
  par(mar=mar)
  #empty plot
  plot(0,0, xlim=c(min(bar.mp)-0.5,max(bar.mp)+0.5), ylim=c(0,ymax),
       type='n', xlab='', ylab=ylab, main='', axes=FALSE)
  #color background of the plotting region
  rect(par('usr')[1], par('usr')[3],par('usr')[2],par('usr')[4], col='grey95')
  #add customized axes with 45 degree labels on the x-axis
  box()
  axis(1, at=bar.mp, labels=rep('',length(bar.mp)), tcl=-0.3)
  axis(2)
  #add labels to the x-axis (45-degree angle)
  y.text <- grconvertY(-0.04, from = 'npc', to = 'user') #y-coordinate for the labels
  text(bar.mp, y.text, labels=labels, xpd=TRUE, srt=-45, cex=cex.lab, adj=c(0,0.5))
  #add bars
  barplot(freq, add=TRUE, axes=FALSE, axisnames=FALSE, border=NA, col=col)
  
  #add labels above the bars
  if(!is.null(labels.ab)){
    text(bar.mp, freq+0.02*max(freq), labels=labels.ab, adj=c(0.5,0), cex=cex.lab.ab)
  }
  
  ##add legend
  if(add.legend==TRUE){
    x.leg <- (par('usr')[1]+par('usr')[2])/2
    y.leg <- grconvertY(0, from='ndc', to='user')
    l <- legend(x.leg, y.leg, legend=leg.lab, fill=leg.col, border='white',
                xjust=0.5, yjust=0, title=leg.title, xpd=TRUE, horiz=FALSE,
                bty='n', cex=0.8)
  }
  
  #add x-label
  if(add.legend==TRUE){
    text((par('usr')[1]+par('usr')[2])/2, l$rect$top, paste(xlab,'\n\n',sep=''), xpd=TRUE, cex=cex.x)
  }else{
    mtext(xlab, cex=cex.x, side=1, line=par('mar')[1]-2)
  }
  
  dev.off()
  
  #return PNG size (useful for importing plot into an RTF document)
  size <- c(png.height,png.width)
  names(size) <- c('height','width')
  return(size)
}


#' Function to create boxplots of a metric variable, grouped by two nominal variables
#' 
#' @author Jens Hooge
#' 
#' @import ggplot2
#'
#' @param data A data frame with group variables and a column with numeric values. (Required)
#' @param group_by1 First grouping variable for which the plot will be facetted by. (Required)
#' @param group_by2 Second grouping variable which defines the subgroup of group1 (Required)
#' @param col.jitter_by The variable name by which the jitter points will be colored by. (Default: NULL)
#' @param shape_jitter_by The variable name by which the jitter points will be shaped by. (Default: NULL)
#' @param title The title of the plot (Default: NULL)
#' @param xlab The x-axis label (Default: NULL)
#' @param ylab The y-axis label (Default: NULL)
#' @param legend Flag indicating whether the legend should be displayed (Default: FALSE)
#' @param rotate Flag indicating whether the x-axis labels should be rotated by 45° (Default: FALSE)
#'
#' @return A ggplot figure
#' @export
#'
#' @examples
#' \dontrun{
#' require(reshape2)
#' n <- 500
#' df <- data.frame(A=sample(c("a1", "a2", "a3"), n, replace=TRUE),
#'                  B=sample(c("b1", "b2"), n, replace=TRUE),
#'                  C=sample(c("c1", "c2", "c3"), n, replace=TRUE),
#'                  D=sample(c("d1", "d2", "d3"), n, replace=TRUE),
#'                  value=rnorm(n, 100, 1))
#'                  
#' general.box_plot_facetted(df, group_by1 = "A", group_by2 = "B",
#'                   col.jitter_by = "A", shape_jitter_by = "C", 
#'                   rotate=TRUE)
#' }
general.box_plot_facetted <- function(data, group_by1, group_by2, 
                              col.jitter_by=NULL, shape_jitter_by=NULL,
                              title=NULL, xlab=NULL, ylab=NULL, legend=FALSE,
                              rotate=FALSE) {
  
  ## Check Inputs
  stopifnot(group_by1 %in% colnames(data),
            group_by2 %in% colnames(data))
  
  facet <- paste0("~", group_by2)
  fig <- ggplot(df, aes_string(x=group_by1, y="value")) + 
    geom_boxplot(outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", width = 0.5, linetype="dashed") +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=3, size=5, show.legend = FALSE) +
    geom_jitter(aes_string(shape=shape_jitter_by, 
                           col=col.jitter_by), alpha=.7) +
    ggtitle(label = title) +
    xlab(xlab) + 
    ylab(ylab) + 
    facet_grid(facet)
    theme_bw()
  
  if(!legend) {
    fig <- fig + theme(legend.position = "none")
  }
  if(rotate) {
    fig <- fig + theme(axis.text.x=element_text(angle=45, hjust=1))
  }
  
  return(fig)
}


#' Function to create a correlation plot.
#' 
#' @author Sebastian Voss, Adam Skubala
#' 
#' @import ellipse
#' 
#' @TODO This function has to be refactored and the number of 
#' input parameters have to be reduced.
#' 
#' @description 
#' The upper diagonale of the correlation matrix is visualized by colored ellipses,
#' with the color, shape and direction of the ellipsis showing the dependence of
#' the respective variable pair. the lower diagonale of the matrix-like plot
#' contains the corresponding correlation values.
#' 
#' @param data data set in wide format (Required)
#' @param cor.type type of correlation coefficient, 'pearson', 'spearman' or 'kendall' (Default: 'spearman')
#' @param col vector of length two with the colors for negative and positive correlations (Default: c('tomato1','green3'))
#' @param plot.type plot type in the upper diagonale, 'ellipse' or 'points', i.e. scatterplots (Default='ellipse')
#' @param size.cor text size of the correlation numbers (Default: 1)
#' @param add.hist add histograms on the left hand side of the plot (Defaul=FALSE)
#' @param disp.sig.lev ellipses of correlations with a p-value above this threshold will not be displayed
#'     and the corresponding correlation number will be displayed in grey (Default: 1)
#' @param disp.thresh ellipses of correlations with an absolute value below this threshold will not be displayed
#      and the corresponding correlation number will be displayed in grey (Default: 0.1) 
#' @param conf.level level for the confidence intervals (Does not work for 'kendall') (Default: Null) 
#' @param n.boot number of bootstrap iterations, if confidence intervals for 'cor.type='spearman'' are calculated (Default: 1000)
#' @param pt.size size of the points, if plot.type='points' (Default: 0.4)
#' @param group.cut number of the last variable in the first group if variables can be split up into two groups,
#      within group correlations will be faded out (Default: Null)
#' @param border color of the border around the ellipses (Default: NA)
#'
#' @param png.out output as PNG? (Default: FALSE)
#' @param png.width width of PNG in inches (Default: 5)
#' @param png.height height of PNG in inches (Default: 5) 
#' @param file.name filename of the PNG without '.png'
#' @param ... Additional parameters for par()
#'
#' @return The correlation plot
#' @export
#'
#' @examples
#' \dontrun{
#' data(mtcars)
#' general.corplot(data=mtcars, cor.type='pearson', disp.sig.lev=0.05, group.cut=3)
#' general.corplot(data=mtcars, cor.type='pearson', disp.sig.lev=0.05, cex=0.95)
#' general.corplot(data=mtcars, cor.type='pearson', disp.sig.lev=0.05, group.cut=4, 
#'                 conf.level=0.95, cex=0.95)
#' general.corplot(data=mtcars, cor.type='pearson', disp.sig.lev=0.05, group.cut=4, 
#'                 plot.type='points', conf.level=0.95, cex=0.95)
#' general.corplot(data=mtcars, cor.type='pearson', disp.thresh=0.5, group.cut=4,
#'                 plot.type='points', conf.level=0.95, cex=0.95, add.hist=TRUE)
#' }
general.corplot <- function(data, cor.type=c('pearson','spearman','kendall'), 
                            col=c('tomato1','green3'), plot.type=c('ellipse','points'), 
                            size.cor=1, add.hist=FALSE, disp.sig.lev=1, disp.thresh=0.1, 
                            conf.level=NULL, n.boot=1000, pt.size=0.4, group.cut=NULL, border=NA,
                            png.out=FALSE, png.width=5, png.height=5, file.name=NULL, ...) {
  #check and recycle parameters
  cor.type <- match.arg(cor.type)
  plot.type <- match.arg(plot.type)
  
  #calculate correlation matrix
  corm <- cor(data, method=cor.type, use='pairwise.complete.obs')
  #also in text format, rounded to 2 digits
  corm.text <- format(round(corm,2), nsmall=2)
  
  #create color gradient function for the specified colors
  col.fun <- colorRamp(c(col[1],'white',col[2]), bias=1)
  
  #plot function
  create.plot <- function(){
    
    #calculate margin size, based on par-settings and variable name lengths 
    par(mar=c(0.5,0.5,0.5,0.5), ...)
    ##string length
    mar.var <- max(strwidth(colnames(data), units='inches'))
    ##transfrom from inches to lines
    mar.var <- mar.var*(par('mar')/par('mai'))[1] 
    ##set margins
    par(mar=c(0.5,1+mar.var+0.5,1+mar.var+0.5,0.5), ...)
    #empty plot
    plot(0,0, xlim=c(0.5,ncol(corm)+ifelse(add.hist,1.5,0)+0.5), ylim=c(-(nrow(corm)+0.5),-0.5), type='n',
         axes=FALSE, xlab='', ylab='')
    #background rectangle and grid lines
    rect(0.5, -(nrow(corm)+0.5), ncol(corm)+ifelse(add.hist,1.5,0)+0.5, -0.5, col='grey95', border=NA)
    abline(v=seq(0.5,ncol(corm)+0.5,1), col='white')
    abline(h=seq(-(nrow(corm)+0.5),-0.5,1), col='white')
    abline(0,-1, col='white', lwd=6)
    #add histograms
    if (add.hist==TRUE){
      #separate histogram column from the rest by white rectangle
      rect(ncol(corm)+0.5, -(nrow(corm)+0.5), ncol(corm)+0.5+0.5, -0.5, col='white', border=NA)
      for (i in 1:ncol(data)){
        #calculate histogram
        hist.i <- hist(data[,i], plot=FALSE)
        x <- hist.i$breaks
        y <- c(0,hist.i$counts)
        #transform into panel coordinates
        tr.i <- utils.trafo_pan(x=x, y=y, pan.row=i, pan.col=ncol(corm)+1.5)
        #add histogram into panel
        rect(xleft=tr.i$xt[-length(tr.i$xt)], ybottom=rep(tr.i$yt[1],length(tr.i$yt)-1),
             xright=tr.i$xt[-1], ytop=tr.i$yt[-1], col='grey60', border='grey95')
      }
    }
    #add row- and column names
    text(x=rep(0, ncol(data)), y=-(1:ncol(data)), labels=colnames(data), adj=c(1,0.5),
         xpd=TRUE, cex=1*par('cex'))
    text(x=1:ncol(data), y=rep(0, ncol(data)), labels=colnames(data), adj=c(0,0.5),
         srt=90, xpd=TRUE, cex=1*par('cex'))
    
    #add ellipses and text for variables i and j
    for (i in 1:(ncol(corm)-1)){
      for (j in (i+1):nrow(corm)){
        
        #calculate p-value 
        p.ij <- cor.test(data[,i], data[,j], method=cor.type, exact=ifelse(cor.type=='spearman',FALSE,TRUE))$p.value
        #if correlation is significant and greater than the specified threshold, add ellipse in the upper triangle
        if (p.ij<=disp.sig.lev & abs(corm[i,j])>=disp.thresh){
          if(plot.type=='ellipse'){
            ell <- ellipse(corm[i,j], t=0.4, npoints=1000)
            ell[,1] <- ell[,1] + j
            ell[,2] <- ell[,2] - i
            polygon(ell, col=rgb(col.fun((corm[i,j]+1)/2), maxColorValue=255), border=border, lwd=0.5)
          }else if(plot.type=='points'){
            tr.ij <- utils.trafo_pan(x=data[,j], y=data[,i], pan.row=i, pan.col=j)
            points(tr.ij$xt, tr.ij$yt, pch=20, cex=pt.size, col=rgb(col.fun((corm[i,j]+1)/2), maxColorValue=255))
          }
        }
        #add correlation number (and confidence interval) in the lower triangle
        
        if(!is.null(conf.level) & cor.type!='kendall'){
          
          cor.height <- strheight(corm.text[i,j], cex=0.7*size.cor*par('cex'))
          
          segments(i, -j-2.25*cor.height/2, i, -j+2.25*cor.height/2, col='grey66', lwd=0.75)
          segments(i-0.1, -j-2.25*cor.height/2, i+0.1, -j-2.25*cor.height/2, col='grey66', lwd=0.75)
          segments(i-0.1, -j+2.25*cor.height/2, i+0.1, -j+2.25*cor.height/2, col='grey66', lwd=0.75)
          
          rect(i-0.45, -j-1.25*cor.height/2, i+0.45, -j+1.25*cor.height/2, col='grey95', border=NA)
          
          if (cor.type=='spearman'){
            
            set.seed(19121909)
            cor.ci <- cor.ci(na.exclude(data[,c(i,j)]), n.iter=n.boot,  p=1-conf.level, method=cor.type, plot=FALSE)
            ci.l <- format(round(cor.ci$ci$low.e, 2), nsmall=2)
            ci.u <- format(round(cor.ci$ci$up.e, 2), nsmall=2)
            
          } else if (cor.type=='pearson'){
            
            cor.ci <- cor.test(data[,i], data[,j], method=cor.type, exact=TRUE, conf.level=conf.level)$conf.int
            ci.l <- format(round(cor.ci[1], 2), nsmall=2)
            ci.u <- format(round(cor.ci[2], 2), nsmall=2)
            
          }
          
          text(i,-j-2*cor.height/2, ci.l, adj=c(0.5, 1.5), cex=0.4*size.cor*par('cex'), col='grey66')
          text(i,-j+2*cor.height/2, ci.u, adj=c(0.5,-0.5), cex=0.4*size.cor*par('cex'), col='grey66')
          
        } else if (!is.null(conf.level) & cor.type=='kendall'){
          
          message('no CI for cor.type=\'kendall\'')
          
        }
        
        text(i,-j, corm.text[i,j], adj=c(0.5,0.5), cex=0.7*size.cor*par('cex'),
             col=ifelse(p.ij<=disp.sig.lev & abs(corm[i,j])>=disp.thresh,'black','grey50'))
        
        
      }
    }
    
    #if two groups are specified by a cutoff, the within group correlation matrix is faded out by a transparent rectangle
    if (!is.null(group.cut)){
      trans <- 0.75
      rect(0.5, -(group.cut+0.5), group.cut+0.5, -0.5, col=adjustcolor('white', alpha.f=trans), border=NA)
      rect(0.5+group.cut, -(nrow(corm)+0.5), ncol(corm)+0.5, -(0.5+group.cut), col=adjustcolor('white', alpha.f=trans), border=NA)
    }
    
  }
  
  if (png.out==TRUE){
    #output as PNG...
    png(paste0(file.name,'.png'), height=png.height, width=png.width, units='in', res=600)
    create.plot()
    dev.off()
    size <- c(png.height, png.width)
    names(size) <- c('height','width')
    invisible(size)
  }else{
    #...or in R device
    create.plot()
  }
}


#' Function to create histograms to check for the necessity of a log-transformation of metric data.
#' 
#' @author Sebastian Voss, Adam Skubala
#' 
#' @description 
#' The function creates two histograms, one of the original data
#' and one of the log-tramsformed data. If desired, a second time point can be added,
#' so that a plot with four histograms is created.
#'
#' @param x Numeric vector (Required)
#' @param file.name begining of the name for the PNG output files (without '.png' ending) (Required)
#' @param by Character or factor vector of the same length than 'x' to 
#'     the time, if two time points should be analyzed. should
#'     contain only two levels (Default: Null)
#' @param nbars number of bars for the histograms (Default: 15)
#' @param title Title of the plot (Default: '')
#' @param xlab Label for the x-axis (Default: '')
#' @param col.orig Color of the histogram bars of the original values (Default: 'lightblue')
#' @param col.log Color of the histogram bars of the log-transformed values (Default: 'lightyellow')
#'
#' @return A .png file with two or four histograms and a numeric vector
#'     with the size of the .png
#' @export
#'
#' @examples
#' \dontrun{
#' x <- exp(rnorm(1000))
#' by <- factor(rep(c('baseline','day 30'), each=500), 
#'              level=c('baseline','day 30'))
#' general.hist_log(x=x, by=by, file.name='test_hist',
#'          title='Test Historgram', xlab='Measurement')
#' }
general.hist_log <- function(x, by=NULL, file.name, nbars=15, title='', xlab='',
                                        col.orig='lightblue', col.log='lightyellow'){
  
  #define some graphical parameters for the case of just one time point
  height.png <- 8
  png.row <- 1
  oma3 <- 5
  mar1 <- 10
  xlab.log <- 'log. values'
  if (xlab==''){
    xlab <- 'orig. values'
  }
  #rename data vector
  x1 <- x
  if (!is.null(by)){
    #change graphical parameters if two time points have to be analyzed
    height.png <- 14
    png.row <- 2
    oma3 <- 7
    mar1 <- 12
    #convert 'by' into factor
    by <- factor(by)
    #split data vector
    x1 <- x[by==levels(by)[1]]
    x2 <- x[by==levels(by)[2]]
  }
  #create PNG
  png(paste0(file.name,'.png'), width=14, height=height.png, units='in', res=300)
  #set graphical parameters
  par(mfrow=c(png.row,2), oma=c(1,0,oma3,0), mar=c(mar1,7,3,2), cex.lab=2, cex.axis=2,
      cex.main=2, mgp=c(4.5,1.5,0))
  
  #create two histograms and add appropriate normal density
  hist.data <- hist(x1, freq=TRUE, col=col.orig, xlab=xlab, ylab='frequency', main='',
                    breaks=utils.exBar(x1,nbars))
  #normal density is scaled to fit the frequency distribution
  curve(sum(hist.data$count)*diff(hist.data$breaks)[1]*dnorm(x, mean(x1), sd(x1)),
        from=min(x1), to=max(x1), add=TRUE, lty=2, lwd=2)
  
  hist.data <- hist(log(x1), freq=TRUE, col=col.log, xlab=xlab.log, ylab='frequency', main='',
                    breaks=utils.exBar(log(x1),nbars))
  curve(sum(hist.data$count)*diff(hist.data$breaks)[1]*dnorm(x, mean(log(x1)), sd(log(x1))),
        from=log(min(x1)), to=log(max(x1)), add=TRUE, lty=2, lwd=2)
  
  #if a second time point is specified by 'by', add two more histograms + time labels
  if (!is.null(by)){
    
    hist.data <- hist(x2, freq=TRUE, col=col.orig, xlab=xlab, ylab='frequency', main='',
                      breaks=utils.exBar(x2,nbars))
    curve(sum(hist.data$count)*diff(hist.data$breaks)[1]*dnorm(x, mean(x2), sd(x2)),
          from=min(x2), to=max(x2), add=TRUE, lty=2, lwd=2)
    
    hist.data <- hist(log(x2), freq=TRUE, col=col.log, xlab=xlab.log, ylab='frequency', main='',
                      breaks=utils.exBar(log(x2),nbars))
    curve(sum(hist.data$count)*diff(hist.data$breaks)[1]*dnorm(x, mean(log(x2)), sd(log(x2))),
          from=log(min(x2)), to=log(max(x2)), add=TRUE, lty=2, lwd=2)
    
    mtext(levels(by)[1], side=3, line=1, outer=TRUE, cex=2)  
    mtext(levels(by)[2], side=3, line=-37, outer=TRUE, cex=2)
    
  }  
  #add title and subtitle
  mtext(title, side=3, line=oma3-3, outer=TRUE, font=2, cex=2.25)
  mtext('----- rescaled density of normal distribution with appropriate mean and sd', side=1,
        outer=TRUE, line=-2, cex=1.5)
  
  dev.off()
  
  #export size information to allow for convenient rescaling
  size <- c(height.png,14)
  names(size) <- c('height','width')
  return(size)
}


#' Function to create scatterplots, color-coded by subgroups
#' 
#' @author Sebastian Voss, Adam Skubala
#'
#' @param x numeric vector (Required)
#' @param y second numeric vector of the same length than 'x' (Required)
#' @param by character or factor vector of the same length than 'x' to specify
#'     by-groups used for color-coding (Default: Null)
#' @param col.by colors for the 'by'-groups(Default: 1:length(unique(by)))
#' @param name.by character, name of the 'by'-variable, used as legend-title
#' @param mar margins of the plot as in 'par()' (Default: c(6,4,1,1))
#' @param xlab character, labels for the x-axis (Default: '')
#' @param ylab character, labels for the y-axis (Default: '')
#' @param xlim Range of the x-axis (Default: Null)
#' @param ylim Range of the y-axis (Default: Null)
#' @param mark TODO: This parameter was not documented in original code (Default: Null)
#' @param xlog TRUE/FALSE is x log-transformed? if yes, axis annotation will account for that (Default: FALSE)
#' @param ylog TRUE/FALSE is y log-transformed? if yes, axis annotation will account for that (Default: FALSE)
#' @param pch symbol for the jitter points (default=1, i.e. the same symbol for all points),
#'     if at least as many symbols as the number of levels of the 'by.col' variable are given,
#'     symbol coding is used in addition to color coding. (Default: 1)
#' @param cex.pch size of the points (Default: 0.75) 
#' @param add.reg TRUE/FALSE if true, adds a regression line (Default: FALSE)
#' @param add.mean TODO: This parameter was not documented in original code (Default: FALSE)
#' @param add.leg TRUE/FALSE if true and 'by' is specified, a legend is added to the plot (Default: TRUE)
#' @param cex.mean TODO: This parameter was not documented in original code (Default: 1)
#' @param file.name name for the PNG output file (without '.png' ending) (Default: NULL)
#' @param png.out TRUE/FALSE if true, a png file is created as output, otherwise
#'     the plot is created in the R graphic device (Default: TRUE)
#' @param png.width width of PNG file in inches (Default: 6)
#' @param png.height height of PNG file in inches (Default: 5)
#' @param ... Additional parameters for par()
#'
#' @return PNG files with the scatterplot ('file.name.png')
#'     numeric vector with size of the PNG-file in inches for rescaling in RTF
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- rnorm(500)
#' y <- x^2+rnorm(500)
#' by <- factor(c(rep('male', 250),
#'                rep('female',250)),
#'              levels=c('male','female'))
#' 
#' general.scatter_by(x, y, by=by, col.by=c('skyblue','tomato'), name.by="Sex",
#'                    xlab="Sex", ylab="Biomarker", file.name="test", add.reg=TRUE,
#'                    mar=c(7,4,4,2), xlim=range(x), ylim=range(y), 
#'                    pch <- c(1, 4), xlog=TRUE, add.mean=TRUE)
#' }
general.scatter_by <- function(x, y, by=NULL, col.by=1:length(unique(by)), name.by='', mar=c(6,4,2,1),
                       xlab='', ylab='', xlim=NULL, ylim=NULL, mark=NULL,
                       xlog=FALSE, ylog=FALSE, pch=1, cex.pch=0.75, add.reg=FALSE, add.mean=FALSE, add.leg=TRUE,
                       cex.mean=1, file.name=NULL, png.out=TRUE, png.width=6, png.height=5, ...){
  
  #delete all observations with at least one NA
  data <- cbind.data.frame(x, y)
  if (!is.null(by)){
    data <- cbind.data.frame(data, factor(by))
  }
  data <- na.exclude(data)
  x <- data[,1]
  y <- data[,2]
  if (!is.null(by)){
    by <- data[,3]
  }
  
  #create scatterplot
  ##plot function
  create.plot <- function(){
    
    #if a by-variable is specified, a legend is added to the plot
    #-> bottom margin has to be bigger
    if (!is.null(by) & add.leg==TRUE){
      mar[1] <- mar[1]+0.5+0.25*length(levels(by))  
    }
    par(mar=mar, ...)
    #empty plot
    plot(0,0, type='n', xlim=xlim, ylim=ylim, main='', ylab='', xlab='', axes=FALSE)
    #color background of the plotting region
    rect(par('usr')[1], par('usr')[3],par('usr')[2],par('usr')[4], col='grey95')
    #add axis and box
    if (xlog==FALSE){
      axis(1)
    }else{
      axis(1, at=axTicks(1), labels=format(round(exp(axTicks(1)),2),nsmall=2))
    }
    if (ylog==FALSE){
      axis(2)
    }else{
      axis(2, at=axTicks(2), labels=format(round(exp(axTicks(2)),2),nsmall=2))
    }
    box()
    #color for points
    if (!is.null(by)){
      col.points <- col.by[as.numeric(by)]
    }else{
      col.points <- col.by[1]
    }
    #if at least as many symbols are given in 'pch' than the number of levels of the color variable,
    #use symbol coding in addition to color coding
    if (!is.null(by) & length(pch)>=length(levels(by))){
      pch.by <- pch[as.numeric(by)]
      #else use just the first symbol in 'pch'
    }else{
      pch.by <- pch[1] #for plot
      pch <- rep(pch[1], length(levels(by))) #for legend
    }
    
    #add vertical and horizontal line if mean for overall mean
    if (add.mean==TRUE & is.null(mark)){
      abline(v=mean(x), lty=2, col='grey66')
      abline(h=mean(y), lty=2, col='grey66')
    }else if(!is.null(mark)){
      abline(v=ifelse(xlog==TRUE, log(mark[1]), mark[1]), lty=2, col='grey66')
      abline(h=ifelse(ylog==TRUE, log(mark[2]), mark[2]), lty=2, col='grey66')
    }
    
    #add points
    points(x,y, col=col.points, cex=cex.pch, pch=pch.by)
    
    #add means (crosses)
    if (add.mean==TRUE){
      #points(mean(x), mean(y), col='grey66', cex=6.5*cex.mean, pch=3, lwd=3)
      if(!is.null(by)){
        mby.x <- as.vector(by(x, INDICES=by, FUN=mean))
        mby.y <- as.vector(by(y, INDICES=by, FUN=mean))
        points(mby.x, mby.y, col=col.by, cex=5*cex.mean, pch=3, lwd=2.5)
      }
    }
    
    #add regression line
    if (add.reg==TRUE){
      fit <- lm(y~x)
      x.str <- ifelse(xlog==TRUE, 'log(x)','x')
      y.str <- ifelse(ylog==TRUE, 'log(y)','y')
      reg.label1 <- paste0('Regression line: ',y.str,' = ',format(round(summary(fit)$coef[1,1],2), nsmall=2))
      reg.label2 <- ifelse(summary(fit)$coef[2,1]>=0,paste0(' + ',format(round(summary(fit)$coef[2,1],2), nsmall=2),' ',x.str),
                           paste0(' - ',format(round(-summary(fit)$coef[2,1],2), nsmall=2),' ',x.str))
      reg.label3 <- paste0('\n(Pearson cor. = ',format(round(cor(x,y,method='p'),2),nsmall=2),', Spearman cor. = ',format(round(cor(x,y, method='s'),2),nsmall=2),', n = ',nrow(data),')')
      reg.label <- paste0(reg.label1, reg.label2)
      
      mtext(reg.label, side=3, line=par('mar')[3]-1.5, cex=0.8*par('cex'))
      mtext(reg.label3, side=3, line=par('mar')[3]-2.5, cex=0.8*par('cex'))
      
      abline(fit, col='grey40', lwd=2.5, lty=3)
      
    }
    
    #add legend, if 'by' is specified
    if (!is.null(by) & add.leg==TRUE){
      x.leg <- (par('usr')[1]+par('usr')[2])/2
      y.leg <- grconvertY(0, from='ndc', to='user')
      l <- legend(x.leg, y.leg, legend=levels(as.factor(by)), pch=pch[1:length(by)],
                  col=col.by, xjust=0.5, yjust=0,
                  title=name.by, xpd=TRUE, horiz=FALSE, bty='n', cex=0.8*par('cex'), ncol=2)
      #add x-label 
      text((par('usr')[1]+par('usr')[2])/2, l$rect$top, paste(xlab,'\n\n',sep=''), xpd=TRUE, cex=1*par('cex'))
    }else{
      #add x-label 
      mtext(xlab, cex=1*par('cex'), side=1, line=par('mar')[1]-2)
    }
    #add y-label 
    mtext(ylab, cex=1*par('cex'), side=2, line=par('mar')[2]-1)
    
  }
  
  #PNG or standard R graphic device
  if (png.out==TRUE){
    
    png(paste0(file.name,'.png'), width=png.width, height=png.height, units='in', res=300)
    create.plot()
    dev.off()
    #return PNG size (useful for importing plot into an RTF document)
    size <- c(png.height, png.width)
    names(size) <- c('height','width')
    return(size)
    
  }else{
    
    create.plot()
    
  }
}

#' Custom geom for general.lineplot
#'
#' @keywords internal
#' 
#' @author Jens Hooge
#'
#' @seealso ggplot2::layer
#'
#' @param mapping Set of aesthetic mappings created by aes or aes_. 
#'                If specified and inherit.aes = TRUE (the default), 
#'                it is combined with the default mapping at the top 
#'                level of the plot. You must supply mapping if there 
#'                is no plot mapping.
#' @param data The data to be displayed in this layer. There are three options: 
#'             If NULL, the default, the data is inherited from the plot data as 
#'             specified in the call to ggplot. 
#'             A data.frame, or other object, will override the plot data. 
#'             All objects will be fortified to produce a data frame. See fortify 
#'             for which variables will be created. 
#'             A function will be called with a single argument, the plot data. 
#'             The return value must be a data.frame., and will be used as the layer data.
#' @param var.type The way the variability region is displayed (e.g. 'errorbars' or 'ribbon')
#' @param stat The statistical transformation to use on the data for this layer, as a string.
#' @param position Position adjustment, either as a string, or the result of a call to a position 
#'                 adjustment function. 
#' @param na.rm If FALSE, the default, missing values are removed with a warning. 
#'              If TRUE, missing values are silently removed.
#' @param show.legend logical. Should this layer be included in the legends? 
#'                    NA, the default, includes if any aesthetics are mapped. 
#'                    FALSE never includes, and TRUE always includes.
#' @param inherit.aes If FALSE, overrides the default aesthetics, rather than combining with them. 
#'                    This is most useful for helper functions that define both data and aesthetics 
#'                    and shouldn't inherit behaviour from the default plot specification, e.g. borders.
#' @param ... other arguments passed on to layer. These are often aesthetics, used to set an 
#'            aesthetic to a fixed value, like color = "red" or size = 3. They may also be parameters to the 
#'            paired geom/stat.
#'
#' @return ggplot2 layer
geom_variability <- function(mapping = NULL, data = NULL, var.type=c("errorbar", "ribbon"), stat = "identity", position = "identity",
                             na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...) {
  
  geom <- geom_errorbar(mapping, data = data,
                        stat = stat, position = position,
                        na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes, 
                        width=0.2, ...)
  if (var.type=="ribbon") {
    geom <- geom_ribbon(mapping, data = data,
                        stat = stat, position = position,
                        na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes, 
                        alpha=0.3, ...)
  }
  return(geom)
}


#' Function generate (non-)facetted lineplot
#' 
#' @author Jens Hooge
#'
#' @param data A data frame in long format (Required)
#' @param by Variable the plot should be facetted by (Default: NULL)
#' @param facet.rows Number of rows in which the facets should be displayed. (Default: 1)
#' @param upper.err Upper bound width on errorbars or ribbon variability region (Default: 1)
#' @param lower.err Lower bound width on errorbars or ribbon variability region (Default: 1)
#' @param var.type Type of how to display variability ("errorbar", "ribbon") (Default: "errorbar")
#' @param title Title of the plot (Default: NULL)
#' @param xlab x-axis label (Default: NULL)
#' @param ylab y-axis label (Default: NULL)
#' @param xlim Upper and lower limit (e.g. c(0, 10)) of x-axis (Default: NULL)
#' @param ylim Upper and lower limit (e.g. c(0, 10)) of y-axis (Default: NULL)
#' @param legend.position Position of the legend ("none", "left", "right", "top", "bottom") (Default: "left")
#' @param size Point Size (Default: 3)
#' @param x.log10trans Logical flag for log10 x-axis trasformation
#' @param y.log10trans Logical flag for log10 y-axis trasformation
#'
#' @return a ggplot2 plot
#' @export
#'
#' @examples
#' \dontrun{
#' library(reshape2)
#' library(ggplot2)
#' n=7
#' df <- data.frame(Dose1=exp((1/2)*1:n),
#'                  Dose2=exp((1/3)*1:n),
#'                  Dose3=exp((1/4)*1:n),
#'                  Dose4=exp((1/5)*1:n),
#'                  Dose5=exp((1/6)*1:n),
#'                  Dose6=exp((1/7)*1:n),
#'                  Visit=1:n)
#' df <- melt(df, id=c("Visit"))
#' colnames(df) <- c("Visit", "Dose", "Measure") 
#' 
#' ## Lineplot with Errorbars
#' general.lineplot(df, upper.err=1, lower.err=1, var.type="errorbar")
#' 
#' ## Lineplot with Errorbars and log transformed y axis
#' general.lineplot(df, upper.err=1, lower.err=1, var.type="errorbar", y.log10trans = T)
#' 
#' ## Facetted Lineplot
#' general.lineplot(df, by="Dose", facet.row=3, upper.err=1, lower.err=1, 
#'                  var.type="errorbar", y.log10trans = T)
#'                  
#' ## Lineplot with Error Ribbon 
#' general.lineplot(df, by="Dose", upper.err=1, lower.err=1, var.type="ribbon")
#' 
#' ## Lineplot with Error Ribbon and log transformed y-axis
#' general.lineplot(df, by="Dose", upper.err=1, lower.err=1, var.type="ribbon", y.log10trans = T)
#' }
general.lineplot <- function(data, by=NULL, facet.rows=1,
                             upper.err, lower.err, 
                             var.type=c("errorbar", "ribbon"), 
                             title=NULL, xlab=NULL, ylab=NULL,
                             xlim=NULL, ylim=NULL,
                             legend.position="right", size=3,
                             x.log10trans=FALSE, y.log10trans=FALSE) {
  
  # Define the top and bottom of the errorbars
  fig <- ggplot(data, aes(data$Visit, data$Measure)) +
    
    geom_point(aes(col=data$Dose), size=size) +
    geom_line(aes(col=data$Dose)) + 
    geom_variability(aes(ymax = data$Measure+upper.err, ymin=data$Measure-lower.err, 
                         col=data$Dose, fill=data$Dose), var.type = var.type) +
    ggtitle(label = title) +
    xlab(xlab) + 
    ylab(ylab) +
    theme_bw() +
    theme(legend.position=legend.position)
  
  
  if(!is.null(by)) {
    stopifnot(by %in% colnames(data))
    facet <- as.formula(paste0("~", by))
    fig <- fig + facet_wrap(facet, scales = "free", nrow=facet.rows)
  }
  if(!is.null(xlim)) {
    fig <- fig + xlim(xlim)
  }
  if(!is.null(ylim)) {
    fig <- fig + ylim(ylim)
  }
  if(x.log10trans) {
    fig <- fig + scale_x_log10()
  }
  if(y.log10trans) {
    fig <- fig + scale_y_log10()
  }
  return(fig)
}



#' Function to generate a facetted scatter plot
#' 
#' @author Jens Hooge
#'
#' @keywords internal
#'
#' @param data A dataframe with at least three columns (x,y and a facetting variable) (Required)
#' @param by Grouping variable (Required)
#' @param fun Function for statistical aggregation function (Default: "mean") 
#' @param smooth.fun Smoothing method (function) to use, eg. "lm", "glm", "gam", "loess", "rlm". (Default: "auto")
#'                   For method = "auto" the smoothing method is chosen based on the size of the largest group 
#'                   (across all panels). loess is used for less than 1,000 observations; otherwise gam is used with 
#'                   formula = y ~ s(x, bs = "cs"). Somewhat anecdotally, loess gives a better appearance, 
#'                   but is O(n^2) in memory, so does not work for larger datasets.
#' @param smooth.formula Formula to use in smoothing function, eg. y ~ x, y ~ poly(x, 2), y ~ log(x) (Default: y ~ x)
#' @param density Display density contour (Default: FALSE)
#' @param title Plot Title (Default: NULL)
#' @param xlab x-axis label (Default: NULL)
#' @param ylab y-axis label (Default: NULL)
#' @param legend.pos (Default: "none")
#' 
#' @return ggplot2 figure
general.scatter.facetted <- function(data, by, fun="mean", smooth.fun="auto", smooth.formula=y ~ x,
                                     density=FALSE,
                                     title=NULL, xlab=NULL, ylab=NULL,
                                     legend.pos="none") {
  xagg <- aggregate(cbind(data$x, data$y), by=list(data[, by]), fun)
  colnames(xagg) <- c(by, paste0("x.", fun), paste0("y.", fun))
  data <- left_join(data, xagg, by)
  
  facet <- paste0("~", by)
  fig <- ggplot(data, aes_string(x=data$x, y=data$y, color=by)) + 
      geom_point() +
      geom_point(aes_string(x=paste0("x.", fun), y=paste0("y.", fun)), 
                 shape=3, size=18, show.legend = FALSE) +
      geom_smooth(method=smooth.fun, formula = smooth.formula) +
      facet_grid(facet) +
      ggtitle(label = title) +
      xlab(xlab) + 
      ylab(ylab) + 
      theme_bw() +
      theme(legend.position = legend.pos)
  
  if(density) {
    fig <- fig + stat_density2d(aes(fill=..level..,alpha=..level..),
                                geom='polygon',colour='black') +
      scale_fill_continuous(low="green", high="red") +
      guides(alpha="none")
  }
  
  return(fig)
}

#' Function to generate a simple scatter plot
#' 
#' @author Jens Hooge
#'
#' @keywords internal
#'
#' @param data A dataframe with at least three columns (x,y and a facetting variable) (Required)
#' @param by Grouping variable (Default: NULL)
#' @param fun Function for statistical aggregation function (Default: "mean") 
#' @param smooth.fun Smoothing method (function) to use, eg. "lm", "glm", "gam", "loess", "rlm". (Default: "auto")
#'                   For method = "auto" the smoothing method is chosen based on the size of the largest group 
#'                   (across all panels). loess is used for less than 1,000 observations; otherwise gam is used with 
#'                   formula = y ~ s(x, bs = "cs"). Somewhat anecdotally, loess gives a better appearance, 
#'                   but is O(n^2) in memory, so does not work for larger datasets.
#' @param smooth.formula Formula to use in smoothing function, eg. y ~ x, y ~ poly(x, 2), y ~ log(x) (Default: y ~ x)
#' @param density Display density contour (Default: FALSE)
#' @param title Plot Title (Default: NULL)
#' @param xlab x-axis label (Default: NULL)
#' @param ylab y-axis label (Default: NULL)
#' @param legend.pos (Default: "none")
#' 
#' @return ggplot2 figure
general.scatter.simple <- function(data, by=NULL, fun="mean", smooth.fun="auto", smooth.formula=y ~ x,
                                   density=FALSE,
                                   title=NULL, xlab=NULL, ylab=NULL,
                                   legend.pos="none") {
  xagg <- apply(data[, c("x", "y")], 2, fun, na.rm=T)
  xagg <- data.frame(rep(xagg[1], nrow(data)),
                     rep(xagg[2], nrow(data)))
  colnames(xagg) <- c(paste0("x.", fun), paste0("y.", fun))
  
  data <- cbind(data, xagg)
  
  fig <- ggplot(data, aes_string(x=data$x, y=data$y)) + 
    geom_point() +
    geom_point(aes_string(x=paste0("x.", fun), y=paste0("y.", fun)), 
               shape=3, size=18, show.legend = FALSE) +
    geom_smooth(method=smooth.fun, formula = smooth.formula) +
    ggtitle(label = title) +
    xlab(xlab) + 
    ylab(ylab) + 
    theme_bw() +
    theme(legend.position = legend.pos)
  
  if(density) {
    fig <- fig + stat_density2d(aes(fill=..level..,alpha=..level..),
                                geom='polygon',colour='black') + 
      scale_fill_continuous(low="green", high="red") +
      guides(alpha="none")
  }
  
  return(fig)
}

#' Function factory to generate a scatter plots
#'
#' @author Jens Hooge
#'
#' @description 
#' Depending on the grouping variable 'by', this function factory generates
#' either an ungrouped (by=NULL) or a grouped scatterplot (by=<grouping_variable>).
#' This function provides the option for different regression fits, based on an inbuilt
#' smooth function or a custom model, defined by smooth.formula.
#'
#' @param data A dataframe with at least three columns (x,y and a facetting variable) (Required)
#' @param by Grouping variable (Required)
#' @param fun Function for statistical aggregation function (Default: "mean") 
#' @param smooth.fun Smoothing method (function) to use, eg. "lm", "glm", "gam", "loess", "rlm". (Default: "auto")
#'                   For method = "auto" the smoothing method is chosen based on the size of the largest group 
#'                   (across all panels). loess is used for less than 1,000 observations; otherwise gam is used with 
#'                   formula = y ~ s(x, bs = "cs"). Somewhat anecdotally, loess gives a better appearance, 
#'                   but is O(n^2) in memory, so does not work for larger datasets.
#' @param smooth.formula Formula to use in smoothing function, eg. y ~ x, y ~ poly(x, 2), y ~ log(x) (Default: y ~ x)
#' @param density Display density contour (Default: FALSE)
#' @param title Plot Title (Default: NULL)
#' @param xlab x-axis label (Default: NULL)
#' @param ylab y-axis label (Default: NULL)
#' @param legend.pos (Default: "none")
#'
#' @return ggplot2 figure
#' @export
#' 
#' @examples
#' \dontrun{
#' x <- rnorm(500)
#' y <- rnorm(500)+x^2
#' df <- data.frame(x=x, y=y, 
#'                  A=sample(c("a1", "a2"), 500, replace = TRUE),
#'                  B=sample(c("b1", "b2"), 500, replace = TRUE))
#' 
#' ## A simple ungrouped scatter plot with a loess fit
#' general.scatter(df)
#' 
#' ## A grouped scatter plot with a loess fit
#' general.scatter(df, by="A")
#' 
#' ## An ungrouped scatter plot and a marker indicating the data variance
#' general.scatter(df, fun="var")
#' 
#' ## Scatterplots with different regression fits
#' # Linear Fit
#' general.scatter(df, by="A", smooth.fun="lm")
#' 
#' ## Robust Linear Fit
#' library(MASS)
#' general.scatter(df, by="A", smooth.fun="rlm")
#' 
#' ## User defined polynomial fit
#' # Degree = 1 resulting in a linear fit
#' general.scatter(df, smooth.fun="glm", smooth.formula=y ~ poly(x, 1)) 
#' # Degree = 2 resulting in a quadratic fit
#' general.scatter(df, smooth.fun="glm", smooth.formula=y ~ poly(x, 2))
#' 
#' ## Equivalently for grouped scatterplots
#' # Degree = 1 resulting in a linear fit
#' general.scatter(df, by="A", smooth.fun="glm", smooth.formula=y ~ poly(x, 1)) 
#' # Degree = 2 resulting in a quadratic fit
#' general.scatter(df, by="B", smooth.fun="glm", smooth.formula=y ~ poly(x, 2))
#' 
#' ## Including a density contour
#' general.scatter(df, by="B", smooth.fun="glm", smooth.formula=y ~ poly(x, 2),
#'                 density=T, legend.pos="bottom")
#' }
general.scatter <- function(data, by=NULL, fun="mean", smooth.fun="auto", smooth.formula=y ~ x,
                            density=FALSE,
                            title=NULL, xlab=NULL, ylab=NULL,
                            legend.pos="none") {
  
  args <- list(data=data, by=by, fun=fun, 
               density=density,
               smooth.fun=smooth.fun, smooth.formula=smooth.formula,
               title=title, xlab=xlab, ylab=ylab,
               legend.pos=legend.pos)
  
  fig <- NULL
  
  if(is.null(by)) {
    fig <- do.call(general.scatter.simple, args)
  } else {
    fig <- do.call(general.scatter.facetted, args)
  }
  
  return(fig)
}