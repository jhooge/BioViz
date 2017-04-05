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
#' general.volcano_plot(iPSC)
#' general.volcano_plot(iPSC, groups=c("A", "B"))
#' general.volcano_plot(iPSC, top_labeled=10, xcutoff=c(-log2(2), log2(4)))
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
#' @author Sebastian Voss
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
#' @param add.legend 
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


# n <- 500
# df <- data.frame(A=sample(c("a1", "a2", "a3"), n, replace=TRUE),
#                  B=sample(c("b1", "b2"), n, replace=TRUE),
#                  C=sample(c("c1", "c2", "c3"), n, replace=TRUE),
#                  D=sample(c("d1", "d2", "d3"), n, replace=TRUE),
#                  value=rnorm(n, 100, 1))
# 
# general.box_plot_facetted(df, group_by1 = "A", group_by2 = "B",
#                   col.jitter_by = "A", shape_jitter_by = "C", 
#                   rotate=TRUE)
