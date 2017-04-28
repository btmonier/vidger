#' @title 
#' MA plot from mean expression and log fold changes from different analytical objects
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to extract necessary results-based data from 
#' different output files to create a MA plot (i.e. a scatter plot) of log2 
#' fold changes versus normalized mean counts while implementing ggplot2 
#' aesthetics. 
#' 
#' @param x treatment `x` for comparison (log2(y/x)).
#' @param y treatment `y` for comparison (log2(y/x)).
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param d.factor a specified factor; for use with DESeq2 objects only. Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the data. Must be either `cuffdiff`, `deseq`, or `edgeR`.
#' @param padj a user defined adjusted p-value cutoff point. Defaults to `0.1`.
#' @param y.lim set manual limits to the y axis. Defaults to `NULL`.
#' @param lfc log fold change level for setting conditonals. If no user input is added (`NULL`), value defaults to `1`.
#' @param title show title of plot. Defaults to `TRUE`.
#' @param legend shows legend of plot. Defaults to `TRUE`.
#' @param grid show major and minor axis lines. Defaults to `TRUE`.
#' @param data.return returns data output of plot if set to `TRUE`. Defaults to `FASLSE`.
#' 
#' @export
#' 
#' @examples
#' # Cuffdiff example
#' load('df.cuff.RData')
#' vsMAPlot(x = 'hESC', y = 'iPS', data = df.cuff, d.factor = NULL, 
#'          type = 'cuffdiff', padj = 0.05, y.lim = NULL, lfc = 2, 
#'          title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE)
#' 
#' # DESeq2 example
#' load('df.deseq.RData')
#' vsMAPlot(x = 'trt', y = 'untrt', data = df.deseq, d.factor = 'dex', 
#'          type = 'deseq', padj = 0.05, y.lim = NULL, lfc = NULL, 
#'          title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE)
#' 
#' # edgeR example
#' load('df.edger.RData')
#' vsMAPlot(x = 'WM', y = 'MM', data = df.edger, d.factor = NULL, 
#'          type = 'edger', padj = 0.1, y.lim = NULL, lfc = 2, 
#'          title = FALSE, legend = TRUE, grid = TRUE, data.return = FALSE)
#'                 
#' # Extract data frame from visualization
#' load('df.cuff.RData')
#' tmp <- vsMAPlot(x = 'hESC', y = 'iPS', data = df.cuff, 
#'                 d.factor = NULL, type = 'cuffdiff', padj = 0.05, 
#'                 y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
#'                 data.return = TRUE)
#' df.ma <- tmp[[1]]
#' head(df.ma)

vsMAPlot <- function(x, y, data, d.factor = NULL, type, padj = 0.1, y.lim = NULL,
                     lfc = NULL, title = TRUE, legend = TRUE, grid = TRUE,
                     data.return = FALSE){
  if (missing(type)) {
    stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
  }
  if(type == 'cuffdiff'){
    dat <- getCuffMA(x, y, data)
  } else if (type == 'deseq') {
    dat <- getDeseqMA(x, y, data, d.factor)
  } else if (type == 'edger') {
    dat <- getEdgeMA(x, y, data)
  } else {
    stop('Please enter correct analysis type.')
  }
  if (!isTRUE(title)) {
    m.lab <- NULL
  } else {
    m.lab  <- ggtitle(paste(y, 'vs.', x))
  }
  if (!isTRUE(legend)) {
    leg <- theme(legend.position = 'none')
  } else {
    leg <- guides(colour = guide_legend(override.aes = list(size = 3)),
                  shape = guide_legend(override.aes = list(size = 3)))
  }
  if (!isTRUE(grid)) {
    grid <- theme_classic()
  } else {
    grid <- theme_bw()
  }
  
  #' Configure data
  dat$isDE <- ifelse(dat$padj <= padj, TRUE, FALSE)
  py <- dat$M
  
  #' Conditional plot components
  if (is.null(y.lim)) {
    y.lim = c(-1, 1) * quantile(abs(py[is.finite(py)]), probs = 0.99) * 0.8
  }
  if (is.null(lfc)) {
    lfc = 1
  }
  
  dat <- ma.ranker(dat, padj, lfc, y.lim)
  
  #' See ranker functions
  tmp.size <- ma.out.ranker(py, y.lim[2])
  tmp.col <- ma.col.ranker(dat$isDE, py, lfc)
  tmp.shp <- ma.shp.ranker(py, y.lim)
  
  #' See counting function
  tmp.cnt <- ma.col.counter(dat, lfc)
  b <- tmp.cnt[[1]]
  g <- tmp.cnt[[2]]
  
  #' See component functions
  comp1 <- ma.comp1(y.lim, padj, lfc, b, g)
  point <- geom_point(
    alpha = 0.7, 
    aes(color = tmp.col, shape = tmp.shp, size = tmp.size)
  )
  comp2 <- ma.comp2(comp1[[4]], comp1[[6]], comp1[[5]], comp1[[1]], comp1[[2]], 
                    comp1[[3]])
  
  #' Plot layers
  tmp.plot <- ggplot(dat, aes(x = A, y = pmax(y.lim[1], pmin(y.lim[2], py)))) +
    point + comp2$color + comp2$shape + comp1$hline1 + comp1$hline2 + 
    comp1$hline3 + comp1$x.lab + comp1$y.lab + m.lab + ylim(y.lim) + 
    comp2$size + grid + leg
  
  if (isTRUE(data.return)) {
    plot.l <- list(data = dat, plot = tmp.plot)
  } else {
    print(tmp.plot)
  }
} 