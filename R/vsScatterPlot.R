#' @title 
#' Scatter plot of log10(FPKM or CPM) values
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to visualize comparisons of log10 values of either 
#' FPKM or CPM measurements of two treatments depending on analytical type.
#' 
#' @param x treatment for the x-axis of the plot
#' @param y treatment for the y-axis of the plot
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param d.factor a specified factor; for use with DESeq2 objects only.
#'  Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the 
#'  data. Must be either `cuffdiff`, `deseq`, or `edgeR`.
#' @param title show title of plot. Defaults to `TRUE`.
#' @param grid show major and minor axis lines. Defaults to `TRUE`.
#' 
#' @return An object created by \code{ggplot}
#' 
#' @export
#' 
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsScatterPlot(
#'  x = 'hESC', y = 'iPS', data = df.cuff, d.factor = NULL, 
#'  type = 'cuffdiff', title = TRUE, grid = TRUE
#' )
#' 
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)               
#' vsScatterPlot(
#'  x = 'treated_paired.end', y = 'untreated_paired.end', 
#'  data = df.deseq, d.factor = 'condition', type = 'deseq', 
#'  title = TRUE, grid = TRUE
#' )
#' 
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsScatterPlot(
#'  x = 'WW', y = 'WM', data = df.edger, d.factor = NULL, 
#'  type = 'edger', title = TRUE, grid = TRUE
#' )

vsScatterPlot <- function(
    x, y, data, d.factor = NULL, type, title = TRUE,
    grid = TRUE
) {
    if (missing(type)) {
        stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
    }
    if(type == 'cuffdiff'){
        dat <- .getCuffScatter(x, y, data)
    } else if (type == 'deseq') {
        dat <- .getDeseqScatter(x, y, data, d.factor)
    } else if (type == 'edger') {
        dat <- .getEdgeScatter(x, y, data)
    } else {
        stop('Please enter correct analysis type.')
    }
    if (!isTRUE(title)) {
        m.title <- NULL
    } else {
        m.title <- ggtitle(paste(y, 'vs.', x))
    }
    if (!isTRUE(grid)) {
        grid <- theme_classic()
    } else {
        grid <- theme_bw()
    }
    if (type == 'edger') {
        aes.xlab <- bquote('log'['10'] ~ '(CPM) -' ~ .(x))
        aes.ylab <- bquote('log'['10'] ~ '(CPM) -' ~ .(y))
    } else {
        aes.xlab <- bquote('log'['10'] ~ '(FPM) -' ~ .(x))
        aes.ylab <- bquote('log'['10'] ~ '(FPM) -' ~ .(y))
    }
    
    tmp.plot <- ggplot(dat, aes(x = log10(x + 1), y = log10(y + 1))) +
        geom_point(size = 1) +
        xlab(aes.xlab) +
        ylab(aes.ylab) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        grid + m.title
    print(tmp.plot)
}
