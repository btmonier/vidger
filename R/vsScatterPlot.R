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
#' @param x treatment `x` for comparison (log2(x/control)). This will be a
#'  factor level in your data.
#' @param y treatment `y` for comparison (log2(y/control)). This will be a 
#'  factor level in your data.
#' @param data output generated from calling the main routines of either
#'  `cuffdiff`, `DESeq2`, or `edgeR` analyses. For `cuffdiff`, this will be a 
#'  `*_exp.diff` file. For `DESeq2`, this will be a generated object of class 
#'  `DESeqDataSet`. For `edgeR`, this will be a generated object of class 
#'  `DGEList`.
#' @param d.factor a specified factor; for use with `DESeq2` objects only.
#'  This input equates to the first parameter for the contrast argument when 
#'  invoking the `results()` function in `DESeq2`. Defaults to `NULL`.
#' @param type an analysis classifier to tell the function how to process the 
#'  data. Must be either `cuffdiff`, `deseq`, or `edger`. `cuffdiff` must be 
#'  used with `cuffdiff` data; `deseq` must be used for `DESeq2` output; 
#'  `edgeR` must be used with `edgeR` data. See the `data` parameter for 
#'  further details.
#' @param title display the main title of plot. Logical; defaults to `TRUE`.
#'  If set to `FALSE`, no title will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to `TRUE`.
#'  If set to `FALSE`, no axis lines will display in plot.
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
    x, y, data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"), 
    title = TRUE, grid = TRUE
) {
    if (missing(type) || !type %in% c("cuffdiff", "deseq", "edger")) {
        stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
    }
    
    type <- match.arg(type)
    if(type == 'cuffdiff'){
        dat <- .getCuffScatter(x, y, data)
    } else if (type == 'deseq') {
        dat <- .getDeseqScatter(x, y, data, d.factor)
    } else if (type == 'edger') {
        dat <- .getEdgeScatter(x, y, data)
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
