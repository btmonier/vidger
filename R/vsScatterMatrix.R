#' @title 
#' Scatter plot matrix of log10(FPKM or CPM) values
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function will generate a matrix of scatterplots for all possible
#' treatment combinations with additional distribution info.
#' 
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
#' @param comp treatments you would like to compare in the form of a vector. 
#'  If no parameter is specified, all possible treatment comparisons will be 
#'  made. Defaults to `NULL`.
#' @param title display the main title of plot. Logical; defaults to `TRUE`.
#'  If set to `FALSE`, no title will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to `TRUE`.
#'  If set to `FALSE`, no axis lines will display in plot.
#' @param man.title a manually specified title at the authors discretion. 
#'  Defaults to `NULL`.
#' 
#' @return An object created by \code{ggplot}
#' 
#' @export
#' 
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsScatterMatrix(
#'  data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
#'  comp = NULL, title = TRUE, grid = TRUE, 
#'  man.title = 'Example title'
#' )
#' 
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsScatterMatrix(
#'  data = df.deseq, d.factor = 'condition', type = 'deseq',
#'  comp = NULL, title = TRUE, grid = FALSE, man.title = NULL
#' )
#' 
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsScatterMatrix(
#'  data = df.edger, d.factor = NULL, type = 'edger', 
#'  comp = c('WM', 'MM'), title = TRUE, grid = TRUE, 
#'  man.title = NULL
#' )

vsScatterMatrix <- function(
    data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"), 
    comp = NULL, title = TRUE, grid = TRUE, man.title = NULL
) {
    if (missing(type) || !type %in% c("cuffdiff", "deseq", "edger")) {
        stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
    }
    
    type <- match.arg(type)
    if(type == 'cuffdiff'){
        dat <- .getCuffScatterMatrix(data)
    } else if (type == 'deseq') {
        dat <- .getDeseqScatterMatrix(data, d.factor)
    } else if (type == 'edger') {
        dat <- .getEdgeScatterMatrix(data)
    }
    
    if(!is.null(comp)) {
        comp <- comp
    } else {
        comp <- 1:ncol(dat)
    }
    
    if (!isTRUE(title)) {
        m.title <- NULL
    } else if (isTRUE(title) & !is.null(man.title)) {
        m.title <- ggtitle(man.title)
    } else if (isTRUE(title) & type == 'edger') {
        m.title <- 'CPM Comparisons'
    } else {
        m.title <- 'FPKM Comparisons'
    }
    if (!isTRUE(grid)) {
        grid <- theme_classic()
    } else {
        grid <- theme_bw()
    }
    if (type == 'edger') {
        aes.xlab <- bquote('log'['10'] ~ '(CPM)')
        aes.ylab <- bquote('log'['10'] ~ '(CPM)')
    } else {
        aes.xlab <- bquote('log'['10'] ~ '(FPM)')
        aes.ylab <- bquote('log'['10'] ~ '(FPM)')
    }
    
    GGally::ggpairs(
        log10(dat + 1), 
        title = m.title, 
        xlab = aes.xlab, 
        ylab = aes.ylab, 
        columns = comp,
        lower = list(continuous = wrap('points', size = 0.5))
    ) + grid
}
