#' @title 
#' Differential gene expression matrix
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to visualize the number of differentially 
#'  expressed genes 
#' (DEG) at a given #' adjusted P-value for each experimental treatment level.
#'  Higher color intensity correlates to a higher number of DEGs. 
#'  
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param padj a user defined adjusted p-value cutoff point. Defaults to 
#'  `0.05`.
#' @param d.factor a specified factor; for use with DESeq2 objects only. 
#'  Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the
#'  data. Must be either `cuffdiff`, `deseq`, or `edgeR`.
#' @param title show title of plot. Defaults to `TRUE`.
#' @param legend show legend of plot. Defaults to `TRUE`.
#' @param grid show major and minor axis lines. Defaults to `TRUE`.
#' 
#' @return An object created by \code{ggplot}
#' 
#' @export
#' 
#' @examples
#' # cuffdiff example
#' data("df.cuff")
#' vsDEGMatrix(
#'  df.cuff, padj = 0.05, d.factor = NULL, type = 'cuffdiff', 
#'  title = TRUE, legend = TRUE, grid = TRUE
#' )
#' 
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsDEGMatrix(
#'  df.deseq, padj = 0.05, d.factor = 'condition', type = 'deseq', 
#'  title = TRUE, legend = TRUE, grid = TRUE
#' )
#' 
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsDEGMatrix(
#'  df.edger, padj = 0.05, d.factor = NULL, type = 'edger', 
#'  title = TRUE, legend = TRUE, grid = TRUE
#' )

vsDEGMatrix <- function(data, padj = 0.05, d.factor = NULL, type, title = TRUE,
                                                legend = TRUE, grid = TRUE) {
    if (missing(type)) {
        stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
    }
    if(type == 'cuffdiff'){
        dat <- .getCuffDEGMat(data, padj)
    } else if (type == 'deseq') {
        dat <- .getDeseqDEGMat(data, d.factor, padj)
    } else if (type == 'edger') {
        dat <- .getEdgeDEGMat(data, padj)
    } else {
        stop('Please enter correct analysis type.')
    }
    if (!isTRUE(legend)) {
        leg <- guides(fill = FALSE)
    } else {
        leg <- NULL
    }
    if (!isTRUE(title)) {
        m.lab <- NULL
    } else {
        m.lab  <- ggtitle(
            bquote('Significant transcripts at ' ~ alpha ~ ' = ' ~ .(padj))
        ) 
    }
    if (!isTRUE(grid)) {
        grid <- theme_classic()
    } else {
        grid <- theme_bw()
    }
    x <- y <- ..n.. <- NULL
    tmp.plot <- ggplot(dat, aes(x = x, y = y)) +
        stat_sum(aes(fill = ..n..), 
                         color = "black", size = 0.3, geom = "tile") + 
        scale_fill_continuous(
            low = "white", high = "royalblue3", 
            name = 'Number of \ntranscripts'
        ) + 
        expand_limits(fill = 0) + leg +
        stat_sum(
            aes(label = ..n..), geom = "text", size = 5, show.legend = FALSE
        ) +
        m.lab + grid + xlab('') + ylab('')
    print(tmp.plot)
}