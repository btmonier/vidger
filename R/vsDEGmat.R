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
#' @param data output generated from calling the main routines of either
#'  `cuffdiff`, `DESeq2`, or `edgeR` analyses. For `cuffdiff`, this will be a 
#'  `*_exp.diff` file. For `DESeq2`, this will be a generated object of class 
#'  `DESeqDataSet`. For `edgeR`, this will be a generated object of class 
#'  `DGEList`. 
#' @param padj a user defined adjusted p-value cutoff point. Defaults to 
#'  `0.05`.
#' @param d.factor a specified factor; for use with `DESeq2` objects only.
#'  This input equates to the first parameter for the contrast argument when 
#'  invoking the `results()` function in `DESeq2`. Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the 
#'  data. Must be either `cuffdiff`, `deseq`, or `edger`. `cuffdiff` must be 
#'  used with `cuffdiff` data; `deseq` must be used for `DESeq2` output; 
#'  `edgeR` must be used with `edgeR` data. See the `data` parameter for 
#'  further details.
#' @param title display the main title of plot. Logical; defaults to `TRUE`.
#'  If set to `FALSE`, no title will display in plot.
#' @param legend display legend of plot. Logical; defaults to `TRUE`.
#'  If set to `FALSE`, no legend will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to `TRUE`.
#'  If set to `FALSE`, no axis lines will display in plot.
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

vsDEGMatrix <- function(
    data, padj = 0.05, d.factor = NULL, 
    type = c("cuffdiff", "deseq", "edger"), title = TRUE, legend = TRUE, 
    grid = TRUE
) {
    if (missing(type) || !type %in% c("cuffdiff", "deseq", "edger")) {
        stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
    }
    
    type <- match.arg(type)
    if(type == 'cuffdiff'){
        dat <- .getCuffDEGMat(data, padj)
    } else if (type == 'deseq') {
        dat <- .getDeseqDEGMat(data, d.factor, padj)
    } else if (type == 'edger') {
        dat <- .getEdgeDEGMat(data, padj)
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
        stat_sum(
            aes(fill = ..n..),
            color = "black", size = 0.3, geom = "tile"
        ) + 
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
