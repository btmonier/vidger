#' @title 
#' Scatter plot matrix of log10(FPKM or CPM) values
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @import ggplot2
#' @import GGally
#' @import tidyr
#' @import dplyr
#'
#' @description
#' This function will generate a matrix of scatterplots for all possible treatment combinations  
#' with additional distribution info.
#' 
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param d.factor a specified factor; for use with DESeq2 objects only. Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the data. Must be either `cuffdiff`, `deseq`, or `edgeR`.
#' @param comp treatments you would like to compare in the form of a vector. If no parameter is specified, all possible treatment comparisons will be made. Defaults to `NULL`.
#' @param title show title of plot. Defaults to `TRUE`.
#' @param grid show major and minor axis lines. Defaults to `TRUE`.
#' @param man.title a manually specified title at the authors discretion. Defaults to `NULL`.
#' 
#' @export
#' 
#' @examples
#' # Cuffdiff example
#' load('df.cuff.RData')
#' vsScatterMatrix(data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
#'                 comp = NULL, title = TRUE, grid = TRUE, 
#'                 man.title = 'Example title')
#' 
#' # DESeq2 example
#' load('df.deseq.RData')
#' vsScatterMatrix(data = df.deseq, d.factor = 'cell', type = 'deseq',
#'                 comp = NULL, title = TRUE, grid = FALSE, man.title = NULL)
#' 
#' # edgeR example
#' load('df.edger.RData')
#' vsScatterMatrix(data = df.edger, d.factor = NULL, type = 'edger', 
#'                 comp = c('WM', 'MM'), title = TRUE, grid = TRUE, 
#'                 man.title = NULL)

vsScatterMatrix <- function(data, d.factor = NULL, type, comp = NULL, 
                            title = TRUE, grid = TRUE, man.title = NULL) {
  if (missing(type)) {
    stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
  }
  
  if(type == 'cuffdiff'){
    dat <- getCuffScatterMatrix(data)
  } else if (type == 'deseq') {
    dat <- getDeseqScatterMatrix(data, d.factor)
  } else if (type == 'edger') {
    dat <- getEdgeScatterMatrix(data)
  } else {
    stop('Please enter correct analysis type.')
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
    aes.xlab <- bquote('log'['10'] ~ '(FPKM)')
    aes.ylab <- bquote('log'['10'] ~ '(FPKM)')
  }
  
  GGally::ggpairs(log10(dat + 1), title = m.title, xlab = aes.xlab, 
                  ylab = aes.ylab, columns = comp,
                  lower = list(continuous = wrap('points', size = 0.5))) + grid
}
