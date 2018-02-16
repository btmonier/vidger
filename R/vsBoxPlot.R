#' @title 
#' Box plot for log10(FPKM or CPM) distributions
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to extract necessary results-based data from analytical objects
#' to create a box plot comparing log10(FPKM or CPM) distributions for experimental treatments.
#' 
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param d.factor a specified factor; for use with DESeq2 objects only. Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the data. Must be either `cuffdiff`, `deseq`, or `edgeR`.
#' @param title show title of plot. Defaults to `TRUE`.
#' @param legend shows legend of plot. Defaults to `TRUE`.
#' @param grid show major and minor axis lines. Defaults to `TRUE`.
#' 
#' @export
#' 
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsBoxPlot(data = df.cuff, d.factor = NULL, type = 'cuffdiff', title = TRUE,
#'           legend = TRUE, grid = TRUE)
#' 
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsBoxPlot(data = df.deseq, d.factor = 'cell', type = 'deseq', title = TRUE,
#'           legend = TRUE, grid = TRUE)
#' 
#' # edgeR example
#' data("df.deseq")
#' require(edgeR)
#' vsBoxPlot(data = df.edger, d.factor = NULL, type = 'edger', title = TRUE,
#'           legend = TRUE, grid = TRUE)

vsBoxPlot <- function(data, d.factor = NULL, type, title = TRUE, legend = TRUE,
                      grid = TRUE){
  if (missing(type)) {
    stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
  }
  if(type == 'cuffdiff'){
    dat <- .getCuffBox(data)
  } else if (type == 'deseq') {
    dat <- .getDeseqBox(data, d.factor)
  } else if (type == 'edger') {
    dat <- .getEdgeBox(data)
  } else {
    stop('Please enter correct analysis type.')
  }
  if (!isTRUE(title)) {
    m.title <- NULL
  } else if (isTRUE(title) & type == 'edger') {
    m.title <- ggtitle('CPM distribution')
  } else {
    m.title <- ggtitle('FPM distribution')
  }
  if (!isTRUE(legend)) {
    leg <- guides(fill = FALSE)
  } else {
    leg <- NULL
  }
  if (!isTRUE(grid)) {
    grid <- theme_classic()
  } else {
    grid <- theme_bw()
  }
  if (type == 'edger') {
    y.lab <- expression(paste('log'['10'], ' (CPM)'))
  } else {
    y.lab <- expression(paste('log'['10'], ' (FPM)'))
  }
  key <- value <- NULL
  tmp.plot <- ggplot(dat, aes(x = key, y = log10(value + 1), fill = key)) +
    geom_boxplot() +
    xlab('Condition') +
    ylab(y.lab) +
    guides(fill = guide_legend(title = 'Condition')) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    leg + m.title + grid
  print(tmp.plot)
}