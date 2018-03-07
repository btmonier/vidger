#' @title 
#' Box plot for log10(FPKM or CPM) distributions
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to extract necessary results-based data from
#'  analytical objects
#' to create a box plot comparing log10(FPKM or CPM) distributions for 
#'  experimental treatments.
#' 
#' @param data output generated from calling the main routines of either
#'  `cuffdiff`, `DESeq2`, or `edgeR` analyses. For `cuffdiff`, this will be a 
#'  `*_exp.diff` file. For `DESeq2`, this will be a generated object of class 
#'  `DESeqDataSet`. For `edgeR`, this will be a generated object of class 
#'  `DGEList`. 
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
#' # Cuffdiff example
#' data("df.cuff")
#' vsBoxPlot(
#'  data = df.cuff, d.factor = NULL, type = 'cuffdiff', title = TRUE,
#'  legend = TRUE, grid = TRUE
#' )
#' 
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsBoxPlot(
#'  data = df.deseq, d.factor = 'condition', type = 'deseq', 
#'  title = TRUE, legend = TRUE, grid = TRUE
#' )
#' 
#' # edgeR example
#' data("df.deseq")
#' require(edgeR)
#' vsBoxPlot(
#'  data = df.edger, d.factor = NULL, type = 'edger', title = TRUE,
#'  legend = TRUE, grid = TRUE
#' )

vsBoxPlot <- function(
    data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"), 
    title = TRUE, legend = TRUE, grid = TRUE
) {
    if (missing(type) || !type %in% c("cuffdiff", "deseq", "edger")) {
        stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
    }
    
    type <- match.arg(type)
    if(type == 'cuffdiff') {
        dat <- .getCuffBox(data)
    } else if (type == 'deseq') {
        dat <- .getDeseqBox(data, d.factor)
    } else if (type == 'edger') {
        dat <- .getEdgeBox(data)
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
