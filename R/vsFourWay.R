#' @title 
#' Four-Way plot for comparison of log fold changes in a multiple factor RNA
#'  seq
#' experiment from different analytical objects
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to extract necessary results-based data from a 
#' DESeq object class to create a .four-way plot to compare log fold changes 
#' in various treatments using ggplot2 aesthetics. 
#' 
#' @details  
#' This function allows the user to extract various elements from a different 
#' analytical object class which in turn, creates a temporary data frame to 
#' plot the necessary ggplot aesthetics. In order for this function to work, 
#' RNA seq experiments must have multiple factors (i.e. two treatments and a 
#' control) and levels including treatments and controls. By having the 
#' recommended criteria, this function will extract the necessary data 
#' dependent on the analysis performed. Data points with 'extreme' values that 
#' exceed the default viewing frame of the plot will change character classes 
#' (i.e. points of interest a substantially large log fold change). 
#'  
#' @param x treatment `x` for comparison (log2(x/control)).
#' @param y treatment `y` for comparison (log2(y/control)).
#' @param control `control` treatment for comparisons of the x and y axes.
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param d.factor a specified factor; for use with DESeq2 objects only. 
#'  Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the 
#'  data. Must be either `cuffdiff`, `deseq`, or `edgeR`.
#' @param padj a user defined adjusted p-value cutoff point. 
#'  Defaults to `0.05`.
#' @param x.lim set manual limits to the x axis. Defaults to `NULL`.
#' @param y.lim set manual limits to the y axis. Defaults to `NULL`.
#' @param lfc log fold change level for setting conditonals. If no user input 
#'  is added (`NULL`), value defaults to `1`.
#' @param legend shows legend of plot. Defaults to `TRUE`.
#' @param title show title of plot. Defaults to `TRUE`.
#' @param grid show major and minor axis lines. Defaults to `TRUE`.
#' @param data.return returns data output of plot if set to `TRUE`. 
#'  Defaults to `FASLSE`.
#' 
#' @return An object created by \code{ggplot}
#' 
#' @export
#' 
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsFourWay(
#'  x = 'hESC', y = 'iPS', control = 'Fibroblasts', data = df.cuff, 
#'  d.factor = NULL, type = 'cuffdiff', padj = 0.05, x.lim = NULL, 
#'  y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
#'  data.return = FALSE
#' )
#' 
#' # DESeq2 example
#' data("df.deseq")
#' vsFourWay(
#'  x = 'treated_paired.end', y = 'untreated_paired.end', 
#'  control = 'untreated_single.read', data = df.deseq, 
#'  d.factor = 'condition', type = 'deseq', padj = 0.05, 
#'  x.lim = NULL, y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
#'  data.return = FALSE
#' )
#' 
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsFourWay(
#'  x = 'WM', y = 'WW', control = 'MM', data = df.edger, 
#'  d.factor = NULL, type = 'edger', padj = 0.05, x.lim = NULL, 
#'  y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
#'  data.return = FALSE
#' )
#'                 
#' # Extract data frame from visualization
#' data("df.cuff")
#' tmp <- vsFourWay(
#'  x = 'WM', y = 'WW', control = 'MM', data = df.edger, 
#'  d.factor = NULL, type = 'edger', padj = 0.05, 
#'  x.lim = NULL, y.lim = NULL, lfc = 2, title = TRUE, 
#'  grid = TRUE, data.return = FALSE
#' )
#' df.four <- tmp[[1]]
#' head(df.four)


vsFourWay <- function(
    x, y, control, data, d.factor = NULL, type, padj = 0.05, 
    x.lim = NULL, y.lim = NULL, lfc = NULL, legend = TRUE, 
    title = TRUE, grid = TRUE, data.return = FALSE
){
    if (missing(type)) {
        stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
    }
    if(type == 'cuffdiff'){
        dat <- .getCuffFourWay(x, y, control, data)
    } else if (type == 'deseq') {
        dat <- .getDeseqFourWay(x, y, control, data, d.factor)
    } else if (type == 'edger') {
        dat <- .getEdgeFourWay(x, y, control, data)
    } else {
        stop('Please enter correct analysis type.')
    }
    if (!isTRUE(title)) {
        m.lab <- NULL
    } else {
        m.lab  <- ggtitle(paste(
            y, 'and', x, 'expression relative to', control))
    }
    if (!isTRUE(legend)) {
        leg <- theme(legend.position = 'none')
    } else {
        leg <- guides(
            colour = guide_legend(override.aes = list(size = 3)),
            shape = guide_legend(override.aes = list(size = 3))
        )
    }
    if (!isTRUE(grid)) {
        grid <- theme_classic()
    } else {
        grid <- theme_bw()
    }
    
    
    dat$isDE_x   <- ifelse(dat$padj_x <= padj, TRUE, FALSE)
    dat$isDE_y   <- ifelse(dat$padj_y <= padj, TRUE, FALSE)
    dat$isDE_all <- ifelse(
        (dat$isDE_x == TRUE) & (dat$isDE_y == TRUE), TRUE, FALSE
    )
    px   <- dat$logFC_x
    py   <- dat$logFC_y
    pall <- dat$isDE_all == TRUE
    all.lab  <- c(x, y, control)
    
    
    if (is.null(x.lim))
        x.lim = c(-1, 1) * quantile(abs(px[is.finite(px)]), probs = 0.99) * 2
    if (is.null(y.lim))
        y.lim = c(-1, 1) * quantile(abs(py[is.finite(py)]), probs = 0.99) * 2
    if (is.null(lfc))
        lfc = 1 
    
    dat <- .four.ranker(dat, padj, lfc, x.lim, y.lim)
    
    
    aes.x <- pmax(x.lim[1], pmin(x.lim[2], px))
    aes.y <- pmax(y.lim[1], pmin(y.lim[2], py))
    
    
    tmp.size <- .four.out.ranker(px, py, x.lim[2], y.lim[2])
    tmp.col <- .four.col.ranker(dat, lfc)
    tmp.shp <- .four.shp.ranker(px, py, x.lim, y.lim)
    
    
    tmp.cnt <- .four.col.counter(dat, lfc)
    b <- tmp.cnt[[1]]
    g <- tmp.cnt[[2]]
    r <- tmp.cnt[[3]]
    
    
    comp1 <- .four.comp1(x.lim, y.lim, padj, all.lab, lfc, b, g, r)
    point <- geom_point(
        alpha = 0.7, 
        aes(color = tmp.col, shape = tmp.shp, size = tmp.size)
    )
    comp2 <- .four.comp2(
        comp1[[3]], comp1[[4]], comp1[[6]], comp1[[5]], comp1[[1]], comp1[[2]]
    )
    
    
    tmp.plot <- ggplot(dat, aes(x = aes.x, y = aes.y)) +
        point + comp1$vline1 + comp1$vline2 + comp1$vline3 + comp1$hline1 + 
        comp1$hline2 + comp1$hline3 + comp2$color + comp2$shape + comp2$size + 
        comp1$x.lab + comp1$y.lab + grid  + m.lab + xlim(x.lim) + ylim(y.lim) +
        leg
    
    if (isTRUE(data.return)) {
        dat2 <- dat[, -ncol(dat)]
        plot.l <- list(data = dat2, plot = tmp.plot)
    } else {
        print(tmp.plot)
    }
}
