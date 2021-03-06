#' @title
#' Four-Way plot for comparison of log fold changes in a multiple factor
#' RNA-seq experiment from different analytical objects
#'
#' @author
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to extract necessary results-based data from a
#' DESeq object class to create a four-way plot to compare log fold changes
#' in various treatments using ggplot2 aesthetics.
#'
#' @details
#' This function allows the user to extract various elements from a different
#' analytical object class which in turn, creates a temporary data frame to
#' plot the necessary ggplot aesthetics. In order for this function to work,
#' RNA seq experiments must have multiple factors (i.e. two treatments and a
#' control) and levels including treatments and controls. By having the
#' recommended criteria, this function will extract the necessary data
#' dependent on the analysis performed. Data points with "extreme" values that
#' exceed the default viewing frame of the plot will change character classes
#' (i.e. points of interest a substantially large log fold change).
#'
#' @param x treatment \code{x} for comparison (\eqn{log_{2}}(x/control)). This
#'  will be a factor level in your data.
#' @param y treatment \code{y} for comparison (\eqn{log_{2}}(y/control)). This
#'  will be a factor level in your data.
#' @param control \code{control} treatment for comparisons of the x and y axes.
#'  This will be a factor level in your data.
#' @param data output generated from calling the main routines of either
#'  \code{cuffdiff}, \code{DESeq2}, or \code{edgeR} analyses. For
#'  \code{cuffdiff}, this will be a \code{*_exp.diff} file. For \code{DESeq2},
#'  this will be a generated object of class \code{DESeqDataSet}. For
#'  \code{edgeR}, this will be a generated object of class \code{DGEList}.
#' @param d.factor a specified factor; for use with \code{DESeq2} objects only.
#'  This input equates to the first parameter for the contrast argument when
#'  invoking the \code{results()} function in \code{DESeq2}. Defaults to
#'  \code{NULL}.
#' @param type an analysis classifier to tell the function how to process the
#'  data. Must be either \code{cuffdiff}, \code{deseq}, or \code{edger}.
#'  \code{cuffdiff} must be used with \code{cuffdiff} data; \code{deseq} must
#'  be used for \code{DESeq2} output; \code{edgeR} must be used with
#'  \code{edgeR} data. See the \code{data} parameter for further details.
#' @param padj a user defined adjusted \eqn{p}-value cutoff point.
#'  Defaults to \code{0.05}.
#' @param x.lim set manual limits (boundaries) to the x axis. Defaults to
#'  \code{NULL}.
#' @param y.lim set manual limits (boundaries) to the y axis. Defaults to
#'  \code{NULL}.
#' @param lfc log fold change level for setting conditonals. If no user input
#'  is added (\code{NULL}), value defaults to \code{1}.
#' @param title display the main title of plot. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no title will display in plot.
#' @param legend display legend of plot. Logical; defaults to \code{TRUE}.
#'  If set to \code{FALSE}, no legend will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no axis lines will display in plot.
#' @param highlight character string of IDs that will be highlighted. Set to
#'  \code{NULL} if you do not want highlighted data.
#' @param data.return returns data output of plot. Logical; defaults to
#'  \code{FALSE}. If set to \code{TRUE}, a data frame will also be called.
#'  Assign to object for reproduction and saving of data frame. See final
#'  example for further details.
#' @param xaxis.title.size change the font size of the \code{x}-axis title 
#'  text. Defaults to \code{10}.
#' @param xaxis.text.size change the font size of the \code{x}-axis text. 
#'  Defaults to \code{10}.
#' @param yaxis.title.size change the font size of the \code{y}-axis title 
#'  text. Defaults to \code{10}.
#' @param yaxis.text.size change the font size of the \code{y}-axis text. 
#'  Defaults to \code{10}.
#' @param main.title.size change the font size of the plot title text. 
#'  Defaults to \code{15}.
#' @param legend.text.size change the font size of the legend body text.
#'  Defaults to \code{9}.
#'
#' @return An object created by \code{ggplot}
#'
#' @export
#'
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsFourWay(
#'      x = "hESC", y = "iPS", control = "Fibroblasts", data = df.cuff,
#'      d.factor = NULL, type = "cuffdiff", padj = 0.05, x.lim = NULL,
#'      y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE,
#'      data.return = FALSE
#' )
#'
#' # DESeq2 example
#' data("df.deseq")
#' vsFourWay(
#'      x = "treated_paired.end", y = "untreated_paired.end",
#'      control = "untreated_single.read", data = df.deseq,
#'      d.factor = "condition", type = "deseq", padj = 0.05,
#'      x.lim = NULL, y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE,
#'      data.return = FALSE
#' )
#'
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsFourWay(
#'      x = "WM", y = "WW", control = "MM", data = df.edger,
#'      d.factor = NULL, type = "edger", padj = 0.05, x.lim = NULL,
#'      y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE,
#'      data.return = FALSE
#' )
#'
#' # Highlight IDs
#' data("df.edger")
#' require(edgeR)
#' hl <- c(
#'     "ID_639",
#'     "ID_518",
#'     "ID_602",
#'     "ID_449",
#'     "ID_076"
#' )
#' vsFourWay(
#'     x = "WM", y = "WW", control = "MM", data = df.edger,
#'     d.factor = NULL, type = "edger", padj = 0.05, x.lim = NULL,
#'     y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE,
#'     data.return = FALSE, highlight = hl
#' )
#'
#' # Extract data frame from visualization
#' data("df.cuff")
#' tmp <- vsFourWay(
#'      x = "WM", y = "WW", control = "MM", data = df.edger,
#'      d.factor = NULL, type = "edger", padj = 0.05,
#'      x.lim = NULL, y.lim = NULL, lfc = 2, title = TRUE,
#'      grid = TRUE, data.return = TRUE
#' )
#' df_four <- tmp[[1]] ## or use tmp$data
#' head(df_four)
#'
#' # Show plot from object (see prior example for more details)
#' tmp[[2]] ## or use tmp$plot


vsFourWay <- function(
    x, y, control, data, d.factor = NULL,
    type = c("cuffdiff", "deseq", "edger"), padj = 0.05, x.lim = NULL,
    y.lim = NULL, lfc = NULL, legend = TRUE, title = TRUE, grid = TRUE,
    highlight = NULL, data.return = FALSE, xaxis.text.size = 10, 
    yaxis.text.size = 10, xaxis.title.size = 10, yaxis.title.size = 10, 
    main.title.size = 15, legend.text.size = 9
) {
    if (missing(type) || !type %in% c("cuffdiff", "deseq", "edger")) {
        stop(
            paste(
                "Please specify analysis type",
                "(\"cuffdiff\", \"deseq\", or \"edger\")"
            )
        )
    }

    type <- match.arg(type)
    if (type == "cuffdiff") {
        dat <- .getCuffFourWay(x, y, control, data)
    } else if (type == "deseq") {
        dat <- .getDeseqFourWay(x, y, control, data, d.factor)
    } else if (type == "edger") {
        dat <- .getEdgeFourWay(x, y, control, data)
    }

    if (!isTRUE(title)) {
        m.lab <- NULL
    } else {
        m.lab  <- ggtitle(paste(
            y, "and", x, "expression relative to", control))
    }

    if (!isTRUE(legend)) {
        leg <- theme(legend.position = "none")
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

    text.size <- theme(
        axis.text.x = element_text(size = xaxis.text.size),
        axis.text.y = element_text(size = yaxis.text.size),
        axis.title.x = element_text(size = xaxis.title.size),
        axis.title.y = element_text(size = yaxis.title.size),
        plot.title = element_text(size = main.title.size),
        legend.text = element_text(size = legend.text.size)
    )   

    id <- logFC_x <- logFC_y <- NULL
    if (is.null(highlight)) {
        tmp.plot <- ggplot(dat, aes(x = aes.x, y = aes.y)) +
            point + comp1$vline1 + comp1$vline2 + comp1$vline3 +
            comp1$hline1 + comp1$hline2 + comp1$hline3 + comp2$color +
            comp2$shape + comp2$size + comp1$x.lab + comp1$y.lab +
            grid  + m.lab + xlim(x.lim) + ylim(y.lim) + leg + text.size
    } else {
        tl <- length(setdiff(highlight, dat$id))
        if (!is.atomic(highlight)) {
            stop("\"highlight\" must be vector.")
        } else if (all(highlight %in% dat$id)) {
            hl <- highlight
        } else if (tl > 0 && tl < length(highlight)) {
            remove <- setdiff(highlight, dat$id)
            message("Some IDs not found in data frame:")
            print(remove)
            message("Plotting the remaining samples...")
            hl <- highlight[!highlight %in% remove]
        } else if (!all(highlight %in% dat$id)) {
            stop("No IDs in highlight vector are present in data frame.")
        }

        tmp.plot <- ggplot(dat, aes(x = aes.x, y = aes.y)) +
            geom_point(
                alpha = 0.4,
                aes(color = tmp.col, shape = tmp.shp, size = tmp.size)
            ) +
            comp1$vline1 + comp1$vline2 + comp1$vline3 +
            comp1$hline1 + comp1$hline2 + comp1$hline3 + comp2$color +
            comp2$shape + comp2$size + comp1$x.lab + comp1$y.lab +
            ggrepel::geom_label_repel(
                data = dat[which(dat$id %in% hl), ],
                aes(label = id, x = logFC_x, y = logFC_y),
                segment.size = 1,
                segment.color = "gray10",
                box.padding = unit(0.4, "lines"),
                point.padding = unit(0.4, "lines")
            ) +
            geom_point(
                data = dat[which(dat$id %in% hl), ],
                aes(x = logFC_x, y = logFC_y),
                color = "grey10",
                size = 3
            ) +
            grid  + m.lab + xlim(x.lim) + ylim(y.lim) + leg + text.size
    }

    if (isTRUE(data.return)) {
        dat2 <- dat[, -ncol(dat)]
        plot.l <- list(data = dat2, plot = tmp.plot)
    } else {
        print(tmp.plot)
    }
}
