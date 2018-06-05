#' @title
#' Scatter plot of \eqn{log_{10}}(FPKM or CPM) values
#'
#' @author
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to visualize comparisons of \eqn{log_{10}} values
#' of either FPKM or CPM measurements of two treatments depending on
#' analytical type.
#'
#' @param x treatment \code{x} for comparison (\eqn{log_{2}}(x/control)). This
#'  will be a factor level in your data.
#' @param y treatment \code{y} for comparison (\eqn{log_{2}}(y/control)). This
#'  will be a factor level in your data.
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
#' @param title display the main title of plot. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no title will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no axis lines will display in plot.
#' @param highlight character string of IDs that will be highlighted. Set to
#'  \code{NULL} if you do not want highlighted data.
#' @param data.return returns data output of plot Logical; defaults to
#'  \code{FALSE}. If set to \code{TRUE}, a data frame will also be called.
#'  Assign to object for reproduction and saving of data frame. See final
#'  example for further details.
#' @param xaxis.title.size change the font size of the \code{x}-axis title 
#'  text. Defaults to \code{12}.
#' @param xaxis.text.size change the font size of the \code{x}-axis text. 
#'  Defaults to \code{10}.
#' @param yaxis.title.size change the font size of the \code{y}-axis title 
#'  text. Defaults to \code{12}.
#' @param yaxis.text.size change the font size of the \code{y}-axis text. 
#'  Defaults to \code{10}.
#' @param main.title.size change the font size of the plot title text. 
#'  Defaults to \code{15}.
#'
#' @return An object created by \code{ggplot}
#'
#' @export
#'
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsScatterPlot(
#'      x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL,
#'      type = "cuffdiff", title = TRUE, grid = TRUE
#' )
#'
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsScatterPlot(
#'      x = "treated_paired.end", y = "untreated_paired.end",
#'      data = df.deseq, d.factor = "condition", type = "deseq",
#'      title = TRUE, grid = TRUE
#' )
#'
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsScatterPlot(
#'      x = "WW", y = "WM", data = df.edger, d.factor = NULL,
#'      type = "edger", title = TRUE, grid = TRUE
#' )
#'
#' # Highlight IDs
#' data("df.cuff")
#' hl <- c(
#'      "XLOC_000033",
#'      "XLOC_000099",
#'      "XLOC_001414",
#'      "XLOC_001409"
#' )
#' vsScatterPlot(
#'      x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL,
#'      type = "cuffdiff", title = TRUE, grid = TRUE, highlight = hl
#' )
#'
#' # Extract data frame from visualization
#' data("df.cuff")
#' tmp <- vsScatterPlot(
#'      x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL,
#'      type = "cuffdiff", title = TRUE, grid = TRUE, data.return = TRUE
#' )
#' df_scatter <- tmp[[1]] ## or use tmp$data
#' head(df_scatter)
#'
#' # Show plot from object (see prior example for more details)
#' tmp[[2]] ## or use tmp$plot

vsScatterPlot <- function(
    x, y, data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"),
    title = TRUE, grid = TRUE, highlight = NULL, data.return = FALSE,
    xaxis.text.size = 10, yaxis.text.size = 10, xaxis.title.size = 12, 
    yaxis.title.size = 12, main.title.size = 15
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
        dat <- .getCuffScatter(x, y, data)
    } else if (type == "deseq") {
        dat <- .getDeseqScatter(x, y, data, d.factor)
    } else if (type == "edger") {
        dat <- .getEdgeScatter(x, y, data)
    }

    if (!isTRUE(title)) {
        m.title <- NULL
    } else {
        m.title <- ggtitle(paste(y, "vs.", x))
    }

    if (!isTRUE(grid)) {
        grid <- theme_classic()
    } else {
        grid <- theme_bw()
    }

    if (type == "edger") {
        aes.xlab <- bquote("log"["10"] ~ "(CPM) -" ~ .(x))
        aes.ylab <- bquote("log"["10"] ~ "(CPM) -" ~ .(y))
    } else {
        aes.xlab <- bquote("log"["10"] ~ "(FPM) -" ~ .(x))
        aes.ylab <- bquote("log"["10"] ~ "(FPM) -" ~ .(y))
    }

    text.size <- theme(
        axis.text.x = element_text(size = xaxis.text.size),
        axis.text.y = element_text(size = yaxis.text.size),
        axis.title.x = element_text(size = xaxis.title.size),
        axis.title.y = element_text(size = yaxis.title.size),
        plot.title = element_text(size = main.title.size)
    )

    id <- NULL
    if (is.null(highlight)) {
        tmp.plot <- ggplot(dat, aes(x = log10(x + 1), y = log10(y + 1))) +
            geom_point(size = 1) +
            xlab(aes.xlab) +
            ylab(aes.ylab) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
            grid + m.title + text.size


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
        tmp.plot <- ggplot(dat, aes(x = log10(x + 1), y = log10(y + 1))) +
            geom_point(size = 1, color = "grey73") +
            xlab(aes.xlab) +
            ylab(aes.ylab) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
            ggrepel::geom_label_repel(
                data = dat[which(dat$id %in% hl), ],
                aes(label = id),
                segment.size = 1,
                segment.color = "royalblue1",
                box.padding = unit(0.4, "lines"),
                point.padding = unit(0.4, "lines")
            ) +
            geom_point(
                data = dat[which(dat$id %in% hl), ],
                aes(x = log10(x), y = log10(y)),
                color = "royalblue1",
                size = 3
            ) +
            theme(text = element_text(size = 30)) +
            grid + m.title + text.size
    }
    if (isTRUE(data.return)) {
        plot.l <- list(data = dat, plot = tmp.plot)
        return(plot.l)
    } else {
        print(tmp.plot)
    }
}
