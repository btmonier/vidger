#' @title
#' MA plot from mean expression and \eqn{log} fold changes from different
#' analytical objects
#'
#' @author
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to extract necessary results-based data from
#' different output files to create a MA plot (i.e. a scatter plot) of
#' \eqn{log_{2}} fold changes versus normalized mean counts while implementing
#' ggplot2 aesthetics.
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
#' @param padj a user defined adjusted \code{p}-value cutoff point.
#'  Defaults to \code{0.05}.
#' @param y.lim set manual limits (boundaries) to the y axis. Defaults to
#'  \code{NULL}.
#' @param lfc \eqn{log} fold change level for setting conditonals. If no user
#'  input is added (\code{NULL}), value defaults to \code{1}.
#' @param title display the main title of plot. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no title will display in plot.
#' @param legend display legend of plot. Logical; defaults to \code{TRUE}.
#'  If set to \code{FALSE}, no legend will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no axis lines will display in plot.
#' @param highlight character string of IDs that will be highlighted. Set to
#'  \code{NULL} if you do not want highlighted data.
#' @param data.return returns data output of plot Logical; defaults to
#'  \code{FALSE}. If set to \code{TRUE}, a data frame will also be called.
#'  Assign to object for reproduction and saving of data frame. See final
#'  example for further details.
#'
#' @return An object created by \code{ggplot}
#'
#' @export
#'
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsMAPlot(
#'  	x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL,
#'  	type = "cuffdiff", padj = 0.05, y.lim = NULL, lfc = 1,
#'  	title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE
#' )
#'
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsMAPlot(
#'  	x = "treated_paired.end", y = "untreated_paired.end",
#'  	data = df.deseq, d.factor = "condition", type = "deseq",
#'  	padj = 0.05, y.lim = NULL, lfc = NULL, title = TRUE,
#'  	legend = TRUE, grid = TRUE, data.return = FALSE
#' )
#'
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsMAPlot(
#'  	x = "WM", y = "MM", data = df.edger, d.factor = NULL,
#'  	type = "edger", padj = 0.1, y.lim = NULL, lfc = 1,
#'  	title = FALSE, legend = TRUE, grid = TRUE, data.return = FALSE
#' )
#'
#' # Highlight IDs
#' data("df.deseq")
#' require(DESeq2)
#' hl <- c(
#'     "FBgn0022201",
#'     "FBgn0003042",
#'     "FBgn0031957",
#'     "FBgn0033853",
#'     "FBgn0003371"
#' )
#' vsMAPlot(
#'  	x = "treated_paired.end", y = "untreated_paired.end",
#'  	data = df.deseq, d.factor = "condition", type = "deseq",
#'  	padj = 0.05, y.lim = NULL, lfc = NULL, title = TRUE,
#'  	legend = TRUE, grid = TRUE, data.return = FALSE, highlight = hl
#' )
#'
#' # Extract data frame from visualization
#' data("df.cuff")
#' tmp <- vsMAPlot(
#'  	x = "hESC", y = "iPS", data = df.cuff,
#'  	d.factor = NULL, type = "cuffdiff", padj = 0.05,
#'  	y.lim = NULL, lfc = 1, title = TRUE, grid = TRUE,
#'  	data.return = TRUE
#' )
#' df_ma <- tmp[[1]] ## or use tmp$data
#' head(df_ma)
#'
#' # Show plot from object (see prior example for more details)
#' tmp[[2]] ## or use tmp$plot

vsMAPlot <- function(
	x, y, data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"),
	padj = 0.05, y.lim = NULL, lfc = NULL, title = TRUE, legend = TRUE,
	grid = TRUE, highlight = NULL, data.return = FALSE
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
		dat <- .getCuffMA(x, y, data)
	} else if (type == "deseq") {
		dat <- .getDeseqMA(x, y, data, d.factor)
	} else if (type == "edger") {
		dat <- .getEdgeMA(x, y, data)
	}

	if (!isTRUE(title)) {
		m.lab <- NULL
	} else {
		m.lab  <- ggtitle(paste(y, "vs.", x))
	}

	if (!isTRUE(legend)) {
		leg <- theme(legend.position = "none")
	} else {
		leg <- guides(
			colour = guide_legend(
				override.aes = list(size = 3)
			),
			shape = guide_legend(
				override.aes = list(size = 3)
			)
		)
	}

	if (!isTRUE(grid)) {
		grid <- theme_classic()
	} else {
		grid <- theme_bw()
	}

	dat$isDE <- ifelse(dat$padj <= padj, TRUE, FALSE)
	py <- dat$M

	if (is.null(y.lim)) {
		y.lim = c(-1.5, 1.5) * quantile(abs(py[is.finite(py)]), probs = 0.99)
	}
	if (is.null(lfc)) {
		lfc = 1
	}

	dat <- .ma.ranker(dat, padj, lfc, y.lim)

	tmp.size <- .ma.out.ranker(py, y.lim[2])
	tmp.col <- .ma.col.ranker(dat$isDE, py, lfc)
	tmp.shp <- .ma.shp.ranker(py, y.lim)

	tmp.cnt <- .ma.col.counter(dat, lfc)
	b <- tmp.cnt[[1]]
	g <- tmp.cnt[[2]]

	comp1 <- .ma.comp1(y.lim, padj, lfc, b, g)
	point <- geom_point(
		alpha = 0.7,
		aes(color = tmp.col, shape = tmp.shp, size = tmp.size)
	)
	comp2 <- .ma.comp2(
		comp1[[4]], comp1[[6]], comp1[[5]], comp1[[1]], comp1[[2]], comp1[[3]]
	)


	A <- id <- M <- NULL
	if (is.null(highlight)) {
		tmp.plot <- ggplot(
			dat, aes(x = A, y = pmax(y.lim[1], pmin(y.lim[2], py)))
		) +
			point +
			comp2$color + comp2$shape + comp1$hline1 + comp1$hline2 +
			comp1$hline3 + comp1$x.lab + comp1$y.lab + m.lab + ylim(y.lim) +
			comp2$size + grid + leg
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

		A <- NULL
		tmp.plot <- ggplot(
			dat, aes(x = A, y = pmax(y.lim[1], pmin(y.lim[2], py)))
		) +
			geom_point(
				alpha = 0.4,
				aes(color = tmp.col, shape = tmp.shp, size = tmp.size)
			) +
			comp2$color + comp2$shape + comp1$hline1 + comp1$hline2 +
			comp1$hline3 + comp1$x.lab + comp1$y.lab + m.lab + ylim(y.lim) +
			comp2$size +
			ggrepel::geom_label_repel(
				data = dat[which(dat$id %in% hl), ],
				aes(label = id, x = A, y = M),
				segment.size = 1,
				segment.color = "gray10",
				box.padding = unit(0.4, "lines"),
				point.padding = unit(0.4, "lines")
			) +
			geom_point(
				data = dat[which(dat$id %in% hl), ],
				aes(x = A, y = M),
				color = "red",
				size = 3
			) +
			grid + leg
	}

	if (isTRUE(data.return)) {
		dat2 <- dat[, -ncol(dat)]
		plot.l <- list(data = dat2, plot = tmp.plot)
	} else {
		print(tmp.plot)
	}
}