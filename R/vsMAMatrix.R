#' @title
#' MA plot matrix from \eqn{log_{2}} fold changes and
#' \eqn{-log_{10}}(\eqn{p}-values)
#'
#' @author
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to generate MA plots for all possible treatment
#' combinations for a given factor in either a cuffdiff, DESeq2, or edgeR
#' data set.
#'
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
#' @param y.lim set manual limits (boundaries) to the y axis. Defaults to
#'  \code{NULL}.
#' @param lfc \eqn{log} fold change level for setting conditonals. If no user
#'  input is added (\code{NULL}), value defaults to \code{1}.
#' @param title display the main title of plot. Logical; defaults to
#' 	\code{TRUE}. If set to \code{FALSE}, no title will display in plot.
#' @param legend display legend of plot. Logical; defaults to \code{TRUE}.
#'  If set to \code{FALSE}, no legend will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no axis lines will display in plot.
#' @param counts displays the number of differentially expressed genes for
#'  each treatment comparison. Defaults to \code{TRUE}.
#' @param data.return returns data output of plot. Logical; defaults to
#'  \code{FALSE}. If set to \code{TRUE}, a data frame will also be called.
#'  Assign to object for reproduction and saving of data frame. See final
#'  example for further details.
#' @param xaxis.title.size change the font size of the \code{x}-axis title 
#'  text. Defaults to \code{10}.
#' @param xaxis.text.size change the font size of the \code{x}-axis text. 
#'  Defaults to \code{9}.
#' @param yaxis.title.size change the font size of the \code{y}-axis title 
#'  text. Defaults to \code{10}.
#' @param yaxis.text.size change the font size of the \code{y}-axis text. 
#'  Defaults to \code{9}.
#' @param main.title.size change the font size of the plot title text. 
#'  Defaults to \code{15}.
#' @param legend.text.size change the font size of the legend body text.
#'  Defaults to \code{9}.
#' @param facet.title.size change the font size of the facet wrap title text.
#'  Defaults to \code{10}.
#'
#' @return An object created by \code{ggplot}
#'
#' @export
#'
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsMAMatrix(
#'  	data = df.cuff, d.factor = NULL, type = "cuffdiff",
#'  	padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE,
#'  	grid = TRUE, counts = TRUE, data.return = FALSE
#' )
#'
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsMAMatrix(
#'  	data = df.deseq, d.factor = "condition", type = "deseq",
#'  	padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE,
#'  	grid = TRUE, counts = TRUE, data.return = FALSE
#' )
#'
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsMAMatrix(
#'  	data = df.edger, d.factor = NULL, type = "edger",
#'  	padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE,
#'  	grid = TRUE, counts = TRUE, data.return = FALSE
#' )
#'
#' # Extract data frame from visualization
#' data("df.cuff")
#' tmp <- vsMAMatrix(
#'  	data = df.cuff, d.factor = NULL,
#'  	type = "cuffdiff", padj = 0.05, y.lim = NULL,
#'  	lfc = 1, title = TRUE, grid = TRUE,
#'  	counts = TRUE, data.return = TRUE
#' )
#' df_mamat <- tmp[[1]] ## or use tmp$data
#' head(df_mamat)
#'
#' # Show plot from object (see prior example for more details)
#' tmp[[2]] ## or use tmp$plot

vsMAMatrix <- function(
	data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"),
	padj = 0.05, y.lim = NULL, lfc = NULL, title = TRUE, legend = TRUE,
	grid = TRUE, counts = TRUE, data.return = FALSE, xaxis.text.size = 9, 
    yaxis.text.size = 9, xaxis.title.size = 10, yaxis.title.size = 10, 
    main.title.size = 15, legend.text.size = 9, facet.title.size = 10
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
		dat <- .getCuffMAMatrix(data)
	} else if (type == "deseq") {
		dat <- .getDeseqMAMatrix(data, d.factor)
	} else if (type == "edger") {
		dat <- .getEdgeMAMatrix(data)
	}

	if (!isTRUE(title)) {
		m.lab <- NULL
	} else {
		m.lab <- ggtitle("MA Matrix")
	}

	if (!isTRUE(legend)) {
		leg <- theme(legend.position = "none")
	} else {
		leg <- guides(colour = guide_legend(override.aes = list(size = 3)))
	}

	if (!isTRUE(grid)) {
		grid <- theme_classic()
	} else {
		grid <- theme_bw()
	}

	dat$isDE <- ifelse(dat$padj <= padj, TRUE, FALSE)
	py <- dat$M

	if (is.null(y.lim)) {
		y.lim = c(-1, 1) * quantile(abs(py[is.finite(py)]), probs = 0.99) * 0.8
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

	tmp.l <- .mamat.col.count(dat)

   if (isTRUE(counts)) {
      b.count <- ggplot2::geom_text(
         ggplot2::aes(label = .data$Freq, x = -Inf, y = Inf),
         data = tmp.l$blue,
         x = 1.5,
         y = -1,
         color = "royalblue1",
         fontface = 2
      )
      g.count <- ggplot2::geom_text(
         ggplot2::aes(label = .data$Freq, x = Inf, y = Inf),
         data = tmp.l$green,
         x = 1.5,
         y = 1.5,
         color = "green",
         fontface = 2
      )
   } else {
      b.count <- NULL
      g.count <- NULL
   }

	text.size <- theme(
        axis.text.x = element_text(size = xaxis.text.size),
        axis.text.y = element_text(size = yaxis.text.size),
        axis.title.x = element_text(size = xaxis.title.size),
        axis.title.y = element_text(size = yaxis.title.size),
        plot.title = element_text(size = main.title.size),
        legend.text = element_text(size = legend.text.size),
        strip.text = element_text(size = facet.title.size)
    )

	M <- A <- pval <- color <- size <- shape <- NULL
	tmp.plot <- ggplot(
		dat, aes(x = A, y = pmax(y.lim[1], pmin(y.lim[2], py)))
	) +
		point +
		comp2$shape +
		comp2$color +
		comp2$size +
		comp1$hline1 +
		comp1$hline2 +
		comp1$hline3 +
		comp1$x.lab +
		comp1$y.lab +
		m.lab +
		grid +
		leg +
		b.count +
		g.count +
		ylim(y.lim) +
		facet_wrap(id_x ~ id_y) +
		text.size

	if (isTRUE(data.return)) {
		dat2 <- dat[, -ncol(dat)]
		plot.l <- list(data = dat2, plot = tmp.plot)
	} else {
		print(tmp.plot)
	}
}
