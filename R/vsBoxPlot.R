#' @title
#' Box plot for \eqn{log_{10}}(FPKM or CPM) distributions
#'
#' @author
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to extract necessary results-based data from
#' analytical objects to create a box plot comparing log10(FPKM or CPM)
#' distributions for experimental treatments.
#'
#' @param data output generated from calling the main routines of either
#'  \code{cuffdiff}, \code{DESeq2}, or \code{edgeR} analyses. For
#'  \code{cuffdiff}, this will be a \code{*_exp.diff} file. For \code{DESeq2},
#'  this will be a generated object of class \code{DESeqDataSet}. For
#'  \code{edgeR}, this will be a generated object of class \code{DGEList}.
#' @param d.factor a specified factor; for use with \code{DESeq2} objects only.
#'  This input equates to the first parameter for the contrast argument when
#'  invoking the \code{results()} function in \code{DESeq2}. Defaults to
#'  \code{NULL}
#' @param type an analysis classifier to tell the function how to process the
#'  data. Must be either \code{cuffdiff}, \code{deseq}, or \code{edger}.
#'  \code{cuffdiff} must be used with \code{cuffdiff} data; \code{deseq} must
#'  be used for \code{DESeq2} output; \code{edgeR} must be used with
#'  \code{edgeR} data. See the \code{data} parameter for further details.
#' @param title display the main title of plot. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no title will display in plot.
#' @param legend display legend of plot. Logical; defaults to \code{TRUE}.
#'  If set to \code{FALSE}, no legend will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to
#'  \code{TRUE}. If set to \code{FALSE}, no axis lines will display in plot.
#' @param aes changes overall layout of the graph. \code{box}: box plot;
#'  \code{violin}: violin plot; \code{boxdot}: box plot with dots;
#'  \code{viodot}: violin plot with dots; \code{viosumm}: violin plot with
#'  summary statistics; \code{notch}: box plots with notches.
#'  Defaults to \code{box}.
#' @param fill.color changes the fill color for the plots. See
#'  \code{RColorBrewer::display.brewer.all()} function for further details.
#'  If \code{NULL}, colors will default to standard \code{ggplot2}
#'  aesthetics.
#' @param data.return returns data output of plot. Logical; defaults to
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
#' vsBoxPlot(
#'  	data = df.cuff, d.factor = NULL, type = "cuffdiff", title = TRUE,
#'  	legend = TRUE, grid = TRUE
#' )
#'
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsBoxPlot(
#'  	data = df.deseq, d.factor = "condition", type = "deseq",
#'  	title = TRUE, legend = TRUE, grid = TRUE
#' )
#'
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsBoxPlot(
#'  	data = df.edger, d.factor = NULL, type = "edger", title = TRUE,
#'  	legend = TRUE, grid = TRUE
#' )
#'
#' # Display different colors for plot
#' data("df.edger")
#' vsBoxPlot(
#'  	data = df.edger, d.factor = NULL, type = "edger", title = TRUE,
#'  	legend = TRUE, grid = TRUE, fill.color = "RdGy",
#' 		data.return = FALSE
#' )
#'
#' # Extract data frame from visualization
#' data("df.edger")
#' require(edgeR)
#' tmp <- vsBoxPlot(
#'  	data = df.edger, d.factor = NULL, type = "edger", title = TRUE,
#'  	legend = TRUE, grid = TRUE, data.return = FALSE
#' )
#' df_box <- tmp[[1]] ## or use tmp$data
#' head(df_box)
#'
#' # Show plot from object (see prior example for more details)
#' tmp[[2]] ## or use tmp$plot

vsBoxPlot <- function(
	data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"),
	title = TRUE, legend = TRUE, grid = TRUE,
	aes = c("box", "violin", "boxdot", "viodot", "viosumm", "notch"),
	fill.color = NULL, data.return = FALSE
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
	if(type == "cuffdiff") {
		dat <- .getCuffBox(data)
	} else if (type == "deseq") {
		dat <- .getDeseqBox(data, d.factor)
	} else if (type == "edger") {
		dat <- .getEdgeBox(data)
	}

	if (!isTRUE(title)) {
		m.title <- NULL
	} else if (isTRUE(title) & type == "edger") {
		m.title <- ggtitle("CPM distribution")
	} else {
		m.title <- ggtitle("FPM distribution")
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

	if (type == "edger") {
		y.lab <- expression(paste("log"["10"], ' (CPM)'))
	} else {
		y.lab <- expression(paste("log"["10"], ' (FPM)'))
	}

	aes_types <- c("box", "violin", "boxdot", "viodot", "viosumm", "notch")
	if (missing(aes)) {
		aes <- "box"
	} else if (!aes %in% aes_types) {
		stop(
			paste0(
				"Please specfiy correct aesthetic:", "\n",
				"\t", "\"box\"....... (box plot)", "\n",
				"\t", "\"violin\".... (violin plot)", "\n",
				"\t", "\"boxdot\".... (box plot with dots)", "\n",
				"\t", "\"viodot\".... (violin plot with dots)", "\n",
				"\t", "\"viosumm\"... (violin plot with summary stats)", "\n",
				"\t", "\"notch\"..... (box plot with notch)", "\n"
			)
		)
	}

	brew_cols <- rownames(RColorBrewer::brewer.pal.info)
	if (is.null(fill.color)) {
		scale_col <- NULL
	} else if (!fill.color %in% brew_cols) {
		stop(
			paste(
				"Please specify correct RColorBrewer pallette.", "\n",
				"See \"RColorBrewer::display.brewer.all()\" for",
				"available options."
			)
		)
	} else {
		scale_col <- scale_fill_brewer(palette = fill.color)
	}

	key <- value <- NULL

	if (aes == "box") {
		aes_plot_1 <- geom_boxplot()
		aes_plot_2 <- NULL
	} else if (aes == "violin") {
		aes_plot_1 <- geom_violin(trim = FALSE)
		aes_plot_2 <- NULL
	} else if (aes == "boxdot") {
		aes_plot_1 <- geom_boxplot()
		aes_plot_2 <- geom_jitter(
			shape = 16,
			color = "grey14",
			alpha = 0.7,
			position = position_jitter(0.2)
		)
	} else if (aes == "viodot") {
		aes_plot_1 <- geom_violin(trim = FALSE)
		aes_plot_2 <- geom_jitter(
			shape = 16,
			color = "grey14",
			alpha = 0.7,
			position = position_jitter(0.2)
		)
	} else if (aes == "viosumm") {
		aes_plot_1 <- geom_violin(trim = FALSE)
		aes_plot_2 <- geom_boxplot(width = 0.1, fill = "white")
	} else if (aes == "notch") {
		aes_plot_1 <- geom_boxplot(notch = TRUE)
		aes_plot_2 <- NULL
	}


	tmp.plot <- ggplot(
			dat, aes(x = key, y = log10(value + 1), fill = key)
		) +
		aes_plot_1 +
		aes_plot_2 +
		xlab("Condition") +
		ylab(y.lab) +
		guides(fill = guide_legend(title = "Condition")) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		leg + m.title + grid + scale_col

	if (isTRUE(data.return)) {
		plot.l <- list(data = dat, plot = tmp.plot)
	} else {
		print(tmp.plot)
	}
}
