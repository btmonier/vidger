#' @title 
#' Scatter plot matrix of \eqn{log_{10}}(FPKM or CPM) values
#'
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function will generate a matrix of scatterplots for all possible
#' treatment combinations with additional distribution info.
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
#' @param comp treatments you would like to compare in the form of a vector. 
#'  If no parameter is specified, all possible treatment comparisons will be 
#'  made. Defaults to \code{NULL}.
#' @param title display the main title of plot. Logical; defaults to 
#'  \code{TRUE}. If set to \code{FALSE}, no title will display in plot.
#' @param grid display major and minor axis lines. Logical; defaults to 
#'  \code{TRUE}. If set to \code{FALSE}, no axis lines will display in plot.
#' @param man.title a manually specified title at the authors discretion. 
#'  Defaults to \code{NULL}.
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
#' vsScatterMatrix(
#'  	data = df.cuff, d.factor = NULL, type = "cuffdiff", 
#'  	comp = NULL, title = TRUE, grid = TRUE, 
#'  	man.title = "Example title"
#' )
#' 
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsScatterMatrix(
#'  	data = df.deseq, d.factor = "condition", type = "deseq",
#'  	comp = NULL, title = TRUE, grid = FALSE, man.title = NULL
#' )
#' 
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsScatterMatrix(
#'  	data = df.edger, d.factor = NULL, type = "edger", 
#'  	comp = c("WM", "MM"), title = TRUE, grid = TRUE, 
#'  	man.title = NULL
#' )
#' 
#' # Extract data frame from visualization
#' data("df.edger")
#' tmp <- vsScatterMatrix(
#'  	data = df.edger, d.factor = NULL, type = "edger", 
#'  	comp = c("WM", "MM"), title = TRUE, grid = TRUE, 
#'  	man.title = NULL, data.return = TRUE
#' )
#' df_scatmat <- tmp[[1]] ## or use tmp$data
#' head(df_scatmat)
#' 
#' # Show plot from object (see prior example for more details)
#' tmp[[2]] ## or use tmp$plot

vsScatterMatrix <- function(
	data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"), 
	comp = NULL, title = TRUE, grid = TRUE, man.title = NULL,
	data.return = FALSE
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
		dat <- .getCuffScatterMatrix(data)
	} else if (type == "deseq") {
		dat <- .getDeseqScatterMatrix(data, d.factor)
	} else if (type == "edger") {
		dat <- .getEdgeScatterMatrix(data)
	}
	
	if(!is.null(comp)) {
		comp <- comp
	} else {
		comp <- seq_len(ncol(dat))
	}
	
	if (!isTRUE(title)) {
		m.title <- NULL
	} else if (isTRUE(title) & !is.null(man.title)) {
		m.title <- ggtitle(man.title)
	} else if (isTRUE(title) & type == "edger") {
		m.title <- "CPM Comparisons"
	} else {
		m.title <- "FPKM Comparisons"
	}

	if (!isTRUE(grid)) {
		grid <- theme_classic()
	} else {
		grid <- theme_bw()
	}
	
	if (type == "edger") {
		aes.xlab <- bquote("log"["10"] ~ "(CPM)")
		aes.ylab <- bquote("log"["10"] ~ "(CPM)")
	} else {
		aes.xlab <- bquote("log"["10"] ~ "(FPM)")
		aes.ylab <- bquote("log"["10"] ~ "(FPM)")
	}
	
	tmp.plot <- GGally::ggpairs(
		log10(dat + 1), 
		title = m.title, 
		xlab = aes.xlab, 
		ylab = aes.ylab, 
		columns = comp,
		lower = list(continuous = wrap("points", size = 0.5))
	) + grid

	if (isTRUE(data.return)) {
		plot.l <- list(data = dat, plot = tmp.plot)
	} else {
		print(tmp.plot)
	}	
}
