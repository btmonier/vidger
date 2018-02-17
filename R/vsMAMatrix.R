#' @title 
#' MA plot matrix from log2 fold changes and -log10(p-values)
#'   
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#' 
#' @description
#' This function allows you to generate MA plots for all possible treatment
#' combinations for a given factor in either a cuffdiff, DESeq2, or edgeR
#' data set.
#'  
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param d.factor a specified factor; for use with DESeq2 objects only. 
#'  Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the 
#'  data. Must be either `cuffdiff`, `deseq`, or `edgeR`.
#' @param padj a user defined adjusted p-value cutoff point. Defaults to `0.1`.
#' @param y.lim set manual limits to the y axis. Defaults to `NULL`.
#' @param lfc log fold change level for setting conditonals. If no user input 
#'  is added (`NULL`), value defaults to `1`.
#' @param title show title of plot. Defaults to `TRUE`.
#' @param legend shows legend of plot. Defaults to `TRUE`.
#' @param grid show major and minor axis lines. Defaults to `TRUE`.
#' @param counts displays the number of differentially expressed genes for 
#'  each treatment comparison. Defaults to `TRUE`.
#' @param data.return returns data output of plot if set to `TRUE`. Defaults 
#'  to `FASLSE`.
#' 
#' @return An object created by \code{ggplot}
#' 
#' @export
#' 
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsMAMatrix(data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
#'                 padj = 0.05, y.lim = NULL, lfc = 2, title = TRUE, 
#'                 grid = TRUE, counts = TRUE, data.return = FALSE)
#' 
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsMAMatrix(data = df.deseq, d.factor = 'condition', type = 'deseq', 
#'                 padj = 0.05, y.lim = NULL, lfc = 2, title = TRUE, 
#'                 grid = TRUE, counts = TRUE, data.return = FALSE)
#' 
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsMAMatrix(data = df.edger, d.factor = NULL, type = 'edger', 
#'                 padj = 0.05, y.lim = NULL, lfc = 2, title = TRUE, 
#'                 grid = TRUE, counts = TRUE, data.return = FALSE)
#'                 
#' # Extract data frame from visualization
#' data("df.cuff")
#' tmp <- vsMAMatrix(data = df.cuff, d.factor = NULL, 
#'                        type = 'cuffdiff', padj = 0.05, y.lim = NULL,
#'                        lfc = 2, title = TRUE, grid = TRUE, 
#'                        counts = TRUE, data.return = TRUE)
#' df.vmat <- tmp[[1]]
#' head(df.vmat)

vsMAMatrix <- function(
    data, d.factor = NULL, type, padj = 0.1, y.lim = NULL, lfc = NULL, 
    title = TRUE, legend = TRUE, grid = TRUE, counts = TRUE,
    data.return = FALSE
) {
    if (missing(type)) {
        stop(
            'Please specify analysis type ("cuffdiff", "deseq", or "edger")'
        )
    }
    if(type == 'cuffdiff') {
        dat <- .getCuffMAMatrix(data)
    } else if (type == 'deseq') {
        dat <- .getDeseqMAMatrix(data, d.factor)
    } else if (type == 'edger') {
        dat <- .getEdgeMAMatrix(data)
    } else {
        stop('Please enter correct analysis type.')
    }
    if (!isTRUE(title)) {
        m.lab <- NULL
    } else {
        m.lab <- ggtitle('MA Matrix')
    }
    if (!isTRUE(legend)) {
        leg <- theme(legend.position = 'none')
    } else {
        leg <- guides(colour = guide_legend(override.aes = list(size = 3)))
    }
    if (!isTRUE(grid)) {
        grid <- theme_classic()
    } else {
        grid <- theme_bw()
    }
    
    py <- dat$M
    p <- padj
    
    if (is.null(y.lim)) {
        y.lim = c(-1, 1) * quantile(abs(py[is.finite(py)]), probs = 0.99) * 0.8
    }
    if (is.null(lfc)) {
        lfc = 1
    }
    
    pc <- .mamat.comp(padj, lfc)
    dat <- .mamat.ranker(dat, padj, lfc, y.lim)
        
    tmp.l <- .mamat.col.count(dat)
    
    if (isTRUE(counts)) {
        b.count <- annotate(
            'text', 
            x = -Inf, 
            y = Inf, 
            vjust = 1.5, 
            hjust = -1, 
            label = tmp.l$blue$Freq, 
            color = 'royalblue1', 
            fontface = 2
        )
        g.count <- annotate(
            'text', 
            x = Inf, 
            y = Inf, 
            vjust = 1.5, 
            hjust = 1.5, 
            label = tmp.l$green$Freq, 
            color = 'green', 
            fontface = 2
        )
    } else {
        b.count <- NULL
        g.count <- NULL
    }
   


    M <- A <- pval <- color <- size <- shape <- NULL
    tmp.plot <- ggplot(
    	dat,
    	if (type != "cuffdiff") {
    		aes(
    			x = A, y = pmax(y.lim[1], pmin(y.lim[2], py))
    		)
    	} else if (type == "cuffdiff") {
	    	aes(
	    		x = A, y = M
	    	)
    	}

    ) +
        geom_point(
        	aes(
            	color = color, shape = shape, size = size
            ), 
            na.rm = TRUE
        ) +
        scale_color_manual(
            name = '',
            values = c(
                'grey' = 'grey73', 
                'blue' = 'royalblue1', 
                'green' = 'green'
            ),
            labels = c(
                'grey' = pc$gry, 
                'blue' = pc$blu,
                'green' = pc$grn
            )
        ) +
        scale_shape_manual(
            name = '',
            values = c('circle' = 16, 'l.triangle' = 2, 'r.triangle' = 6),
            labels = c(
            	"circle" = paste(
            		round(y.lim[1], 2), "< lfc <", round(y.lim[2]), 2
            	),
            	"l.triangle" = paste("lfc >", round(y.lim[2]), 2),
            	"r.triangle" = paste("lfc <", round(y.lim[1]), 2)
            )
            # guide = 'none'
        ) +
        scale_size_manual(
            name = '',
            values = c(
                'sub' = 1, 
                't1' = 1.5,
                't2' = 2,
                't3' = 3,
                't4' = 5
            ),
            labels = c(
            	"sub" = "SUB",
            	"t1" = "T-1",
            	"t2" = "T-2",
            	"t3" = "T-3",
            	"t4" = "T-4"
            )
            # guide = 'none'
        ) +
        theme_bw() +
        facet_wrap(id_x ~ id_y) +
        pc$vline1 + pc$vline2 + pc$vline3 + pc$y.lab + pc$x.lab +
        b.count + g.count + ylim(y.lim) + m.lab + leg
    
    if (isTRUE(data.return)) {
        plot.l <- list(data = dat, plot = tmp.plot)
    } else {
        print(tmp.plot)
    }
}