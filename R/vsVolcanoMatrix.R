#' @title 
#' Volcano plot matrix from log2 fold changes and -log10(p-values)
#'   
#' @author 
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#' 
#' @description
#' This function allows you to extract necessary results-based data from a 
#' DESEq object class to create a volcano plot (i.e. a scatter plot) of the
#' negative log of the p-value versus the log of the fold change while
#' implementing ggplot2 aesthetics for all possible combinations of treatments. 
#'  
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param d.factor a specified factor; for use with DESeq2 objects only. Defaults to `NULL`
#' @param type an analysis classifier to tell the function how to process the data. Must be either `cuffdiff`, `deseq`, or `edgeR`.
#' @param padj a user defined adjusted p-value cutoff point. Defaults to `0.1`.
#' @param x.lim set manual limits to the x axis. Defaults to `NULL`.
#' @param lfc log fold change level for setting conditonals. If no user input is added (`NULL`), value defaults to `1`.
#' @param title show title of plot. Defaults to `TRUE`.
#' @param legend shows legend of plot. Defaults to `TRUE`.
#' @param grid show major and minor axis lines. Defaults to `TRUE`.
#' @param counts displays the number of differentially expressed genes for each treatment comparison. Defaults to `TRUE`.
#' @param data.return returns data output of plot if set to `TRUE`. Defaults to `FASLSE`.
#' 
#' @export
#' 
#' @examples
#' # Cuffdiff example
#' data("df.cuff")
#' vsVolcanoMatrix(data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
#'                 padj = 0.05, x.lim = NULL, lfc = 2, title = TRUE, 
#'                 grid = TRUE, counts = TRUE, data.return = FALSE)
#' 
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsVolcanoMatrix(data = df.deseq, d.factor = 'condition', type = 'deseq', 
#'                 padj = 0.05, x.lim = NULL, lfc = 2, title = TRUE, 
#'                 grid = TRUE, counts = TRUE, data.return = FALSE)
#' 
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsVolcanoMatrix(data = df.edger, d.factor = NULL, type = 'edger', 
#'                 padj = 0.05, x.lim = NULL, lfc = 2, title = TRUE, 
#'                 grid = TRUE, counts = TRUE, data.return = FALSE)
#'                 
#' # Extract data frame from visualization
#' data("df.cuff")
#' tmp <- vsVolcanoMatrix(data = df.cuff, d.factor = NULL, 
#'                        type = 'cuffdiff', padj = 0.05, x.lim = NULL,
#'                        lfc = 2, title = TRUE, grid = TRUE, 
#'                        counts = TRUE, data.return = TRUE)
#' df.vmat <- tmp[[1]]
#' head(df.vmat)

vsVolcanoMatrix <- function(data, d.factor = NULL, type, padj = 0.1, 
                            x.lim = NULL, lfc = NULL, title = TRUE, 
                            legend = TRUE, grid = TRUE, counts = TRUE,
                            data.return = FALSE) {
  if (missing(type)) {
    stop('Please specify analysis type ("cuffdiff", "deseq", or "edger")')
  }
  if(type == 'cuffdiff'){
    dat <- .getCuffVolcanoMatrix(data)
  } else if (type == 'deseq') {
    dat <- .getDeseqVolcanoMatrix(data, d.factor)
  } else if (type == 'edger') {
    dat <- .getEdgeVolcanoMatrix(data)
  } else {
    stop('Please enter correct analysis type.')
  }
  if (!isTRUE(title)) {
    m.lab <- NULL
  } else {
    m.lab  <- ggtitle('Volcano Matrix')
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
  
  px <- dat$logFC
  p <- padj
  
  if (is.null(x.lim)) {
    x.lim = c(-1, 1) * quantile(abs(px[is.finite(px)]), probs = 0.99) * 0.8
  }
  if (is.null(lfc)) {
    lfc = 1
  }
  
  dat <- .vomat.ranker(dat, padj, lfc, x.lim)
  
  pc <- .vomat.comp(padj, lfc)
  
  tmp.l <- .vomat.col.count(dat)
  
  if (isTRUE(counts)) {
    b.count <- annotate('text', x = -Inf, y = Inf, vjust = 1.5, hjust = -1, 
                        label = tmp.l$blue$Freq, color = 'royalblue1', 
                        fontface = 2)
    g.count <- annotate(
      'text', x = Inf, y = Inf, vjust = 1.5, hjust = 1.5, 
      label = tmp.l$green$Freq, color = 'green', fontface = 2
    )
  } else {
    b.count <- NULL
    g.count <- NULL
  }
  
  logFC <- pval <- color <- size <- shape <- NULL
  tmp.plot <- ggplot(dat, aes(x = logFC, y = -log10(pval))) +
    geom_point(aes(color = color, size = size, shape = shape), 
               na.rm = TRUE) +
    scale_color_manual(name = '',
                       values = c('grey' = 'grey73', 
                                  'blue' = 'royalblue1', 
                                  'green' = 'green'),
                       labels = c('grey' = pc$gry, 'blue' = pc$blu, 
                                  'green' = pc$grn)) +
    scale_shape_manual(name = '',
                       values = c('circle' = 16, 
                                  'l.triangle' = 60,
                                  'r.triangle' = 62),
                       guide = 'none') +
    scale_size_manual(name = '',
                      values = c('sub' = 1, 
                                 't1' = 1.5,
                                 't2' = 2,
                                 't3' = 3,
                                 't4' = 5),
                      guide = 'none') +
    theme_bw() +
    facet_wrap(id_x ~ id_y) +
    pc$vline1 + pc$vline2 + pc$vline3 + pc$hline1 + pc$x.lab + pc$y.lab +
    b.count + g.count + xlim(x.lim) + m.lab + leg
  
  if (isTRUE(data.return)) {
    plot.l <- list(data = dat, plot = tmp.plot)
  } else {
    print(tmp.plot)
  }
}
