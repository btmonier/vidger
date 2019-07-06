#' @title
#' Volcano plot matrix from \eqn{log_{2}} fold changes and
#' \eqn{-log_{10}}(\eqn{p}-values)
#'
#' @author
#' Brandon Monier, \email{brandon.monier@sdstate.edu}
#'
#' @description
#' This function allows you to extract necessary results-based data from
#' either a \code{DESeq2} object, \code{edgeR} object, or \code{cuffdiff} data
#' frame to create a volcano plot (i.e. a scatter plot) of the negative
#' \eqn{log} of the \eqn{p}-value versus the \eqn{log} of the fold change
#' while implementing ggplot2 aesthetics for all possible combinations of
#' treatments.
#'
#' @param data a cuffdiff, DESeq2, or edgeR object.
#' @param d.factor a specified factor; for use with DESeq2 objects only.
#'  Defaults to \code{NULL}
#' @param type an analysis classifier to tell the function how to process the
#'  data. Must be either \code{cuffdiff}, \code{deseq}, or \code{edgeR}.
#' @param padj a user defined adjusted \eqn{p}-value cutoff point.
#'  Defaults to \code{0.05}.
#' @param x.lim set manual limits to the x axis. Defaults to \code{NULL}.
#' @param lfc \eqn{log} fold change level for setting conditonals. If no user
#'  input is added (\code{NULL}), value defaults to \code{1}.
#' @param title show title of plot. Defaults to \code{TRUE}.
#' @param legend shows legend of plot. Defaults to \code{TRUE}.
#' @param grid show major and minor axis lines. Defaults to \code{TRUE}.
#' @param counts displays the number of differentially expressed genes for
#'  each treatment comparison. Defaults to \code{TRUE}.
#' @param data.return returns data output of plot if set to \code{TRUE}.
#'  Defaults to \code{FALSE}.
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
#' vsVolcanoMatrix(
#'      data = df.cuff, d.factor = NULL, type = "cuffdiff",
#'      padj = 0.05, x.lim = NULL, lfc = 2, title = TRUE,
#'      grid = TRUE, counts = TRUE, data.return = FALSE
#' )
#'
#' # DESeq2 example
#' data("df.deseq")
#' require(DESeq2)
#' vsVolcanoMatrix(
#'      data = df.deseq, d.factor = "condition", type = "deseq",
#'      padj = 0.05, x.lim = NULL, lfc = 2, title = TRUE,
#'      grid = TRUE, counts = TRUE, data.return = FALSE
#' )
#'
#' # edgeR example
#' data("df.edger")
#' require(edgeR)
#' vsVolcanoMatrix(
#'      data = df.edger, d.factor = NULL, type = "edger",
#'      padj = 0.05, x.lim = NULL, lfc = 2, title = TRUE,
#'      grid = TRUE, counts = TRUE, data.return = FALSE
#' )
#'
#' # Extract data frame from visualization
#' data("df.cuff")
#' tmp <- vsVolcanoMatrix(
#'      data = df.cuff, d.factor = NULL,
#'      type = "cuffdiff", padj = 0.05, x.lim = NULL,
#'      lfc = 2, title = TRUE, grid = TRUE,
#'      counts = TRUE, data.return = TRUE
#' )
#' df_vmat <- tmp[[1]]
#' head(df_vmat)
#'
#' # Show plot from object (see prior example for more details)
#' tmp[[2]] ## or use tmp$plot

vsVolcanoMatrix <- function(
    data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"),
    padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, legend = TRUE,
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
        dat <- .getCuffVolcanoMatrix(data)
    } else if (type == "deseq") {
        dat <- .getDeseqVolcanoMatrix(data, d.factor)
    } else if (type == "edger") {
        dat <- .getEdgeVolcanoMatrix(data)
    }

    if (!isTRUE(title)) {
        m.lab <- NULL
    } else {
        m.lab <- ggtitle("Volcano Matrix")
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
    tmp.l_blue <- data.frame(
        do.call('rbind', strsplit(as.character(tmp.l$blue$Var1)," ",fixed=TRUE))
    )
    tmp.l_green <- data.frame(
        do.call('rbind', strsplit(as.character(tmp.l$green$Var1)," ",fixed=TRUE))
    )
    tmp.l$blue$id_x <- tmp.l_blue$X1
    tmp.l$blue$id_y <- tmp.l_blue$X2
    tmp.l$green$id_x <- tmp.l_green$X1
    tmp.l$green$id_y <- tmp.l_green$X2

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
        # b.count <- annotate(
        #     "text",
        #     x = -Inf,
        #     y = Inf,
        #     vjust = 1.5,
        #     hjust = -1,
        #     label = tmp.l$blue$Freq,
        #     color = "royalblue1",
        #     fontface = 2
        # )
        # g.count <- annotate(
        #     "text",
        #     x = Inf,
        #     y = Inf,
        #     vjust = 1.5,
        #     hjust = 1.5,
        #     label = tmp.l$green$Freq,
        #     color = "green",
        #     fontface = 2
        # )
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

    logFC <- pval <- color <- size <- shape <- NULL
    tmp.plot <- ggplot(dat, aes(x = logFC, y = -log10(pval))) +
        geom_point(aes(
            color = color, size = size, shape = shape), na.rm = TRUE
        ) +
        scale_color_manual(
            name = "",
            values = c(
                "grey" = "grey73",
                "blue" = "royalblue1",
                "green" = "green"
            ),
            labels = c(
                "grey" = pc$gry,
                "blue" = pc$blu,
                "green" = pc$grn)
        ) +
        scale_shape_manual(
            name = "",
            values = c("circle" = 16, "l.triangle" = 60, "r.triangle" = 62),
            guide = "none"
        ) +
        scale_size_manual(
            name = "",
            values = c(
                "sub" = 1,
                "t1" = 1.5,
                "t2" = 2,
                "t3" = 3,
                "t4" = 5
            ),
            guide = "none"
        ) +
        theme_bw() +
        facet_wrap(id_x ~ id_y) +
        pc$vline1 + pc$vline2 + pc$vline3 + pc$hline1 + pc$x.lab + pc$y.lab +
        b.count + g.count + xlim(x.lim) + m.lab + leg + text.size

    if (isTRUE(data.return)) {
        dat2 <- dat[, -ncol(dat)]
        plot.l <- list(data = dat2, plot = tmp.plot)
    } else {
        print(tmp.plot)
    }
}
