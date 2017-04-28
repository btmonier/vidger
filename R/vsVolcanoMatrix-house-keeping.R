#'-----------------------------------------------------#
#' Title:  ggviseq - House Keeping - volcano matrix    #
#' Author: Brandon Monier (brandon.monier@sdstate.edu) #
#' Date:   04.27.17                                    #
#'-----------------------------------------------------#

#'---------
#' Preamble
#'---------

#'...



#'--------------------------
#' Volcano matrix components
#'--------------------------

#' @export
vomat.comp <- function(padj, lfc) {
  # Color
  gry <- paste0('padj > ', padj)
  blu <- paste0('padj < ', padj, ' & |lfc| > ', lfc)
  grn <- paste0('padj < ', padj, ' & |lfc| < ', lfc)
  
  # Lines and labels
  vline1 <- geom_vline(xintercept = 0, color = 'red3', size = 0.5, 
                       alpha = 0.8, linetype = 'longdash')
  vline2 <- geom_vline(xintercept = -lfc, color = 'grey32', size = 0.5, 
                       alpha = 0.8, linetype = 'dashed')
  vline3 <- geom_vline(xintercept = lfc, color = 'grey32', size = 0.5,
                       alpha = 0.8, linetype = 'dashed')
  hline1 <- geom_hline(yintercept = -log10(padj), color = 'grey32', size = 0.5, 
                       alpha = 0.8, linetype = 'dashed')
  x.lab <- xlab(expression(paste('log'['2'], ' fold change')))
  y.lab <- ylab(expression(paste('-log'['10'], '(p-value)')))
  comp.l <- list(gry = gry, blu = blu, grn = grn,
                 vline1 = vline1, vline2 = vline2, vline3 = vline3, 
                 hline1 = hline1, x.lab = x.lab, y.lab = y.lab)
}

#' @export
vomat.ranker <- function(data, padj, lfc, x.lim) {
  dat <- data
  # Color
  dat$color <- 'grey'
  dat$color[dat$padj <= padj & abs(dat$logFC) > lfc] <- 'blue'
  dat$color[dat$padj <= padj & abs(dat$logFC) < lfc] <- 'green'
  # Size
  dat$size <- vo.out.ranker(dat$logFC, x.lim[2])
  # Shape
  dat$shape <- 'circle'
  dat$shape[dat$logFC < x.lim[1]] <- 'l.triangle'
  dat$shape[dat$logFC > x.lim[2]] <- 'r.triangle'
  return(dat)
}

#' @export
vomat.col.count <- function(data) {
  tab <- as.data.frame(t(table(data$color, paste(data$id_x, data$id_y))))
  b.l <- tab[which(tab$Var2 == 'blue'), ]
  g.l <- tab[which(tab$Var2 == 'green'), ]
  col.l <- list(blue = b.l, green = g.l)
  return(col.l)
}



#'------------------------
#' Volcano plot extraction
#'------------------------

#' edgeR - REQUIRES `getEdgeScatter()`
#' @export
getEdgeVolcano <- function(x, y, data) {
  dat <- getEdgeScatter(x, y, data)
  deg <- exactTest(data, pair = c(x, y))
  deg <- topTags(deg, n = nrow(deg$table))
  deg <- deg$table[order(as.numeric(rownames(deg$table))),]
  dat$logFC <- deg$logFC
  dat$pval <- deg$PValue
  dat$padj <- deg$FDR
  return(dat)
}


#' cuffdiff - NO PREREQUISITES
#' @export
getCuffVolcano <- function(x, y, data) {
  deg <- data
  deg <- subset(deg, (sample_1 == x & sample_2 == y) | 
                  (sample_1 == y & sample_2 == x))
  dat <- data.frame(test_id = deg$test_id)
  if (x %in% deg$sample_1 & y %in% deg$sample_2) {
    dat$x <- deg$value_1
    dat$y <- deg$value_2
  } else if (y %in% deg$sample_1 & x %in% deg$sample_2) {
    dat$x <- deg$value_2
    dat$y <- deg$value_1
  }
  dat$logFC <- log2(dat$y / dat$x)
  dat$pval <- deg$p_value
  dat$padj <- deg$q_value
  dat[is.na(dat)] <- 0
  return(dat)
}


#' DESeq2 - NO PREREQUISITES
#' @export
getDeseqVolcano <- function(x, y, data, d.factor) {
  if(missing(d.factor)) {
    stop('This appears to be a DESeq object. Please state d.factor variable.')
  }
  dat1 <- as.data.frame(colData(data))
  dat2 <- fpkm(data)
  dat3 <- results(data, contrast = c(d.factor, y, x))
  nam_x <- row.names(dat1[which(dat1[d.factor] == x),])
  nam_y <- row.names(dat1[which(dat1[d.factor] == y),])
  x <- rowMeans(dat2[, nam_x])
  y <- rowMeans(dat2[, nam_y])
  dat4 <- data.frame(x, y)
  dat4$logFC <- log2(dat4$y / dat4$x)
  dat4$pval <- dat3$pvalue
  dat4$padj <- dat3$padj
  dat4 <- dat4[complete.cases(dat4),]
  return(dat4)
}