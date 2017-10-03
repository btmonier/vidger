#'-----------------------------------------------------#
#' Title:  ggviseq - House Keeping - scatter plots     #
#' Author: Brandon Monier (brandon.monier@sdstate.edu) #
#' Date:   04.04.17                                    #
#'-----------------------------------------------------#

#'---------
#' Preamble
#'---------

#'...


#'------------------------
#' Scatter plot extraction
#'------------------------

#' edgeR
#' @export
getEdgeScatter <- function(x, y, data) {
  dat.cpm <- cpm(data$counts)
  tmp.x <- row.names(data$samples[which(data$samples$group == x), ])
  tmp.y <- row.names(data$samples[which(data$samples$group == y), ])
  x <- rowMeans(dat.cpm[, tmp.x])
  y <- rowMeans(dat.cpm[, tmp.y])
  dat.plt <- data.frame(x, y)
  return(dat.plt)
}


#' cuffdiff
#' @export
getCuffScatter <- function(x, y, data) {
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
  return(dat)
}


#' DESeq2
#' @export
getDeseqScatter <- function(x, y, data, d.factor) {
  if(is.null(factor)) {
    stop('This appears to be a DESeq object. Please state factor variable.')
  }
  dat1 <- as.data.frame(colData(data))
  dat2 <- fpm(data)
  nam_x <- row.names(dat1[which(dat1[d.factor] == x),])
  nam_y <- row.names(dat1[which(dat1[d.factor] == y),])
  x <- rowMeans(dat2[, nam_x])
  y <- rowMeans(dat2[, nam_y])
  dat3 <- data.frame(x, y)
  return(dat3)
}















