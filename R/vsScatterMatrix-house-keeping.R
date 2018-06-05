#--------------------------------------------------------#
# Title:  ggviseq - House Keeping - scatterplot matrix   #
# Author: Brandon Monier (brandon.monier@sdstate.edu)    #
# Date:   04.10.17                                       #
#--------------------------------------------------------#

.getCuffScatterMatrix <- function(data) {
    key <- value <- NULL
    dat1 <- data[, c('test_id', 'sample_1', 'value_1')]
    dat2 <- data[, c('test_id', 'sample_2', 'value_2')]
    names(dat1) <- c('id', 'key', 'value')
    names(dat2) <- c('id', 'key', 'value')
    dat3 <- rbind(dat1, dat2)
    dat3 <- unique(dat3)
    dat3$key <- as.factor(dat3$key)
    dat3 <- tidyr::spread(dat3, key, value)
    dat4 <- dat3[, -1]
    row.names(dat4) <- dat3[, 1]
    return(dat4)
}



.getDeseqScatterMatrix <- function(data, d.factor = NULL) {
    if(is.null(d.factor)) {
        stop(
            paste(
                "This appears to be a DESeq object.",
                "Please state \"d.factor\" variable."
            )
        )
    }
    dat1 <- as.data.frame(colData(data))
    dat2 <- fpm(data)
    nam <- as.vector(unique(dat1[[d.factor]]))
    ls.nam <- lapply(nam, function(i) {
        row.names(dat1[which(dat1[d.factor] == i), ])
    })
    ls.mean <- lapply(seq_along(ls.nam), function(i) {
        rowMeans(dat2[, ls.nam[[i]]])
    })
    names(ls.mean) <- sapply(nam, paste)
    dat3 <- as.data.frame(ls.mean)
    return(dat3)
}



.getEdgeScatterMatrix <- function(data) {
    dat.cpm <- cpm(data$counts)
    nam <- as.vector(unique(data$sample$group))

    ls.nam <- lapply(nam, function(i) {
        row.names(data$samples[which(data$samples$group == i), ])
    })
    ls.mean <- lapply(seq_along(ls.nam), function(i) {
        rowMeans(dat.cpm[, ls.nam[[i]]])
    })

    names(ls.mean) <- sapply(nam, paste)
    dat <- as.data.frame(ls.mean)
    return(dat)
}
