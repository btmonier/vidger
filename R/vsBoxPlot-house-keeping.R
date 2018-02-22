#-----------------------------------------------------#
# Title:  ggviseq - Data Extractors - Box plots       #
# Author: Brandon Monier (brandon.monier@sdstate.edu) #
# Date:   04.04.17                                    #
#-----------------------------------------------------#

.getEdgeBox <- function(data) {
    dat.cpm <- cpm(data$counts)
    tmp <- as.vector(unique(data$sample$group))
    ls1 <- list()
    ls2 <- list()
    for (i in tmp) {
        ls1[[i]] <- row.names(data$samples[which(data$samples$group == i), ])
        for (j in 1:length(ls1)) {
            ls2[[j]] <- rowMeans(dat.cpm[, ls1[[j]]])
        }
    }
    names(ls2) <- sapply(tmp, paste)
    dat <- as.data.frame(ls2)
    dat <- tidyr::gather(as.data.frame(dat))
    dat$key <- as.factor(dat$key)
    return(dat)
}



.getCuffBox <- function(data) {
    dat1 <- data[, c('test_id', 'sample_1', 'value_1')]
    dat2 <- data[, c('test_id', 'sample_2', 'value_2')]
    names(dat1) <- c('id', 'key', 'value')
    names(dat2) <- c('id', 'key', 'value')
    dat3 <- rbind(dat1, dat2)
    dat3 <- unique(dat3)
    dat3$key <- as.factor(dat3$key)
    return(dat3)
}



.getDeseqBox <- function(data, d.factor) {
    if(is.null(d.factor)) {
        stop(
            'This appears to be a DESeq object. 
            Please state d.factor variable.'
        )
    }
    dat1 <- as.data.frame(colData(data))
    dat2 <- fpm(data)
    nam <- as.vector(unique(dat1[[d.factor]]))
    ls.nam <- list()
    ls.mean <- list()
    for (i in nam) {
        ls.nam[[i]] <- row.names(dat1[which(dat1[d.factor] == i), ])
        for (j in 1:length(ls.nam)){
            ls.mean[[j]] <- rowMeans(dat2[, ls.nam[[j]]])
        }
    }
    names(ls.mean) <- sapply(nam, paste)
    dat3 <- as.data.frame(ls.mean)
    dat3 <- tidyr::gather(as.data.frame(dat3))
    dat3$key <- as.factor(dat3$key)
    return(dat3)
}
