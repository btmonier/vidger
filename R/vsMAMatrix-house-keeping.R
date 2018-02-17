#-----------------------------------------------------#
# Title:  ViDGER - House Keeping - MA Plot Matrix     #
# Author: Brandon Monier (brandon.monier@sdstate.edu) #
# Date:   02.16.2018                                  #
#-----------------------------------------------------#

.mamat.comp <- function(padj, lfc) {
    # Color
    gry <- paste0('padj > ', padj)
    blu <- paste0('padj < ', padj, ' & |lfc| > ', lfc)
    grn <- paste0('padj < ', padj, ' & |lfc| < ', lfc)
    
    # Lines and labels
    vline1 <- geom_hline(
        yintercept = 0, 
        color = 'red3', 
        size = 0.5, 
        alpha = 0.8, 
        linetype = 'longdash'
    )
    vline2 <- geom_hline(
        yintercept = -lfc, 
        color = 'grey32', 
        size = 0.5, 
        alpha = 0.8, 
        linetype = 'dashed'
    )
    vline3 <- geom_hline(
        yintercept = lfc, 
        color = 'grey32', 
        size = 0.5,
        alpha = 0.8, 
        linetype = 'dashed'
    )
    x.lab <- xlab('mean expression [log(x)]')
    y.lab <- ylab(expression(paste('log'['2'], ' fold change')))
    comp.l <- list(
        gry = gry, 
        blu = blu, 
        grn = grn,
        vline1 = vline1, 
        vline2 = vline2, 
        vline3 = vline3,
        x.lab = x.lab, 
        y.lab = y.lab
    )
}



.mamat.ranker <- function(data, padj, lfc, y.lim) {
    dat <- data
    # Color
    dat$color <- 'grey'
    dat$color[dat$padj <= padj & abs(dat$M) > lfc] <- 'blue'
    dat$color[dat$padj <= padj & abs(dat$M) < lfc] <- 'green'
    # Size
    dat$size <- .ma.out.ranker(dat$M, y.lim[2])
    # Shape
    dat$shape <- 'circle'
    dat$shape[dat$M < y.lim[1]] <- 'l.triangle'
    dat$shape[dat$M > y.lim[2]] <- 'r.triangle'
    return(dat)
}



.mamat.col.count <- function(data) {
    tab <- as.data.frame(t(table(data$color, paste(data$id_x, data$id_y))))
    b.l <- tab[which(tab$Var2 == 'blue'), ]
    g.l <- tab[which(tab$Var2 == 'green'), ]
    col.l <- list(blue = b.l, green = g.l)
    return(col.l)
}



.getEdgeMAMatrix <- function(data) {
    v_1 <- as.vector(unique(data$sample$group))
    m_a <- expand.grid(v_1, v_1)
    m_a <- as.matrix(m_a[which(m_a$Var1 != m_a$Var2), ])
    l_a <- split(m_a, row(m_a))
    
    l1 <- list()
    for(i in 1:length(l_a)){
        l1[[i]] <- .getEdgeMA(l_a[[i]][1], l_a[[i]][2], data)
        l1[[i]]$id_x <- l_a[[i]][1]
        l1[[i]]$id_y <- l_a[[i]][2]
    }
    dat1 <- do.call('rbind', l1)
    # dat1 <- dat1[, c(6:7, 1:2, 3:5)]
    # dat1$id_x <- as.factor(dat1$id_x)
    # dat1$id_y <- as.factor(dat1$id_y)
    return(dat1)
}



.getCuffMAMatrix <- function(data) {
    dat <- data
    v_1 <- union(dat$sample_1, dat$sample_2)
    m_a <- expand.grid(v_1, v_1)
    m_a <- as.matrix(m_a[which(m_a$Var1 != m_a$Var2), ])
    l_a <- split(m_a, row(m_a))
    
    l1 <- list()
    for(i in 1:length(l_a)){
        l1[[i]] <- .getCuffMA(l_a[[i]][1], l_a[[i]][2], data)
        l1[[i]]$id_x <- l_a[[i]][1]
        l1[[i]]$id_y <- l_a[[i]][2]
    }
    dat1 <- do.call('rbind', l1)



    # dat1 <- dat1[, c(1, 7:8, 2:6)]
    dat1$id_x <- as.factor(dat1$id_x)
    dat1$id_y <- as.factor(dat1$id_y)
    return(dat1)
}



.getDeseqMAMatrix <- function(data, d.factor) {
    if(is.null(d.factor)) {
        stop(
            'This appears to be a DESeq object. 
            Please state d.factor variable.'
        )
    }
    
    dat <- as.data.frame(colData(data))
    v_1 <- as.vector(unique(dat[[d.factor]]))
    m_a <- expand.grid(v_1, v_1)
    m_a <- as.matrix(m_a[which(m_a$Var1 != m_a$Var2), ])
    l_a <- split(m_a, row(m_a))
    
    l1 <- list()
    for(i in 1:length(l_a)) {
        l1[[i]] <- .getDeseqMA(l_a[[i]][1], l_a[[i]][2], data, d.factor)
        l1[[i]]$id_x <- l_a[[i]][1]
        l1[[i]]$id_y <- l_a[[i]][2]
    }
    dat1 <- do.call('rbind', l1)
    # dat1 <- dat1[, c(1:2, 7:8, 3:6)]
    # dat1$id_x <- as.factor(dat1$id_x)
    # dat1$id_y <- as.factor(dat1$id_y)
    return(dat1)
}
