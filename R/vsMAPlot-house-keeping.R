#-----------------------------------------------------#
# Title:  ggviseq - House Keeping - MA plots          #
# Author: Brandon Monier (brandon.monier@sdstate.edu) #
# Date:   04.04.17                                    #
#-----------------------------------------------------#

.getEdgeMA <- function(x, y, data) {
    dat <- .getEdgeScatter(x, y, data)
    deg <- exactTest(data, pair = c(x, y))
    deg <- topTags(deg, n = nrow(deg$table))
    deg <- deg$table[order(as.numeric(rownames(deg$table))),]
    dat$x <- log10(dat$x + 1)
    dat$y <- log10(dat$y + 1)
    dat$A <- rowMeans(data$counts)
    dat$A <- log10(dat$A)
    # dat$A <- 0.5 * log2(dat$x * dat$y)
    dat$M <- deg$logFC
    dat$pval <- deg$PValue
    dat$padj <- deg$FDR
    dat <- dat[complete.cases(dat), ]
    return(dat)
}



.getCuffMA <- function(x, y, data) {
    sample_1 <- sample_2 <- NULL
    deg <- data
    deg <- subset(deg, (sample_1 == x & sample_2 == y) | 
                                    (sample_1 == y & sample_2 == x))
    dat <- data.frame(test_id = deg$test_id)
    
    if (x %in% deg$sample_1 && y %in% deg$sample_2) {
        dat$x <- deg$value_1
        dat$y <- deg$value_2
    } else if (y %in% deg$sample_1 && x %in% deg$sample_2) {
        dat$x <- deg$value_2
        dat$y <- deg$value_1
    }
    dat$x <- log10(dat$x + 1)
    dat$y <- log10(dat$y + 1)
    # dat$A <- 0.5 * log2(dat$x * dat$y)
    dat$A <- 0.5 * (dat$x + dat$y)
    # dat$M <- deg[, "log2.fold_change."]
    dat$M <- log2(dat$x / dat$y)
    dat$pval <- deg$p_value
    dat$padj <- deg$q_value
    dat <- do.call(
        data.frame, lapply(dat, function(x) replace(x, is.infinite(x), NA))
    )

    dat <- dat[complete.cases(dat), ]
    return(dat)
}



.getDeseqMA <- function(x, y, data, d.factor) {
    if(missing(d.factor)) {
        stop(
            'This appears to be a DESeq object. 
            Please state d.factor variable.'
        )
    }
    dat1 <- as.data.frame(colData(data))
    dat2 <- fpm(data)
    dat3 <- results(data, contrast = c(d.factor, y, x))
    nam_x <- row.names(dat1[which(dat1[d.factor] == x),])
    nam_y <- row.names(dat1[which(dat1[d.factor] == y),])
    x <- log10(rowMeans(dat2[, nam_x]) + 1)
    y <- log10(rowMeans(dat2[, nam_y]) + 1)
    dat4 <- data.frame(x, y)
    dat4$A <- log10(dat3$baseMean)
    # dat4$A <- 0.5 * log2(dat4$x * dat4$y)
    dat4$M <- dat3$log2FoldChange
    dat4$pval <- dat3$pvalue
    dat4$padj <- dat3$padj
    dat4 <- dat4[complete.cases(dat4), ]
    return(dat4)  
}



.ma.comp1 <- function(y.lim, padj, lfc, b, g) {
    list(
        sh1 = paste('lfc < ', round(y.lim[1], 2)),
        sh2 = paste('lfc > ', round(y.lim[2], 2)),
        sh3 = paste(round(y.lim[1], 2), '< lfc <', round(y.lim[2], 2)),
        col1 = paste('padj >', padj),
        col2 = paste('padj <', padj, '; lfc >' , lfc, ' (', b ,')'),
        col3 = paste('padj <', padj, '; lfc <', lfc, ' (', g ,')'),
        hline1 = geom_hline(
            yintercept = 0, color = 'red3', size = 0.5,
            alpha = 0.8, linetype = 'longdash'
        ), 
        hline2 = geom_hline(
            yintercept = -lfc, 
            color = 'grey32', 
            size = 1, 
            alpha = 0.8, 
            linetype = 'dashed'
        ),
        hline3 = geom_hline(
            yintercept = lfc, 
            color = 'grey32', 
            size = 1,
            alpha = 0.8, 
            linetype = 'dashed'
        ),
        x.lab = xlab('mean expression [log(x)]'),
        y.lab = ylab(expression(paste('log'['2'], ' fold change')))
             
    )
}



.ma.comp2 <- function(a, b, c, d, e, f) {
    list(
        color = scale_color_manual(
            name = '', 
            values = c('gry' = 'grey73', 'grn' ='green', 'blu' ='royalblue1'), 
            labels = c('gry' = a, 'grn' = b, 'blu' = c)
        ),
        shape = scale_shape_manual(
            name = '', 
            values = c('tri1' = 6, 'tri2' = 2, 'circ' = 16),
            labels = c('tri1' = d, 'tri2' = e, 'circ' = f)
        ),
        size = scale_size_manual(
            name = '', 
            values = c(
                'sub' = 1, 
                't1' = 1.5, 
                't2' = 2, 
                't3' = 3, 
                't4' = 4
            ),
            labels = c(
                'sub' = 'SUB', 
                't1' = 'T-1', 
                't2' = 'T-2', 
                't3' = 'T-3', 
                't4' = 'T-4'
            )
        )
    )
}



.ma.out.ranker <- function(log2fc, lim) {
    vec2 <- log2fc[which(abs(log2fc) >= lim)]
    tmp <- quantile(abs(vec2))
    ifelse(
        abs(log2fc) < lim, 'sub',
        ifelse(
            abs(vec2) >= tmp[[4]], 't4',
            ifelse(
                abs(vec2) >= tmp[[3]] & abs(vec2) < tmp[[4]], 't3',
                ifelse(
                    abs(vec2) >= tmp[[2]] & abs(vec2) < tmp[[3]], 't2','t1'
                )
            )
        )
    ) 
}



.ma.col.ranker <- function(isDE, log2fc, lfc) {
    ifelse(
        isDE == TRUE & abs(log2fc) < lfc, 'grn',
        ifelse(isDE == TRUE & abs(log2fc) > lfc, 'blu', 'gry')
    )
}



.ma.shp.ranker <- function(log2fc, lim){
    ifelse(log2fc < lim[1], 'tri1', ifelse(log2fc > lim[2], 'tri2', 'circ'))
}



.ma.ranker <- function(data, padj, lfc, y.lim) {
    dat <- data
    dat$color <- 'grey'
    dat$color[dat$isDE == TRUE & abs(dat$M) >= lfc] <- 'blue'
    dat$color[dat$isDE == TRUE & abs(dat$M) < lfc] <- 'green'
    dat$size <- .ma.out.ranker(dat$A, y.lim[2])
    dat$shape <- 'circle'
    dat$shape[dat$logFC < y.lim[1]] <- 'down.triangle'
    dat$shape[dat$logFC > y.lim[2]] <- 'up.triangle'
    return(dat)
}



.ma.col.counter <- function(dat, lfc) {
    de <- dat$isDE
    py <- abs(dat$M)
    
    blu.c <- nrow(dat[which(py >= lfc & de == TRUE), ])
    grn.c <- nrow(dat[which(py < lfc & de == TRUE), ])
    
    l.count <- list(blu.c, grn.c)
    return(l.count)
}
