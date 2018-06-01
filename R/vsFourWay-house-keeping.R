#-----------------------------------------------------#
# Title:  ggviseq - House Keeping - four way plot     #
# Author: Brandon Monier (brandon.monier@sdstate.edu) #
# Date:   03.28.17                                    #
#-----------------------------------------------------#


.getEdgeFourWay <- function(x, y, control, data) {
    dat_x <- .getEdgeScatter(control, x, data)
    deg_x <- exactTest(data, pair = c(control, x))
    deg_x <- topTags(deg_x, n = nrow(deg_x$table))
    deg_x <- deg_x$table[order(as.numeric(rownames(deg_x$table))),]
    
    dat_y <- .getEdgeScatter(control, y, data)
    deg_y <- exactTest(data, pair = c(control, y))
    deg_y <- topTags(deg_y, n = nrow(deg_y$table))
    deg_y <- deg_y$table[order(as.numeric(rownames(deg_y$table))),]
    
    dat_all <- data.frame(
        id = dat_x$id, x = dat_x$y, y = dat_y$y, control = dat_x$x,
        logFC_x = deg_x$logFC, pval_x = deg_x$PValue, padj_x = deg_x$FDR,
        logFC_y = deg_y$logFC, pval_y = deg_y$PValue, padj_y = deg_y$FDR
    )
    dat_all <- dat_all[complete.cases(dat_all), ]
    return(dat_all)
}



.getCuffFourWay <- function(x, y, control, data) {
    sample_1 <- sample_2 <- NULL
    deg <- data
    deg_x <- subset(deg, (sample_1 == x & sample_2 == control) | 
                                        (sample_1 == control & sample_2 == x))
    dat_x <- data.frame(id = deg_x$test_id)
    if (x %in% deg_x$sample_1 && control %in% deg_x$sample_2) {
        dat_x$x <- deg_x$value_1
        dat_x$control <- deg_x$value_2
    } else if (control %in% deg_x$sample_1 && x %in% deg_x$sample_2) {
        dat_x$x <- deg_x$value_2
        dat_x$control <- deg_x$value_1
    }
    dat_x$logFC_x <- log2(dat_x$x / dat_x$control)
    dat_x$pval_x <- deg_x$p_value
    dat_x$padj_x <- deg_x$q_value
    
    deg_y <- subset(deg, (sample_1 == y & sample_2 == control) | 
                                        (sample_1 == control & sample_2 == y))
    dat_y <- data.frame(id = deg_y$test_id)
    if (y %in% deg_y$sample_1 && control %in% deg_y$sample_2) {
        dat_y$y <- deg_y$value_1
        dat_y$control <- deg_y$value_2
    } else if (control %in% deg_y$sample_1 && y %in% deg_y$sample_2) {
        dat_y$y <- deg_y$value_2
        dat_y$control <- deg_y$value_1
    }
    dat_y$logFC_y <- log2(dat_y$y / dat_y$control)
    dat_y$pval_y <- deg_y$p_value
    dat_y$padj_y <- deg_y$q_value
    
    dat_all <- data.frame(dat_x, dat_y[, c(2, 4:6)])
    dat_all <- dat_all[, c(1, 2, 7, 3:6, 8:10)]
    dat_all <- dat_all[complete.cases(dat_all), ]
    return(dat_all)
}



.getDeseqFourWay <- function(x, y, control, data, d.factor) {
    if(is.null(d.factor)) {
        stop(
            'This appears to be a DESeq object. 
            Please state d.factor variable.'
        )
    }
    dat1 <- as.data.frame(colData(data))
    dat2 <- fpm(data)
    dat3_x <- results(data, contrast = c(d.factor, x, control))
    dat3_y <- results(data, contrast = c(d.factor, y, control))
    
    nam_x <- row.names(dat1[which(dat1[d.factor] == x),])
    nam_y <- row.names(dat1[which(dat1[d.factor] == y),])
    nam_c <- row.names(dat1[which(dat1[d.factor] == control),])
    
    x <- rowMeans(dat2[, nam_x])
    y <- rowMeans(dat2[, nam_y])
    id <- rownames(dat2)
    control <- rowMeans(dat2[, nam_c])
    dat_all <- data.frame(id, x, y, control)
    dat_all$logFC_x <- log2(dat_all$x / dat_all$control)
    dat_all$pval_x <- dat3_x$pvalue
    dat_all$padj_x <- dat3_x$padj
    dat_all$logFC_y <- log2(dat_all$y / dat_all$control)
    dat_all$pval_y <- dat3_y$pvalue
    dat_all$padj_y <- dat3_y$padj
    dat_all <- dat_all[complete.cases(dat_all), ]
    return(dat_all)
}



.four.comp1 <- function(x.lim, y.lim, padj, all.lab, lfc, b, g, r) {
    list(
        sh1 = paste(
            '|lfc(x)| > ', round(x.lim[2], 1), 'OR  |lfc(y)| > ',
            round(y.lim[2], 1)
        ),
        sh2 = paste(
            round(x.lim[1], 1), '< lfc(x) <', round(x.lim[2], 1), 
            ' AND ', round(y.lim[1], 1), '< lfc(y) <', round(y.lim[2]), 1
        ),
        col1 = paste0('padj(x, y) > ', padj),
        col2 = paste0(
            'padj(x) < ', padj, ', |lfc(x)| > ', lfc, 
            ' |lfc(y)| < ', lfc, ' (', g ,')'
        ),
        col3 = paste0(
            'padj(y) < ', padj, ', |lfc(x)| < ', lfc, ' |lfc(y)| > ',
            lfc, ' (', r ,')'
        ),
        col4 = paste0(
            'padj(x, y) < ', padj, ', |lfc(x, y)| > ', lfc, ' (', b ,')'
        ),
        vline1 = geom_vline(
            xintercept = 0, 
            color = 'red3', 
            size = 0.5, 
            alpha = 0.8, 
            linetype = 'longdash'
        ), 
        vline2 = geom_vline(
            xintercept = -lfc, 
            color = 'grey32', 
            size = 1, 
            alpha = 0.8, 
            linetype = 'dashed'),
        vline3 = geom_vline(
            xintercept = lfc, 
            color = 'grey32', 
            size = 1,
            alpha = 0.8, 
            linetype = 'dashed'
        ),
        hline1 = geom_hline(
            yintercept = 0, 
            color = 'red3', 
            size = 0.5, 
            alpha = 0.8, 
            linetype = 'longdash'
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
        x.lab = xlab(
            bquote(
                .(all.lab[1]) ~ 'expression relative to' ~ 
                .(all.lab[3]) ~ 'log'['2'] ~ 'fold change'
            )
        ),
        y.lab = ylab(
            bquote(
                .(all.lab[2]) ~ 'expression relative to' ~ 
                .(all.lab[3]) ~ 'log'['2'] ~ 'fold change'
            )
        )    
    )
}



.four.comp2 <- function(a, b, c, d, e, f) {
    list(
        color = scale_color_manual(
            name = '', 
            values = c(
                'gry' = 'grey73', 
                'grn' ='green', 
                'blu' ='royalblue1', 
                'red' = 'red'
            ), 
            labels = c(
                'gry' = a, 
                'grn' = b, 
                'blu' = c, 
                'red' = d
            )
        ),
        shape = scale_shape_manual(
            name = '', 
            values = c(
                'tri1' = 17, 
                'circ' = 16
            ),
            labels = c(
                'tri1' = e, 
                'circ' = f
            )
        ),
        size = scale_size_manual(
            name = '', 
            values = c(
                'sub' = 1,
                't1' = 2, 
                't2' = 3, 
                't3' = 4, 
                't4' = 5
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



.four.out.ranker <- function(px, py, x.lim, y.lim) {
    vec2_x <- px[which(abs(px) >= x.lim)]
    vec2_y <- py[which(abs(py) >= y.lim)]
    vec3 <- c(vec2_x, vec2_y)
    tmp <- quantile(abs(vec3))
    ifelse(
        abs(px) < x.lim & abs(py) < y.lim, 'sub',
        ifelse(
            abs(vec3) >= tmp[[4]], 't4',
            ifelse(
                abs(vec3) >= tmp[[3]] & abs(vec3) < tmp[[4]], 't3',
                ifelse(
                    abs(vec3) >= tmp[[2]] & abs(vec3) < tmp[[3]], 't2','t1'
                )
            )
        )
    ) 
}



.four.col.ranker <- function(dat, lfc) {
    de_x <- dat$isDE_x
    de_y <- dat$isDE_y
    de_a <- dat$isDE_all
    px <- dat$logFC_x
    py <- dat$logFC_y
    
    ifelse(
        de_a == TRUE & abs(px) >= lfc & abs(py) >= lfc, 'blu',
        ifelse(
            abs(px) >= lfc & abs(py) < lfc & 
            de_x == TRUE & de_y == FALSE, 'grn',
            ifelse(
                abs(px) < lfc & abs(py) >= lfc & 
                de_x == FALSE & de_y == TRUE, 'red', 'gry'
            )
        )
    )
}



.four.shp.ranker <- function(px, py, x.lim, y.lim) {
    ifelse(
        (px < x.lim[1]) | 
        (px > x.lim[2]) | 
        (py < y.lim[1]) | 
        (py > y.lim[2]), 
        'tri1', 'circ'
    )
}



.four.ranker <- function(data, padj, lfc, x.lim, y.lim) {
    dat <- data
    dat$color <- 'grey'
    dat$color[dat$padj_x <= padj & dat$padj_y <= padj & 
    abs(dat$logFC_x) > lfc & abs(dat$logFC_y) > lfc] <- 'blue'
    dat$color[dat$padj_x <= padj & dat$padj_y > padj & 
    abs(dat$logFC_x) > lfc & abs(dat$logFC_y) <= lfc] <- 'green'
    dat$color[dat$padj_y <= padj & dat$padj_x > padj & 
    abs(dat$logFC_x) <= lfc & abs(dat$logFC_y) > lfc] <- 'red'
    dat$size <- .four.out.ranker(dat$logFC_x, dat$logFC_y, x.lim[2], y.lim[2])
    dat$shape <- 'circle'
    dat$shape[dat$logFC_x < x.lim[1] | 
    dat$logFC_x > x.lim[2] | 
    dat$logFC_y < y.lim[1] | 
    dat$logFC_y > y.lim[2]] <- 'triangle'
    return(dat)
}



.four.col.counter <- function(dat, lfc) {
    de_x <- dat$isDE_x
    de_y <- dat$isDE_y
    de_a <- dat$isDE_all
    px <- abs(dat$logFC_x)
    py <- abs(dat$logFC_y)
    
    blu.c <- nrow(dat[which(px >= lfc & py > lfc & de_a == TRUE), ])
    grn.c <- nrow(dat[which(px >= lfc & py < lfc & de_x == TRUE & 
        de_y == FALSE), ])
    red.c <- nrow(dat[which(px < lfc & py >= lfc & de_x == FALSE & 
        de_y == TRUE), ])
    
    l.count <- list(blu.c, grn.c, red.c)
    return(l.count)
}
