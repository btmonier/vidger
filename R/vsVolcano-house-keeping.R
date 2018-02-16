#-----------------------------------------------------#
# Title:  ggviseq - Data Extractors - Volcano Plots   #
# Author: Brandon Monier (brandon.monier@sdstate.edu) #
# Date:   04.04.17                                    #
#-----------------------------------------------------#

.vo.comp1 <- function(x.lim, padj, lfc, b, g) {
    list(
        sh1 = paste('lfc < ', round(x.lim[1], 2)),
        sh2 = paste('lfc > ', round(x.lim[2], 2)),
        sh3 = paste(round(x.lim[1], 2), '< lfc <', round(x.lim[2], 2)),
        col1 = paste('padj >', padj),
        col2 = paste('padj <', padj, '; |lfc| >', lfc, ' (', b ,')'),
        col3 = paste('padj <', padj, '; |lfc| <', lfc, ' (', g ,')'),
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
            linetype = 'dashed'
        ),
        vline3 = geom_vline(
            xintercept = lfc, 
            color = 'grey32', 
            size = 1,
            alpha = 0.8, 
            linetype = 'dashed'
        ),
        hline1 = geom_hline(
            yintercept = -log10(padj), 
            color = 'grey32', 
            size = 1, 
            alpha = 0.8, 
            linetype = 'dashed'
        ),
        x.lab = xlab(expression(paste('log'['2'], ' fold change'))),
        y.lab = ylab(expression(paste('-log'['10'], '(p-value)')))
             
    )
}



.vo.comp2 <- function(a, b, c, d, e, f){
    list(
        color = scale_color_manual(
            name = '', 
            values = c('gry' = 'grey73', 'grn' ='green', 'blu' ='royalblue1'), 
            labels = c('gry' = a, 'grn' = b, 'blu' = c)
        ),
        shape = scale_shape_manual(
            name = '', 
            values = c('tri1' = 60, 'tri2' = 62, 'circ' = 16),
            labels = c('tri1' = d, 'tri2' = e, 'circ' = f)
        ),
        size = scale_size_manual(
            name = '', 
            values = c(
                'sub' = 1,
                't1' = 2.5, 
                't2' = 3.5, 
                't3' = 4.5, 
                't4' = 5.5
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



.vo.out.ranker <- function(log2fc, lim) {
    vec2 <- log2fc[which(abs(log2fc) >= lim)]
    tmp <- quantile(abs(vec2))
    ifelse(
        abs(log2fc) < lim, 'sub',
        ifelse(
            abs(vec2) >= tmp[[4]], 't4',
            ifelse(
                abs(vec2) >= tmp[[3]] & abs(vec2) < tmp[[4]], 't3',
                ifelse(
                    abs(vec2) >= tmp[[2]] & abs(vec2) < tmp[[3]], 't2', 't1')
            )
        )
    ) 
}



.vo.col.ranker <- function(isDE, log2fc, lfc) {
    ifelse(isDE == TRUE & abs(log2fc) < lfc, 'grn', 
                 ifelse(isDE == TRUE & abs(log2fc) > lfc, 'blu', 'gry'))
}



.vo.shp.ranker <- function(log2fc, lim){
    ifelse(log2fc < lim[1], 'tri1', ifelse(log2fc > lim[2], 'tri2', 'circ'))
}



.vo.ranker <- function(data, padj, lfc, x.lim) {
    dat <- data
    dat$color <- 'grey'
    dat$color[dat$padj <= padj & abs(dat$logFC) > lfc] <- 'blue'
    dat$color[dat$padj <= padj & abs(dat$logFC) < lfc] <- 'green'
    dat$size <- .vo.out.ranker(dat$logFC, x.lim[2])
    dat$shape <- 'circle'
    dat$shape[dat$logFC < x.lim[1]] <- 'l.triangle'
    dat$shape[dat$logFC > x.lim[2]] <- 'r.triangle'
    return(dat)
}



.vo.col.counter <- function(dat, lfc) {
    de <- dat$isDE
    px <- abs(dat$logFC)
    
    blu.c <- nrow(dat[which(px >= lfc & de == TRUE), ])
    grn.c <- nrow(dat[which(px < lfc & de == TRUE), ])
    
    l.count <- list(blu.c, grn.c)
    return(l.count)
}



.getEdgeVolcano <- function(x, y, data) {
    dat <- .getEdgeScatter(x, y, data)
    deg <- exactTest(data, pair = c(x, y))
    deg <- topTags(deg, n = nrow(deg$table))
    deg <- deg$table[order(as.numeric(rownames(deg$table))),]
    dat$logFC <- deg$logFC
    dat$pval <- deg$PValue
    dat$padj <- deg$FDR
    return(dat)
}



.getCuffVolcano <- function(x, y, data) {
    sample_1 <- sample_2 <- NULL
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
    dat <- dat[complete.cases(dat),]  
    return(dat)
}



.getDeseqVolcano <- function(x, y, data, d.factor) {
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
    x <- rowMeans(dat2[, nam_x])
    y <- rowMeans(dat2[, nam_y])
    dat4 <- data.frame(x, y)
    dat4$logFC <- log2(dat4$y / dat4$x)
    dat4$pval <- dat3$pvalue
    dat4$padj <- dat3$padj
    dat4 <- dat4[complete.cases(dat4),]
    return(dat4)
}
