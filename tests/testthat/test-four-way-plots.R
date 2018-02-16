context("four way plots")

test_that("four way plots work", {
    # Load data
    data(df.cuffdiff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsFourWay(x = 'iPS', y = 'hESC', control = 'Fibroblasts', data = df.cuff,
              d.factor = NULL, type = 'cuffdiff', padj = 0.05, x.lim = NULL,
              y.lim = NULL, lfc = NULL, legend = TRUE, title = TRUE, 
              grid = TRUE)

    # DESeq2

    # edgeR
    vsFourWay(x = 'WW', y = 'WM', control = 'MM', data = df.edger,
              d.factor = NULL, type = 'edger', padj = 0.05, x.lim = NULL,
              y.lim = NULL, lfc = NULL, legend = TRUE, title = TRUE, 
              grid = TRUE)
})