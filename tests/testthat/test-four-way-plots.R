context("four way plots")

test_that("four way plots work", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsFourWay(x = 'iPS', y = 'hESC', control = 'Fibroblasts', data = df.cuff,
              d.factor = NULL, type = 'cuffdiff', padj = 0.05, x.lim = NULL,
              y.lim = NULL, lfc = NULL, legend = TRUE, title = TRUE, 
              grid = TRUE)

    # DESeq2
    vsFourWay(x = 'treated_paired.end', y = 'untreated_paired.end', 
              control = 'untreated_single.read', data = df.deseq, 
              d.factor = 'condition', type = 'deseq', padj = 0.05, 
              x.lim = NULL, y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
              data.return = FALSE)

    # edgeR
    vsFourWay(x = 'WW', y = 'WM', control = 'MM', data = df.edger,
              d.factor = NULL, type = 'edger', padj = 0.05, x.lim = NULL,
              y.lim = NULL, lfc = NULL, legend = TRUE, title = TRUE, 
              grid = TRUE)
})