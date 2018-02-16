context("scatter plots")

test_that("scatter plots work", {
    # Load data
    data(df.cuffdiff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsScatterPlot(x = 'hESC', y = 'iPS', data = df.cuff, type = 'cuffdiff',
                  d.factor = NULL, title = TRUE, grid = TRUE)

    # DESeq2
    vsScatterPlot(x = 'treated', y = 'untreated', data = df.deseq, 
                  type = 'deseq', d.factor = 'condition', title = TRUE, 
                  grid = TRUE)

    # edgeR
    vsScatterPlot(x = 'WM', y = 'MM', data = df.edger, type = 'edger',
                  d.factor = NULL, title = TRUE, grid = TRUE)
})