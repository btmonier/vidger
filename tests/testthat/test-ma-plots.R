context("MA plots")

test_that("MA plots work", {
    # Load data
    data(df.cuffdiff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsMAPlot(x = 'iPS', y = 'hESC', data = df.cuff, d.factor = NULL, 
             type = 'cuffdiff', padj = 0.05, y.lim = NULL, lfc = NULL, 
             title = TRUE, legend = TRUE, grid = TRUE)

    # DESeq2
    vsMAPlot(x = 'treated', y = 'untreated', data = df.deseq, 
             d.factor = 'condition', type = 'deseq', padj = 0.05, y.lim = NULL,
             lfc = NULL, title = TRUE, legend = TRUE, grid = TRUE)

    # edgeR
    vsMAPlot(x = 'treated', y = 'untreated', data = df.deseq, 
             d.factor = 'condition', type = 'deseq', padj = 0.05, y.lim = NULL,
             lfc = NULL, title = TRUE, legend = TRUE, grid = TRUE)
})