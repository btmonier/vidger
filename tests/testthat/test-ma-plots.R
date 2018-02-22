context("MA plots")

test_that("MA plots work", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsMAPlot(x = 'iPS', y = 'hESC', data = df.cuff, d.factor = NULL, 
             type = 'cuffdiff', padj = 0.05, y.lim = NULL, lfc = NULL, 
             title = TRUE, legend = TRUE, grid = TRUE)

    # DESeq2
    vsMAPlot(x = 'treated_paired.end', y = 'untreated_paired.end', 
             data = df.deseq, d.factor = 'condition', type = 'deseq', 
             padj = 0.05, y.lim = NULL, lfc = NULL, title = TRUE, 
             legend = TRUE, grid = TRUE)

    # edgeR
    vsMAPlot(x = 'WM', y = 'MM', data = df.edger, d.factor = NULL, 
             type = 'edger', padj = 0.1, y.lim = NULL, lfc = 1, 
             title = FALSE, legend = TRUE, grid = TRUE, data.return = FALSE)
})