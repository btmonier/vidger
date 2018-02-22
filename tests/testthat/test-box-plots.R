context("box plots")

test_that("box plots work", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsBoxPlot(data = df.cuff, d.factor = NULL, type = 'cuffdiff', title = TRUE,
              legend = TRUE, grid = TRUE)

    # DESeq2
    vsBoxPlot(data = df.deseq, d.factor = 'condition', type = 'deseq', 
              title = TRUE, legend = TRUE, grid = TRUE)

    # edgeR
    vsBoxPlot(data = df.edger, d.factor = NULL, type = 'edger', title = TRUE,
              legend = TRUE, grid = TRUE)
})