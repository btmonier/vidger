context("DEG matrices")

test_that("DEG matrices work", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsDEGMatrix(data = df.cuff, padj = 0.05, d.factor = NULL, 
                type = 'cuffdiff', title = TRUE, legend = TRUE, grid = TRUE)

    # DESeq2
    vsDEGMatrix(data = df.deseq, padj = 0.05, d.factor = 'condition', 
                type = 'deseq', title = TRUE, legend = TRUE, grid = TRUE)

    # edgeR
    vsDEGMatrix(data = df.edger, padj = 0.05, d.factor = NULL, type = 'edger', 
                title = TRUE, legend = TRUE, grid = TRUE)
})