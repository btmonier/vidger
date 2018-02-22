context("scatter plot matrices")

test_that("scatter plot matrices work", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsScatterMatrix(data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
                    comp = NULL, title = TRUE, grid = TRUE, man.title = NULL)

    # DESeq2
    vsScatterMatrix(data = df.deseq, d.factor = 'condition', type = 'deseq',
                    comp = NULL, title = TRUE, grid = TRUE, man.title = NULL)

    # edgeR
    vsScatterMatrix(data = df.edger, d.factor = NULL, type = 'edger', 
                    comp = NULL, title = TRUE, grid = TRUE, man.title = NULL)
})