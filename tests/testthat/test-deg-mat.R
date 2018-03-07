context("vsDEGMatrix")

test_that("vsDEGMatrix() gives proper errors", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff with wrong `type` parameter
    expect_error(
        vsDEGMatrix(
            df.cuff, padj = 0.05, d.factor = NULL, type = 'cuff', 
            title = TRUE, legend = TRUE, grid = TRUE
        )
    )

    # DESeq2 with no `d.factor` parameter
    expect_error(
        vsDEGMatrix(
            df.deseq, padj = 0.05, d.factor = NULL, type = 'deseq', 
            title = TRUE, legend = TRUE, grid = TRUE
        )
    )

    # edgeR with missing `type` parameter
    expect_error(
        vsDEGMatrix(
            df.edger, padj = 0.05, d.factor = NULL, type = "", 
            title = TRUE, legend = TRUE, grid = TRUE
        )
    )

    # cuffdiff with completely wrong `type` parameter
    expect_error(
        vsDEGMatrix(
            df.cuff, padj = 0.05, d.factor = NULL, type = "cufdif", 
            title = TRUE, legend = TRUE, grid = TRUE
        )
    ) 
})