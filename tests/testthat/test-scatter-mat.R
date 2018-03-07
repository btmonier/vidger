context("vsScatterMatrix")

test_that("vsScatterMatrix() gives proper errors", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff with wrong `type` parameter
    expect_error(
        vsScatterMatrix(
            data = df.cuff, d.factor = NULL, type = "edger", 
            comp = NULL, title = TRUE, grid = TRUE, 
            man.title = "Example title"
        )
    )

    # DESeq2 with no `d.factor` parameter
    expect_error(
        vsScatterMatrix(
            data = df.deseq, d.factor = NULL, type = "deseq",
            comp = NULL, title = TRUE, grid = FALSE, man.title = NULL
        )
    )

    # edgeR with missing `type` parameter
    expect_error(
        vsScatterMatrix(
            data = df.edger, d.factor = NULL, type = , 
            comp = c("WM", "MM"), title = TRUE, grid = TRUE, 
            man.title = NULL
        )
    )

    # cuffdiff with completely wrong `type` parameter
    expect_error(
        vsScatterMatrix(
            data = df.cuff, d.factor = NULL, type = "edger", 
            comp = NULL, title = TRUE, grid = TRUE, 
            man.title = "Example title"
        )
    ) 
})