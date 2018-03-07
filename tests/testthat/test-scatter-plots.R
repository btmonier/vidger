context("vsScatterPlot")

test_that("vsScatterPlot() give proper errors", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff with wrong `type` parameter
    expect_error(
        vsScatterPlot(
            x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL, 
            type = "deseq", title = TRUE, grid = TRUE
        )
    )

    # DESeq2 with no `d.factor` parameter
    expect_error(
        vsScatterPlot(
            x = "treated_paired.end", y = "untreated_paired.end", 
            data = df.deseq, d.factor = NULL, type = "deseq", 
            title = TRUE, grid = TRUE
        )
    )

    # edgeR with missing `type` parameter
    expect_error(
        vsScatterPlot(
            x = "WW", y = "WM", data = df.edger, d.factor = NULL, 
            type = , title = TRUE, grid = TRUE
        )
    )

    # cuffdiff with completely wrong `type` parameter
    expect_error(
        vsScatterPlot(
            x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL, 
            type = "cufdif", title = TRUE, grid = TRUE
        )
    )

})