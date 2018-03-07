context("vsBoxPlot")

test_that("vsBoxPlot() gives proper errors", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff with wrong `type` parameter
    expect_error(
      vsBoxPlot(
        data = df.cuff, d.factor = NULL, type = "deseq", title = TRUE,
        legend = TRUE, grid = TRUE
      )
    )

    # DESeq2 with no `d.factor` parameter
    expect_error(
      vsBoxPlot(
        data = df.deseq, d.factor = NULL, type = "deseq", title = TRUE,
        legend = TRUE, grid = TRUE
      )
    )

    # edgeR with missing `type` parameter
    expect_error(
      vsBoxPlot(
        data = df.edger, d.factor = NULL, type = , title = TRUE, 
        legend = TRUE, grid = TRUE
      )
    )

    # cuffdiff with completely wrong `type` parameter
    expect_error(
      vsBoxPlot(
        data = df.cuff, d.factor = NULL, type = "godzilla", title = TRUE, 
        legend = TRUE, grid = TRUE
      )
    )    
})