context("vsMAMatrix")

test_that("vsMAMatrix() gives proper errors", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff with wrong `type` parameter
    expect_error(
        vsMAMatrix(
            data = df.cuff, d.factor = NULL, type = "deseq", 
            padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE, 
            grid = TRUE, counts = TRUE, data.return = FALSE
        )        
    )

    # DESeq2 with no `d.factor` parameter
    expect_error(
        vsMAMatrix(
            data = df.deseq, d.factor = NULL, type = "deseq", 
            padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE, 
            grid = TRUE, counts = TRUE, data.return = FALSE
        )
    )

    # edgeR with missing `type` parameter
    expect_error(
        vsMAMatrix(
            data = df.edger, d.factor = NULL, type = "", 
            padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE, 
            grid = TRUE, counts = TRUE, data.return = FALSE
        )        
    )

    # cuffdiff with completely wrong `type` parameter
    expect_error(
        vsMAMatrix(
            data = df.cuff, d.factor = NULL, type = "cufdif", 
            padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE, 
            grid = TRUE, counts = TRUE, data.return = FALSE
        )        
    )

})