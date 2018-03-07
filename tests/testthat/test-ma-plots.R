context("vsMAPlot")

test_that("vsMAPlot() give proper errors", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff with wrong `type` parameter
    expect_error(
        vsMAPlot(
            x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL, 
            type = "edger", padj = 0.05, y.lim = NULL, lfc = 1, 
            title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE
        )    
    )

    # DESeq2 with no `d.factor` parameter
    expect_error(
        vsMAPlot(
            x = "treated_paired.end", y = "untreated_paired.end", 
            data = df.deseq, d.factor = NULL, type = "deseq", 
            padj = 0.05, y.lim = NULL, lfc = NULL, title = TRUE, 
            legend = TRUE, grid = TRUE, data.return = FALSE
        )
    )

    # edgeR with missing `type` parameter
    expect_error(
        vsMAPlot(
            x = "WM", y = "MM", data = df.edger, d.factor = NULL, 
            type = , padj = 0.1, y.lim = NULL, lfc = 1, 
            title = FALSE, legend = TRUE, grid = TRUE, data.return = FALSE
        )        
    )

    # cuffdiff with completely wrong `type` parameter
    expect_error(
        vsMAPlot(
            x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL, 
            type = "cuffdif", padj = 0.05, y.lim = NULL, lfc = 1, 
            title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE
        )    
    )

    # DESeq2 with wrong `x` parameter
    expect_error(
        vsMAPlot(
            x = "treat.end", y = "untreated_paired.end", 
            data = df.deseq, d.factor = "condition", type = "deseq", 
            padj = 0.05, y.lim = NULL, lfc = NULL, title = TRUE, 
            legend = TRUE, grid = TRUE, data.return = FALSE
        )
    )

    # edgeR with wrong `y` parameter
    expect_error(
        vsMAPlot(
            x = "WM", y = "condition", data = df.edger, d.factor = NULL, 
            type = "edger", padj = 0.1, y.lim = NULL, lfc = 1, 
            title = FALSE, legend = TRUE, grid = TRUE, data.return = FALSE
        )
    )

})