context("vsVolcano")

test_that("vsVolcano() give proper errors", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff with wrong `type` parameter
    expect_error(
        vsVolcano(
            x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL, 
            type = "edger", padj = 0.05, x.lim = NULL, lfc = 2, 
            title = TRUE, grid = TRUE, data.return = FALSE
        )
    )

    # DESeq2 with no `d.factor` parameter
    expect_error(
        vsVolcano(
            x = "treated_paired.end", y = "untreated_paired.end", 
            data = df.deseq, d.factor = NULL, 
            type = "deseq", padj = 0.05, x.lim = NULL, lfc = NULL, 
            title = TRUE, grid = TRUE, data.return = FALSE
        )
    )

    # edgeR with missing `type` parameter
    expect_error(
        vsVolcano(
            x = "WM", y = "MM", data = df.edger, d.factor = NULL, 
            type = , padj = 0.1, x.lim = NULL, lfc = 2, 
            title = FALSE, grid = TRUE, data.return = FALSE
        )
    )

    # cuffdiff with completely wrong `type` parameter
    expect_error(
        vsVolcano(
            x = "hESC", y = "iPS", data = df.cuff, d.factor = NULL, 
            type = "edg", padj = 0.05, x.lim = NULL, lfc = 2, 
            title = TRUE, grid = TRUE, data.return = FALSE
        )
    )

    # DESeq2 with wrong `x` parameter
    expect_error(
        vsVolcano(
            x = "treated_p.end", y = "untreated_paired.end", 
            data = df.deseq, d.factor = "condition", 
            type = "deseq", padj = 0.05, x.lim = NULL, lfc = NULL, 
            title = TRUE, grid = TRUE, data.return = FALSE
        )
    )

    # edgeR with wrong `y` parameter
    expect_error(
        vsVolcano(
            x = "WM", y = "MMM", data = df.edger, d.factor = NULL, 
            type = "edger", padj = 0.1, x.lim = NULL, lfc = 2, 
            title = FALSE, grid = TRUE, data.return = FALSE
        )
    )
})