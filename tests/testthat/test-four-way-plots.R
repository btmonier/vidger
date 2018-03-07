context("vsFourWay")

test_that("vsFourWay() gives proper errores", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff with wrong `type` parameter
    expect_error(
        vsFourWay(
            x = "hESC", y = "iPS", control = "Fibroblasts", data = df.cuff, 
            d.factor = NULL, type = "edger", padj = 0.05, x.lim = NULL, 
            y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
            data.return = FALSE
        )
    )


    # DESeq2 with no `d.factor` parameter
    expect_error(
        vsFourWay(
            x = "treated_paired.end", y = "untreated_paired.end", 
            control = "untreated_single.read", data = df.deseq, 
            d.factor = NULL, type = "deseq", padj = 0.05, 
            x.lim = NULL, y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
            data.return = FALSE
        )
    )

    # edgeR with missing `type` parameter
    expect_error(
        vsFourWay(
            x = "WM", y = "WW", control = "MM", data = df.edger, 
            d.factor = NULL, type = "", padj = 0.05, x.lim = NULL, 
            y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
            data.return = FALSE
        )
    )

    # cuffdiff with completely wrong `type` parameter
    expect_error(
        vsFourWay(
            x = "hESC", y = "iPS", control = "Fibroblasts", data = df.cuff, 
            d.factor = NULL, type = "cufdif", padj = 0.05, x.lim = NULL, 
            y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
            data.return = FALSE
        )
    )

    # DESeq2 with wrong `x` parameter
    expect_error(
        vsFourWay(
            x = "treated", y = "untreated_paired.end", 
            control = "untreated_single.read", data = df.deseq, 
            d.factor = "condition", type = "deseq", padj = 0.05, 
            x.lim = NULL, y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
            data.return = FALSE
        )
    )

    # edgeR with wrong `y` parameter
    expect_error(
        vsFourWay(
            x = "WM", y = "ww", control = "MM", data = df.edger, 
            d.factor = NULL, type = "edger", padj = 0.05, x.lim = NULL, 
            y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
            data.return = FALSE
        )
    )

    # cuffdiff with wrong `control` parameter
    expect_error(
        vsFourWay(
            x = "hESC", y = "iPS", control = "fibro", data = df.cuff, 
            d.factor = NULL, type = "cuffdiff", padj = 0.05, x.lim = NULL, 
            y.lim = NULL, lfc = 2, title = TRUE, grid = TRUE, 
            data.return = FALSE
        )
    )
})