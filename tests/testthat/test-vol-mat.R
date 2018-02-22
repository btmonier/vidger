context("volcano plot matrices")

test_that("volcano plot matrices work", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsVolcanoMatrix(data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
                    padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
                    legend = TRUE, grid = TRUE, counts = TRUE)

    # DESeq2
    vsVolcanoMatrix(data = df.deseq, d.factor = 'condition', type = 'deseq', 
                    padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
                    legend = TRUE, grid = TRUE, counts = TRUE)

    # edgeR
    vsVolcanoMatrix(data = df.edger, d.factor = NULL, type = 'edger', 
                    padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
                    legend = TRUE, grid = TRUE, counts = TRUE)
})