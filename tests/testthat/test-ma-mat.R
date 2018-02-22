context("MA Matrices")

test_that("MA matrices work", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # Cuffdiff test
    vsMAMatrix(data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
                    padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE, 
                    grid = TRUE, counts = TRUE, data.return = FALSE)
    
    # DESeq2 test
    require(DESeq2)
    vsMAMatrix(data = df.deseq, d.factor = 'condition', type = 'deseq', 
                    padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE, 
                    grid = TRUE, counts = TRUE, data.return = FALSE)
    
    # edgeR test
    require(edgeR)
    vsMAMatrix(data = df.edger, d.factor = NULL, type = 'edger', 
                    padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE, 
                    grid = TRUE, counts = TRUE, data.return = FALSE)
})