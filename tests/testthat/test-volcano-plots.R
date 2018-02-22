context("volcano plots")

test_that("volcano plots work", {
    # Load data
    data(df.cuff)
    data(df.deseq)
    data(df.edger)

    # cuffdiff
    vsVolcano(x = 'iPS', y = 'hESC', data = df.cuff, d.factor = NULL, 
              type = 'cuffdiff', padj = 0.05, x.lim = NULL, lfc = NULL, 
              title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE)

    # DESeq2
    vsVolcano(x = 'treated_paired.end', y = 'untreated_paired.end', 
              data = df.deseq, d.factor = 'condition', type = 'deseq', 
              padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
              legend = TRUE, grid = TRUE, data.return = FALSE)

    # edgeR
    vsVolcano(x = 'WW', y = 'MM', data = df.edger, d.factor = NULL, 
              type = 'edger', padj = 0.05, x.lim = NULL, lfc = NULL, 
              title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE)
})