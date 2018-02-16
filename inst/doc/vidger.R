## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.pos = 'h')

## ---- echo=FALSE, message=FALSE------------------------------------------
require(vidger)
require(DESeq2)
require(edgeR)
data("df.cuff")
data("df.deseq")
data("df.edger")

## ---- message=FALSE------------------------------------------------------
# Install from GitHub
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("btmonier/vidger")

# Load vidger
require(vidger)

## ---- eval=FALSE---------------------------------------------------------
#  data(<data_set>)

## ----  message=FALSE, fig.align='center', fig.cap='A boxplot example using the `vsBoxPlot()` function with `cuffdiff} data. In this example, FPKM distributions for each treatment within an experiment are shown in the form of a "box and whisker plot". '----
vsBoxPlot(data = df.cuff, d.factor = NULL, type = 'cuffdiff', title = TRUE,
          legend = TRUE, grid = TRUE)

## ----  message=FALSE, fig.align='center', fig.cap='A boxplot example using the `vsBoxPlot()` function with `DESeq2} data. In this example, FPKM distributions for each treatment within an experiment are shown in the form of a "box and whisker plot". '----
vsBoxPlot(data = df.deseq, d.factor = 'condition', type = 'deseq', 
          title = TRUE, legend = TRUE, grid = TRUE)

## ----  message=FALSE, fig.align='center', fig.cap='A boxplot example using the `vsBoxPlot()` function with `edgeR} data. In this example, CPM distributions for each treatment within an experiment are shown in the form of a "box and whisker plot". '----
vsBoxPlot(data = df.edger, d.factor = NULL, type = 'edger', title = TRUE,
          legend = TRUE, grid = TRUE)

## ---- message=FALSE, fig.align='center', fig.cap='A scatterplot example using the `vsScatterPlot()` function with `Cuffdiff} data. In this visualization, $log_{10}$ comparisons are made of fragments per kilobase of transcript per million mapped reads (FPKM) measurments. The dashed line represents regression line for the comparison.'----
vsScatterPlot(x = 'hESC', y = 'iPS', data = df.cuff, type = 'cuffdiff',
              d.factor = NULL, title = TRUE, grid = TRUE)

## ---- message=FALSE, fig.align='center', fig.cap='A scatterplot example using the `vsScatterPlot()` function with `DESeq2} data. In this visualization, $log_{10}$ comparisons are made of fragments per kilobase of transcript per million mapped reads (FPKM) measurments. The dashed line represents regression line for the comparison.'----
vsScatterPlot(x = 'treated', y = 'untreated', data = df.deseq, type = 'deseq',
              d.factor = 'condition', title = TRUE, grid = TRUE)

## ---- message=FALSE, fig.align='center', fig.cap='A scatterplot example using the `vsScatterPlot()` function with `edgeR} data. In this visualization, $log_{10}$ comparisons are made of counts per million mapped reads (CPM) measurments. The dashed line represents regression line for the comparison.'----
vsScatterPlot(x = 'WM', y = 'MM', data = df.edger, type = 'edger',
              d.factor = NULL, title = TRUE, grid = TRUE)

## ---- message=FALSE, fig.align='center', fig.cap='A scatterplot matrix example using the `vsScatterMatrix()` function with `Cuffdiff} data. Similar to the scatterplot function, this visualization allows for all comparisons to be made within an experiment. In addition to the scatterplot visuals, FPKM distributions (histograms) and correlation (Corr) values are generated.'----
vsScatterMatrix(data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
                comp = NULL, title = TRUE, grid = TRUE, man.title = NULL)

## ---- message=FALSE, fig.align='center', fig.cap='A scatterplot matrix example using the `vsScatterMatrix()` function with `DESeq2} data. Similar to the scatterplot function, this visualization allows for all comparisons to be made within an experiment. In addition to the scatterplot visuals, FPKM distributions (histograms) and correlation (Corr) values are generated.'----
vsScatterMatrix(data = df.deseq, d.factor = 'condition', type = 'deseq',
                comp = NULL, title = TRUE, grid = TRUE, man.title = NULL)

## ---- message=FALSE, fig.align='center', fig.cap='A scatterplot matrix example using the `vsScatterMatrix()` function with `edgeR} data. Similar to the scatterplot function, this visualization allows for all comparisons to be made within an experiment. In addition to the scatterplot visuals, CPM distributions (histograms) and correlation (Corr) values are generated.'----
vsScatterMatrix(data = df.edger, d.factor = NULL, type = 'edger', comp = NULL,
                title = TRUE, grid = TRUE, man.title = NULL)

## ----  message=FALSE, fig.align='center', fig.cap='A matrix of differentially expressed genes (DEGs) at a given *p*-value using the `vsDEGMatrix()` function with `Cuffdiff} data. With this function, the user is able to visualize the number of DEGs ata given adjusted *p*-value for each experimental treatment level. Higher color intensity correlates to a higher number of DEGs.'----
vsDEGMatrix(data = df.cuff, padj = 0.05, d.factor = NULL, type = 'cuffdiff', 
            title = TRUE, legend = TRUE, grid = TRUE)

## ----  message=FALSE, fig.align='center', fig.cap='A matrix of differentially expressed genes (DEGs) at a given *p*-value using the `vsDEGMatrix()` function with `DESeq2} data. With this function, the user is able to visualize the number of DEGs ata given adjusted *p*-value for each experimental treatment level. Higher color intensity correlates to a higher number of DEGs.'----
vsDEGMatrix(data = df.deseq, padj = 0.05, d.factor = 'condition', 
            type = 'deseq', title = TRUE, legend = TRUE, grid = TRUE)

## ----  message=FALSE, fig.align='center', fig.cap='A matrix of differentially expressed genes (DEGs) at a given *p*-value using the `vsDEGMatrix()` function with `edgeR} data. With this function, the user is able to visualize the number of DEGs ata given adjusted *p*-value for each experimental treatment level. Higher color intensity correlates to a higher number of DEGs.'----
vsDEGMatrix(data = df.edger, padj = 0.05, d.factor = NULL, type = 'edger', 
            title = TRUE, legend = TRUE, grid = TRUE)

## ----  message=FALSE, fig.align='center', fig.cap='MA plot visualization using the `vsMAPLot()` function with `Cuffdiff} data. LFCs are plotted mean counts to determine the variance between two treatments in terms of gene expression. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. Numerical values in parantheses for each legend color indicate the number of transcripts that meet the prior conditions. Triangular shapes represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Dashed lines indicate user-defined LFC values.'----
vsMAPlot(x = 'iPS', y = 'hESC', data = df.cuff, d.factor = NULL, 
         type = 'cuffdiff', padj = 0.05, y.lim = NULL, lfc = NULL, 
         title = TRUE, legend = TRUE, grid = TRUE)

## ----  message=FALSE, fig.align='center', fig.cap='MA plot visualization using the `vsMAPLot()` function with `DESeq2` data. LFCs are plotted mean counts to determine the variance between two treatments in terms of gene expression. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. Numerical values in parantheses for each legend color indicate the number of transcripts that meet the prior conditions. Triangular shapes represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Dashed lines indicate user-defined LFC values.'----
vsMAPlot(x = 'treated', y = 'untreated', data = df.deseq, 
         d.factor = 'condition', type = 'deseq', padj = 0.05, y.lim = NULL, 
         lfc = NULL, title = TRUE, legend = TRUE, grid = TRUE)

## ----  message=FALSE, fig.align='center', fig.cap='MA plot visualization using the `vsMAPLot()` function with `edgeR` data. LFCs are plotted mean counts to determine the variance between two treatments in terms of gene expression. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. Numerical values in parantheses for each legend color indicate the number of transcripts that meet the prior conditions. Triangular shapes represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Dashed lines indicate user-defined LFC values.'----
vsMAPlot(x = 'WW', y = 'MM', data = df.edger, d.factor = NULL, 
         type = 'edger', padj = 0.05, y.lim = NULL, lfc = NULL, 
         title = TRUE, legend = TRUE, grid = TRUE)

## ----  message=FALSE, fig.align='center', fig.cap='A volcano plot example using the `vsVolcano()` function with `Cuffdiff} data. In this visualization, comparisons are made between the $-log_{10}$ *p*-value versus the $log_2$ fold change (LFC) between two treatments. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. Numerical values in parantheses for each legend color indicate the number of transcripts that meet the prior conditions. Left and right brackets (< and >) represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Vertical and horizontal lines indicate user-defined LFC and adjusted *p*-values, respectively.'----
vsVolcano(x = 'iPS', y = 'hESC', data = df.cuff, d.factor = NULL, 
          type = 'cuffdiff', padj = 0.05, x.lim = NULL, lfc = NULL, 
          title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE)

## ----  message=FALSE, fig.align='center', fig.cap='A volcano plot example using the `vsVolcano()` function with `DESeq2} data. In this visualization, comparisons are made between the $-log_{10}$ *p*-value versus the $log_2$ fold change (LFC) between two treatments. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. Numerical values in parantheses for each legend color indicate the number of transcripts that meet the prior conditions. Left and right brackets (< and >) represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Vertical and horizontal lines indicate user-defined LFC and adjusted *p*-values, respectively.'----
vsVolcano(x = 'treated', y = 'untreated', data = df.deseq, 
          d.factor = 'condition', type = 'deseq', padj = 0.05, 
          x.lim = NULL, lfc = NULL, title = TRUE, legend = TRUE, 
          grid = TRUE, data.return = FALSE)

## ----  message=FALSE, fig.align='center', fig.cap='A volcano plot example using the `vsVolcano()` function with `edgeR} data. In this visualization, comparisons are made between the $-log_{10}$ *p*-value versus the $log_2$ fold change (LFC) between two treatments. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. Numerical values in parantheses for each legend color indicate the number of transcripts that meet the prior conditions. Left and right brackets (< and >) represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Vertical and horizontal lines indicate user-defined LFC and adjusted *p*-values, respectively.'----
vsVolcano(x = 'WW', y = 'MM', data = df.edger, d.factor = NULL, 
          type = 'edger', padj = 0.05, x.lim = NULL, lfc = NULL, 
          title = TRUE, legend = TRUE, grid = TRUE, data.return = FALSE)

## ----  message=FALSE, fig.width=8, fig.height=5, fig.align='center', fig.cap='A volcano plot matrix using the `vsVolcanoMatrix()` function with `Cuffdiff} data. Similar to the `vsVolcano()` function, `vsVolcanoMatrix()` will generate a matrix of volcano plots for all comparisons within an experiment. Comparisons are made between the $-log_{10}$ *p*-value versus the $log_2$ fold change (LFC) between two treatments. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. The blue and green numbers in each facet represent the number of transcripts that meet the criteria for blue and green nodes in each comparison. Left and right brackets (< and >) represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Vertical and horizontal lines indicate user-defined LFC and adjusted *p*-values, respectively.'----
vsVolcanoMatrix(data = df.cuff, d.factor = NULL, type = 'cuffdiff', 
                padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
                legend = TRUE, grid = TRUE, counts = TRUE)

## ----  message=FALSE, fig.width=8, fig.height=5, fig.align='center', fig.cap='A volcano plot matrix using the `vsVolcanoMatrix()` function with `DESeq2} data. Similar to the `vsVolcano()` function, `vsVolcanoMatrix()` will generate a matrix of volcano plots for all comparisons within an experiment. Comparisons are made between the $-log_{10}$ *p*-value versus the $log_2$ fold change (LFC) between two treatments. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. The blue and green numbers in each facet represent the number of transcripts that meet the criteria for blue and green nodes in each comparison. Left and right brackets (< and >) represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Vertical and horizontal lines indicate user-defined LFC and adjusted *p*-values, respectively.'----
vsVolcanoMatrix(data = df.deseq, d.factor = 'condition', type = 'deseq', 
                padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
                legend = TRUE, grid = TRUE, counts = TRUE)

## ----  message=FALSE, fig.width=8, fig.height=5, fig.align='center', fig.cap='A volcano plot matrix using the `vsVolcanoMatrix()` function with `edgeR} data. Similar to the `vsVolcano()` function, `vsVolcanoMatrix()` will generate a matrix of volcano plots for all comparisons within an experiment. Comparisons are made between the $-log_{10}$ *p*-value versus the $log_2$ fold change (LFC) between two treatments. Blue nodes on the graph represent statistically significant LFCs which are greater than a given value than a user-defined LFC parameter. Green nodes indicate statistically significant LFCs which are less than the user-defined LFC parameter. Gray nodes are data points that are not statistically significant. The blue and green numbers in each facet represent the number of transcripts that meet the criteria for blue and green nodes in each comparison. Left and right brackets (< and >) represent values that exceed the viewing area of the graph. Node size changes represent the magnitude of the LFC values (i.e. larger shapes reflect larger LFC values). Vertical and horizontal lines indicate user-defined LFC and adjusted *p*-values, respectively.'----
vsVolcanoMatrix(data = df.edger, d.factor = NULL, type = 'edger', padj = 0.05, 
                x.lim = NULL, lfc = NULL, title = TRUE, legend = TRUE, 
                grid = TRUE, counts = TRUE)

## ----  message=FALSE, fig.width=8, fig.height=5, fig.align='center', fig.cap='A four way plot visualization using the `vsFourWay()` function with `Cuffdiff} data. In this example, LFCs comparisons between two treatments and a control are made. Blue nodes indicate statistically significant LFCs which are greater than a given user-defined value for both x and y-axes. Green nodes reflect statistically significant LFCs which are less than a user-defined value for treatment y and greater than said value for treatment x. Similar to green nodes, red nodes reflect statistically significant LFCs which are greater than a user-defined vlaue treatment y and less than said value for treatment x. Gray nodes are data points that are not statistically significant for both x and y-axes. Triangular shapes indicate values which exceed the viewing are for the graph. Size change reflects the magnitude of LFC values (i.e. larger shapes reflect larger LFC values). Vertical and horizontal dashed lines indicate user-defined LFC values.'----
vsFourWay(x = 'iPS', y = 'hESC', control = 'Fibroblasts', data = df.cuff,
          d.factor = NULL, type = 'cuffdiff', padj = 0.05, x.lim = NULL,
          y.lim = NULL, lfc = NULL, legend = TRUE, title = TRUE, grid = TRUE)

## ----  message=FALSE, fig.width=8, fig.height=5, fig.align='center', fig.cap='A four way plot visualization using the `vsFourWay()` function with `edgeR} data. In this example, LFCs comparisons between two treatments and a control are made. Blue nodes indicate statistically significant LFCs which are greater than a given user-defined value for both x and y-axes. Green nodes reflect statistically significant LFCs which are less than a user-defined value for treatment y and greater than said value for treatment x. Similar to green nodes, red nodes reflect statistically significant LFCs which are greater than a user-defined vlaue treatment y and less than said value for treatment x. Gray nodes are data points that are not statistically significant for both x and y-axes.'----
vsFourWay(x = 'WW', y = 'WM', control = 'MM', data = df.edger,
          d.factor = NULL, type = 'edger', padj = 0.05, x.lim = NULL,
          y.lim = NULL, lfc = NULL, legend = TRUE, title = TRUE, grid = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='An illustration detailing the principles behind the node size for the differntial gene expression functions. In this figure, the data points increase in size depending on which quartile they reside as the absolute LFC increases (top bar). Data points that fall within the viewing area classified as "SUB" while data points that exceed this area are classified as "T-1" through "T-4". '----
knitr::include_graphics("lfc-shape.png", auto_pdf = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='Benchmarks for the `vsScatterPlot()` function. Time (ms) distributions were generated for this function using 100 trials for each of the three RNAseq data objects. Cuffdiff, DESeq2, and edgeR example data sets contained 1200, 724, and 29391 transcripts, respectively. '----
knitr::include_graphics("eff-scatter.png", auto_pdf = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='Benchmarks for the `vsScatterMatrix()` function. Time (ms) distributions were generated for this function using 100 trials for each of the three RNAseq data objects. Cuffdiff, DESeq2, and edgeR example data sets contained 1200, 724, and 29391 transcripts, respectively.'----
knitr::include_graphics("eff-smatrix.png", auto_pdf = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='Benchmarks for the `vsBoxPlot()` function. Time (ms) distributions were generated for this function using 100 trials for each of the three RNAseq data objects. Cuffdiff, DESeq2, and edgeR example data sets contained 1200, 724, and 29391 transcripts, respectively. '----
knitr::include_graphics("eff-box.png", auto_pdf = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='Benchmarks for the `vsDEGMatrix()` function. Time (ms) distributions were generated for this function using 100 trials for each of the three RNAseq data objects. Cuffdiff, DESeq2, and edgeR example data sets contained 1200, 724, and 29391 transcripts, respectively.'----
knitr::include_graphics("eff-deg.png", auto_pdf = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='Benchmarks for the `vsVolcano()` function. Time (ms) distributions were generated for this function using 100 trials for each of the three RNAseq data objects. Cuffdiff, DESeq2, and edgeR example data sets contained 1200, 724, and 29391 transcripts, respectively.'----
knitr::include_graphics("eff-volcano.png", auto_pdf = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='Benchmarks for the `vsVolcanoMatrix()` function. Time (s) distributions were generated for this function using 100 trials for each of the three RNAseq data objects. Cuffdiff, DESeq2, and edgeR example data sets contained 1200, 724, and 29391 transcripts, respectively.'----
knitr::include_graphics("eff-vmatrix.png", auto_pdf = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='Benchmarks for the `vsMAPlot()` function. Time (ms) distributions were generated for this function using 100 trials for each of the three RNAseq data objects. Cuffdiff, DESeq2, and edgeR example data sets contained 1200, 724, and 29391 transcripts, respectively. '----
knitr::include_graphics("eff-maplot.png", auto_pdf = TRUE)

## ---- echo=FALSE, fig.align='center', fig.cap='Benchmarks for the `vsFourWay()` function. Time (ms) distributions were generated for this function using 100 trials for each of the three RNAseq data objects. Cuffdiff, DESeq2, and edgeR example data sets contained 1200, 724, and 29391 transcripts, respectively. '----
knitr::include_graphics("eff-four.png", auto_pdf = TRUE)

