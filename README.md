# ggviseq
Make rapid visualizations of RNAseq data in R

Installation
------------

The easiest way to obtain this package is to install `devtools` and pull the package contents directly from GitHub.

``` r
# Development version from GitHub
# install.packages('devtools')

devtools::install_github('btmonier/ggviseq')
```

New Functions
-------------
* `vsBoxplot()`
* `vsDEGmat()`
* `vsFourWay()`
* `vsMAPlot()`
* `vsScatterMatrix()`
* `vsScatterPlot()`
* `vsVolcano()`
* `vsVolcanoMatrix()`

Package Dependencies
--------------------
This package requires a few additional packages to work properly. Make sure the following packages are installed and loaded into your `R` directory:

* `ggplot2` for main visualizations
* `dplyr` for data munging
* `tidyr` for additional data munging
* `GGally` for matrix visualizations

Loading test data
-----------------
I have loaded three objects from three seperate analysis types:

* `df.cuff.RData` A `cuffdiff` output file.
* `df.deseq.RData` A `DESeq2` object class.
* `df.edger.RData` An `edgeR` object class.

To load these data sets, navigate to the working directory to where they are located (`setwd()`), and load by:

``` r
load("<object-type>.RData")
```
...where `<object-type>` is one of the previously mentioned data sets.

Getting help
------------
For additional information on these functions, please see the given documentation in the `ggviseq` package by adding the `?` help operator before any of the given functions in this package or by using the `help()` function. 
