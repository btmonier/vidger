
ViDGER <img src="man/figures/logo-02.png" align="right" />
==========================================================

[![Build Status](https://travis-ci.org/btmonier/vidger.svg?branch=master)](https://travis-ci.org/btmonier/vidger)

Overview
--------

`vidger` (**Vi**sualization of **D**ifferential **G**ene **E**xpression using **R**), is an `R` package that can rapidly generate information-rich visualizations for the interpretation of differential gene expression results from three widely-used tools: `Cuffdiff`, `DESeq2`, and `edgeR`.

Installation
------------

The stable version of this package is available on [Bioconductor](http://bioconductor.org/). You can install it by:

``` r
source("http://bioconductor.org/biocLite.R")
biocLite("vidger")
```

If you want the latest version, install it directly from this GitHub repo:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("btmonier/vidger", ref = "devel")
```

Functions
---------

The stable release of `vidger` has 9 visualization functions:

-   `vsScatterPlot()`
-   `vsScatterMatrix()`
-   `vsBoxplot()`
-   `vsDEGMatrix()`
-   `vsVolcano()`
-   `vsVolcanoMatrix()`
-   `vsMAPlot()`
-   `vsMAMatrix()`
-   `vsFourWay()`

Loading test data
-----------------

To simulate the usage of the three aformentioned tools, "toy" data sets have been implemented in this package. Each of these data sets represents their respective `R` class:

-   `df.cuff` A `cuffdiff` output file.
-   `df.deseq` A `DESeq2` object class.
-   `df.edger` An `edgeR` object class.

To load these data sets, use the following command:

``` r
data("<object-type>")
```

...where `"<object-type>"` is one of the previously mentioned data sets.

Getting help
------------

For additional information on these functions, please see the given documentation in the `vidger` package by adding the `?` help operator before any of the given functions in this package or by using the `help()` function.

For a more in-depth analysis, consider reading the vignette provided with this package:

``` r
vignette("vidger")
```
