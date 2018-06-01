# ViDGER 1.1.0

* Added labelling parameters for user to highlight/label IDs of interest.
  This was added to the following functions:

    + `vsScatterPlot()`
    + `vsMAPlot()`
    + `vsVolcano()`
    + `vsFourWay()`

  This feature can be used by adding a vector of IDs found in their respective
  data class using the `highlight` parameter.

* Added `data.return` parameter to the `vsScatterPlot()` function.



# ViDGER 1.0.0

* Added a `NEWS.md` file to track changes to the package.

* Added new logo to commemorate its acceptance into Bioconductor.

* Deployment of 9 functions:
    + `vsScatterPlot()`
    + `vsScatterMatrix()`
    + `vsBoxplot()`
    + `vsDEGMatrix()`
    + `vsVolcano()`
    + `vsVolcanoMatrix()`
    + `vsMAPlot()`
    + `vsMAMatrix()`
    + `vsFourWay()`