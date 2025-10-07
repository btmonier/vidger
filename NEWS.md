# ViDGER 1.4.5
* Fixed data return examples


# ViDGER 1.1.3

* Added `grey.scale` parameter to `vsDEGMatrix()` function.

* Added text size parameters for *all* functions. These parameters are as
  follows:

    + `xaxis.text.size`
    + `yaxis.text.size`
    + `xaxis.title.size`
    + `yaxis.title.size`
    + `main.title.size`
    + `legend.text.size`
    + `legend.title.size`
    + `facet.title.size`

  Each of these parameters controls font size in terms of typographic points.
  


# ViDGER 1.1.2

* Added `data.return` parameters to *all* functions. 

* Added new aesthetic design parameter (`aes`) for the `vsBoxPlot()` function:

    + `box`: standard box plot
    + `violin`: violin plot
    + `boxdot`: box plot with dot plot overlay
    + `viodot`: violin plot with dot plot overlay
    + `viosumm`: violin plot with summary stats overlay
    + `notch`: box plot with notch

* Added new fill color parameter (`fill.color`) for the `vsBoxplot()` function:

    + Based on palettes found in the `RColorBrewer` 
      [package.](https://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf)
    + A visual list of all the palettes can be found
      [here.](http://www.r-graph-gallery.com/38-rcolorbrewers-palettes/)



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
