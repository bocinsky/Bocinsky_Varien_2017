### Calibrating Maize Paleoproduction Models using Experimental Data

[![DOI](https://zenodo.org/badge/67527309.svg)](https://zenodo.org/badge/latestdoi/67527309)

This is the script for the analyses in Bocinsky and Varien 2017:

R. Kyle Bocinsky and Mark D. Varien. 2017. Calibrating Maize Paleoproduction Models using Experimental Data. *Ethnobiology* 37(2).

To run, clone this repository, then open `Bocinsky_Varien_2017.Rproj` in RStudio and then the `Bocinsky_Varien_2017.R` script. Run it all (⌘-A then ⌘-return).

The analysis consists of three components. First, VEP-style paleoproduction estimates are generated for the locations of the Pueblo Farming Project gardens. Then, PFP data are analysed to estimate standard productivity across gardens. Finally, these data are brought together in a series of figures to be used in the manuscript.

This script makes heavy use of piping from the [magrittr](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html) package and the data manipulation functions in the [dplyr](https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html) package. External soils data from the USGS and weather data from NOAA are downloaded using the [FedData](https://github.com/bocinsky/FedData) package.

For figure output, [Ghostscript](http://www.ghostscript.com/) and [ImageMagick](http://www.imagemagick.org/script/index.php) command-line tools are also required.

Throughout the VEP productivity reconstruction, references are made to the canonical VEP reconstruction reference:

Kohler, Timothy A. 2012. Modeling agricultural productivity and farming effort. In Kohler, T. A. and Varien, M. D., editors, *Emergence and Collapse of Early Villages: Models of Central Mesa Verde Archaeology*, chapter 6, pages 85–112. University of California Press, Berkeley, California.

This code has been developed and tested in Mac OS X 10.12.