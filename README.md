### Calibrating Maize Paleoproduction Models using Experimental Data

This is the script for the analyses in Varien and Bocinsky, under review:

Varien, Mark D. and R. Kyle Bocinsky. Under review. Calibrating Maize Paleoproduction Models using Experimental Data. *Ethnobiology*.

To run, clone this repository, then open `Ethnobiology.Rproj` in RStudio and then the `Ethnobiology.R` script. Run it all (⌘-A then ⌘-return).

The analysis consists of three components. First, VEP-style paleoproduction estimates are generated for the locations of the Pueblo Farming Project gardens. Then, PFP data are analysed to estimate standard productivity across gardens. Finally, these data are brought together in a series of figures to be used in the manuscript.

This script makes heavy use of piping from the [magrittr](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html) package and the data manipulation functions in the [dplyr](https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html) package.

Throughout the VEP productivity reconstruction, references are made to the canonical VEP reconstruction reference:

Kohler, Timothy A. 2012. Modeling agricultural productivity and farming effort. In Kohler, T. A. and Varien, M. D., editors, *Emergence and Collapse of Early Villages: Models of Central Mesa Verde Archaeology*, chapter 6, pages 85–112. University of California Press, Berkeley, California.

This code has been developed and tested in Mac OS X 10.10.