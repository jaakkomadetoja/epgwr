# epgwr
Error Propagation in Geographically Weighted Regression

This R package applies error propagation in geographically weighted regression using Monte Carlo simulation.

To run this package, you need:
- The statistical software R
- (Recommended) RStudio
- The package has only been tested using Windows OS; it might not work in Linux or Mac

In R (or RStudio), run the following commands to install the package:
install.packages("devtools")
library("devtools")
devtools::install_github("jaakkomadetoja/epgwr")
library(epgwr)

For examples and how to use the package, run the following command in R
?epgwr_mc
