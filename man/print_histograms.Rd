% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{print_histograms}
\alias{print_histograms}
\title{Print histograms}
\usage{
print_histograms(data, print_file = TRUE, breaks = 20)
}
\arguments{
\item{data}{output from the function \code{\link{epgwr_mc}}}

\item{print_file}{TRUE (default) for creating .png images and saving them to
disk, FALSE for just showing histograms in R plot device one at a time}

\item{breaks}{the number of breaks in the histograms}
}
\value{

}
\description{
This function prints histograms for single-value metrics
  resulting from error propagation in geographically weighted regression
}
\examples{

}
\references{
Madetoja, J. (2018). Error propagation in geographically weighted
  regression. (Doctoral dissertation, Aalto University). Manuscript in
  preparation.
}
\seealso{
\code{\link{epgwr_mc}}, \code{\link{print_maps}}
}
\author{
Jaakko Madetoja
}
