# Visualization tools for the results of epgwr_mc


#' @title Print histograms
#'
#' @description This function prints histograms for single-value metrics
#'   resulting from error propagation in geographically weighted regression
#'
#' @param data output from the function \code{\link{epgwr_mc}}
#' @param print_file TRUE (default) for creating .png images and saving them to
#'   disk, FALSE for just showing histograms in R plot device one at a time
#' @param breaks the number of breaks in the histograms
#'
#' @return
#'
#' @examples
#'
#' @export
#'
#' @seealso \code{\link{epgwr_mc}}, \code{\link{print_maps}}
#'
#' @author Jaakko Madetoja
#'
#' @references Madetoja, J. (2018). Error propagation in geographically weighted
#'   regression. (Doctoral dissertation, Aalto University). Manuscript in
#'   preparation.

# Function for printing the histograms for single-value metrics and highlighting the original value in them, as well as calculating the number of NA's
# data is the output from epgwr_mc
# print_file TRUE (default) for creating .png images and saving them to disk, FALSE for just showing histograms in R plot device one at a time
# breaks is the number of breaks in the histograms
print_histograms <- function(data, print_file=TRUE, breaks=20) {
  opar <- par() # Need to change parameters, so storing the old here
  if(print_file) {
    output_folder <- choose.dir(caption="Select folder for the outputs")
  }
  for(i in 1:ncol(data$simul_metrics)) {
    if(print_file){
      png(filename=paste0(output_folder, "/histogram_", colnames(data$simul_metrics)[i], ".png"), width=720, height=480)
    }
    par(mar=c(7, 4, 4, 2) + 0.1) # New parameters for margins
    par(mgp=c(5,1,0)) # More lines for text in margin
    if(sum(!is.na(data$simul_metrics[, i]))>0) { # Check that at least one value exists
      if(sum(!is.na(data$simul_metrics[,i])) != 1) { # There's more than one value
        if(sd(data$simul_metrics[,i], na.rm=TRUE)!=0) { # There's more than one unique value
          hist(data$simul_metrics[,i], main=paste0("Histogram of ", colnames(data$simul_metrics)[i], "\n", nrow(data$simul_metrics), " simulations; ", signif(sum(is.na(data$simul_metrics[,i])) / nrow(data$simul_metrics) * 100, 3), "% N/A's", "\n","original value:", signif(data$original_simul_metrics[,i], 3)), xlab=paste(strwrap(strsplit(data$info_simul_metrics[i], split=";")[[1]], width=100), collapse="\n"), breaks=seq(min(data$simul_metrics[,i], na.rm=TRUE), max(data$simul_metrics[,i], na.rm=TRUE), length=breaks), right=FALSE, xlim=range(c(data$simul_metrics[,i], data$original_simul_metrics[,i]), na.rm=TRUE))
          abline(v=data$original_simul_metrics[,i], col="red", lwd=3)
        }
        else { # All values are the same
          hist(data$simul_metrics[,i], main=paste0("Histogram of ", colnames(data$simul_metrics)[i], "\n", nrow(data$simul_metrics), " simulations; ", signif(sum(is.na(data$simul_metrics[,i])) / nrow(data$simul_metrics) * 100, 3), "% N/A's", "\n","original value:", signif(data$original_simul_metrics[,i], 3)), xlab=paste(strwrap(strsplit(data$info_simul_metrics[i], split=";")[[1]], width=100), collapse="\n"), right=FALSE, xlim=range(c(data$simul_metrics[,i], data$original_simul_metrics[,i]), na.rm=TRUE))
          abline(v=data$original_simul_metrics[,i], col="red", lwd=3)
        }
      }
      else { # There's only one value
        hist(data$simul_metrics[,i], main=paste0("Histogram of ", colnames(data$simul_metrics)[i], "\n", nrow(data$simul_metrics), " simulations; ", signif(sum(is.na(data$simul_metrics[,i])) / nrow(data$simul_metrics) * 100, 3), "% N/A's", "\n","original value:", signif(data$original_simul_metrics[,i], 3)), xlab=paste(strwrap(strsplit(data$info_simul_metrics[i], split=";")[[1]], width=100), collapse="\n"), right=FALSE, xlim=range(c(data$simul_metrics[,i], data$original_simul_metrics[,i]), na.rm=TRUE))
        abline(v=data$original_simul_metrics[,i], col="red", lwd=3)
      }
    }
    else { # If all are NA's, output an empty histogram
      hist(0, main=paste0("Histogram of ", colnames(data$simul_metrics)[i], "\n", nrow(data$simul_metrics), " simulations; ", signif(sum(is.na(data$simul_metrics[,i])) / nrow(data$simul_metrics) * 100, 3), "% N/A's", "\n","original value:", signif(data$original_simul_metrics[,i], 3)), xlab=paste(strwrap(strsplit(data$info_simul_metrics[i], split=";")[[1]], width=100), collapse="\n"), border="white")
      abline(v=data$original_simul_metrics[,i], col="red", lwd=3)
    }
    if(print_file) {
      invisible(dev.off())
    }
    else {
      readline(prompt="Press [Enter] to continue")
    }
  }
  suppressWarnings(par(opar)) # Restore old parameters
  invisible()
}

#' @title Print maps
#'
#' @description This function creates maps for point metrics resulting from
#'   error propagation in geographically weighted regression
#'
#' @param data output from the function \code{\link{epgwr_mc}}
#' @param spatialdataframe a SpatialPointsDataFrame or SpatialPolygonsDataFrame
#'   object as defined in package \code{\pkg{sp}} that will be used for plotting
#' @param print_file TRUE (default) for creating .png images and saving them to
#'   disk, FALSE for just showing maps in R plot device one at a time
#' @param point_size the size of the point when plotting SpatialPointDataFrame
#' @param inner_point_size the size for the inner points used for original
#'   values
#' @param fill_colors vector of colors used to fill the points/polygons; default
#'   is from white ot black (rev(grey(0:20 / 20))); first and last will be used
#'   for original points; needs at least 10 colors (otherwise they will be
#'   recycled)
#' @param plot_border TRUE (default) for plotting black borders for points;
#'   FALSE is useful when there are lots of points
#'
#' @return
#'
#' @examples
#'
#' @export
#'
#' @seealso \code{\link{epgwr_mc}}, \code{\link{print_maps}}
#'
#' @author Jaakko Madetoja
#'
#' @references Madetoja, J. (2018). Error propagation in geographically weighted
#'   regression. (Doctoral dissertation, Aalto University). Manuscript in
#'   preparation.


# Function to create maps for point metrics
# data is the output from epgwr_mc
# spatialdataframe is a Spatial*DataFrame object (Point or Polygons) that will be used for plotting; it needs @data and @polygons or @points
# print_file TRUE (default) for creating .png images and saving them to disk, FALSE for just showing the maps in R plot device one at a time
# point_size is the size (or actually scaling value) of the point when plotting SpatialPointDataFrame
# inner_point_size is the size (or actually scaling value) for the inner point used for original values
# fill_colors is the vector of colors used to fill the points/polygons; default is from white ot black (rev(grey(0:20 / 20))); first and last will be used for original points; needs at least 10 colors (otherwise they will be recycled)
# plot_border TRUE (default) for plotting black borders for points, FALSE for not; FALSE is useful when there are lots of points
print_maps <- function(data, spatialdataframe, print_file=TRUE, plot_original=TRUE, point_size=2, inner_point_size=0.8, fill_colors=rev(grey(0:20 / 20)), plot_border=TRUE) {
  require(sp)
  if(print_file) {
    output_folder <- choose.dir(caption="Select folder for the outputs")
  }
  for(i in 1:ncol(data$point_metrics)) {
    if(print_file){
      png(filename=paste0(output_folder, "/map_", colnames(data$point_metrics)[i], ".png"), width=720, height=720)
    }
    spatialdataframe@data$new_attribute <- data$point_metrics[,i] # Add new attribute the be plotted
    if(plot_original) { # This code will plot original values on top of the map
      original_colours <- data$original_point_metrics[,i] # These are created for the fill colours
      original_colours[original_colours==0] <- fill_colors[1]
      original_colours[original_colours==1] <- rev(fill_colors)[1]
      original_points.spdf <- SpatialPointsDataFrame(coordinates(spatialdataframe), data.frame(1:nrow(data$point_metrics))) # Empty SpatialPointsDataFrame object for plotting original values
      if(class(spatialdataframe) == "SpatialPolygonsDataFrame") { # For polygons
        print(spplot(obj=spatialdataframe, zcol="new_attribute", at=c(0:9 / 10, 1.01), main=paste0("Map of ", colnames(data$point_metrics)[i], "\nwith original values on top", "\n", paste(strwrap(strsplit(data$info_point_metrics[i], split=";")[[1]], width=100), collapse="\n")), col.regions=fill_colors, sp.layout=list(list("sp.points", original_points.spdf, pch=16, cex=inner_point_size, col=original_colours))))
      }
      else if(class(spatialdataframe) == "SpatialPointsDataFrame") { # For points
        if(plot_border) { # Plot borders
          print(spplot(obj=spatialdataframe, zcol="new_attribute", cex=point_size, cuts=c(0:10 / 10), key.space="right", main=paste0("Map of ", colnames(data$point_metrics)[i], "\nwith original values on top", "\n", paste(strwrap(strsplit(data$info_point_metrics[i], split=";")[[1]], width=100), collapse="\n")), col.regions=fill_colors, sp.layout=list(list("sp.points", original_points.spdf, pch=1, cex=point_size, col="black"), list("sp.points", original_points.spdf, pch=16, cex=inner_point_size, col=original_colours))))
        }
        else { # Don't plot borders
          print(spplot(obj=spatialdataframe, zcol="new_attribute", cex=point_size, cuts=c(0:10 / 10), key.space="right", main=paste0("Map of ", colnames(data$point_metrics)[i], "\nwith original values on top", "\n", paste(strwrap(strsplit(data$info_point_metrics[i], split=";")[[1]], width=100), collapse="\n")), col.regions=fill_colors, sp.layout=list(list("sp.points", original_points.spdf, pch=16, cex=inner_point_size, col=original_colours))))
        }
      }
      else stop("Need Spatial*DataFrame")
    }
    else { # Print only the simulation values
      if(class(spatialdataframe) == "SpatialPolygonsDataFrame") {
        print(spplot(obj=spatialdataframe, zcol="new_attribute", at=c(0:9 / 10, 1.01), main=paste0("Map of ", colnames(data$point_metrics)[i]), col.regions=fill_colors))
      }
      else if(class(spatialdataframe) == "SpatialPointsDataFrame") {
        if(plot_border) {
          print(spplot(obj=spatialdataframe, zcol="new_attribute", cex=point_size, cuts=c(0:10 / 10), key.space="right", main=paste0("Map of ", colnames(data$point_metrics)[i], "\n", paste(strwrap(strsplit(data$info_point_metrics[i], split=";")[[1]], width=100), collapse="\n")), col.regions=fill_colors, sp.layout=list(list("sp.points", spatialdataframe, pch=1, cex=point_size, col="black"))))
        }
        else {
          print(spplot(obj=spatialdataframe, zcol="new_attribute", cex=point_size, cuts=c(0:10 / 10), key.space="right", main=paste0("Map of ", colnames(data$point_metrics)[i], "\n", paste(strwrap(strsplit(data$info_point_metrics[i], split=";")[[1]], width=100), collapse="\n")), col.regions=fill_colors))
        }
      }
      else stop("Need Spatial*DataFrame")
    }
    if(print_file) {
      invisible(dev.off())
    }
    else {
      readline(prompt="Press [Enter] to continue")
    }
  }
  invisible()
}

#' @title Plot confidence intervals
#'
#' @description This function creates a plot for values and their 95 %
#'   confidence intervals ordered from left to right according to the value
#'
#' @param data a vector of data values
#' @param SE a vector of standard deviations
#' @param max maximum value
#' @param min minimum value
#' @param point_size the size of the points
#' @param title the title of the plot
#'
#' @return
#'
#' @examples
#'
#' @export
#'
#' @seealso
#'
#' @author Jaakko Madetoja
#'
#' @references Madetoja, J. (2018). Error propagation in geographically weighted
#'   regression. (Doctoral dissertation, Aalto University). Manuscript in
#'   preparation.


plot_confidence_intervals <- function(data, SE, max, min, point_size=0.8, title) {
  sorted_data <- data[order(data)]
  sorted_SE <- SE[order(data)]
  plot(x=1:length(data), y=sorted_data, pch=20, cex=point_size, col="brown", main=title)
  upper_values <- sorted_data+2*sorted_SE
  upper_values[upper_values>max] <- 1
  lower_values <- sorted_data-2*sorted_SE
  lower_values[lower_values<min] <- 0
  points(x=1:length(data), y=upper_values, pch=20, cex=point_size, col="chocolate3")
  points(x=1:length(data), y=lower_values, pch=20, cex=point_size, col="chocolate3")
  invisible()
}

#' @title Plot boxplots
#'
#' @description This function creates boxplots for single-value metrics
#'   resulting from error propagation in geographically weighted regression. It
#'   is useful for comparing the results of different simulations from the same
#'   GWR model, for example, using different errors.
#'
#' @param data output from the function \code{\link{epgwr_mc}}; for multiple
#'   results, use list(data1, data2, data3...) where data1 etc. are the results
#'   from the function \code{\link{epgwr_mc}}
#' @param names a vector of names for different simulations, for example
#'   names=c("large error", "small error", "spatial autocorrelation")
#' @param range this determines how far the plot whiskers extend out from the
#'   box; ff range is positive, the whiskers extend to the most extreme data
#'   point which is no more than range times the interquartile range from the
#'   box; a value of zero causes the whiskers to extend to the data extremes
#' @param plot_original TRUE (default) for plotting original value to the box
#'   plot as a red line; if using multiple simulation results, this will be
#'   taken from the first result, so make sure the original GWR cases are all
#'   the same
#' @param horizontal logical indicating if the boxplots should be horizontal;
#'   default FALSE means vertical boxes
#' @param print_file TRUE (default) for creating .png images and saving them to
#'   disk, FALSE for just showing boxplots in R plot device one at a time
#' @param margins the parameter mar which is a numerical vector of the form
#'   c(bottom, left, top, right) which gives the number of lines of margin to be
#'   specified on the four sides of the plot; the default is c(5, 4, 4, 2) +
#'   0.1.
#'
#'
#' @return
#'
#' @examples
#'
#' @export
#'
#' @seealso \code{\link{epgwr_mc}}
#'
#' @author Jaakko Madetoja
#'
#' @references Madetoja, J. (2018). Error propagation in geographically weighted
#'   regression. (Doctoral dissertation, Aalto University). Manuscript in
#'   preparation.


# For plotting boxplots of metrics; useful for comparing multiple simulations from the same GWR model (for example using different errors)
# data is the result from simulations; for multiple results, use list(data1, data2, data3...) where data1 etc. are the results from the GWR simulation function
# names is a vector of names for different simulations, for example names=c("large error", "small error", "spatial autocorrelation")
# range is straight from boxplot command: this determines how far the plot whiskers extend out from the box. If range is positive, the whiskers extend to the most extreme data point which is no more than range times the interquartile range from the box. A value of zero causes the whiskers to extend to the data extremes
# plot_original TRUE if you want to plot original value to the box plot as a red line; if using multiple simulation results, this will be taken from the first result, so make sure your original GWR cases are all the same if using this one
# horizontal is logical indicating if the boxplots should be horizontal; default FALSE means vertical boxes
# print_file TRUE creates .png files to given folder
# margins is the parameter mar which is "A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1."
plot_boxplots <- function(data, names, range=1.5, plot_original=TRUE, horizontal=FALSE, print_file=TRUE, margins=c(5,4,4,2)+0.1) {
  opar <- par() # Need to change parameters, so storing the old here
  if(print_file) {
    output_folder <- choose.dir(caption="Select folder for the outputs")
  }
  for(i in 1:ncol(data[[1]]$simul_metrics)) { # i runs through different metrics
    if(print_file){
      png(filename=paste0(output_folder, "/boxplot_", colnames(data[[1]]$simul_metrics)[i], ".png"), width=720, height=480)
    }

    par(mar=margins) # New parameters for margins

    data_temp <- cbind()
    names_temp <- c()
    empty_columns <- c() # For possible later use (inform user about which simulations had only N/A's)

    # Remove data that has only N/A's (boxplot can't plot if all values are NA)
    for (k in 1:length(data)) { # k runs through different simulations
      if(sum(!is.na(data[[k]]$simul_metrics[, i]))>0) { # If all values are NA (boxplot can't plot if all values are NA)
        data_temp <- cbind(data_temp, data[[k]]$simul_metrics[, i])
        names_temp <- c(names_temp, names[k])
      }
    }
    if(length(data_temp)>0) { # There's at least one boxplot that can be made
      data_temp <- data_temp[,c(ncol(data_temp):1)] # This mirrors the matrix so that first row in boxplot will be the uppest one
      names_temp <- names_temp[c(length(names_temp):1)] # Same for names
      boxplot(data_temp, names=names_temp, horizontal=horizontal, main=paste0("Boxplot of ", colnames(data[[1]]$simul_metrics)[i], "\n", nrow(data[[1]]$simul_metrics), " simulations", "\n","original value:", signif(data[[1]]$original_simul_metrics[,i], 3)), las=1)
      if(plot_original) { # Plot original value; different code for horizontal and vertical case
        if(horizontal) {
          abline(v=data[[1]]$original_simul_metrics[,i], col="red", lwd=3)
        }
        else {
          abline(h=data[[1]]$original_simul_metrics[,i], col="red", lwd=3)
        }
      }
    }
    else { # Plot empty plot
      plot.new()
      text(x=0.5, y=0.5, labels=paste0("Boxplot of ", colnames(data[[1]]$simul_metrics)[i], "\n", nrow(data[[1]]$simul_metrics), " simulations; ", "\n","original value:", signif(data[[1]]$original_simul_metrics[,i], 3), "\n", "all values N/A"), cex=2)
    }

    if(print_file) {
      invisible(dev.off())
    }
    else {
      readline(prompt="Press [Enter] to continue")
    }
  }
  suppressWarnings(par(opar)) # Restore old parameters
  invisible()
}
