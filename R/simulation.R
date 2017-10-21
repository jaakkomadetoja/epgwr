# Functions for running error propagation in geographically weighted regression
# NOTE: This has been tested only using Windows OS. It's possible that this will
# not work for Linux or Mac


#' @keywords internal

# data_in is the data
# formula_in is the formula for regression
# n_obs is the number of observations
# y_name is the name of the independent variable
# x_names is a vector of names of the depedent variables (needs to be the same size as sd_x_vec)
# sd_y is the standard deviation of dependent variable; can also be a vector of different deviations for different points; the realizations are taken from normal distribution
# sd_x_vec is a vector of standard deviations of independent variables; can also be a list of vectors of different deviations for different points; the realizations are taken from normal distribution
# sd_coord is the standard deviation of coordinates; can also be a vector of different deviations for different points; the realizations are taken from multivariate normal distribution
# y_min and y_max are minimum and maximum values for the dependent variable
# x_min and x_max are vectors of minimum and maximum values for independent variables
# adapt_in TRUE for adaptive or FALSE for fixed bandwidth
# bw in FALSE for letting software calibrate the bandwidth or a constant value
# gweight_in is geographical weighting function
# rho_in and neigh_dist are the rho value and the distance treshold for autoregressive random values
# normalize TRUE for making sure that autocorrelated random values have the same mean and sd as original random values
# skip_F123 FALSE for calculating Leung's F statistics and TRUE for skipping them
# sar_operator_in is the SAR operator for the input data
# coefficient_vectors and R2_vector are RC class vectors created in .middle_function; they will be used for calculating point metrics
.gwr_simulation <- function(data_in, formula_in, n_obs, y_name, x_names, coord_names, sd_y, sd_x_vec, sd_coord, y_min, y_max, x_min, x_max, adapt_in, bw_in, gweight_in, rho_in, neigh_dist, normalize, skip_F123, print_progress, sar_operator_in, coefficient_vectors, R2_vector) {

  if(length(sd_coord)==1) {
    rvalues <- rmvnorm(n_obs, mean=c(0,0), sigma=diag(2)*sd_coord*sd_coord) # Creates random values for coordinates when standard deviation is the same for all points
  }
  else {
    rvalues <- t(sapply(sd_coord, function(x) rmvnorm(1, mean=c(0,0), sigma=diag(2)*x*x))) # Creates random values for coordinates when standard deviation varies for different points
  }

  if(rho_in) { # Calculate autocorrelated random values for the coordinates
    rvalues1 <- sar_operator_in %*% rvalues # Spatially autocorrelated error values
    if(normalize) { # Normalize the autocorrelated errors
      if(sd(rvalues1)!=0) {
        rvalues1 <- (rvalues1/sd(rvalues1))*sd(rvalues) # sd the same as in rvalues
      }
      rvalues1 <- rvalues1 - mean(rvalues1) + mean(rvalues) # mean the same as in rvalues
    }
    rvalues <- rvalues1
  }

  data_in[[coord_names[1]]] <- data_in[[coord_names[1]]] + rvalues[,1] # Add the random values to the X coordinate
  data_in[[coord_names[2]]] <- data_in[[coord_names[2]]] + rvalues[,2] # Add the random values to the Y coordinate

  if(rho_in) { # Calculate the new sar_operator to match the new coordinates
    coord_vectors <- cbind(data_in[[coord_names[1]]], data_in[[coord_names[2]]])
    neighbours <- dnearneigh(coord_vectors, 0, neigh_dist)
    distances <- nbdists(neighbours, coord_vectors)
    inv_distances <- lapply(distances, function(x) 1/x)
    sar_operator <- invIrM(neighbours=neighbours, rho=rho_in, glist=inv_distances, style="W", feasible=TRUE)
  }

  for (i in 1:length(x_names)) {
    rvalues <- rnorm(n_obs, mean=0, sd=sd_x_vec[[i]]) # Random values for independent variables with the correct standard deviation
    if(rho_in) {
      rvalues1 <- sar_operator %*% rvalues # Spatially autocorrelated error values
      if(normalize) { # Normalize the autocorrelated errors
        if(sd(rvalues1)!=0) {
          rvalues1 <- (rvalues1/sd(rvalues1))*sd(rvalues)
        }
        rvalues1 <- rvalues1 - mean(rvalues1) + mean(rvalues)
      }
      rvalues <- rvalues1
    }
    data_in[[x_names[i]]] <- data_in[[x_names[i]]] + rvalues # Add the random values to the correct independent variable
    if(!identical(x_min, FALSE)) { # This is the right check for if x_min is not FALSE; this works with value 0 correctly as well
      data_in[[x_names[i]]][data_in[[x_names[i]]]<x_min[i]] <- x_min[i] # Values smaller than defined in x_min will be x_min
    }
    if(!identical(x_max, FALSE)) {
      data_in[[x_names[i]]][data_in[[x_names[i]]]>x_max[i]] <- x_max[i] # Same for x_max
    }
  }
  rvalues <- rnorm(n_obs, mean=0, sd=sd_y) # Random values for the dependent variable
  if(rho_in) {
    rvalues1 <- sar_operator %*% rvalues # Spatially autocorrelated error values
    if(normalize) { # Normalize the autocorrelated errors
      if(sd(rvalues1)!=0) {
        rvalues1 <- (rvalues1/sd(rvalues1))*sd(rvalues)
      }
      rvalues1 <- rvalues1 - mean(rvalues1) + mean(rvalues)
    }
    rvalues <- rvalues1
  }
  data_in[[y_name]] <- data_in[[y_name]] + rvalues # Add the random values to the dependent variable
  if(!identical(y_min, FALSE)) { # Same for y
    data_in[[y_name]][data_in[[y_name]]<y_min[i]] <- y_min[i]
  }
  if(!identical(y_max, FALSE)) {
    data_in[[y_name]][data_in[[y_name]]>y_max[i]] <- y_max[i]
  }


  lm_model <- lm(formula=formula_in, data=data_in) # Calculate global model; gwr also does this but doesn't store all necessary values (like R squared)
  if(bw_in) { # If the bandwidth has been given by the user, skip the calibration
    bw <- bw_in
  }
  else {
    bw <- gwr.sel(formula=formula_in, data=data_in, coords=cbind(data_in[[coord_names[1]]], data_in[[coord_names[2]]]), adapt=adapt_in, gweight=gweight_in, verbose=FALSE) # Bandwidth calculation
  }

  if(adapt_in) { # GWR model created here. Need two different codes for adaptive and fixed bandwidth, because gwr-function takes first with "adapt" (and ignores "bandwidth") and the second with bandwidth (and adapt=NULL)
    gwr_model <- gwr(formula=formula_in, data=data_in, coords=cbind(data_in[[coord_names[1]]], data_in[[coord_names[2]]]), adapt=bw, gweight=gweight_in, hatmatrix=TRUE)
  }
  else {
    gwr_model <- gwr(formula=formula_in, data=data_in, coords=cbind(data_in[[coord_names[1]]], data_in[[coord_names[2]]]), bandwidth=bw, gweight=gweight_in, adapt=NULL, hatmatrix=TRUE)
  }

  if(adapt_in) {
    bw <- round(gwr_model$adapt*n_obs) # The code will give out the number of points used in adaptive bandwidth; "floor" rounding is used in the results of gwr; in reality weights are calculated using weighted average of floor() and floor()+1 as can be seen in gwe.R file
  }

  # Metrics-----------------------------------------------------------

  # Simulation metrics describing OLS:

  # How good is the model
  lm_R2 <- summary(lm_model)$r.squared # R squared for global model
  lm_AIC <- AIC(lm_model) # AIC for the global model
  # Calculate Moran's I for residuals
  inv_distance_matrix <- 1/as.matrix(dist(gwr_model$SDF@coords)) # Inverses of distance values
  diag(inv_distance_matrix) <- 0 # Diagonals as zero
  lm_residuals_moran <- Moran.I(gwr_model$lm$residuals, inv_distance_matrix) # This calculates Moran's index for the residuals with inverse distance as weighting matrix
  lm_moran <- lm_residuals_moran$observed # The index value
  lm_moran_p <- lm_residuals_moran$p.value # p value

  # Simulation metrics describing GWR:

  # First Calculate local T values for all independent variables
  T_values <- list()
  for (i in 1:length(x_names)) {
    T_values[[i]] <- gwr_model$SDF@data[[x_names[i]]] / gwr_model$SDF@data[[paste0(x_names[i], "_se_EDF")]] # Yes, _EDF is the right column (it gives the (almost) same values as the ones used in GWmodel and uses the right effective degrees of freedom; I don't know what the other even is)
  }

  # Stationarity
  F1 <- NA
  F2 <- NA
  F3 <- NA
  gwr_F1_p <- NA
  gwr_F2_p <- NA
  gwr_F3_sum <- NA
  if(!skip_F123) {
    capture.output(F3 <- LMZ.F3GWR.test(gwr_model), file='NUL') # This creates F3 statistics from Leung et al. (2000) silently; works for Windows
    F1 <- LMZ.F1GWR.test(gwr_model) # F1 statistics from Leung et al. (2000)
    F2 <- LMZ.F2GWR.test(gwr_model) # F2 statistics from Leung et al. (2000)
    gwr_F1_p <- F1$p.value # p value for the whole GWR model stationarity
    gwr_F2_p <- F2$p.value # p value for the whole GWR model stationarity
    gwr_F3_sum <- sum(F3[2:nrow(F3),"Pr(>)"]>0.05, na.rm=TRUE) # Number of stationary variables
  }
  # How good is the model
  gwr_R2 <- 1 - (gwr_model$results$rss/gwr_model$gTSS) # This is the Quasi-global R2 that can be seen in results of GWR-model
  gwr_AIC <- gwr_model$results$AICh # This is the AIC; AICc would be $AICb
  gwr_R2_sd <- sd(gwr_model$SDF@data$localR2, na.rm=TRUE) # Standard deviation of local R squared values; describes how much the explanatory power varies
  # How local are the relationships; bandwidth describes this
  # Moran's I for residuals from GWR
  gwr_residuals_moran <- Moran.I(gwr_model$SDF@data$gwr.e, inv_distance_matrix) # Using inv_distance_matrix calculated before
  gwr_moran <- gwr_residuals_moran$observed # The index value
  gwr_moran_p <- gwr_residuals_moran$p.value # The p value

  # Simulation metrics comparing OLS and GWR

  # Which model is better
  dif_R2 <- gwr_R2 - lm_R2 # Difference in R squared values
  dif_AIC <- gwr_AIC - lm_AIC # Difference in AIC values
  # Model improvement with F1 and F2 tests
  # Which model has more random residuals
  dif_moran <- gwr_moran - lm_moran
  dif_moran_p <- gwr_moran_p - lm_moran_p
  # Where is GWR better
  dif_R2_percentage <-  sum((gwr_model$SDF@data$localR2 - lm_R2) > 0, na.rm=TRUE) / n_obs # Percentage of points where local R squared is better than global model R squared


  # This loop will be used for all simulation metrics that utilize separate coefficients
  lm_vif <- vector(length=length(x_names)) # VIF for global model
  lm_coef <- vector(length=length(x_names)) # Coefficient values
  lm_coef_t <- vector(length=length(x_names)) # T values for coefficients

  gwr_coef_F3_p <- vector(length=length(x_names))*NA # F3 statistics (Leung et al., 2000) for analyzing the stationarity of each independent variable
  gwr_coef_significant <- vector(length=length(x_names)) # Percentage of significant coefficient values for each independent variable
  gwr_coef_pos_significant <- vector(length=length(x_names)) # Percentage of significant positive coefficient values for each independent variable
  gwr_coef_neg_significant <- vector(length=length(x_names)) # Percentage of significant negative coefficient values for each independent variable
  gwr_coef_sd_significant <- vector(length=length(x_names)) # Standard deviation of significant T values
  gwr_coef_sd_pos_significant <- vector(length=length(x_names)) # Standard deviation of significant positive T values
  gwr_coef_sd_neg_significant <- vector(length=length(x_names)) # Standard deviation of significant negative T values
  for (i in 1:length(x_names)) {
    lm_vif[i] <- car::vif(lm_model)[[x_names[i]]]
    lm_coef[i] <- lm_model$coefficients[[x_names[i]]] # These could also be retrieved from GWR results
    lm_coef_t[i] <- summary(lm_model)$coefficients[[x_names[i],"t value"]]

    if(!skip_F123) {
      gwr_coef_F3_p[i] <- F3[[x_names[i], "Pr(>)"]] # p value for stationarity
    }
    gwr_coef_significant[i] <-  sum(abs(T_values[[i]]) > 1.96, na.rm=TRUE) / n_obs
    gwr_coef_pos_significant[i] <- sum(T_values[[i]] > 1.96, na.rm=TRUE) / n_obs
    gwr_coef_neg_significant[i] <- sum(T_values[[i]] < (-1.96), na.rm=TRUE) / n_obs
    gwr_coef_sd_significant[i] <- sd(T_values[[i]][abs(T_values[[i]]) > 1.96], na.rm=TRUE)
    gwr_coef_sd_pos_significant[i] <- sd(T_values[[i]][T_values[[i]] > 1.96], na.rm=TRUE)
    gwr_coef_sd_neg_significant[i] <- sd(T_values[[i]][T_values[[i]] < (-1.96)], na.rm=TRUE)
  }

  # Point metrics
  if(!is.null(coefficient_vectors)) { # In case point metrics need to be skipped
    for (i in 1:length(x_names)) { # This adds significant coefficients (as TRUE or FALSE) to the coefficient_vectors
      coefficient_vectors[[i]]$positive_coef_sum <- coefficient_vectors[[i]]$positive_coef_sum + (T_values[[i]] > 1.96)
      coefficient_vectors[[i]]$negative_coef_sum <- coefficient_vectors[[i]]$negative_coef_sum + (T_values[[i]] < (-1.96))
    }
    R2_vector$R2_sum <- R2_vector$R2_sum + (gwr_model$SDF@data$localR2 > lm_R2) # Points where local R squared is bigger than global
  }

  if(print_progress) {
    print(paste0(Sys.time(), ": A simulation done")) # Print this line if print_progress=TRUE
  }

  # Output from the function; note that point metrics will not be given as output, but they are stored in RC object created in .middle_function
  c(lm_R2, lm_AIC, lm_moran, lm_moran_p, gwr_F1_p, gwr_F2_p, gwr_F3_sum, gwr_R2, gwr_AIC, gwr_R2_sd, bw, gwr_moran, gwr_moran_p, dif_R2, dif_AIC, dif_moran, dif_moran_p, dif_R2_percentage, lm_vif, lm_coef, lm_coef_t, gwr_coef_F3_p, gwr_coef_significant, gwr_coef_pos_significant, gwr_coef_neg_significant, gwr_coef_sd_significant, gwr_coef_sd_pos_significant, gwr_coef_sd_neg_significant)
}


#' @keywords internal

# This is the middle function called by the calling function. Here the actual
# gwr function is called. This function is needed for the parallel computing to
# work because parSapply-command can't utilize Reference Classes (RC) correctly
# (i.e. one RC-object can't be modified by different cores at the same time).
# Here one RC-object is given to one core, so it can be modified in each
# simulation.
.middle_function <- function(n_sim, data_in, formula_in, n_obs, y_name, x_names, coord_names, sd_y, sd_x_vec, sd_coord, y_min, y_max, x_min, x_max, adapt_in, bw_in, gweight_in, rho_in, neigh_dist, normalize, skip_F123, print_progress){
  my_env <- new.env() # Required to set reference classes in a package

  # Create Reference Class -objects
  coefficient_RC <- setRefClass("coefficient_RC", fields=list(positive_coef_sum = "numeric", negative_coef_sum = "numeric"), where=my_env)
  coefficient_vectors = list()
  for (i in 1:length(x_names)) { # Create positive_ and negative_coef_sum for each independent variable
    coefficient_vectors[i] <- coefficient_RC$new(positive_coef_sum=numeric(n_obs), negative_coef_sum=numeric(n_obs))
  }
  R2_RC <- setRefClass("R2_RC", fields=list(R2_sum = "numeric"), where=my_env)
  R2_vector <- R2_RC$new(R2_sum=numeric(n_obs)) # Vector for comparing local and global R squared

  if(rho_in) { # Calculate SAR operator for the input data; this will be used to generating spatially autocorrelated random values for coordinates; the other random values utilize a SAR operator calculated using new coordinates
    coord_vectors <- cbind(data_in[[coord_names[1]]], data_in[[coord_names[2]]]) # Original coordinates
    neighbours <- dnearneigh(coord_vectors, 0, neigh_dist) # List of neighbors within given neigh_dist
    distances <- nbdists(neighbours, coord_vectors) # Distances to those neighbors
    inv_distances <- lapply(distances, function(x) 1/x) # Inverse of those distances
    sar_operator_in <- invIrM(neighbours=neighbours, rho=rho_in, glist=inv_distances, style="W", feasible=TRUE) # Matrix that can be "used for generating simultaneous autoregressive random variables"
  }
  else sar_operator_in <- NULL

  out <- t(replicate(n_sim, .gwr_simulation(data_in=data_in, formula_in=formula_in, n_obs=n_obs, y_name=y_name, x_names=x_names, coord_names=coord_names, sd_y=sd_y, sd_x_vec=sd_x_vec, sd_coord=sd_coord, y_min=y_min, y_max=y_max, x_min=x_min, x_max=x_max, adapt_in=adapt_in, bw_in=bw_in, gweight_in=gweight_in, rho_in=rho_in, neigh_dist=neigh_dist, normalize=normalize, skip_F123=skip_F123, print_progress=print_progress, sar_operator_in=sar_operator_in, coefficient_vectors=coefficient_vectors, R2_vector=R2_vector)))

  # Process the RC objects into a matrix for the output
  point_metrics <- cbind(as.vector(coefficient_vectors[[1]]$positive_coef_sum), as.vector(coefficient_vectors[[1]]$negative_coef_sum))
  for (i in 2:length(x_names)) {
    point_metrics <- cbind(point_metrics, as.vector(coefficient_vectors[[i]]$positive_coef_sum), as.vector(coefficient_vectors[[i]]$negative_coef_sum))
  }
  point_metrics <- cbind(point_metrics, as.vector(R2_vector$R2_sum))

  output <- list(simul_metrics=out, point_metrics=point_metrics)
  output
}

#' @title Error propagation in geographically weighted regression using Monte
#'   Carlo simulation
#'
#' @description This function applies error propagation in a given
#'   geographically weighted regression model using Monte Carlo simulation
#'
#' @param data_in data for the GWR model as a (data.frame) table, a
#'   SpatialPointsDataFrame or a SpatialPolygonsDataFrame object with columns as
#'   metric coordinates; for Spatial*DataFrame objects, only the data table is
#'   used, everything else is ignored
#' @param coord_names vector of X and Y coordinate names in the data_in table
#'   (this is also required for the Spatial*DataFrame data table!); defaults are
#'   "X" and "Y"; this is case-sensitive!; for example coord_names=c("X", "Y")
#' @param y_name name of the independent variable; for example y_name="PctBach"
#' @param x_names vector of names of the depedent variables; needs to be the
#'   same size as sd_x_vec; for example x_names=c("TotPop90", "PctRural",
#'   "PctEld", "PctFB", "PctPov", "PctBlack")
#' @param sd_y standard deviation of dependent variable; can also be a vector of
#'   different deviations for different points; for example sd_y=0.5 or
#'   sd_y=c(0.5, 0.7, 0.5, 0.2, 0.5, 0.8, 0.2, 0.5, 0.4, 0.9)
#' @param sd_x_vec vector of standard deviations of independent variables; can
#'   also be a list of vectors of different deviations for different points; for
#'   example sd_x_vec=c(500, 8, 1, 0.05, 1.5, 2) or sd_x_vec=list(500, 8, 1,
#'   0.05, c(0.5, 0.7, 0.5, 0.2, 0.5, 0.8, 0.2, 0.5, 0.4, 0.9), 2)
#' @param sd_coord standard deviation of coordinates; can also be a vector of
#'   different deviations for different points; for example sd_coord=5000 or
#'   sd_coord=c(5000, 7000, 5000, 2000, 5000, 8000, 2000, 5000, 4000, 9000)
#' @param y_min minimum value for the dependent variable; use default FALSE if
#'   no such values are valid
#' @param y_max maximum value for the dependent variable; use default FALSE if
#'   no such values are valid
#' @param x_min vector of minimum values for independent variables; use default
#'   FALSE if no such values are valid; for example x_min=C(0,0,0,0,0,0)
#' @param x_max vector of maximum values for independent variables; use default
#'   FALSE if no such values are valid; for example x_max=c(999999,1,1,1,1,1)
#' @param n_sim number of simulations
#' @param multicore number of cores (or actually processes that are run on
#'   different cores; default is 4) used in calculation; set to FALSE (or 1) to
#'   use just one
#' @param adapt_in TRUE for adaptive or FALSE (default) for fixed bandwidth
#' @param bw_in FALSE (default) for letting the software calibrate the
#'   bandwidth, or a set value to skip bandwidth calibration and use a constant;
#'   NOTE: if using adaptive bandwidth, give a value between ]0,1] as the
#'   percentage of points used
#' @param gweight_in geographical weighting function gwr.Gauss, gwr.gauss or
#'   gwr.bisquare (default)
#' @param rho_in rho value for autoregressive random values (more than 0 and
#'   less than 1) or FALSE (default) if spatial autocorrelation is not used;
#'   NOTE: high rho values will distort the distribution (mean values will
#'   increase/decrease and sd will increase if neigh_dist is small) of error
#'   values for each realization; use normalize=TRUE to fix this problem
#' @param neigh_dist distance treshold for autoregressive random values, d_max;
#'   increasing this will decrease level of spatial autocorrelation of error
#'   values; minimum value for this is the largest nearest neighbor distance
#'   after the realizations of errors have been applied to coordinates, so use a
#'   small distance for which around each point there are a couple of points
#' @param normalize TRUE (default) for normalizing the autocorrelated errors to
#'   match the same mean and sd as the non-correlated errors; heavily
#'   recommended
#' @param skip_F123 FALSE (default) for calculating Leung's F123 statistics;
#'   TRUE for skipping them as they can take quite a long time to calculate
#' @param seed seed for the RNG simulations; default 42; NULL for no seed
#'   (hasn't been tested)
#' @param print_progress TRUE for printing a message each time a simulation is
#'   done, FALSE (default) for skipping it; if multicore=FALSE, printing happens
#'   via the normal R console, otherwise the user is asked a directory where a
#'   text file will be stored and the messages are printed there; this is
#'   because normal printing is not possible with parallel computing in Windows
#'
#' @note The tool has been tested using Windows OS and might not work for Linux
#'   or Mac
#'
#'@return A list including
#'  \item{simul_metrics}{a matrix for single-value metrics}
#'  \item{point_metrics}{a matrix for point metrics}
#'  \item{original_simul_metrics}{a vector for original values for single-value metrics}
#'  \item{original_point_metrics}{a vector for original values for point metrics}
#'  \item{info_simul_metrics}{a vector for explanations for single-value metrics}
#'  \item{info_point_metrics}{a vector for explanations for point metrics}
#'  \item{function_call}{function call}
#'
#'
#' @examples
#' \dontrun{
#' ## An example with artificial data
#' size <- 20 # size of the square; there will be (size+1)^2 points
#' u <- rep(0:size, times=size+1) # x coordinate, runs like 0, 1, 2, 0, 1, 2, 0, 1, 2...
#' v <- rep(0:size, each=size+1) # y coordinate, runs like 0, 0, 0, 1, 1, 1, 2, 2, 2...
#'
#' # Beta's:
#' b0 <- 0 # Skip intercept
#' b1 <- 0.5 # Constant coefficient
#' b2 <- (u + v) / (max(u) + max(v)) # Linear coefficient
#' b3 <- sin((1/max(u))*pi*u) # y-coordinate constant, sin curve in x-direction
#'
#' # x's and e
#' RNGkind("L'Ecuyer-CMRG") # Using the same RNG as in simulations
#' set.seed(42) # Set the seed
#' x1 <- runif((size+1)^2) # Uniform distribution [0,1]
#' x2 <- runif((size+1)^2)
#' x3 <- runif((size+1)^2)
#' e0 <- rnorm((size+1)^2, mean=0, sd=0.25) # Random residuals
#'
#' # y
#' y <- b0 + b1*x1 + b2*x2 + b3*x3 + e0
#'
#' simul_data <- data.frame(y, x1, x2, x3, b0, b1, b2, b3, e0, u, v)
#'
#' artificial_results <- epgwr_mc(data_in=simul_data, y_name="y",
#' x_names=c("x1", "x2", "x3"), coord_names=c("u", "v"), sd_y=0.4,
#' sd_x_vec=c(0.2, 0.2, 0.2), sd_coord=0, n_sim=10, multicore=2, adapt_in=FALSE,
#' gweight_in=gwr.bisquare, rho_in = FALSE, neigh_dist = FALSE)
#'
#' # Create histograms for the single-value metrics
#' print_histograms(data=artificial_results)
#' # Create maps for the point metrics
#' simul_points <- SpatialPointsDataFrame(cbind(u, v), simul_data)
#' print_maps(data=artificial_results, spatialdataframe=simul_points)
#'
#' rm(size, u, v, b0, b1, b2, b3, x1, x2, x3, e0, y, simul_data,
#' artificial_results) # Remove results from the environment
#'
#' ## An example with Georgia data
#' data(georgia) # Georgia data set from the package spgwr
#' georgia_table <- gSRDF@data
#' data(georgia_sd) # Standard deviations for the Georgia data
#'
#' georgia_results <- epgwr_mc(data_in=georgia_table, y_name="PctBach",
#' x_names=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov", "PctBlack"),
#' coord_names=c("X", "Y"), sd_y=georgia_sd[["PctBach_sd"]],
#' sd_x_vec=list(georgia_sd[["TotPop90_sd"]], georgia_sd[["PctRural_sd"]],
#' georgia_sd[["PctEld_sd"]], georgia_sd[["PctFB_sd"]],
#' georgia_sd[["PctPov_sd"]], georgia_sd[["PctBlack_sd"]]), sd_coord=0, y_min=0,
#' y_max=100, x_min=c(0,0,0,0,0,0), x_max=c(9999999, 100, 100, 100, 100, 100),
#' n_sim=100, multicore=4, adapt_in=FALSE, gweight_in=gwr.bisquare, rho_in =
#' FALSE, neigh_dist = FALSE)
#'
#' # Visualize results
#' print_histograms(data=georgia_results)
#' print_maps(data=georgia_results, spatialdataframe=gSRDF, print_file=FALSE)
#'
#' }
#'
#' @export
#'
#' @seealso \code{\link{print_histograms}}, \code{\link{print_maps}},
#'   \code{\link{plot_boxplots}}
#'
#' @author Jaakko Madetoja
#'
#' @references Madetoja, J. (2018). Error propagation in geographically weighted
#'   regression. (Doctoral dissertation, Aalto University). Manuscript in
#'   preparation.


# data_in is the data; should be a (data.frame) table, a SpatialPointsDataFrame or a SpatialPolygonsDataFrame object with columns as coordinates (metric); for Spatial*DataFrame objects, only @data table is used, everything else is ignored
# coord_names is a vector of X and Y coordinate names in the data_in table (this is also required for the Spatial*DataFrame @data table!); defaults are "X" and "Y"; this is case-sensitive!; for example coord_names=c("X", "Y")
# y_name is the name of the independent variable; for example y_name="PctBach"
# x_names is a vector of names of the depedent variables; needs to be the same size as sd_x_vec; for example x_names=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov", "PctBlack")
# sd_y is the standard deviation of dependent variable; can also be a vector of different deviations for different points; for example sd_y=0.5 or sd_y=c(0.5, 0.7, 0.5, 0.2, 0.5, 0.8, 0.2, 0.5, 0.4, 0.9)
# sd_x_vec is a vector of standard deviations of independent variables; can also be a list of vectors of different deviations for different points; for example sd_x_vec=c(500, 8, 1, 0.05, 1.5, 2) or sd_x_vec=list(500, 8, 1, 0.05, c(0.5, 0.7, 0.5, 0.2, 0.5, 0.8, 0.2, 0.5, 0.4, 0.9), 2)
# sd_coord is the standard deviation of coordinates; can also be a vector of different deviations for different points; for example sd_coord=5000 or sd_coord=c(5000, 7000, 5000, 2000, 5000, 8000, 2000, 5000, 4000, 9000)
# y_min and y_max are minimum and maximum values for the dependent variable; use default FALSE if no such values are valid; for example y_min=0 y_max=1
# x_min and x_max are vectors of minimum and maximum values for independent variables; use default FALSE if no such values are valid; for example x_min=C(0,0,0,0,0,0) m_max=c(999999,1,1,1,1,1)
# n_sim is the number of simulations
# multicore is the number of cores (or actually processes that are run on different cores; default is 4) used in calculation; set to FALSE (or 1) to use just one
# adapt_in TRUE for adaptive or FALSE (default) for fixed bandwidth
# bw_in FALSE (default) for letting the software calibrate the bandwidth, or a set value to skip bandwidth calibration and use a constant; NOTE: if using adaptive bandwidth, give a value between ]0,1] as the percentage of points used
# gweight_in is geographical weighting function gwr.Gauss, gwr.gauss or gwr.bisquare (I changed bisquare to default)
# rho_in is the rho value for autoregressive random values (more than 0 and less than 1) or FALSE (default) if spatial autocorrelation is not used; NOTE: high rho values will distort the distribution (mean values will increase/decrease and sd will increase if neigh_dist is small) of error values for each realization; use normalize=TRUE to fix this problem
# neigh_dist is the distance treshold for autoregressive random values; increasing this will descrease level of spatial autocorrelation of error values; minimum value for this is the largest nearest neighbor distance after the realizations of errors have been applied to coordinates; TLDR: use a small distance for which around each point there are a couple of points
# normalize TRUE (default) for normalizing the autocorrelated errors to match the same mean and sd as the non-correlated errors; heavily recommended
# skip_F123 FALSE (default) for calculating Leung's F123 statistics; TRUE for skipping them (as they can take quite a long time to calculate)
# seed is the seed for the RNG simulations; default 42; NULL for no seed (hasn't been tested)
# print_progress TRUE for printing a message each time a simulation is done, FALSE (default) for skipping it; if multicore=FALSE, this happens via the normal R console, otherwise the user is asked a directory where a text file will be stored and the message are printed there; this is because normal printing is not possible with parallel computing in Windows
epgwr_mc <- function(data_in, y_name, x_names, coord_names = c("X", "Y"), sd_y, sd_x_vec, sd_coord, y_min = FALSE, y_max = FALSE, x_min = FALSE, x_max = FALSE, n_sim, multicore = 4, adapt_in = FALSE, bw_in = FALSE, gweight_in = gwr.bisquare, rho_in = FALSE, neigh_dist = FALSE, normalize = TRUE, skip_F123 = FALSE, seed = 42, print_progress = FALSE) {

  # The required packages
  #require(GWmodel) # not using this package
  require(spgwr)
  require(mvtnorm)
  require(car)
  require(ape)
  require(sp)

  function_call <- match.call()

  if(multicore>1 && (n_sim/multicore)%%1!=0) { # If using multicore and n_sim/multicore is not an integer number
    n_sim <- (round(n_sim/multicore))*multicore
    print(paste("Each core needs to have the same number of simulations; rounding the number of simulations to", n_sim))
  }

  if(length(x_names)!=length(sd_x_vec)) stop("Wrong number of standard deviations or names of independent variables")

  if(class(data_in) == "SpatialPolygonsDataFrame" || class(data_in) == "SpatialPointsDataFrame") {
    data_in <- data_in@data # This takes the actual data (as a data.frame) from Spatial*DataFrame object
  }

  if(rho_in) {
    require(spdep)
    try(if(!neigh_dist) stop("Distance treshold for autoregressive random variables missing"))
    if(rho_in<=0 || rho_in>=1) { # rho_in is out of feasible range ]0,1[
      stop(paste(rho_in, "is not a feasible parameter for rho"))
    }
  }

  f <- formula(paste(y_name, "~", paste(x_names, collapse="+"))) # This should work for creating a formula from text input

  n_obs <- length(data_in[[y_name]]) # This calculates the number of observations (i.e. points)

  assign("gweight_in", gweight_in, envir=globalenv()) # Assigning the weighting function to global environment variable gweight_in; this seems to be the easiest way to make sure the F test calculations work; without this, when those functions access gwr_model$gweight and get "gweight_in", the functions don't understand what it is

  if(multicore>1) { # using multiple cores
    require(parallel)
    if(print_progress) {
      output_folder <- choose.dir(caption="Select folder for the log.txt file")
    }
    tic <- proc.time()[3] # Used for measuring how long it takes to run the simulations; adding it here so that choosing the folder above won't distort the time to run the simulations
    n_sim_per_core <- n_sim/multicore # Number of simulations for each core
    if(print_progress) {
      cl <- makeCluster(multicore, outfile=paste0(output_folder, "/log.txt")) # Initiate cluster with logging
    }
    else {
      cl <- makeCluster(multicore) # Initiate cluster
    }
    data_in <- data_in; y_name <- y_name; x_names <- x_names; coord_names <- coord_names; sd_y <- sd_y; sd_x_vec <- sd_x_vec; sd_coord <- sd_coord; y_min <- y_min; y_max <- y_max; x_min <- x_min; x_max <- x_max; adapt_in <- adapt_in; bw_in <- bw_in; gweight_in <- gweight_in; rho_in <- rho_in; neigh_dist <- neigh_dist; normalize <- normalize; skip_F123 <- skip_F123; print_progress <- print_progress # These are needed or the clusterExport won't work; it's because clusterExport exports variables from some environment, so we'll have to put these in there
    RNGkind("L'Ecuyer-CMRG") # This is suggested type for RNG with parallel (it might come automatically with the next command)
    clusterSetRNGStream(cl, iseed=seed) # This makes sure that different processes run independent and reproducible random-number streams
    clusterExport(cl=cl, varlist=c("n_sim_per_core", ".gwr_simulation", ".middle_function", "data_in", "coord_names", "y_name", "x_names", "sd_y", "sd_x_vec", "sd_coord", "y_min", "y_max", "x_min", "x_max", "n_sim", "adapt_in", "bw_in", "gweight_in", "f", "n_obs", "rho_in", "neigh_dist", "normalize", "skip_F123", "print_progress"), envir=environment())
    # NOTE: I also need to export special packages if needed with GWR with clusterEvalQ(cl, library(GWmodel))
    clusterEvalQ(cl, library(spgwr))
    clusterEvalQ(cl, library(mvtnorm))
    clusterEvalQ(cl, library(car))
    clusterEvalQ(cl, library(ape))
    clusterEvalQ(cl, library(spdep))
    out <-
      t(parSapply(cl, integer(multicore), eval.parent(substitute(
        function(...)
          .middle_function(
            n_sim = n_sim_per_core, data_in = data_in, formula_in = f, n_obs = n_obs, y_name = y_name, x_names =
              x_names, coord_names = coord_names, sd_y = sd_y, sd_x_vec = sd_x_vec, sd_coord = sd_coord, y_min = y_min, y_max = y_max, x_min = x_min, x_max = x_max, adapt_in = adapt_in, bw_in = bw_in, gweight_in = gweight_in, rho_in = rho_in, neigh_dist = neigh_dist, normalize = normalize, skip_F123 = skip_F123, print_progress = print_progress
          )
      )), simplify = FALSE)) # This is the same code as used in replicate command (using this because there's no replicate for parallel); using simplify=FALSE to make the processing of output easier
    stopCluster(cl)

    # The following is used for formatting the data into simple matrices
    simul_metrics <- out[[1]]$simul_metrics # This will be a matrix containing metrics for each simulation
    point_metrics <- out[[1]]$point_metrics # This will be a matrix containing  metrics for each point
    for (i in 2:multicore) {
      simul_metrics <- rbind(simul_metrics, out[[i]]$simul_metrics)
      point_metrics <- point_metrics + out[[i]]$point_metrics
    }
  }
  else { # using just one core
    tic <- proc.time()[3] # Used for measuring how long it takes to run the simulations
    RNGkind("L'Ecuyer-CMRG") # Using the same RNG as with multicore (although the results are not exactly the same with and without parallel)
    set.seed(seed=seed) # This can't reproduce results with parallel computing
    out <-
      .middle_function(
        n_sim = n_sim, data_in = data_in, formula_in = f, n_obs = n_obs, y_name = y_name, x_names =
          x_names, coord_names = coord_names, sd_y = sd_y, sd_x_vec = sd_x_vec, sd_coord = sd_coord, y_min = y_min, y_max = y_max, x_min = x_min, x_max = x_max, adapt_in = adapt_in, bw_in = bw_in, gweight_in = gweight_in, rho_in = rho_in, neigh_dist = neigh_dist, normalize = normalize, skip_F123 = skip_F123, print_progress = print_progress
      )
    simul_metrics <- out$simul_metrics
    point_metrics <- out$point_metrics
  }

  point_metrics <- point_metrics / n_sim # Percentage of simulations that give TRUE value for each location  (for example the percentage of simulations in which a certain coefficient is significantly positive (or local R2 is bigger than global); one value for each location)

  # Calculate the metrics using the original data
  original <- .middle_function(n_sim=1, data_in=data_in, formula_in=f, n_obs=n_obs, y_name=y_name, x_names=x_names, coord_names=coord_names, sd_y=0, sd_x_vec=numeric(length=length(x_names)), sd_coord=0, y_min = FALSE, y_max = FALSE, x_min = FALSE, x_max = FALSE, adapt_in=adapt_in, gweight_in=gweight_in, bw_in=bw_in, rho_in=FALSE, neigh_dist=FALSE, skip_F123=skip_F123, print_progress=FALSE)
  original_simul_metrics <- original$simul_metrics # These will be the simulation metrics for the original data; one value per metric (as there is only one simulation)
  original_point_metrics <- original$point_metrics # These will be the point metrics for the original data; TRUE or FALSE for each location


  # Create names for the columns in the results
  simul_names <- c("lm_R2", "lm_AIC", "lm_moran", "lm_moran_p", "gwr_F1_p", "gwr_F2_p", "gwr_F3_sum", "gwr_R2", "gwr_AIC", "gwr_R2_sd", "bw", "gwr_moran", "gwr_moran_p", "dif_R2", "dif_AIC", "dif_moran", "dif_moran_p", "dif_R2_percentage")

  # Create explanations for different metrics:
  info_simul_metrics <- c("lm_R2: R squared value for global (OLS) model; The bigger the value the better the global model", "lm_AIC: AIC value for global (OLS) model; The smaller the value the better the global model", "lm_moran: Moran's Index for residuals from global (OLS) model; The smaller the index the less autocorrelation is present in the residuals and the better the global model", "lm_moran_p: p value for Moran's Index for residuals from global (OLS) model; The bigger the value the less autocorrelation is present in the residuals and the better the global model", "gwr_F1_p: p value for the F1 statistics for the GWR model stationarity; The smaller the value the better the GWR model is compared to OLS", "gwr_F2_p: p value for the F2 statistics for the GWR model stationarity; The smaller the value the better the GWR model is compared to OLS", "gwr_F3_sum: The number of variables where the p value for F3 statistics > 0.05; The number of stationary variables", "gwr_R2: R squared value for local (GWR) model; The bigger the value the better the local model", "gwr_AIC: AIC value for local (GWR) model; The smaller the value the better the local model", "gwr_R2_sd: Standard deviation of local (GWR) R squared values; The bigger the value the more the explanatory power varies in the local (GWR) model", "bw: The size of the bandwidth in local (GWR) model; Distance for fixed bandwidth and number of points for adaptive", "gwr_moran: Moran's Index for residuals from local (GWR) model; The smaller the index the less autocorrelation is present in the residuals and the better the local model", "gwr_moran_p: p value for Moran's Index for residuals from local (GWR) model; The bigger the value the less autocorrelation is present in the residuals and the better the local model", "dif_R2: The difference between R squared value for local (GWR) and global (OLS) model; The bigger the value the better the local model is compared to the global", "dif_AIC: The difference between AIC value for local (GWR) and global (OLS) model; The more negative the value the better the local model is compared to the global", "dif_moran: The difference between Moran Index values for residuals for local (GWR) and global (OLS) model; The more negative the value the better the local model is compared to the global", "dif_moran_p: The difference between p values Moran Index for residuals for local (GWR) and global (OLS) model; The bigger the value the better the local model is compared to the global", "dif_R2_percentage: Percentage of points where R squared value is better for local (GWR) model than global (OLS) model; The bigger the value the better the local model is compared to the global")


  # Vectors for names and explanations
  lm_vif <- c()
  info_lm_vif <- c()
  lm_coef <- c()
  info_lm_coef <- c()
  lm_coef_t <- c()
  info_lm_coef_t <- c()
  gwr_coef_F3_p <- c()
  info_gwr_coef_F3_p <- c()
  gwr_coef_significant <- c()
  info_gwr_coef_significant <- c()
  gwr_coef_pos_significant <- c()
  info_gwr_coef_pos_significant <- c()
  gwr_coef_neg_significant <- c()
  info_gwr_coef_neg_significant <- c()
  gwr_coef_sd_significant <- c()
  info_gwr_coef_sd_significant <- c()
  gwr_coef_sd_pos_significant <- c()
  info_gwr_coef_sd_pos_significant <- c()
  gwr_coef_sd_neg_significant <- c()
  info_gwr_coef_sd_neg_significant <- c()

  point_names <- c()
  info_point_metrics <- c()

  for(i in 1:length(x_names)){
    lm_vif <- c(lm_vif, paste0("lm_vif_", x_names[i]))
    info_lm_vif <- c(info_lm_vif, paste0("lm_vif_", x_names[i], ": VIF value for the variable ", x_names[i], "; The bigger the value the more multicollinearity exists in global (OLS) model, i.e. the more correlated the variable is with the other independent variables"))
    lm_coef <- c(lm_coef, paste0("lm_coef_", x_names[i]))
    info_lm_coef <- c(info_lm_coef, paste0("lm_coef_", x_names[i], ": Coefficient value for the variable ", x_names[i], " in the global (OLS) model"))
    lm_coef_t <- c(lm_coef_t, paste0("lm_coef_t_", x_names[i]))
    info_lm_coef_t <- c(info_lm_coef_t, paste0("lm_coef_t_", x_names[i], ": T value for the coefficient for the variable ", x_names[i], " in the global (OLS) model; The more positive or negative the value the more significantly different from zero the coefficient"))
    gwr_coef_F3_p <- c(gwr_coef_F3_p, paste0("gwr_coef_F3_p_", x_names[i]))
    info_gwr_coef_F3_p <- c(info_gwr_coef_F3_p, paste0("gwr_coef_F3_p_", x_names[i], ": p value for the F3 statistics for the stationarity of the coefficient for the variable ", x_names[i], "; The smaller the value the less stationary the coefficient and more locally varying the relationship between the independent and dependent variable"))
    gwr_coef_significant <- c(gwr_coef_significant, paste0("gwr_coef_significant_", x_names[i]))
    info_gwr_coef_significant <- c(info_gwr_coef_significant, paste0("gwr_coef_significant_", x_names[i], ": Percentage of significant local (GWR) coefficient values for the variable ", x_names[i], "; The bigger the value the more significant the local relationship between the independent and dependent variable in the area"))
    gwr_coef_pos_significant <- c(gwr_coef_pos_significant, paste0("gwr_coef_pos_significant_", x_names[i]))
    info_gwr_coef_pos_significant <- c(info_gwr_coef_pos_significant, paste0("gwr_coef_pos_significant_", x_names[i], ": Percentage of positive significant local (GWR) coefficient values for the variable ", x_names[i], "; The bigger the value the more areas with significantly positive relationship between the independent and dependent variable"))
    gwr_coef_neg_significant <- c(gwr_coef_neg_significant, paste0("gwr_coef_neg_significant_", x_names[i]))
    info_gwr_coef_neg_significant <- c(info_gwr_coef_neg_significant, paste0("gwr_coef_neg_significant_", x_names[i], ": Percentage of negative significant local (GWR) coefficient values for the variable ", x_names[i], "; The bigger the value the more areas with significantly negative relationship between the independent and dependent variable"))
    gwr_coef_sd_significant <- c(gwr_coef_sd_significant, paste0("gwr_coef_sd_significant_", x_names[i]))
    info_gwr_coef_sd_significant <- c(info_gwr_coef_sd_significant, paste0("gwr_coef_sd_significant_", x_names[i], ": Standard deviation of significant local (GWR) coefficient values for the variable ", x_names[i], "; The bigger the value the more local the relationship between the independent and dependent variable"))
    gwr_coef_sd_pos_significant <- c(gwr_coef_sd_pos_significant, paste0("gwr_coef_sd_pos_significant_", x_names[i]))
    info_gwr_coef_sd_pos_significant <- c(info_gwr_coef_sd_pos_significant, paste0("gwr_coef_sd_pos_significant_", x_names[i], "; Standard deviation of significant positive local (GWR) coefficient values for the variable ", x_names[i]))
    gwr_coef_sd_neg_significant <- c(gwr_coef_sd_neg_significant, paste0("gwr_coef_sd_neg_significant_", x_names[i]))
    info_gwr_coef_sd_neg_significant <- c(info_gwr_coef_sd_neg_significant, paste0("gwr_coef_sd_neg_significant_", x_names[i], "; Standard deviation of significant negative local (GWR) coefficient values for the variable ", x_names[i]))

    point_names <- c(point_names, paste0("positive_coef_sum_", x_names[i]), paste0("negative_coef_sum_", x_names[i]))
    info_point_metrics <- c(info_point_metrics, paste0("positive_coef_sum_", x_names[i], ": Percentage of simulations that result in a positive significant coefficient value for the variable ", x_names[i]), paste0("negative_coef_sum_", x_names[i], ": Percentage of simulations that result in a negative significant coefficient value for the variable ", x_names[i]))
  }

  simul_names <- c(simul_names, lm_vif, lm_coef, lm_coef_t, gwr_coef_F3_p, gwr_coef_significant, gwr_coef_pos_significant, gwr_coef_neg_significant, gwr_coef_sd_significant, gwr_coef_sd_pos_significant, gwr_coef_sd_neg_significant)

  info_simul_metrics <- c(info_simul_metrics, info_lm_vif, info_lm_coef, info_lm_coef_t, info_gwr_coef_F3_p, info_gwr_coef_significant, info_gwr_coef_pos_significant, info_gwr_coef_neg_significant, info_gwr_coef_sd_significant, info_gwr_coef_sd_pos_significant, info_gwr_coef_sd_neg_significant)
  info_point_metrics <- c(info_point_metrics, "R2_sum: Percentage of simulations that result in a bigger local (GWR) than global (OLS) R squared value; The bigger the value the more constantly better the local (GWR) than local (OLS) model in a given location")

  output <- list(simul_metrics=simul_metrics, point_metrics=point_metrics, original_simul_metrics=original_simul_metrics, original_point_metrics=original_point_metrics, info_simul_metrics=info_simul_metrics, info_point_metrics=info_point_metrics, function_call=function_call)

  point_names <- c(point_names, "R2_sum")
  colnames(output$simul_metrics) <- simul_names  # Set the column names
  colnames(output$point_metrics) <- point_names
  colnames(output$original_simul_metrics) <- simul_names # Metrics for the original data are the same as for simulations
  colnames(output$original_point_metrics) <- point_names

  rm(gweight_in, envir=globalenv()) # Remove the created object from before

  toc <- proc.time()[3] - tic
  cat("Total run time:", round(toc), "seconds")

  output
}

#' @title Standard deviations for the Georgia data
#'
#' @description Standard deviations for the attributes for the Georgia census
#'   data set
#'
#' @docType data
#'
#' @keywords data
#'
#' @usage data(georgia_sd)
#'
#' @format A data frame object with 159 observations for the 7 variables
#'   TotPop90_sd, PctRural_sd, PctBach_sd, PctEld_sd, PctFB_sd, PctPov_sd,
#'   PctBlack_sd
#'
#' @source Created by Madetoja (2018) based on the 1990 U.S. Census of
#'   Population documents
#'
#' @references Madetoja, J. (2018). Error propagation in geographically weighted
#'   regression. (Doctoral dissertation, Aalto University). Manuscript in
#'   preparation.
"georgia_sd"
