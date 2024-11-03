library(sf)
library(raster)
library(ggplot2)
library(tigris)


CalSim <- function(scale, a, b, c) {
  #' thsi funcion does stuff
  #' Inputs:
  #' are scale for spatial correlation
  
  USstates = states(cb = FALSE, resolution = "500k", year = NULL)
  California = USstates[USstates$NAME == "California",]

  # load some spatial data. Administrative Boundary

  shp <- st_as_sf(California)
  grid <- shp %>%
    st_make_grid(cellsize = 0.1, what = "centers") %>% # grid of points
  st_intersection(shp)

  coordinates = st_coordinates(grid)

  library("CompRandFld")
  library(RandomFields)
  #and scale 0.1 and sill = 1
  A = RFsim(coordx = coordinates[, "X"], coordy = coordinates[, "Y"], corrmodel = "exponential", model = "BinaryGauss", threshold = 0, param = list(nugget = 0, mean = 0, scale = scale, sill = 0.6))
  #A = matrix(A$data, nrow = 50, ncol = 50, byrow = T)
  #MEdiation simulation
  
  #correlation isn't so high 0.2, maybe collinearity shouldn't be a huge problem in that scenario
  model <- "exponential"
  parameters <- list(mean = 0.5, sill = 1, nugget = 0.5, scale = scale)
  Con = RFsim(coordx = coordinates[, "X"], coordy = coordinates[,"Y"], corrmodel = model, model = "Gaussian", param = parameters)
  #Simulate the A dependent part 
  M_A <- a * A$data + 1.5*Con$data
  #Simulation the covariance function
  #correlation options, cauchy, exponential, gauss, gencauchy, spherical, stable, matern, wave
  model <- "exponential"

  CorrelationParam("exponential")
  parameters <- list(mean = 0, sill = 1, nugget = 0, scale = scale)

  M_cov = RFsim(coordx = coordinates[, "X"], coordy = coordinates[, "Y"], corrmodel = model, model = "Gaussian", param = parameters)
  #M_cov = matrix(M_cov$data, nrow = 50,ncol = 50, byrow = T)

  #Final M 

  M <- M_A + M_cov$data


  #Simulate part of Y dependent on M and A 

  Y_AM = b * A$data + c * M + 1.0*Con$data

  #Covariance part of Y 
  model = "exponential"
  CorrelationParam(model)
  parameters <- list(mean = 0, sill = 1, nugget = 0, scale = scale)

  Y_cov = RFsim(coordx = coordinates[, "X"], coordy = coordinates[, "Y"], corrmodel = model, model = "Gaussian", param = parameters)
  #Y_cov = matrix(Y_cov$data, nrow = 50,ncol = 50, byrow = T)

  Y = Y_AM + Y_cov$data

  Ydf <- as.data.frame(cbind(coordinates, Y))
  names(Ydf) <- c("X", "Y", "data_resp")
  library(raster)
  # create spatial points data frame
  spg <- Ydf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Yraster <- raster(spg)

  Mdf <- as.data.frame(cbind(coordinates, M))
  names(Mdf) <- c("X", "Y", "data_M")
  library(raster)
  # create spatial points data frame
  spg <- Mdf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Mraster <- raster(spg)

  Adf <- as.data.frame(cbind(coordinates, A$data))
  names(Adf) <- c("X", "Y", "data_A")
  library(raster)
  # create spatial points data frame
  spg <- Adf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Araster <- raster(spg)
  
  Cdf <- as.data.frame(cbind(coordinates, Con$data))
  names(Cdf) <- c("X", "Y", "data_C")
  library(raster)
  # create spatial points data frame
  spg <- Cdf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Craster <- raster(spg)

  raster_stack <- stack(Yraster, Araster, Mraster, Craster)
  return(raster_stack)
}



CalSim_inter <- function(scale, a, b, c, l = TRUE) {
  #' thsi funcion does stuff
  #' Inputs:
  #' are scale for spatial correlation
  
  USstates = states(cb = FALSE, resolution = "500k", year = NULL)
  California = USstates[USstates$NAME == "California",]
  
  # load some spatial data. Administrative Boundary
  
  shp <- st_as_sf(California)
  grid <- shp %>%
    st_make_grid(cellsize = 0.1, what = "centers") %>% # grid of points
    st_intersection(shp)
  
  coordinates = st_coordinates(grid)
  
  library("CompRandFld")
  library(RandomFields)
  #and scale 0.1 and sill = 1
  A = RFsim(coordx = coordinates[, "X"], coordy = coordinates[, "Y"], corrmodel = "exponential", model = "BinaryGauss", threshold = 0, param = list(nugget = 0, mean = 0, scale = scale, sill = 0.6))
  #A = matrix(A$data, nrow = 50, ncol = 50, byrow = T)
  #MEdiation simulation
  model <- "exponential"
  parameters <- list(mean = 0.5, sill = 1, nugget = 0.5, scale = scale)
  Con = RFsim(coordx = coordinates[, "X"], coordy = coordinates[,"Y"], corrmodel = model, model = "Gaussian", param = parameters)
  #Simulate the A dependent part 
  L_A <- ((a+b)/2) * A$data 
  #Simulation the covariance function
  #correlation options, cauchy, exponential, gauss, gencauchy, spherical, stable, matern, wave
  model <- "exponential"
  
  CorrelationParam("exponential")
  parameters <- list(mean = 0, sill = 1, nugget = 0, scale = scale)
  
  L_cov = RFsim(coordx = coordinates[, "X"], coordy = coordinates[, "Y"], corrmodel = model, model = "Gaussian", param = parameters)
  #M_cov = matrix(M_cov$data, nrow = 50,ncol = 50, byrow = T)
  
  #Final M 
  
  L <- L_A + L_cov$data
  #correlation isn't so high 0.2, maybe collinearity shouldn't be a huge problem in that scenario
  model <- "exponential"
  parameters <- list(mean = 0.5, sill = 1, nugget = 0.5, scale = scale)
  Con = RFsim(coordx = coordinates[, "X"], coordy = coordinates[,"Y"], corrmodel = model, model = "Gaussian", param = parameters)
  #Simulate the A dependent part 
  M_A <- a * A$data + 1.5*Con$data + (a+b)/2*L
  #Simulation the covariance function
  #correlation options, cauchy, exponential, gauss, gencauchy, spherical, stable, matern, wave
  model <- "exponential"
  
  CorrelationParam("exponential")
  parameters <- list(mean = 0, sill = 1, nugget = 0, scale = scale)
  
  M_cov = RFsim(coordx = coordinates[, "X"], coordy = coordinates[, "Y"], corrmodel = model, model = "Gaussian", param = parameters)
  #M_cov = matrix(M_cov$data, nrow = 50,ncol = 50, byrow = T)
  
  #Final M 
  
  M <- M_A + M_cov$data
  
  
  #Simulate part of Y dependent on M and A 
  
  Y_AM = b * A$data + c * M + 1.0*Con$data + (a+b)/2*L
  
  #Covariance part of Y 
  model = "exponential"
  CorrelationParam(model)
  parameters <- list(mean = 0, sill = 1, nugget = 0, scale = scale)
  
  Y_cov = RFsim(coordx = coordinates[, "X"], coordy = coordinates[, "Y"], corrmodel = model, model = "Gaussian", param = parameters)
  #Y_cov = matrix(Y_cov$data, nrow = 50,ncol = 50, byrow = T)
  
  Y = Y_AM + Y_cov$data
  
  Ydf <- as.data.frame(cbind(coordinates, Y))
  names(Ydf) <- c("X", "Y", "data_resp")
  library(raster)
  # create spatial points data frame
  spg <- Ydf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Yraster <- raster(spg)
  
  Mdf <- as.data.frame(cbind(coordinates, M))
  names(Mdf) <- c("X", "Y", "data_M")
  library(raster)
  # create spatial points data frame
  spg <- Mdf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Mraster <- raster(spg)
  
  Adf <- as.data.frame(cbind(coordinates, A$data))
  names(Adf) <- c("X", "Y", "data_A")
  library(raster)
  # create spatial points data frame
  spg <- Adf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Araster <- raster(spg)
  
  Cdf <- as.data.frame(cbind(coordinates, Con$data))
  names(Cdf) <- c("X", "Y", "data_C")
  library(raster)
  # create spatial points data frame
  spg <- Cdf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Craster <- raster(spg)
  
  Ldf <- as.data.frame(cbind(coordinates, L))
  names(Ldf) <- c("X", "Y", "data_L")
  library(raster)
  # create spatial points data frame
  spg <- Ldf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Lraster <- raster(spg)
  
  raster_stack <- stack(Yraster, Araster, Mraster, Craster, Lraster)
  return(raster_stack)
}
#Simulation Frameworks that we sould 
#High total effect/ low total effect
#high proportion mediated, low proportion mediated

#Large range and small range 
#Aggregated by zip code or not


#SHould I add interaction and a intermediary confounder 

#GWR and varying coefficient models 


#CalSim(scale = 5000, )


