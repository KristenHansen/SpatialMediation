#High total effect, mid proportion mediated #36%
raster_stack = CalSim(5, 0.94, 4.0, 2.35)
mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
col.bw <- gwr.sel(data_resp ~ data_A + data_M, data = med_spdf)

gwr_outcome <- gwr(data_resp ~ data_A + data_M, data = med_spdf, bandwidth = col.bw)
#predicted_surface = spplot(gwr_outcome$SDF, "pred")
#true_sim_outcome = spplot(med_spdf, "data_resp")
#grid.arrange(predicted_surface,true_sim_outcome, nrow = 1, ncol = 2)

#RMSE
rmse = sqrt(mean((gwr_outcome$SDF$pred - med_spdf$data_resp) ^ 2, na.rm = T))
rmse/mean(med_spdf$data_resp)
#0.210
hist(rnorm(100))
med.bw <- gwr.sel(data_M ~ data_A, data = med_spdf)
#med.bw <- 5
gwr_mediator <- gwr(data_M ~ data_A, data = med_spdf, bandwidth = 0.1297171)
rmse_M = sqrt(mean((gwr_mediator$SDF$pred - med_spdf$data_M) ^ 2, na.rm = T))
rmse_M/abs(mean(med_spdf$data_M))
#0.262

#tot.bw <- gwr.sel(data_resp ~ data_A , data = med_spdf)
#tot.bw <- 5
#gwr_tot <- gwr(data_resp ~ data_A , data = med_spdf, bandwidth = tot.bw)


#spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model
gwr_outcome$SDF$NIE = gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A
#spplot(gwr_outcome$SDF, "NIE")

library("raster")

library("akima")
r <- rasterFromXYZ(as.data.frame(gwr_outcome$SDF)[, c("x", "y", "pred", "NIE", "data_A")])


my_col <- terrain.colors(10)
r_df <- as.data.frame(r, xy = TRUE)
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = data_A)) +
  scale_fill_viridis_c(name = "Spatial NDE prediction") +
  coord_quickmap()

ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = NIE)) +
  scale_fill_viridis_c(name = "Spatial NIE") +
  coord_quickmap()

ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = pred)) +
  scale_fill_viridis_c(name = "Y prediction") +
  coord_quickmap()

raster_truedf = as.data.frame(raster_stack, xy = TRUE)
ggplot() +
  geom_raster(data = raster_truedf, aes(x = x, y = y, fill = data_resp)) +
  scale_fill_viridis_c(name = "Y truth") +
  coord_quickmap()


#High total effect low proportion mediated (4%)
raster_stack = CalSim(5, .32, 6.0, 0.8)
mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
col.bw <- gwr.sel(data_resp ~ data_A + data_M, data = med_spdf)

gwr_outcome <- gwr(data_resp ~ data_A + data_M, data = med_spdf, bandwidth = col.bw)
#predicted_surface = spplot(gwr_outcome$SDF, "pred")
#true_sim_outcome = spplot(med_spdf, "data_resp")
#grid.arrange(predicted_surface,true_sim_outcome, nrow = 1, ncol = 2)
# Call:
#   gwr(formula = data_resp ~ data_A + data_M, data = med_spdf, bandwidth = col.bw)
# Kernel function: gwr.Gauss 
# Fixed bandwidth: 0.3421006 
# Summary of GWR coefficient estimates at data points:
#   Min.  1st Qu.   Median  3rd Qu.     Max.  Global
# X.Intercept. -2.61761 -0.67717 -0.42641 -0.17554  0.23807 -0.7624
# data_A        5.65143  5.97389  6.09990  6.29139  7.29337  6.4517
# data_M       -0.40744  0.63299  0.79533  0.92341  1.57921  0.7024

#RMSE
rmse = sqrt(mean((gwr_outcome$SDF$pred - med_spdf$data_resp) ^ 2, na.rm = T))
rmse/mean(med_spdf$data_resp)
#0.219

med.bw <- gwr.sel(data_M ~ data_A, data = med_spdf)
#med.bw <- 5
gwr_mediator <- gwr(data_M ~ data_A, data = med_spdf, bandwidth = 0.3814709 )

# Fixed bandwidth: 0.3421006 
# Summary of GWR coefficient estimates at data points:
#   Min.   1st Qu.    Median   3rd Qu.      Max. Global
# X.Intercept. -1.219611 -0.444849  0.190050  0.605955  1.458223 0.1042
# data_A       -1.672915  0.066254  0.281158  0.447131  1.182027 0.1700

rmse_M = sqrt(mean((gwr_mediator$SDF$pred - med_spdf$data_M) ^ 2, na.rm = T))
rmse_M/abs(mean(med_spdf$data_M))
#0.285


#tot.bw <- gwr.sel(data_resp ~ data_A, data = med_spdf)
#tot.bw <- 5
#gwr_tot <- gwr(data_resp ~ data_A, data = med_spdf, bandwidth = tot.bw)


#spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model
gwr_outcome$SDF$NIE = gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A
#spplot(gwr_outcome$SDF, "NIE")

library("raster")

library("akima")
r <- rasterFromXYZ(as.data.frame(gwr_outcome$SDF)[, c("x", "y", "pred", "NIE", "data_A")])


my_col <- terrain.colors(10)
r_df <- as.data.frame(r, xy = TRUE)
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = data_A)) +
  scale_fill_viridis_c(name = "Spatial NDE prediction") +
  coord_quickmap()

ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = NIE)) +
  scale_fill_viridis_c(name = "Spatial NIE") +
  coord_quickmap()

ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = pred)) +
  scale_fill_viridis_c(name = "Y prediction") +
  coord_quickmap()

raster_truedf = as.data.frame(raster_stack, xy = TRUE)
ggplot() +
  geom_raster(data = raster_truedf, aes(x = x, y = y, fill = data_resp)) +
  scale_fill_viridis_c(name = "Y truth") +
  coord_quickmap()


#Low total effect low proportion mediated
raster_stack = CalSim(5, 0.1012, .6, .252)
mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
col.bw <- gwr.sel(data_resp ~ data_A + data_M, data = med_spdf)

gwr_outcome <- gwr(data_resp ~ data_A + data_M, data = med_spdf, bandwidth = col.bw)
#predicted_surface = spplot(gwr_outcome$SDF, "pred")
#true_sim_outcome = spplot(med_spdf, "data_resp")
#grid.arrange(predicted_surface,true_sim_outcome, nrow = 1, ncol = 2)

# Fixed bandwidth: 0.5701801 
# Summary of GWR coefficient estimates at data points:
#   Min.   1st Qu.    Median   3rd Qu.      Max.  Global
# X.Intercept. -1.160432 -0.725387  0.247810  0.619420  1.506446 -0.5279
# data_A       -0.761044  0.572697  0.734777  0.875215  2.246528  1.4486
# data_M       -0.746126  0.048609  0.220966  0.423923  1.010132  0.0873
#RMSE
rmse = sqrt(mean((gwr_outcome$SDF$pred - med_spdf$data_resp) ^ 2, na.rm = T))
rmse
#0.296 #Higher RMSE for low total low mediated effect

med.bw <- gwr.sel(data_M ~ data_A, data = med_spdf)
#med.bw <- 5
gwr_mediator <- gwr(data_M ~ data_A, data = med_spdf, bandwidth = med.bw) #fit is not good here
# Call:
#   gwr(formula = data_M ~ data_A, data = med_spdf, bandwidth = med.bw)
# Kernel function: gwr.Gauss 
# Fixed bandwidth: 0.5701801 
# Summary of GWR coefficient estimates at data points:
#   Min.   1st Qu.    Median   3rd Qu.      Max.  Global
# X.Intercept. -0.570520 -0.012920  0.467209  0.756817  0.899212  0.1524
# data_A       -1.945858 -0.660287 -0.186638  0.052372  0.906959 -0.1803

rmse_M = sqrt(mean((gwr_mediator$SDF$pred - med_spdf$data_M) ^ 2, na.rm = T))
rmse_M
#0.353

tot.bw <- gwr.sel(data_resp ~ data_A, data = med_spdf)
#tot.bw <- 5
gwr_tot <- gwr(data_resp ~ data_A, data = med_spdf, bandwidth = tot.bw)


#spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model
gwr_outcome$SDF$NIE = gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A
#spplot(gwr_outcome$SDF, "NIE")

library("raster")

library("akima")
r <- rasterFromXYZ(as.data.frame(gwr_outcome$SDF)[, c("x", "y", "pred", "NIE", "data_A")])


my_col <- terrain.colors(10)
r_df <- as.data.frame(r, xy = TRUE)
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = data_A)) +
  scale_fill_viridis_c(name = "Spatial NDE prediction") +
  coord_quickmap()

#on average not well estimated (A effect on M switch to negative)
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = NIE)) +
  scale_fill_viridis_c(name = "Spatial NIE") +
  coord_quickmap()

ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = pred)) +
  scale_fill_viridis_c(name = "Y prediction") +
  coord_quickmap()

raster_truedf = as.data.frame(raster_stack, xy = TRUE)
ggplot() +
  geom_raster(data = raster_truedf, aes(x = x, y = y, fill = data_resp)) +
  scale_fill_viridis_c(name = "Y truth") +
  coord_quickmap()

