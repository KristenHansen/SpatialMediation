# California grid spatial mediation simulation

source("CalSim.R")
# Simulation settings needed
# First we must edit the range of the spatial autocorrelation
# High, what we already have 5
raster_stack <- CalSim(5, 0.3, 0.4, 0.75)
mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
col.bw <- gwr.sel(data_resp ~ data_A + data_M, data = med_spdf)

gwr_outcome <- gwr(data_resp ~ data_A + data_M, data = med_spdf, bandwidth = col.bw)

# RMSE
rmse <- sqrt(mean((gwr_outcome$SDF$pred - med_spdf$data_resp) ^ 2, na.rm = T))
rmse/abs(mean(med_spdf$data_resp))
# 0.1535605
# 0.224
med.bw <- gwr.sel(data_M ~ data_A, data = med_spdf)
gwr_mediator <- gwr(data_M ~ data_A, data = med_spdf, bandwidth = med.bw)

#tot.bw <- gwr.sel(data_resp ~ data_A, data = med_spdf)

#gwr_tot <- gwr(data_resp ~ data_A, data = med_spdf, bandwidth = tot.bw)


# NDE is data_A from this model
gwr_outcome$SDF$NIE <- gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A


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

raster_truedf <- as.data.frame(raster_stack, xy = TRUE)
ggplot() +
  geom_raster(data = raster_truedf, aes(x = x, y = y, fill = data_resp)) +
  scale_fill_viridis_c(name = "Y truth") +
  coord_quickmap()
# RMSE 0.21


# Medium, 2
raster_stack <- CalSim(2, 0.3, 0.4, 0.75)
mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
col.bw <- gwr.sel(data_resp ~ data_A + data_M, data = med_spdf)

gwr_outcome <- gwr(data_resp ~ data_A + data_M, data = med_spdf, bandwidth = col.bw)
# predicted_surface = spplot(gwr_outcome$SDF, "pred")
# true_sim_outcome = spplot(med_spdf, "data_resp")

# grid.arrange(predicted_surface,true_sim_outcome, nrow = 1, ncol = 2)

# RMSE
rmse <- sqrt(mean((gwr_outcome$SDF$pred - med_spdf$data_resp) ^ 2, na.rm = T))

rmse/mean(med_spdf$data_resp)
# 0.224
med.bw <- gwr.sel(data_M ~ data_A, data = med_spdf)
# med.bw <- 5
gwr_mediator <- gwr(data_M ~ data_A, data = med_spdf, bandwidth = med.bw)
rmse <- sqrt(mean((gwr_mediator$SDF$pred - med_spdf$data_M) ^ 2, na.rm = T))

rmse/mean(med_spdf$data_M)
tot.bw <- gwr.sel(data_resp ~ data_A, data = med_spdf)
# tot.bw <- 5
gwr_tot <- gwr(data_resp ~ data_A, data = med_spdf, bandwidth = tot.bw)


# spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model
gwr_outcome$SDF$NIE <- gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A
# spplot(gwr_outcome$SDF, "NIE")

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

raster_truedf <- as.data.frame(raster_stack, xy = TRUE)
ggplot() +
  geom_raster(data = raster_truedf, aes(x = x, y = y, fill = data_resp)) +
  scale_fill_viridis_c(name = "Y truth") +
  coord_quickmap()
# Low, 0.5 More biased as scale becomes smaller

raster_stack <- CalSim(0.5, 0.3, 0.4, 0.75)
mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
col.bw <- gwr.sel(data_resp ~ data_A + data_M, data = med_spdf)

gwr_outcome <- gwr(data_resp ~ data_A + data_M, data = med_spdf, bandwidth = col.bw)
# predicted_surface = spplot(gwr_outcome$SDF, "pred")
# true_sim_outcome = spplot(med_spdf, "data_resp")
# grid.arrange(predicted_surface,true_sim_outcome, nrow = 1, ncol = 2)

# RMSE
rmse <- sqrt(mean((gwr_outcome$SDF$pred - med_spdf$data_resp) ^ 2, na.rm = T))
rmse/mean(med_spdf$data_resp)
# 0.356
med.bw <- gwr.sel(data_M ~ data_A, data = med_spdf)
# med.bw <- 5
gwr_mediator <- gwr(data_M ~ data_A, data = med_spdf, bandwidth = 0.09727422)
rmse <- sqrt(mean((gwr_mediator$SDF$pred - med_spdf$data_M) ^ 2, na.rm = T))

rmse/mean(med_spdf$data_M)

tot.bw <- gwr.sel(data_resp ~ data_A, data = med_spdf)
# tot.bw <- 5
gwr_tot <- gwr(data_resp ~ data_A, data = med_spdf, bandwidth = tot.bw)


# spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model
gwr_outcome$SDF$NIE <- gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A
# spplot(gwr_outcome$SDF, "NIE")

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

raster_truedf <- as.data.frame(raster_stack, xy = TRUE)
ggplot() +
  geom_raster(data = raster_truedf, aes(x = x, y = y, fill = data_resp)) +
  scale_fill_viridis_c(name = "Y truth") +
  coord_quickmap()