#Aggregate raster data into zip codes
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Mediation\ Analysis")
source("CalSim.R")
AggSim <- function(raster_stack){
  # set.seed(360)
  #raster_stack = CalSim(scale, a, b, c)
  #hfp_meso <- st_crop(humanFp, mesoam)
  library("spmoran")
  library(rgdal)
  polygons <- readOGR(dsn = "~/Downloads/ShapeFileZCTA/", layer = "cb_2016_us_zcta510_500k") 
  #plot(polygons,add=TRUE)
  polygons$ZCTA5CE10 = as.numeric(as.character(polygons$ZCTA5CE10))
  calZCTA = polygons[polygons$ZCTA5CE10 >90000,]
  calZCTA = calZCTA[calZCTA$ZCTA5CE10 < 96162,]
  Cal<- map_data("state", region = "california")
  #california <- filter(MainStates,region ==  "california" )
  ## extract mean value for each polygon
  library(raster)
  v1 <- extract(raster_stack, calZCTA, fun=mean, na.rm=TRUE,sp = T)
  v1$data_A = ifelse(v1$data_A>0.5, 1,0)
  
  library(spatialEco)
  
  library(spgwr)
  v1 = sp.na.omit(v1, margin = 1)
  
  start_gwr = Sys.time()
  col.bw <- bw.gwr(data_resp ~ data_A + data_M + data_C, data = med_spdf, kernel = "gaussian", longlat = T)
  
  gwr_outcome <- gwr.basic(data_resp ~ data_A + data_M + data_C, data = med_spdf, bw = col.bw, kernel = "gaussian", longlat = T)
  #predicted_surface = spplot(gwr_outcome$SDF, "pred")
  #true_sim_outcome = spplot(med_spdf, "data_resp")
  #grid.arrange(predicted_surface,true_sim_outcome, nrow = 1, ncol = 2)
  
  #RMSE
  rmse = sqrt(mean((gwr_outcome$SDF$yhat - med_spdf$data_resp) ^ 2, na.rm = T))
  #rmse/mean(med_spdf$data_resp)
  #0.210
  
  med.bw <- bw.gwr(data_M ~ data_A + data_C, data = v1, kernel = "gaussian", longlat = T)
  #med.bw <- 5
  gwr_mediator <- gwr.basic(data_M ~ data_A + data_C, data = v1, bw = med.bw, kernel = "gaussian", longlat = T)
  rmse_M = sqrt(mean((gwr_mediator$SDF$yhat - med_spdf$data_M) ^ 2, na.rm = T))
  #rmse_M/abs(mean(med_spdf$data_M))
  #0.262
  #rmse_Med = rmse_M
  
  # #tot.bw <- gwr.sel(data_resp ~ data_A , data = med_spdf)
  # #tot.bw <- 5
  # #gwr_tot <- gwr(data_resp ~ data_A , data = med_spdf, bandwidth = tot.bw)
  # 
  # 
  # #spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model
  gwr_outcome$SDF$NIE = gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A
  # #gwr_outcome$SDF$Y = med_spdf$data_resp
  gwr_outcome$SDF$nie_SE = sqrt(gwr_mediator$SDF$data_A_SE^2 + gwr_outcome$SDF$data_M_SE^2)
  # #spplot(gwr_outcome$SDF, "NIE")
  # 
  # library("raster")
  # 
  # library("akima")
  # #r <- rasterFromXYZ(as.data.frame(gwr_outcome$SDF)[, c("x", "y", "pred", "NIE", "data_A")])
  end1 = Sys.time() - start_gwr
  
  #SVC
  v2 <- as.data.frame(v1)
  zip_cent <- read.csv("~/Downloads/zip_gps.csv", header = T)

  library(dplyr)
  df_data <- left_join(v2, zip_cent, by = c("ZCTA5CE10" = "ZCTA5"))

  #Now do the SVC stuff
  start <- Sys.time()

  df_data = na.omit(df_data)
  coords <- df_data[,c("FinalLon","FinalLat")]
  #y_mod <- besf_vc(df_data$data_resp, df_data[,c("data_A", "data_C")], coords = coords, x_sel = FALSE)
  
  y_mod <- besf_vc(df_data$data_resp, df_data[,c("data_A", "data_M", "data_C")], coords = coords, x_sel = FALSE)
  rmse_svc = sqrt(mean((y_mod$pred - df_data$data_resp) ^ 2, na.rm = T))
  rmse_svc
  #y_mod <- besf_vc(df_data$data_resp, df_data$data_A, coords = coords)

  Avar <- df_data[,c(4,6)]
  Mvar <- df_data$data_M
  m_mod = besf_vc(df_data$data_M, df_data[,c("data_A", "data_C")], coords = coords)

  rmse_m_svc = sqrt(mean((m_mod$pred - df_data$data_M) ^ 2, na.rm = T))

  df_data$NDE = y_mod$b_vc[,2]

  df_data$NIE = y_mod$b_vc[,3] * m_mod$b_vc[,2]

  df_data$Y_pred = y_mod$pred

  df_data$se_NDE = y_mod$bse_vc$data_A

  df_data$se_NIE = sqrt(y_mod$bse_vc$data_M^2 + m_mod$bse_vc$data_A^2)

  end2 <- Sys.time() - start # 1.73 minutes
  
  #return(list("data_svc" = df_data, "data_gwr" = gwr_outcome$SDF, "rmse_gwr" = rmse, "rmse_m_gwr" = rmse_M, "rmse_svc" = rmse_svc, "rmse_m_svc" = rmse_m_svc, "time_svc" = end2, "time_gwr" = end1))
  return(list( "data_gwr" = gwr_outcome$SDF, "rmse_gwr" = rmse, "rmse_m_gwr" = rmse_M, "time_gwr" = end1))
  
}

S1 = AggSim(raster_S1) #seed 350
S2 = AggSim(raster_S2)
median(S2$data_gwr$data_A_SE)
median(S2$data_gwr$nie_SE)

Ci_low = S2$data_gwr$data_A - 2*S2$data_gwr$data_A_SE
ci_high = S2$data_gwr$data_A + 2*S2$data_gwr$data_A_SE
bool_s2_nde = ifelse(Ci_low < 0.6 & ci_high > 0.6, TRUE, FALSE)
sum(bool_s2_nde)/length(bool_s2_nde)

Ci_low = S2$data_gwr$NIE - 2*S2$data_gwr$nie_SE
ci_high = S2$data_gwr$NIE + 2*S2$data_gwr$nie_SE
bool_s2_nie = ifelse(Ci_low < 0.025 & ci_high > 0.025, TRUE, FALSE)
sum(bool_s2_nie)/length(bool_s2_nde)
S3 = AggSim(raster_S3) #This setting is very difficult

median(S3$data_gwr$data_A_SE)
median(S3$data_gwr$nie_SE)

Ci_low = S3$data_gwr$data_A - 2*S3$data_gwr$data_A_SE
ci_high = S3$data_gwr$data_A + 2*S3$data_gwr$data_A_SE
bool_S3_nde = ifelse(Ci_low < 0.6 & ci_high > 0.6, TRUE, FALSE)
sum(bool_S3_nde)/length(bool_S3_nde)

Ci_low = S3$data_gwr$NIE - 2*S3$data_gwr$nie_SE
ci_high = S3$data_gwr$NIE + 2*S3$data_gwr$nie_SE
bool_s3_nie = ifelse(Ci_low < 0.025 & ci_high > 0.025, TRUE, FALSE)
sum(bool_s3_nie)/length(bool_s3_nie)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s3_nde = ifelse(Ci_low < 0.6 & ci_high > 0.6, TRUE, FALSE)
sum(bool_s3_nde)/length(bool_s1_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s3_nie = ifelse(Ci_low < 0.025 & ci_high > 0.025, TRUE, FALSE)
sum(bool_s3_nie)/length(bool_s1_nde)


S4 = AggSim(raster_S4)

median(S4$data_gwr$data_A_SE)
median(S4$data_gwr$nie_SE)

Ci_low = S4$data_gwr$data_A - 2*S4$data_gwr$data_A_SE
ci_high = S4$data_gwr$data_A + 2*S4$data_gwr$data_A_SE
bool_S4_nde = ifelse(Ci_low < 0.4 & ci_high > 0.4, TRUE, FALSE)
sum(bool_S4_nde)/length(bool_S4_nde)

Ci_low = S4$data_gwr$NIE - 2*S4$data_gwr$nie_SE
ci_high = S4$data_gwr$NIE + 2*S4$data_gwr$nie_SE
bool_s4_nie = ifelse(Ci_low < 0.225 & ci_high > 0.225, TRUE, FALSE)
sum(bool_s4_nie)/length(bool_s4_nie)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s4_nde = ifelse(Ci_low < 0.4 & ci_high > 0.4, TRUE, FALSE)
sum(bool_s4_nde)/length(bool_s1_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s4_nie = ifelse(Ci_low < 0.225 & ci_high > 0.225, TRUE, FALSE)
sum(bool_s4_nie)/length(bool_s1_nde)

S5 = AggSim(raster_S5)

median(S5$data_gwr$data_A_SE)
median(S5$data_gwr$nie_SE)

Ci_low = S5$data_gwr$data_A - 2*S5$data_gwr$data_A_SE
ci_high = S5$data_gwr$data_A + 2*S5$data_gwr$data_A_SE
bool_S5_nde = ifelse(Ci_low < 0.4 & ci_high > 0.4, TRUE, FALSE)
sum(bool_S5_nde)/length(bool_S5_nde)

Ci_low = S5$data_gwr$NIE - 2*S5$data_gwr$nie_SE
ci_high = S5$data_gwr$NIE + 2*S5$data_gwr$nie_SE
bool_S5_nie = ifelse(Ci_low < 0.225 & ci_high > 0.225, TRUE, FALSE)
sum(bool_S5_nie)/length(bool_S5_nie)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s5_nde = ifelse(Ci_low < 0.4 & ci_high > 0.4, TRUE, FALSE)
sum(bool_s5_nde)/length(bool_s1_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s5_nie = ifelse(Ci_low < 0.225 & ci_high > 0.225, TRUE, FALSE)
sum(bool_s5_nie)/length(bool_s1_nde)

S6 = AggSim(raster_S6)

median(S6$data_gwr$data_A_SE)
median(S6$data_gwr$nie_SE)

Ci_low = S6$data_gwr$data_A - 2*S6$data_gwr$data_A_SE
ci_high = S6$data_gwr$data_A + 2*S6$data_gwr$data_A_SE
bool_S6_nde = ifelse(Ci_low < 0.4 & ci_high > 0.4, TRUE, FALSE)
sum(bool_S6_nde)/length(bool_S6_nde)

Ci_low = S6$data_gwr$NIE - 2*S6$data_gwr$nie_SE
ci_high = S6$data_gwr$NIE + 2*S6$data_gwr$nie_SE
bool_S6_nie = ifelse(Ci_low < 0.225 & ci_high > 0.225, TRUE, FALSE)
sum(bool_S6_nie)/length(bool_S6_nie)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s6_nde = ifelse(Ci_low < 0.4 & ci_high > 0.4, TRUE, FALSE)
sum(bool_s6_nde)/length(bool_s1_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s6_nie = ifelse(Ci_low < 0.225 & ci_high > 0.225, TRUE, FALSE)
sum(bool_s6_nie)/length(bool_s6_nde)

S7 = AggSim(raster_S7)

median(S7$data_gwr$data_A_SE)
median(S7$data_gwr$nie_SE)

Ci_low = S7$data_gwr$data_A - 2*S7$data_gwr$data_A_SE
ci_high = S7$data_gwr$data_A + 2*S7$data_gwr$data_A_SE
bool_S7_nde = ifelse(Ci_low < 6 & ci_high > 6, TRUE, FALSE)
sum(bool_S7_nde)/length(bool_S7_nde)

Ci_low = S7$data_gwr$NIE - 2*S7$data_gwr$nie_SE
ci_high = S7$data_gwr$NIE + 2*S7$data_gwr$nie_SE
bool_S7_nie = ifelse(Ci_low < 0.25 & ci_high > 0.25, TRUE, FALSE)
sum(bool_S7_nie)/length(bool_S7_nie)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s7_nde = ifelse(Ci_low < 6.0 & ci_high > 6.0, TRUE, FALSE)
sum(bool_s7_nde)/length(bool_s7_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s7_nie = ifelse(Ci_low < 0.25 & ci_high > 0.25, TRUE, FALSE)
sum(bool_s7_nie)/length(bool_s7_nde)
S8 = AggSim(raster_S8)

median(S8$data_gwr$data_A_SE)
median(S8$data_gwr$nie_SE)

Ci_low = S8$data_gwr$data_A - 2*S8$data_gwr$data_A_SE
ci_high = S8$data_gwr$data_A + 2*S8$data_gwr$data_A_SE
bool_S8_nde = ifelse(Ci_low < 6 & ci_high > 6, TRUE, FALSE)
sum(bool_S8_nde)/length(bool_S8_nde)

Ci_low = S8$data_gwr$NIE - 2*S8$data_gwr$nie_SE
ci_high = S8$data_gwr$NIE + 2*S8$data_gwr$nie_SE
bool_S8_nie = ifelse(Ci_low < 0.25 & ci_high > 0.25, TRUE, FALSE)
sum(bool_S8_nie)/length(bool_S8_nie)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s8_nde = ifelse(Ci_low < 6.0 & ci_high > 6.0, TRUE, FALSE)
sum(bool_s8_nde)/length(bool_s7_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s8_nie = ifelse(Ci_low < 0.25 & ci_high > 0.25, TRUE, FALSE)
sum(bool_s8_nie)/length(bool_s7_nde)

S9 = AggSim(raster_S9)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s9_nde = ifelse(Ci_low < 6.0 & ci_high > 6.0, TRUE, FALSE)
sum(bool_s9_nde)/length(bool_s7_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s9_nie = ifelse(Ci_low < 0.25 & ci_high > 0.25, TRUE, FALSE)
sum(bool_s9_nie)/length(bool_s7_nde)

S10 = AggSim(raster_S10)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s10_nde = ifelse(Ci_low < 4.0 & ci_high > 4.0, TRUE, FALSE)
sum(bool_s10_nde)/length(bool_s7_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s10_nie = ifelse(Ci_low < 2.25 & ci_high > 2.25, TRUE, FALSE)
sum(bool_s10_nie)/length(bool_s7_nde)

S11 = AggSim(raster_S11)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s11_nde = ifelse(Ci_low < 4.0 & ci_high > 4.0, TRUE, FALSE)
sum(bool_s11_nde)/length(bool_s7_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s11_nie = ifelse(Ci_low < 2.25 & ci_high > 2.25, TRUE, FALSE)
sum(bool_s11_nie)/length(bool_s7_nde)

S12 = AggSim(raster_S12)

median(df_data$se_NDE)
median(df_data$se_NIE)

Ci_low = df_data$NDE - 2*df_data$se_NDE
ci_high = df_data$NDE + 2*df_data$se_NDE
bool_s12_nde = ifelse(Ci_low < 4.0 & ci_high > 4.0, TRUE, FALSE)
sum(bool_s12_nde)/length(bool_s7_nde)

Ci_low = df_data$NIE - 2*df_data$se_NIE
ci_high = df_data$NIE + 2*df_data$se_NIE
bool_s12_nie = ifelse(Ci_low < 2.25 & ci_high > 2.25, TRUE, FALSE)
sum(bool_s12_nie)/length(bool_s7_nde)

time_svc <- c(S1$time_svc, S2$time_svc, S3$time_svc, S4$time_svc, S5$time_svc, S6$time_svc, S7$time_svc, S8$time_svc, S9$time_svc, S10$time_svc, S11$time_svc, S12$time_svc)
time_gwr <- c(S1$time_gwr, S2$time_gwr, S3$time_gwr, S4$time_gwr, S5$time_gwr, S6$time_gwr, S7$time_gwr, S8$time_gwr, S9$time_gwr, S10$time_gwr, S11$time_gwr, S12$time_gwr)


rmse_svc <- c(S1$rmse_svc, S2$rmse_svc, S3$rmse_svc, S4$rmse_svc, S5$rmse_svc, S6$rmse_svc, S7$rmse_svc, S8$rmse_svc, S9$rmse_svc, S10$rmse_svc, S11$rmse_svc, S12$rmse_svc)
rmse_gwr <- c(S1$rmse_gwr, S2$rmse_gwr, S3$rmse_gwr, S4$rmse_gwr, S5$rmse_gwr, S6$rmse_gwr, S7$rmse_gwr, S8$rmse_gwr, S9$rmse_gwr, S10$rmse_gwr, S11$rmse_gwr, S12$rmse_gwr)

rmse_m_svc <- c(S1$rmse_m_svc, S2$rmse_m_svc, S3$rmse_m_svc, S4$rmse_m_svc, S5$rmse_m_svc, S6$rmse_m_svc, S7$rmse_m_svc, S8$rmse_m_svc, S9$rmse_m_svc, S10$rmse_m_svc, S11$rmse_m_svc, S12$rmse_m_svc)
rmse_m_gwr <- c(S1$rmse_m_gwr, S2$rmse_m_gwr, S3$rmse_m_gwr, S4$rmse_m_gwr, S5$rmse_m_gwr, S6$rmse_m_gwr, S7$rmse_m_gwr, S8$rmse_m_gwr, S9$rmse_m_gwr, S10$rmse_m_gwr, S11$rmse_m_gwr, S12$rmse_m_gwr)

results = cbind(time_svc, time_gwr, rmse_svc, rmse_m_svc, rmse_gwr, rmse_m_gwr)
write.csv(results, "Aggregated_Results.csv")
shape <- st_read(dsn = "~/Downloads/ShapefileZCTA/cb_2016_us_zcta510_500k.shp")
#m <- merge(shape, S12$data_svc, by = c("ZCTA5CE10"))
#spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model


#shape <- st_read(dsn = "~/Downloads/ShapefileZCTA/cb_2016_us_zcta510_500k.shp")
#m <- S1$data_gwr

#spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model



plotting <- function(S, Sname, NDE, NIE){
    
    m <- st_as_sf(S$data_gwr)
    m$Diff <- m$Y - m$pred
    Cal<- map_data("state", region = "california")
    main = paste(Sname, "_GWR_TruthPred_agg.png", sep = "")
    png(main, width = 6.5, height = 6, units = 'in', res = 300)
    ggplot()  +
      geom_sf(data = m, aes(fill = Diff), colour = NA)  + geom_polygon(data=Cal, aes(x=long, y=lat, group=group),
                                                                       colour="black", fill = NA)+ scale_fill_gradient2(low = "blue",mid = "white",
                                                                                                                        high="red", midpoint = 0) + theme(legend.key.height = unit(1.0, "cm"), legend.key.width = unit(0.3, "cm")) +
      scale_x_continuous(label = I) +
      guides(fill = guide_colourbar(barheight = unit( 2.25 , "in" ),
                                    ticks.colour = "black",
                                    ticks.linewidth = 1, frame.colour = "black",
                                    frame.linewidth = 1)) +
      labs(fill = "Truth-Pred") + theme( panel.grid = element_blank(),text = element_text(size=10), axis.title.x=element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    
   
    dev.off()
    
    main = paste(Sname, "_GWR_NDE_agg.png", sep = "")
    png(main, width = 6.5, height = 6, units = 'in', res = 300)
    ggplot()  +
      geom_sf(data = m, aes(fill = data_A), colour = NA)  + geom_polygon(data=Cal, aes(x=long, y=lat, group=group),
                                                                         colour="black", fill = NA)+ scale_fill_gradient2(low = "yellow",mid = "orange",
                                                                                                                          high="red", midpoint = NDE) + theme(legend.key.height = unit(1.0, "cm"), legend.key.width = unit(0.3, "cm")) +
      scale_x_continuous(label = I) +
      guides(fill = guide_colourbar(barheight = unit( 2.25 , "in" ),
                                    ticks.colour = "black",
                                    ticks.linewidth = 1, frame.colour = "black",
                                    frame.linewidth = 1)) +
      labs(fill = "NDE") + theme( panel.grid = element_blank(),text = element_text(size=10), axis.title.x=element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) + annotate("label", x = -116, y = 40, label = paste0("Overall NDE = ", NDE))
    
    
    
    dev.off()
    
    main = paste(Sname, "_GWR_NIE_agg.png", sep = "")
    png(main, width = 6.5, height = 6, units = 'in', res = 300)
    ggplot()  +
      geom_sf(data = m, aes(fill = NIE), colour = NA)  + geom_polygon(data=Cal, aes(x=long, y=lat, group=group),
                                                                      colour="black", fill = NA)+ scale_fill_gradient2(low = "purple",mid = "blue",
                                                                                                                       high="cyan2", midpoint = NIE) + theme(legend.key.height = unit(1.0, "cm"), legend.key.width = unit(0.3, "cm")) +
      scale_x_continuous(label = I) +
      guides(fill = guide_colourbar(barheight = unit( 2.25 , "in" ),
                                    ticks.colour = "black",
                                    ticks.linewidth = 1, frame.colour = "black",
                                    frame.linewidth = 1)) +
      labs(fill = "NIE") + theme( panel.grid = element_blank(),text = element_text(size=10), axis.title.x=element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) + annotate("label", x = -116, y = 40, label = paste0("Overall NIE = ", NIE))
    
    
    dev.off()
    
   
    m <- merge(shape, S$data_svc, by = c("ZCTA5CE10"))
    #m <- st_as_sf(S$data_svc)
    #SVC plots 
    m$Diff = m$data_resp - m$Y_pred
    
    
    
    
    
    main = paste(Sname, "_SVC_TruthPred_agg.png", sep = "")
    png(main, width = 6.5, height = 6, units = 'in', res = 300)
    ggplot()  +
      geom_sf(data = m, aes(fill = Diff), colour = NA)  + geom_polygon(data=Cal, aes(x=long, y=lat, group=group),
                                                                       colour="black", fill = NA)+ scale_fill_gradient2(low = "blue",mid = "white",
                                                                                                                        high="red", midpoint = 0) + theme(legend.key.height = unit(1.0, "cm"), legend.key.width = unit(0.3, "cm")) +
      scale_x_continuous(label = I) +
      guides(fill = guide_colourbar(barheight = unit( 2.25 , "in" ),
                                    ticks.colour = "black",
                                    ticks.linewidth = 1, frame.colour = "black",
                                    frame.linewidth = 1)) +
      labs(fill = "Truth-Pred") + theme( panel.grid = element_blank(),text = element_text(size=10), axis.title.x=element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    dev.off()
    
    main = paste(Sname, "_SVC_NDE_agg.png", sep = "")
    png(main, width = 6.5, height = 6, units = 'in', res = 300)
    ggplot()  +
      geom_sf(data = m, aes(fill = NDE), colour = NA)  + geom_polygon(data=Cal, aes(x=long, y=lat, group=group),
                                                                         colour="black", fill = NA)+ scale_fill_gradient2(low = "yellow",mid = "orange",
                                                                                                                          high="red", midpoint = NDE) + theme(legend.key.height = unit(1.0, "cm"), legend.key.width = unit(0.3, "cm")) +
      scale_x_continuous(label = I) +
      guides(fill = guide_colourbar(barheight = unit( 2.25 , "in" ),
                                    ticks.colour = "black",
                                    ticks.linewidth = 1, frame.colour = "black",
                                    frame.linewidth = 1)) +
      labs(fill = "NDE") + theme( panel.grid = element_blank(),text = element_text(size=10), axis.title.x=element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) + annotate("label", x = -116, y = 40, label = paste0("Overall NDE = ", NDE))
    
    
    dev.off()
    
    main = paste(Sname, "_SVC_NIE_Agg.png", sep = "")
    png(main, width = 6.5, height = 6, units = 'in', res = 300)
    ggplot()  +
      geom_sf(data = m, aes(fill = NIE), colour = NA)  + geom_polygon(data=Cal, aes(x=long, y=lat, group=group),
                                                                      colour="black", fill = NA)+ scale_fill_gradient2(low = "purple",mid = "blue",
                                                                                                                       high="cyan2", midpoint = NIE) + theme(legend.key.height = unit(1.0, "cm"), legend.key.width = unit(0.3, "cm")) +
      scale_x_continuous(label = I) +
      guides(fill = guide_colourbar(barheight = unit( 2.25 , "in" ),
                                    ticks.colour = "black",
                                    ticks.linewidth = 1, frame.colour = "black",
                                    frame.linewidth = 1)) +
      labs(fill = "NIE") + theme( panel.grid = element_blank(),text = element_text(size=10), axis.title.x=element_blank(), axis.title.y = element_blank(), aspect.ratio = 1, axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) + annotate("label", x = -116, y = 40, label = paste0("Overall NIE = ", NIE))
    dev.off()
}

plotting(S1, "S1",.6,0.025)
plotting(S2, "S2", 0.6, 0.025)
plotting(S3, "S3", 0.6, 0.025)
plotting(S4, "S4", 0.4,.225)
plotting(S5, "S5", 0.4, .225)
plotting(S6, "S6", 0.4, .225)
plotting(S7, "S7", 6.0, .256)
plotting(S8, "S8", 6.0,.256)
plotting(S9, "S9", 6.0, .256)
plotting(S10, "S10", 4.0, 2.21)
plotting(S11, "S11", 4.0, 2.21)
plotting(S12, "S12", 4.0, 2.21)
