hpi = read.csv("~/Downloads/HPI2.csv")
mediator = hpi$pm25
dat$exposure = ifelse(hpi$income_pctile > 50, 0, 1)
outcome = hpi$LEB_pctile


CES_data <- read.csv("~/Downloads/CES_noDemo.csv")
library(dplyr)
dat = left_join(hpi, CES_data, by = c("CensusTract" = "Census.Tract"))
outcome = dat$Low.Birth.Weight 
library(sf)
CES <- st_read(dsn = "~/Library/Mobile Documents/com~apple~CloudDocs/Mediation Analysis/CESJune2018Update_SHP/CES3June2018Update.shp")
dat$tract = dat$CensusTract
CES = merge(CES, dat, by = c("tract"))
censustracts = read.csv("~/Downloads/Census_Tracts_2020.csv")
CES = CES[which(CES$California.County == "Los Angeles" & CES$Latitude.x < 34.2& CES$Longitude.x < 118.6 ),]
library(sp)
CES_test = as_Spatial(CES)
data = CES_test
data$exposure = ifelse(data$income_pctile > 50,0, 1)
library(spmoran)
df_data <- as.data.frame(data)
df_data$exposure = ifelse(df_data$income_pctile >50, 0, 1)
#df_data$exposure = ifelse(df_data$economic_pctile >50, 0, 1)
coords <- cbind(df_data$Longitude.x, df_data$Latitude.x)
df_data = df_data[complete.cases(df_data),]

#Add confounder white_pct, black_pct, latino_pct, asian_pct, twoparents, abovepoverty, emploted, Education.Pctl

#df_data$SB.535.Disadvantaged.Community = ifelse(df_data$SB.535.Disadvantaged.Community == "Yes", 1, 0)
y_mod <- besf_vc(df_data$cvdP, df_data[,c("exposure", "pm25")], coords = coords, x_sel = FALSE)
rmse_svc = sqrt(mean((y_mod$pred - df_data$cvdP) ^ 2, na.rm = T))
#0.36

#Mvar <- df_data$data_M
m_mod = besf_vc(df_data$pm25, df_data[,c("exposure")], coords = coords)

rmse_m_svc = sqrt(mean((m_mod$pred - df_data$pm25) ^ 2, na.rm = T))
rmse_m_svc/mean(df_data$pm25)
#0.1340217
df_data$NDE = y_mod$b_vc[,2]

df_data$NIE = y_mod$b_vc[,3] * m_mod$b_vc[,2]

df_data$Y_pred = y_mod$pred

end2 <- Sys.time() - start # 1.73 minutes
#y_test = lm(Low.Birth.Weight~ pm25 + exposure, data = df_data)
#m_test = lm(pm25 ~ exposure, data = df_data)

m <- merge(CES, df_data, "tract")
m <- st_as_sf(m)
p2 <- ggplot()  +
  geom_sf(data = m, aes(fill = NDE), colour = NA) +scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0.113)

p2
p3 <- ggplot()  +
  geom_sf(data = m, aes(fill = NIE), colour = NA)  + scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0)
p3

#The mean of the errors is 0, this is a good sign
#variance of errors is 1.88
#Varaince of outcome is 2.4 and variance of prediction is only .4 though
m$diff_trunc = ifelse(m$cvdP.x - m$Y_pred > 51, 51, ifelse(m$cvdP.x - m$Y_pred < -51, -51, m$cvdP.x - m$Y_pred))
m$error_over_sd = m$diff_trunc/sd(m$cvdP.x)
p4 <- ggplot()  +
  geom_sf(data = m, aes(fill = error_over_sd), colour = NA) + scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0) + labs(fill = "(True - Pred)/sd(True)") 
p4

library(spgwr)
library("GWmodel")
data = data[complete.cases(data@data),]
col.bw = bw.gwr(cvdP ~ exposure + pm25 + white_pct + employed + Education.Pctl, data = data, kernel = "exponential", longlat = T)
#col.bw <- gwr.sel(Low.Birth.Weight ~ exposure + pm25 + white_pct + employed + Education.Pctl, data = data, longlat = T)

gwr_outcome <- gwr.basic(cvdP ~ exposure + pm25  + white_pct + employed + Education.Pctl, data = data, bw = col.bw, longlat = T, kernel = "exponential")
med.bw <- bw.gwr(pm25 ~ exposure + white_pct + employed + Education.Pctl, data = data, kernel = "exponential", longlat = T)
gwr_mediator <- gwr.basic(pm25 ~ exposure + white_pct + employed + Education.Pctl, data = data,bw = med.bw, longlat = T, kernel = "exponential")


gwr_outcome$SDF$NIE = gwr_outcome$SDF$pm25 * gwr_mediator$SDF$exposure
#test = left_join(as.data.frame(data), gwr_outcome$SDF)
gwr <- st_as_sf(gwr_outcome$SDF)
p2 <- ggplot()  +
  geom_sf(data = gwr, aes(fill = exposure), colour = NA) +scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0.113) +labs(fill = "NDE")
p2

p3 <- ggplot()  +
  geom_sf(data = gwr, aes(fill = NIE), colour = NA)  + scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0)
p3

p4 <- ggplot()  +
  geom_sf(data = gwr, aes(fill = y - yhat), colour = NA) + scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0)+labs(fill = "True Y - Predicted Y")
p4

lm(Low.Birth.Weight ~ exposure + ozone.x + white_pct + employed + Education.Pctl, data = df_data)
lm(ozone.x ~ exposure + white_pct + employed + Education.Pctl, data = df_data)
