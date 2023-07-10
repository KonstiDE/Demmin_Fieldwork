library(sf)
library(terra)
library(dplyr)
library(ranger)
library(basemaps)
library(ggplot2)
library(ggspatial)
library(ggthemes)

esus <- st_read("data/esus.gpkg")
sentinel2 <- rast("data/sentinel2.tif")
gt <- read.csv("data/demmin_gt.csv")

ggplot() +
  annotation_map_tile(zoomin = -1) +
  geom_sf(data = esus)


gt$short <- paste0(gt$field, gt$esu, gt$ssu)
gt$weight_wet <- as.numeric(gt$weight_wet)
gt$weight_bowl <- as.numeric(gt$weight_bowl)

ssus <- cbind(esus, extract(sentinel2, esus))

reg_df <- as.data.frame(ssus)
reg_df <- merge(reg_df, gt, by = "short", all = FALSE)

reg_df$weight_wet <- reg_df$weight_wet - reg_df$weight_bowl

reg_df$evi <- (2.5 * (reg_df$B8 - reg_df$B4)) / (reg_df$B8 + (2.4 * reg_df$B4) + 10000)
reg_df$ndvi <- ((reg_df$B8 - reg_df$B4) / (reg_df$B8 + reg_df$B4))

trainIndex <- nrow(reg_df) * 0.8
train_df <- reg_df[0:trainIndex,]
valid_df <- reg_df[(trainIndex):nrow(reg_df) + 1,]

train_df <- train_df[,c("ndvi", "evi", "weight_wet")]
valid_df <- valid_df[,c("ndvi", "evi", "weight_wet")]

fit_ndvi <- svm(formula = weight_wet ~ ndvi, data = train_df)
fit_evi <- svm(formula = weight_wet ~ evi, data = train_df)

preds_ndvi <- predict(fit_ndvi, valid_df)
preds_evi <- predict(fit_evi, valid_df)

mae_ndvi <- abs(sum(valid_df$weight_wet - preds_ndvi) / nrow(valid_df))
mae_evi <- abs(sum(valid_df$weight_wet - preds_evi) / nrow(valid_df))


sen2pixels <- extract(sentinel2, seq_along(0:length(values(sentinel2))))
sen2pixels$evi <- (2.5 * (sen2pixels$B8 - sen2pixels$B4)) / (sen2pixels$B8 + (2.4 * sen2pixels$B4) + 10000)
sen2pixels$ndvi <- ((sen2pixels$B8 - sen2pixels$B4) / (sen2pixels$B8 + sen2pixels$B4))

preds_fields <- predict(fit_ndvi, sen2pixels)

sentinel2_pred_ndvi <- sentinel2[[1]]
sentinel2_pred_ndvi[!is.na(sentinel2_pred_ndvi)] <- as.numeric(preds_fields)[!is.na(preds_fields)]
sentinel2_pred_ndvi[sentinel2_pred_ndvi < 0] <- 0

ggplot() +
  layer_spatial(data = sentinel2_pred_ndvi) +
  scale_fill_gradientn(
    colors = c("black", "green"),
    guide_colorbar(title = "AGB")
  ) +
  xlab("Above Ground Biomass (Wet) in gram per 5 stems") +
  theme_tufte()






