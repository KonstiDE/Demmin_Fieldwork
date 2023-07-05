library(sf)
library(terra)
library(dplyr)

esus <- vect("data/esus.gpkg")

gt <- read.csv("data/demmin_gt.csv")
gt$short <- paste0(gt$field, gt$esu, gt$ssu)
gt$weight_wet <- as.numeric(gt$weight_wet)
gt$weight_bowl <- as.numeric(gt$weight_bowl)

sentinel2 <- rast("data/sentinel2.tif")

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

fit_ndvi <- lm(formula = weight_wet ~ ndvi, data = train_df)
fit_evi <- lm(formula = weight_wet ~ evi, data = train_df)

preds_ndvi <- predict(fit_ndvi, valid_df)
preds_evi <- predict(fit_evi, valid_df)

mae_ndvi <- sum(valid_df$weight_wet - preds_ndvi) / nrow(valid_df)
mae_evi <- sum(valid_df$weight_wet - preds_evi) / nrow(valid_df)








