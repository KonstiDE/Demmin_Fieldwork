library(sf)
library(terra)
library(dplyr)
library(ranger)
library(basemaps)
library(ggplot2)
library(ggspatial)
library(ggthemes)
library(e1071)
library(ggpubr)

esus <- st_read("data/esus.gpkg")
fields <- st_read("data/fields.gpkg")

sentinel2 <- rast("data/sentinel2_fields.tif")
sentinel2_fields <- mask(sentinel2, fields)

gt <- read.csv("data/demmin_gt.csv")

ggplot() +
  layer_spatial(data = sentinel2_fields) +
  geom_sf(data = esus)


gt$short <- paste0(gt$field, gt$esu, gt$ssu)
gt$weight_wet <- as.numeric(gt$weight_wet)
gt$weight_bowl <- as.numeric(gt$weight_bowl)

ssus <- cbind(esus, extract(sentinel2_fields, esus))

reg_df <- as.data.frame(ssus)
reg_df <- merge(reg_df, gt, by = "short", all = FALSE)

reg_df$weight_wet <- reg_df$weight_wet - reg_df$weight_bowl

reg_df$evi <- (2.5 * (reg_df$B8 - reg_df$B4)) / (reg_df$B8 + (2.4 * reg_df$B4) + 10000)
reg_df$ndvi <- ((reg_df$B8 - reg_df$B4) / (reg_df$B8 + reg_df$B4))

reg_df <- reg_df[sample(1:nrow(reg_df)), ]

trainIndex <- nrow(reg_df) * 0.8
train_df <- reg_df[0:trainIndex,]
valid_df <- reg_df[(trainIndex):nrow(reg_df) + 1,]

train_df <- train_df[,c("short", "ndvi", "evi", "weight_wet")]
valid_df <- valid_df[,c("short", "ndvi", "evi", "weight_wet")]

fit_ndvi <- lm(formula = weight_wet ~ ndvi, data = train_df)
fit_evi <- lm(formula = weight_wet ~ evi, data = train_df)

preds_ndvi <- predict(fit_ndvi, valid_df)
preds_evi <- predict(fit_evi, valid_df)

preds_ndvi <- as.numeric(preds_ndvi)
preds_evi <- as.numeric(preds_evi)
preds_df <- data.frame(x = seq_along(preds_ndvi), y_ndvi = preds_ndvi, weight_wet = valid_df$weight_wet, y_evi = preds_evi)

df_long <- tidyr::pivot_longer(preds_df, cols = c(y_ndvi, weight_wet, y_evi))
df_long$name <- factor(df_long$name, levels = c("y_ndvi", "weight_wet", "y_evi"))

ggplot(df_long, aes(x = as.factor(x), y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  xlab("SSUS") +
  ylab("AGB in g / 5 stems") +
  scale_fill_manual(
    values = c("#0073C2", "#017825", "salmon"),
    labels = c("Est. NDVI", "True AGB", "Est. EVI"),
    guide_colorbar(title = "")
  ) +
  scale_x_discrete(breaks = seq_along(valid_df$short), labels = valid_df$short) +
  coord_flip()

med_ndvi <- mean(valid_df$weight_wet - preds_ndvi)
med_evi <- mean(valid_df$weight_wet - preds_evi)


sen2pixels <- extract(sentinel2_fields, seq_along(0:length(values(sentinel2_fields))))
sen2pixels$evi <- (2.5 * (sen2pixels$B8 - sen2pixels$B4)) / (sen2pixels$B8 + (2.4 * sen2pixels$B4) + 10000)
sen2pixels$ndvi <- ((sen2pixels$B8 - sen2pixels$B4) / (sen2pixels$B8 + sen2pixels$B4))

# Predict with ndvi
preds_fields_ndvi <- predict(fit_ndvi, sen2pixels)

sentinel2_pred_ndvi <- sentinel2_fields[[1]]
sentinel2_pred_ndvi[!is.na(sentinel2_pred_ndvi)] <- as.numeric(preds_fields_ndvi)[!is.na(preds_fields_ndvi)]
sentinel2_pred_ndvi[sentinel2_pred_ndvi < 0] <- 0

plot_ndvi <- ggplot() +
  layer_spatial(data = sentinel2_pred_ndvi) +
  scale_fill_gradientn(
    colors = c("black", "green"),
    na.value = "transparent",
    guide_colorbar(title = "AGB"),
  ) +
  xlab("Above Ground Biomass (Wet) in gram per 5 stems")


# Predict with evi
preds_fields_evi <- predict(fit_evi, sen2pixels)

sentinel2_pred_evi <- sentinel2_fields[[1]]
sentinel2_pred_evi[!is.na(sentinel2_pred_evi)] <- as.numeric(preds_fields_evi)[!is.na(preds_fields_evi)]
sentinel2_pred_evi[sentinel2_pred_evi < 0] <- 0


plot_evi <- ggplot() +
  layer_spatial(data = sentinel2_pred_evi) +
  scale_fill_gradientn(
    colors = c("black", "green"),
    na.value = "transparent",
    guide_colorbar(title = "AGB"),
  ) +
  xlab("Above Ground Biomass (Wet) in gram per 5 stems")

# Plot together
ggarrange(
  plot_ndvi, plot_evi,
  labels = c("NDVI", "EVI")
)

# Plot difference
diff <- abs(sentinel2_pred_ndvi - sentinel2_pred_evi)
plot(diff)




