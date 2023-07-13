library(sf)
library(terra)
library(dplyr)
library(ranger)
library(basemaps)
library(ggplot2)
library(ggspatial)
library(ggthemes)
library(ggpubr)
library(mvtnorm)

esus <- st_read("data/esus.gpkg")
fields <- st_read("data/fields.gpkg")

sentinel2 <- rast("data/sentinel2_fields.tif")
sentinel2_fields <- mask(sentinel2, fields)

gt <- read.csv("data/demmin_gt.csv")

ggplot() +
  layer_spatial(sentinel2) +
  geom_sf(data = fields) +
  geom_sf(data = esus) +
  theme_void()


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



# Sythetic data generation via bootstrapping with randomness
reg_df <- reg_df[,c("short", "ndvi", "evi", "weight_wet")]

additional_indices <- sample(nrow(reg_df), size = 48, replace = TRUE)
synthetic_data <- reg_df[,c("ndvi", "evi", "weight_wet")][additional_indices, ]
randomness_factor <- 0.05

for (col in c("ndvi", "evi", "weight_wet")) {
  additional_data[[col]] <- additional_data[[col]] * (1 + randomness_factor * runif(48, min = -1, max = 1))
}

synthetic_data$short <- "synthetic"
rownames(additional_data) <- NULL

reg_df <- rbind(reg_df, synthetic_data)
rownames(reg_df) <- NULL

reg_df[sample(1:nrow(reg_df)), ]

trainIndex <- nrow(reg_df) * 0.9
train_df <- reg_df[0:trainIndex,]
valid_df <- reg_df[(trainIndex + 1):nrow(reg_df),]

train_df <- train_df[,c("short", "ndvi", "evi", "weight_wet")]
valid_df <- valid_df[,c("short", "ndvi", "evi", "weight_wet")]

fit_ndvi <- lm(formula = weight_wet ~ ndvi, data = train_df)
fit_evi <- lm(formula = weight_wet ~ evi, data = train_df)
fit_both <- lm(formula = weight_wet ~ ndvi + evi, data = train_df)

preds_ndvi <- predict(fit_ndvi, valid_df)
preds_evi <- predict(fit_evi, valid_df)
preds_both <- predict(fit_both, valid_df)

mean_ndvi <- mean(valid_df$weight_wet - preds_ndvi)
mean_evi <- mean(valid_df$weight_wet - preds_evi)
mean_both <- mean(valid_df$weight_wet - preds_both)

preds_ndvi <- as.numeric(preds_ndvi)
preds_evi <- as.numeric(preds_evi)
preds_both <- as.numeric(preds_both)
preds_df <- data.frame(
  x = seq_along(preds_ndvi),
  weight_wet = valid_df$weight_wet,
  y_ndvi = preds_ndvi,
  y_evi = preds_evi,
  y_both = preds_both
)

df_long <- tidyr::pivot_longer(preds_df, cols = c(weight_wet, y_ndvi, y_evi, y_both))
df_long$name <- factor(df_long$name, levels = c("weight_wet","y_ndvi", "y_evi", "y_both"))

ggplot(df_long, aes(x = as.factor(x), y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  xlab("SSUS") +
  ylab("AGB in g / 5 stems") +
  scale_fill_manual(
    values = c("#017825","salmon", "#0073C2", "#f5c242"),
    labels = c("True AGB","Est. NDVI", "Est. EVI", "Est. Both"),
    guide_colorbar(title = "")
  ) +
  scale_x_discrete(breaks = seq_along(valid_df$short), labels = valid_df$short) +
  coord_flip()


sen2pixels <- extract(sentinel2_fields, seq_along(0:length(values(sentinel2_fields))))
sen2pixels$evi <- (2.5 * (sen2pixels$B8 - sen2pixels$B4)) / (sen2pixels$B8 + (2.4 * sen2pixels$B4) + 10000)
sen2pixels$ndvi <- ((sen2pixels$B8 - sen2pixels$B4) / (sen2pixels$B8 + sen2pixels$B4))

# Predict with ndvi
preds_fields_ndvi <- predict(fit_ndvi, sen2pixels)

sentinel2_pred_ndvi <- sentinel2_fields[[1]]
sentinel2_pred_ndvi[!is.na(sentinel2_pred_ndvi)] <- as.numeric(preds_fields_ndvi)[!is.na(preds_fields_ndvi)]
sentinel2_pred_ndvi[sentinel2_pred_ndvi < 0] <- 0

ggplot() +
  layer_spatial(data = sentinel2_pred_ndvi) +
  scale_fill_gradientn(
    colors = c("black", "green"),
    na.value = "transparent",
    guide_colorbar(title = "AGB"),
  ) +
  xlab("Above Ground Biomass (Wet) in gram per 5 stems")


# Predict with evi
preds_fields_evi <- predict(sen2pixels, fit_ndvi)

sentinel2_pred_evi <- sentinel2_fields[[1]]
sentinel2_pred_evi[!is.na(sentinel2_pred_evi)] <- as.numeric(preds_fields_evi)[!is.na(preds_fields_evi)]
sentinel2_pred_evi[sentinel2_pred_evi < 0] <- 0


ggplot() +
  layer_spatial(data = sentinel2_pred_evi) +
  scale_fill_gradientn(
    colors = c("black", "green"),
    na.value = "transparent",
    guide_colorbar(title = "AGB"),
  ) +
  xlab("Above Ground Biomass (Wet) in gram per 5 stems")

# Plot difference
diff <- sentinel2_pred_ndvi - sentinel2_pred_evi

ggplot() +
  layer_spatial(data = diff) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "transparent")


# Large Extent
sentinel2_geom <- rast("data/sentinel2_larger.tif")
fields_geom <- vect("data/fields_larger.gpkg")

sentinel2_geom <- mask(sentinel2_geom, fields_geom)

sen2large <- extract(sentinel2_geom, seq_along(0:length(values(sentinel2_geom))))
colnames(sen2large) <- c("B4", "B3", "B2", "B8")

sen2large$evi <- (2.5 * (sen2large$B8 - sen2large$B4)) / (sen2large$B8 + (2.4 * sen2large$B4) + 10000)
sen2large$ndvi <- ((sen2large$B8 - sen2large$B4) / (sen2large$B8 + sen2large$B4))

pred_large <- predict(fit_ndvi, sen2large)

sentinel2_pred_l <- sentinel2_geom[[1]]
sentinel2_pred_l[!is.na(sentinel2_pred_l)] <- as.numeric(pred_large)[!is.na(pred_large)]
sentinel2_pred_l[sentinel2_pred_l < 0] <- 0

ggplot() +
  layer_spatial(data = sentinel2_pred_l) +
  scale_fill_gradientn(
    colors = c("black", "green"),
    na.value = "transparent",
    guide_colorbar(title = "AGB"),
  ) +
  xlab("Above Ground Biomass (Wet) in gram per 5 stems")



