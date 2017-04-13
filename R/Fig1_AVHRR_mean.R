rm(list = ls())
# load("data/AVHRR_EOS_TP.rda")
load("data/phenology_list20170105.rda")
load("../data/00basemap.rda")
source('R/mainfunc/main.R', encoding = 'UTF-8')
# vars <- vars[c(13:15, 1:12)]
## ------------------ EOS MEAN ---------------------------

# ncol <- length(brks) - 1
# cols <- colorRampPalette(c("firebrick1","orange3", "darkgoldenrod2", "grey90", 
#                            brewer.pal(9, "YlGnBu")[c(4, 6, 7)]))(ncol)#,"firebrick1" 
# cols <- cols[-5]
# cols[ncol] <- "green4"; cols <- rev(cols)
EOS.mean  <- llply(EOS_list, function(x_df) apply(x_df, 1, mean, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)
EOS.sd  <- llply(EOS_list, function(x_df) apply(x_df, 1, sd, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)

SOS.mean  <- llply(SOS_list, function(x_df) apply(x_df, 1, mean, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)
SOS.sd  <- llply(SOS_list, function(x_df) apply(x_df, 1, sd, na.rm = TRUE), .progress = "text") %>% do.call(cbind.data.frame, .)

## --------------- 绘图 --------------------

brks_eos <- c(10, seq(240, 300, 10), 366)
EOS_TP.spatial(EOS.mean, brks=brks_eos, legendName = "DOY", saveFIG = TRUE, file="figures/Fig1_EOS_mean.png")

brks_sos <- c(10, seq(100, 160, 10), 366)
EOS_TP.spatial(SOS.mean, brks=brks_sos, legendName = "DOY", saveFIG = TRUE, file="figures/Fig1_SOS_mean.png")

## ------------------ EOS SD ----------------------------
brks <- c(0, c(5, 10, 20, 30), 100)
EOS_TP.spatial(EOS.sd, brks=brks, legendName = "SD", saveFIG = TRUE, file="figures/Fig1_EOS_sd.png")
EOS_TP.spatial(SOS.sd, brks=brks, legendName = "SD", saveFIG = TRUE, file="figures/Fig1_SOS_sd.png")

## ------------------- agricultural phenological data ----------------------------
## lattice 绘图需要预先设定好shape, size, color,
# spplot(gridClip[Id_clip, ], col = "red")