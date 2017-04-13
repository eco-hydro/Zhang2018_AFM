

VIP_df <- llply(VIP_list, `[[`, "Mean")
VIP <- llply(VIP_df, VIP_tidy)

time_mode <- seq(ymd(sprintf("2010%02d01", beginMonth)), ymd("2011-12-31"), by = "month") 
time_fmt2 <- time_mode %>% format("%b"); time_fmt2[1:6] %<>% paste0("_0")
time_fmt <- month(time_mode)
clusterEN <- c("EBF", "Alpine Meadows", "Alpine Steppe",
               "Alpine Shrub Meadows", "Desert", "DBF",
               "Temperate Meadows", "Temperate Steppe", "TP")
# clusterEN <- c("Evergreen Broadleaf Forests", "Alpine Meadows", "Alpine Steppe", 
#                "Alpine Shrub Meadows", "Desert", "Deciduous Broadleaf Forest", 
#                "Temperate Meadows", "Temperate Steppe", "TP")
## 绘图主函数，先以分区为循环


# for (i in seq_along(clusters)){
#   cluster <- clusters[i]
#   for (j in seq_along(methods)){
#     x <- pls_list[[i]][[j]]
#     pp <- vipCoef_plot(x, i, clusterEN[i], method)
#   }
# }
# pp <- foreach(cluster = clusters, i = icount(), .final = set_fnames(clusters)) %:%
#   foreach(method = methods, j = icount(), .final = set_fnames(methods)) %do% {
#     x <- pls_list[[i]][[j]]
#     vipCoef_plot(x, i, clusterEN[i], method)
#   }

pp <- foreach(cluster = clusters, vip = VIP, i = icount()) %do% {
  vipCoef_plot(vip, i, clusterEN[i], critical = 0.8)
}
library(gridExtra)
handle <- arrangeGrob(grobs = pp, nrow = 2, ncol = 5)

CairoPDF(file = sprintf("Fig7_pls%s.pdf", "mean"), width = 13.5, height = 8)
grid::grid.draw(handle)
dev.off()
