rm(list = ls())
agrValid_plot <- function(stations.sp, zcol=vars, 
                          .x = 0.72, .y = 0.05, 
                          main="bias", legend.cols=1){
  spplot(stations.sp, zcol= zcol, 
         col.regions = cols, #edge.col = "black",
         as.table = TRUE, strip = F, layout = c(4, 4), 
         xlim = xlim + c(1, -1)*0.1, ylim = ylim, 
         cuts = brks, pch = pch, cex = cex, 
         key=list(points=list(pch=pch, cex=cex, fill=cols, col = "transparent"),
                  text = list(levels(cut(1, brks))),
                  # between = 4,
                  title = main, 
                  font=2, cex=1.1, fontfamily='Times', columns=legend.cols),
         key.space = list(x = .x, y = .y, corner = c(0, 0)),
         sp.layout = list(sp_cluster), 
         panel = function(sp.layout, x, y, subscripts, col, ...){
           panel.pointsplot(sp.layout, x, y, subscripts, col, ...)
           i <- panel.number()
           panel.text(76, 40.4, paste0("(",letters[i], ") ", 
                                       vars[zcol][i]), #english name: New_names[i])
                      fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
         }, 
         par.settings = list(axis.line = list(col = "transparent")))
}
# library and load data
source('R/mainfunc/main.R', encoding = 'UTF-8')
load("data/AVHRR_EOS_TP.rda")
load("data/Fig1_agrdata2.rda")
load("data/00basemap.rda")
load("data/Fig1_infoValid.rda")
# vars <- vars[c(13:15, 1:12)]

# initial spplot parameters
brks <- c(-50, -30,  -15, 0, 15, 30, 120)
cex = seq(1.3, 0.8, length.out = 3) %>% c(., rev(.))
pch <- rep(c(25, 24), c(3, 3))
# pch = c(25, 25, 6, 2, 24, 24)
cols <- cols.fun(brks)
cols <- c("firebrick1" ,"orange3", "darkgoldenrod2", 
          brewer.pal(9, "YlGnBu")[c(6, 7)], "green4")
zcol = 1:ncol(infoValid[[1]])

stations.sp@data <- infoValid[[1]]
p <- list()
p[[1]] <- agrValid_plot(stations.sp, zcol, .x = 0.75, .y = 0.1, main="bias", legend.cols = 2)

brks <- c(0, 15, 30, 120)
cex = seq(0.8, 1.3, length.out = 3)
pch = rep(24, 3)
cols <- cols.fun(brks)
cols <- c(brewer.pal(9, "YlGnBu")[c(6, 7)], "green4")

stations.sp@data <- infoValid[[2]]
p[[2]] <- agrValid_plot(stations.sp, zcol, .x = 0.8, .y = 0.1, main="MAE")

stations.sp@data <- infoValid[[3]]
p[[3]] <- agrValid_plot(stations.sp, zcol, .x = 0.8, .y = 0.1, main="RMSE")

fnames <- paste0("figures/Fig2_validity_", c("bias", "MAE", "RMSE"), ".pdf")#[1]
for (i in seq_along(fnames)){
  CairoPDF(fnames[i], width = 11.56, height = 7.35)
  cat(fnames[i], sep = "\n")
  print(p[[i]])
  dev.off()
}

# width = 900; height = 520; ratio = 3.5
# fnames <- paste0("figures/Fig2_validity_", c("bias", "MAE", "RMSE"), ".png")
# for (i in seq_along(fnames)){
#   CairoPNG(filename = fnames[i], width*ratio, height*ratio, dpi = 250, pointsize = 12)
#   print(p[[i]])
#   dev.off()
# }

# a <- melt(infoValid)
# # linedf <- data.frame(yintercept = c(15, -15, 15, 15), variable = )
# ggplot(a, aes(x = variable, y=value, color = variable)) + 
#   geom_boxplot() + 
#   geom_jitter(width = 0.15) +
#   facet_grid(L1~., scales = "free_y") + 
#   geom_vline(xintercept = c(3.5, 9.5), color = "grey40", linetype=1) + 
#   geom_hline(yintercept = 15, color = "red", linetype=2) + 
#   geom_hline(data=data.frame(L1 = "bias"), yintercept = -15, color = "red", linetype=2) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         axis.title.y = element_blank(), 
#         strip.text = element_text(face="bold")) + 
#   guides(color = FALSE)

