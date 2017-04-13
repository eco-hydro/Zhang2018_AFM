nbrk <- 5
cols <- colorRampPalette(c("firebrick1" ,"orange3", "darkgoldenrod2", "grey90", 
                           brewer.pal(9, "YlGnBu")[c(4, 6, 7)]))(nbrk)
cols <- cols[-5]
cols[nbrk] <- "green4"; cols <- rev(cols)
cex = seq(1, 2.5, length.out = nbrk)#[x$val]

# coordinates(x) <- ~ long + lat
# proj4string(x) <- prj
# 
# spplot(x[, "val"], cex = cex, col.regions = cols, 
#        sp.layout = list(sp_text, sp_cluster), 
#        xlim = xlim, ylim = ylim, 
#        # lattice.options = ggplot2like.opts(),
#        par.settings = list(axis.line = list(col = NA),
#                            reference.line = list(col = "grey")))
## 采用arcgis进行绘图
if (class(df) != "SpatialPointsDataFrame"){
  coordinates(df) <- ~ long + lat
  proj4string(df) <- prj
}

x_mean <- subset(df, var == "mean")
x_sd <- subset(df, var == "sd")

p <- list()
p[[1]] <-  spplot(x_mean[, "val"], cex = cex, col.regions = cols, 
       sp.layout = list(sp_text, sp_cluster), 
       xlim = xlim, ylim = ylim, 
       # lattice.options = ggplot2like.opts(),
       par.settings = list(axis.line = list(col = NA),
                           reference.line = list(col = "grey")))
p[[2]] <- spplot(x_sd[, "val"], cex = cex, col.regions = cols, 
                 sp.layout = list(sp_text, sp_cluster), 
                 xlim = xlim, ylim = ylim, 
                 # lattice.options = ggplot2like.opts(),
                 par.settings = list(axis.line = list(col = NA),
                                     reference.line = list(col = "grey")))

print(c(p[[1]], p[[2]]))

library(gridExtra)
grid.arrange(p[[1]], p[[2]], ncol =2)

