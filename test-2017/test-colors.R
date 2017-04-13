colors <- hue_pal()(5)
colors[1] <- "firebrick1"
colors[4] <- "blue"#brewer.pal(9, "YlGnBu")[7]#
colors[5] <- "purple"
colors[2] <- "orange3"
colors[3] <- "green"

colors <- c("firebrick1", "orange2", "chartreuse4", "blue", "purple")#colors()[641]
meteVars <- c("Tmax", "Tavg", "Tmin", "Prec", "Srad")
opt <- apply(x, 1, function(x) which.max(abs(x))) 
opt[sapply(opt, length) == 0] <- NA; opt <- factor(names(unlist(opt)), meteVars)

gridClip@data <- data.frame(opt)

A = 8; by=0.9; ntick = 2
spplot(gridClip, col.regions = colors, sp.layout = list(sp_cluster),
       xlim = xlim, ylim = ylim, strip = F, 
       panel = function (x, y, z, subscripts, ...,  sp.layout) 
       {
         sppanel(list(sp.layout), panel.number(), first = TRUE)
         panel.levelplot.raster(x, y, z, subscripts, ..., interpolate = F)
         # panel.contourplot(x, y, z, subscripts, ..., contour = TRUE, labels = F)
         sppanel(list(sp.layout), panel.number(), first = FALSE)
         
         i <- panel.number()
         panel.text(76, 40.4, paste0("(",letters[i], ") ", colors[i]), #english name: New_names[i])
                    fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
         # panel.text(76.6, 40.8 +heights[i], metricsEN[i], fontfamily = "Times", cex = 1.2, font = 2, adj = 0)
         panel.addbarchart(z, subscripts, colors, showCluster = F, A = A, by = 0.9, ntick = ntick, ...)
       })