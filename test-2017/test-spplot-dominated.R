X = PCOR.list

colors <- c("firebrick1", "orange3", "chartreuse4", "blue", "purple")
meteVars <- c("Tmax", "Tavg", "Tmin", "Prec", "Srad")

x <- llply(X[4:9], function(x) {
  meteVars <- c("Tmax", "Tavg", "Tmin", "Prec", "Srad")
  Rc <- get_critical(0.1, 32 - 6)
  
  # x[x < Rc & x > -Rc] <- NA
  opt <- apply(x, 1, function(x) which.max(abs(x))) 
  opt[sapply(opt, length) == 0] <- NA
  opt <- factor(names(unlist(opt)), meteVars)
  return(opt)
}, .progress = "text") %>% do.call(cbind.data.frame, .) %>% set_colnames(paste0("M", colnames(.)))


gridClip@data <- x
zcol <- colnames(x)
A = 8; by=0.9; ntick = 2
p <- spplot(gridClip, col.regions = colors, 
            sp.layout = list(sp_cluster), as.table = T, layout = c(3, 2), 
            xlim = xlim, ylim = ylim, strip = F, 
       panel = function (x, y, z, subscripts, ...,  sp.layout) 
       {
         sppanel(list(sp.layout), panel.number(), first = TRUE)
         panel.levelplot.raster(x, y, z, subscripts, ..., interpolate = F)
         # panel.contourplot(x, y, z, subscripts, ..., contour = TRUE, labels = F)
         sppanel(list(sp.layout), panel.number(), first = FALSE)
         
         i <- panel.number()
         panel.text(76, 40.4, zcol[i], #english name: New_names[i])
                    fontfamily = "Times", cex = 1.3, font = 2, adj = 0)
         # panel.text(76.6, 40.8 +heights[i], metricsEN[i], fontfamily = "Times", cex = 1.2, font = 2, adj = 0)
         panel.addbarchart(z, subscripts, colors, showCluster = F, A = A, by = 0.9, ntick = ntick, ...)
       }, 
       par.settings = list(axis.line = list(col = "transparent")))

CairoPNG(filename = "test02.png", 10, 4, dpi = 300, pointsize = 12, units = "in")
# spplot(gridClip, as.table = T, layout = c(5, 12))
print(p)
dev.off()